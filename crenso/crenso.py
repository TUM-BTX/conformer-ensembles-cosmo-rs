#!/usr/bin/env python3
# coding: utf-8
"""CRENSO: automated conformer ensemble generation with CREST and CENSO.

This module drives a full conformer-generation workflow:

  1. CREST metadynamics sampling (optionally augmented with artificial
     dispersion and artificial charge runs to broaden the search),
  2. CREST re-ranking, sorting, GFN2 optimization and clustering,
  3. four sequential CENSO calls (prescreening, screening, per-solvent
     optimization, final screening) driven by the ``censo2rc*`` files,
  4. an optional xTB Hessian for the vibrational spectrum.

The result is ``CRENSOconf_final.xyz`` in the molecule's working directory.

External programs
----------------
``crest``, ``censo`` and ``xtb`` must be installed and on ``PATH``; ``molbar``
is optional and only needed for topology filtering (``--noreftopo``). Each can
be overridden with an environment variable:

    CRENSO_CREST_BIN, CRENSO_CENSO_BIN, CRENSO_XTB_BIN, CRENSO_MOLBAR_BIN

CENSO configuration
-------------------
The four ``censo2rc`` files are located via ``--censorc`` or the
``CRENSO_CENSORC`` environment variable; by default the ``config/censo2rc``
next to this file is used. The ``[paths]`` section of those files points at the
local ORCA and xTB executables and must be edited for your machine.

Usage
-----
    python crenso.py --smiles "CC(O)CO" --crest-solvent hexane --censo-solvents h2o
    python crenso.py --xyz molecule.xyz --gas-phase

Optional: set ``CHEMEO_API_KEY`` to enable Chemeo lookups when RDKit fails to
embed a SMILES string; without it the other MoleculeResolver services are used.
"""

import argparse
import os
from pathlib import Path
import sys
import json
import re
import contextlib
import subprocess
import time
import shutil
import tempfile
from dataclasses import dataclass, field, asdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


def _resolve_binary(env_var: str, name: str) -> str | None:
    """Locate an external program: ``$env_var`` first, then ``PATH``."""
    override = os.environ.get(env_var)
    if override:
        if os.path.isfile(override) and os.access(override, os.X_OK):
            return override
        print(f"Warning: {env_var}={override} is not an executable file.")
    return shutil.which(name)


# MolBar is optional: it is only used for topology filtering (--noreftopo).
MOLBAR_BIN = _resolve_binary("CRENSO_MOLBAR_BIN", "molbar")
MOLBAR_AVAILABLE = MOLBAR_BIN is not None
if not MOLBAR_AVAILABLE:
    print("Note: MolBar not found. Topology filtering (--noreftopo) will be disabled. "
          "Set CRENSO_MOLBAR_BIN to enable it.")

# Default location of the CENSO configuration files (config/ next to this file).
DEFAULT_CENSORC = Path(__file__).resolve().parent / "config" / "censo2rc"


def get_molbar_from_file(xyz_path: str, return_data: bool = False) -> str:
    """Call the molbar CLI to get the molecular barcode for an XYZ file."""
    result = subprocess.run(
        [MOLBAR_BIN, str(xyz_path)],
        capture_output=True, text=True, timeout=120
    )
    if result.returncode != 0:
        raise RuntimeError(f"MolBar failed for {xyz_path}: {result.stderr.strip()}")
    return result.stdout.strip()

# Global timing tracking
timings = {}
_current_context = None

def set_context(name):
    """Set the current top-level context (e.g. 'initial_sample_crest')."""
    global _current_context
    _current_context = name
    if name not in timings:
        timings[name] = {}

# ==================== Utility Functions ====================

def flatten_dict(d: dict) -> dict:
    """Flatten nested dictionaries into a single-level dict."""
    flat = {}
    for key, value in d.items():
        if isinstance(value, dict):
            flat.update(flatten_dict(value))
        else:
            flat[key] = value
    return flat

@contextlib.contextmanager
def cd(target: Path):
    prev = Path.cwd()
    target.mkdir(parents=True, exist_ok=True)
    os.chdir(target)
    try:
        yield
    finally:
        os.chdir(prev)

def run(cmd: list[str], output: Path, step: str = None):
    """
    Run a subprocess and measure its runtime.
    Results are stored in timings[_current_context][step].
    """
    env = os.environ.copy()
    start = time.time()
    with output.open("a") as f:
        subprocess.run(
            cmd,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
            env=env
        )
    elapsed = time.time() - start
    if _current_context:
        step_name = step or "unnamed_step"
        timings[_current_context][step_name] = timings[_current_context].get(step_name, 0.0) + elapsed
        print(f"[TIMER] {_current_context}:{step_name} took {elapsed:.2f} s")
    return elapsed

def cat(src: Path, dest: Path):
    with src.open("r") as fsrc, dest.open("a") as fdest:
        for line in fsrc:
            fdest.write(line)

def count_conf(file: Path) -> int:
    if not file.exists():
        return -1
    lines = file.read_text().splitlines()
    try:
        n_atoms = int(lines[0])
    except (ValueError, IndexError):
        return -1
    block = n_atoms + 2
    return len(lines) // block if len(lines) % block == 0 else -1

def get_atom_count(xyz_file: Path) -> int:
    """Read the number of atoms from the first line of an XYZ file."""
    try:
        return int(xyz_file.read_text().splitlines()[0].strip())
    except (ValueError, IndexError):
        return -1


def split_ensemble_xyz(ensemble_file: Path):
    """
    Split a multi-structure XYZ file into individual structures.
    Yields (header_comment, xyz_block_str) for each structure.
    """
    lines = ensemble_file.read_text().splitlines()
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        comment = lines[i + 1] if (i + 1) < len(lines) else ""
        block_lines = lines[i:i + n_atoms + 2]
        yield comment, "\n".join(block_lines) + "\n"
        i += n_atoms + 2


def parse_molbar_output(molbar_string: str) -> dict:
    """
    Parse MolBar output string.
    Format: MolBar | version | formula | charge | topology | heavy_atom | topography | chirality
    """
    parts = [p.strip() for p in molbar_string.split('|')]
    if len(parts) < 7:
        return {}
    return {
        'topology_spectrum': parts[4],
        'heavy_atom_topology_spectrum': parts[5],
    }


def filter_ensemble_by_molbar(ensemble_file: Path, reference_xyz: Path) -> int:
    """
    Filter an ensemble XYZ file in-place, keeping only conformers whose
    topology_spectrum AND heavy_atom_topology_spectrum match the reference.

    Returns the number of conformers kept.
    """
    if not MOLBAR_AVAILABLE:
        print("WARNING: MolBar not available, skipping topology filtering.")
        return count_conf(ensemble_file)

    # Get reference MolBar
    try:
        ref_molbar = get_molbar_from_file(str(reference_xyz), return_data=False)
        ref_data = parse_molbar_output(ref_molbar)
        ref_topo = ref_data.get('topology_spectrum', '')
        ref_heavy = ref_data.get('heavy_atom_topology_spectrum', '')
        print(f"  Reference topology: {ref_topo[:60]}...")
        print(f"  Reference heavy atom topology: {ref_heavy[:60]}...")
    except Exception as e:
        print(f"WARNING: Failed to get MolBar for reference {reference_xyz}: {e}")
        print("Skipping topology filtering.")
        return count_conf(ensemble_file)

    if not ref_topo or not ref_heavy:
        print("WARNING: Empty reference MolBar spectra. Skipping filtering.")
        return count_conf(ensemble_file)

    # Filter each conformer
    kept = []
    removed = 0
    for comment, block in split_ensemble_xyz(ensemble_file):
        # Write conformer to temp file for MolBar
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as tmp:
            tmp.write(block)
            tmp_path = tmp.name
        try:
            conf_molbar = get_molbar_from_file(tmp_path, return_data=False)
            conf_data = parse_molbar_output(conf_molbar)
            conf_topo = conf_data.get('topology_spectrum', '')
            conf_heavy = conf_data.get('heavy_atom_topology_spectrum', '')

            if conf_topo == ref_topo and conf_heavy == ref_heavy:
                kept.append(block)
            else:
                removed += 1
        except Exception as e:
            print(f"  WARNING: MolBar failed for a conformer, removing it: {e}")
            removed += 1
        finally:
            os.unlink(tmp_path)

    print(f"  MolBar filtering: kept {len(kept)}, removed {removed} (topology changed)")

    # Overwrite ensemble file with filtered conformers
    ensemble_file.write_text("".join(kept))
    return len(kept)


def clean_folder(folder: Path):
    """
    Remove everything in `folder` except the original .xyz file,
    CRENSOconf_final.xyz, censo_finalsampling.out, and censo.out.
    """
    keep = {
        f"{folder.name}.xyz",
        f"{folder.name}_log.txt",
        f"{folder.name}_info.json",    # molecule info
        "CRENSOconf_final.xyz",
        "censo_finalsampling.out",
        "censo.out",
        "timings.json",
        "timings_solvated.json",       # Add this
        "timings_gas.json",            # Add this
        "vibspectrum",
        "gas_phase"
    }

    print(f"Cleaning folder: {folder}")

    for item in folder.iterdir():
        if item.name not in keep:
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
    

@dataclass
class CrensoConfig:
    """Settings for a CRENSO run.

    The defaults are the ones used in the accompanying paper. ``crest_solvent``
    is the ALPB solvent driving the CREST search; ``censo_solvents`` is the list
    of CPCM solvents for the CENSO optimization step, each run separately and
    then merged.
    """

    # CREST sampling
    crest_solvent: str = "hexane"
    ewin: float = 6
    mdlen: str = "x1.0"
    mdlen_art: str = "x2.0"
    nclust: str = "500"
    artificial_disp: list = field(default_factory=lambda: [0.5, 0.5, 1.5, 2.0])
    artificial_chrg: list = field(default_factory=lambda: [1, 2, -1, -2])

    # CENSO refinement
    censo_solvents: list = field(default_factory=lambda: ["h2o"])
    censorc_path: str = None
    uhf: int = 0

    # Parallelisation: P tasks x O threads each
    P: int = 1
    O: int = 1
    maxcores: int = 1

    def __post_init__(self):
        if self.censorc_path is None:
            self.censorc_path = os.environ.get("CRENSO_CENSORC", str(DEFAULT_CENSORC))
        self.censorc_path = str(Path(self.censorc_path).expanduser().resolve())
        missing = [
            f"{self.censorc_path}{suffix}"
            for suffix in ("", "_part1", "_part2", "_part1_sp")
            if not Path(f"{self.censorc_path}{suffix}").is_file()
        ]
        if missing:
            raise FileNotFoundError(
                "Missing CENSO configuration file(s): " + ", ".join(missing)
                + "\nPoint --censorc / CRENSO_CENSORC at the 'censo2rc' file; the "
                  "_part1, _part2 and _part1_sp variants must sit next to it."
            )
        if isinstance(self.censo_solvents, str):
            self.censo_solvents = [self.censo_solvents]
        for name, value in (("P", self.P), ("O", self.O), ("maxcores", self.maxcores)):
            if int(value) < 1:
                raise ValueError(f"{name} must be >= 1, got {value}")

    @classmethod
    def from_dict(cls, configuration: dict) -> "CrensoConfig":
        """Build a config from a (possibly nested) dictionary.

        Accepts the legacy nested layout, e.g.
        ``{"crest_config": {...}, "O": 4, "censorc_path": "..."}``.
        Unknown keys are ignored with a warning so that older scripts keep working.
        """
        flat = flatten_dict(configuration)
        known = {f for f in cls.__dataclass_fields__}
        unknown = set(flat) - known
        if unknown:
            print(f"Note: ignoring unrecognised configuration key(s): {sorted(unknown)}")
        return cls(**{k: v for k, v in flat.items() if k in known})

    def to_dict(self) -> dict:
        return asdict(self)


class CRENSOgen:
    crest_bin = _resolve_binary("CRENSO_CREST_BIN", "crest")
    censo_bin = _resolve_binary("CRENSO_CENSO_BIN", "censo")
    xtb_bin = _resolve_binary("CRENSO_XTB_BIN", "xtb")

    def __init__(self, smiles: str, configuration,
                base_dir: Path = None,
                    xyz: str = None, nci_mode: bool = False, folder_name: str = None,
                    gas_phase: bool = False, noreftopo: bool = False,
                    xtb_preopt: bool = False, gentle_sampling: int = 0):
        self.smiles = smiles
        cfg = (configuration if isinstance(configuration, CrensoConfig)
               else CrensoConfig.from_dict(configuration))
        self.config = cfg
        self.__dict__.update(cfg.to_dict())
        for tool, path in (("crest", self.crest_bin), ("censo", self.censo_bin),
                           ("xtb", self.xtb_bin)):
            if path is None:
                raise RuntimeError(
                    f"Required program '{tool}' not found on PATH. Install it or set "
                    f"CRENSO_{tool.upper()}_BIN to its location."
                )
        self.nclust2 = 2 * int(self.nclust)
        self.Threads = self.P * self.O
        self.charge = 0
        self.xyz = xyz
        self.nci_mode = nci_mode
        self.gas_phase = gas_phase
        self.noreftopo = noreftopo
        self.xtb_preopt = xtb_preopt
        self.gentle_sampling = gentle_sampling

        # Prepare molecule formula and folders
        if folder_name:
            # Use provided folder name
            self.folder_name = folder_name
        elif xyz is None:
            # Generate from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string provided.")
            # Naming the folder after the IUPAC name is a convenience that needs
            # network access; fall back to the formula alone when unavailable.
            iupac_name = "unknown"
            try:
                from pubchempy import get_compounds
                compounds = get_compounds(smiles, "smiles")
                if compounds and compounds[0].iupac_name:
                    iupac_name = compounds[0].iupac_name
            except Exception as exc:
                print(f"Note: could not look up an IUPAC name ({exc}); using formula only.")

            self.mol = Chem.AddHs(mol)
            self.folder_name = CalcMolFormula(self.mol) + "_" + iupac_name
        else:
            # Use XYZ filename stem
            self.folder_name = Path(xyz).stem

        # Ensure base output directory exists and is ABSOLUTE
        self.base_dir = Path(base_dir if base_dir is not None else Path.cwd()).resolve()
        self.base_dir.mkdir(parents=True, exist_ok=True)

        # Create and set molecule-specific folder (ABSOLUTE PATH)
        self.folder_path = self.base_dir / self.folder_name
        self.folder_path.mkdir(parents=True, exist_ok=True)
        print(f"Working in {self.folder_path}")

    def get_xyz_from_smiles(self) -> Path:
        mol = Chem.AddHs(Chem.MolFromSmiles(self.smiles))
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            print(f"WARNING: Embedding failed for SMILES: {self.smiles}")
            print("Attempting to standardize SMILES via MoleculeResolver...")
            try:
                from moleculeresolver import MoleculeResolver
                # Chemeo needs a personal API token; supply it via the environment
                # rather than hardcoding it. Without it the remaining services are used.
                api_keys = {}
                if os.environ.get("CHEMEO_API_KEY"):
                    api_keys["chemeo"] = os.environ["CHEMEO_API_KEY"]
                with MoleculeResolver(available_service_API_keys=api_keys) as mr:
                    std_smiles = mr.standardize_SMILES(self.smiles)
                print(f"Standardized SMILES: {std_smiles}")
                mol = Chem.AddHs(Chem.MolFromSmiles(std_smiles))
                result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                if result != 0:
                    raise RuntimeError(f"Embedding failed even with standardized SMILES: {std_smiles}")
            except ImportError:
                raise RuntimeError("Embedding failed and moleculeresolver is not installed.")
        AllChem.UFFOptimizeMolecule(mol)
        conf = mol.GetConformer()

        lines = [str(mol.GetNumAtoms()), ""]
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            lines.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
        
        xyz_file = self.folder_path / "initial_rdkit_struct.xyz"
        xyz_file.write_text("\n".join(lines))
        print(f"XYZ file written to {xyz_file}")
        
        return xyz_file
    
    def is_monoatomic(self) -> bool:
        """Check if the molecule is monoatomic (single atom like Rn, Ar, He)."""
        if self.xyz:
            xyz_path = Path(self.xyz)
            if xyz_path.exists():
                return get_atom_count(xyz_path) == 1
        if self.smiles:
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is not None:
                mol_h = Chem.AddHs(mol)
                return mol_h.GetNumAtoms() == 1
        return False

    def resolve_input_xyz(self) -> str:
        """
        Resolve the input XYZ file, returning the local filename.
        If SMILES is provided, generates XYZ via RDKit.
        This is shared logic used by both the monoatomic path and initial_sample_crest.
        """
        with cd(self.folder_path):
            if self.xyz is not None:
                xyz_path = Path(self.xyz)
                if xyz_path.is_absolute() and xyz_path.exists():
                    local_xyz = Path("input_structure.xyz")
                    shutil.copy2(xyz_path, local_xyz)
                    return str(local_xyz)
                elif Path(self.xyz).exists():
                    return str(self.xyz)
                else:
                    print(f"WARNING: Provided xyz path does not exist: {self.xyz}")
                    print("Generating XYZ from SMILES instead...")
                    self.get_xyz_from_smiles()
                    return "initial_rdkit_struct.xyz"
            else:
                self.get_xyz_from_smiles()
                return "initial_rdkit_struct.xyz"

    def run_monoatomic_sp(self):
        """
        Fast-path for monoatomic molecules: just run a GFN2 single-point
        and create CRENSOconf_final.xyz from the input structure.
        """
        set_context("monoatomic_sp")
        with cd(self.folder_path):
            input_xyz = self.resolve_input_xyz()
            print(f"Monoatomic molecule detected. Running single-point on {input_xyz}")

            xtb_args = [
                self.xtb_bin,
                input_xyz,
                "--gfn", "2",
                *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                "--chrg", str(self.charge),
                "--json"
            ]
            run(xtb_args, Path("xtb_sp.out"), step="xtb_sp")

            # Create CRENSOconf_final.xyz from the input
            shutil.copy2(input_xyz, "CRENSOconf_final.xyz")
            print("Created CRENSOconf_final.xyz from input structure.")

            sub_total = sum(v for k, v in timings["monoatomic_sp"].items())
            timings["monoatomic_sp"]["total"] = sub_total

    def initial_sample_crest(self):
        set_context("CREST_initial_sample_crest")
        with cd(self.folder_path):
            if self.xyz != None:
                # If an XYZ file is provided, use it instead of generating one from SMILES
                xyz_path = Path(self.xyz)

                # Check if the provided path is absolute and exists
                if xyz_path.is_absolute() and xyz_path.exists():
                    # Copy the XYZ file to the working directory with a standard name
                    local_xyz = Path("input_structure.xyz")
                    shutil.copy2(xyz_path, local_xyz)
                    input_xyz = str(local_xyz)
                    print(f"Using provided XYZ file: {xyz_path} (copied to {input_xyz})")
                elif Path(self.xyz).exists():
                    # Relative path that exists
                    input_xyz = str(self.xyz)
                    print(f"Using provided XYZ file: {input_xyz}")
                else:
                    # Path doesn't exist - this might be a SMILES string mistakenly passed as xyz
                    print(f"WARNING: Provided xyz path does not exist: {self.xyz}")
                    print("This appears to be a SMILES string. Generating XYZ from SMILES instead...")
                    xyz = self.get_xyz_from_smiles()
                    print(f"XYZ generated from SMILES at {xyz}")
                    input_xyz = "initial_rdkit_struct.xyz"
            else:
                print("Generating XYZ from SMILES...")
                xyz = self.get_xyz_from_smiles()
                print(f"XYZ generated from SMILES at {xyz}")
                input_xyz = "initial_rdkit_struct.xyz"

            # Optional GFN-FF loose preoptimization of the initial structure
            if self.xtb_preopt:
                print("Running xTB GFN-FF loose preoptimization...")
                xtb_args = [
                    self.xtb_bin,
                    input_xyz,
                    "--gfnff",
                    *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                    "--opt", "loose",
                    "--chrg", str(self.charge),
                ]
                run(xtb_args, Path("xtb_pre_opt.out"), step="xtb_preopt")
                if Path("xtbopt.xyz").exists():
                    input_xyz = "xtbopt.xyz"
                    print(f"Preoptimization done. Using {input_xyz} for CREST.")
                else:
                    print("WARNING: xtbopt.xyz not found after preopt. Using original structure.")

            if self.noreftopo:
                print("NOTE: --noreftopo mode enabled. Will post-filter conformers with MolBar.")

            if self.gentle_sampling == 1:
                print("NOTE: --gentle 1: Constraining bonds (--cbonds) and reducing timestep.")
                gentle_flags = ["--cbonds", "--tstep", "3"]
            elif self.gentle_sampling == 2:
                print("NOTE: --gentle 2: Strong heavy-atom constraints, shorter MD, gentler metadynamics.")
                gentle_flags = ["--cheavy", "0.1", "--tstep", "3", "--len", "x0.5", "--vbdump", "2.0"]
            elif self.gentle_sampling == 3:
                print("NOTE: --gentle 3: Maximum rigidity + reduced search thoroughness.")
                gentle_flags = ["--cheavy", "0.25", "--tstep", "2", "--len", "x0.5", "--vbdump", "2.0", "--mrest", "2", "--squick"]
            elif self.gentle_sampling >= 4:
                print("NOTE: --gentle 4: Ultra-rigid — max force constant, minimal search, no restarts.")
                gentle_flags = ["--cheavy", "0.5", "--tstep", "1.5", "--len", "x0.3", "--vbdump", "3.0", "--mrest", "1", "--mquick"]
            else:
                gentle_flags = []

            if self.nci_mode == False:
                # base CREST args including mdlen
                crest_args = [
                    self.crest_bin,
                    input_xyz,
                    *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                    "-ewin", str(self.ewin),
                    "--gfnff", "-v4", "-norotmd", "-norestart", "-nocross",
                    "--chrg", str(self.charge),
                    "-T", str(self.Threads),
                    "-mdlen", self.mdlen,
                    *(["--noreftopo"] if self.noreftopo else []),
                    *gentle_flags
                ]

            else:
                # CREST args for large complexes using NCI mode
                print("Setting up CREST arguments for NCI mode...")
                crest_args = [
                    self.crest_bin,
                    input_xyz,
                    "--nci",
                    *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                    "--chrg", str(self.charge),
                    "-gfnff",
                    *(["--noreftopo"] if self.noreftopo else [])
                    #"--quick"          ##uncomment if it takes too much
                ]

                # clustering args. Cluster the ensemble for most representative
                #structures. Based on dihedral angles.
            cluster_args = [
                self.crest_bin,
                "--for", "crest_conformers.xyz",
                #"--cluster", str(self.nclust2)
                "--cluster", "normal"
            ]

            # 1. normal sampling
            if self.nci_mode == False:
                print("Running CREST (normal):", crest_args)
            elif self.nci_mode == True:
                print("Running CREST (NCI mode):", crest_args)

            run(crest_args, Path("crest_normal.out"), step="crest")

            # MolBar filtering: remove conformers with changed topology
            if self.noreftopo and Path("crest_conformers.xyz").exists():
                print("Running MolBar topology filter on crest_conformers.xyz...")
                filter_ensemble_by_molbar(Path("crest_conformers.xyz"), Path(input_xyz))

            run(cluster_args, Path("crest_normal_clust.out"), step="clustering")
            cat(Path("crest_clustered.xyz"), Path("search_gfnff.xyz"))
            if count_conf(Path("search_gfnff.xyz")) >= int(self.nclust2):       ## ALL GOOD TIL HERE
                return
            
            # 2. artificial dispersion sampling
            for disp in self.artificial_disp:
                try:
                    crest_args = [
                        self.crest_bin,
                        input_xyz,
                        "-mdlen", self.mdlen_art,           ## longer MD for artificial dispersion
                        "-dispscal", str(disp),
                        *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                        "-ewin", str(self.ewin),
                        "--gfnff", "-v4", "-norotmd", "-norestart", "-nocross",
                        "--chrg", str(self.charge),
                        "-T", str(self.Threads),
                        *(["--noreftopo"] if self.noreftopo else []),
                        *gentle_flags
                    ]
                    print(f"Running CREST (with artificial dispersion): {disp}")
                    run(crest_args,  Path("crest_disp.out"), step=f"crest_artdisp_{disp}")

                    if self.noreftopo and Path("crest_conformers.xyz").exists():
                        print(f"Running MolBar topology filter (disp={disp})...")
                        filter_ensemble_by_molbar(Path("crest_conformers.xyz"), Path(input_xyz))

                    run(cluster_args, Path("crest_disp_clust.out"), step=f"clustering_disp_{disp}")
                    cat(Path("crest_clustered.xyz"), Path("search_gfnff.xyz"))
                    if count_conf(Path("search_gfnff.xyz")) >= int(self.nclust2):
                        print(f"Reached target number of conformers ({self.nclust2}). Stopping further sampling.")
                        return
                except subprocess.CalledProcessError as e:
                    print(f"WARNING: CREST with artificial dispersion {disp} failed (possibly due to topology changes).")
                    print(f"Error: {e}")
                    print("Continuing with structures already collected in search_gfnff.xyz...")
                    continue
                except Exception as e:
                    print(f"WARNING: Unexpected error during artificial dispersion sampling with {disp}: {e}")
                    print("Continuing with structures already collected in search_gfnff.xyz...")
                    continue 

            # 3. artificial charge sampling
            for chrg in self.artificial_chrg:
                try:
                    crest_args = [
                        self.crest_bin,
                        input_xyz,
                        "-mdlen", self.mdlen_art,        ## longer MD for artificial charge
                        *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                        "-ewin", str(self.ewin),
                        "--gfnff", "-v4", "-norotmd", "-norestart", "-nocross",
                        "--chrg", str(chrg),
                        "-T", str(self.Threads),
                        *(["--noreftopo"] if self.noreftopo else []),
                        *gentle_flags
                    ]
                    print(f"Running CREST (with artificial charge): {chrg}")
                    run(crest_args,  Path("crest_chrg.out"), step=f"crest_artchrg_{chrg}")

                    if self.noreftopo and Path("crest_conformers.xyz").exists():
                        print(f"Running MolBar topology filter (chrg={chrg})...")
                        filter_ensemble_by_molbar(Path("crest_conformers.xyz"), Path(input_xyz))

                    run(cluster_args, Path("crest_chrg_clust.out"), step=f"clustering_chrg_{chrg}")
                    cat(Path("crest_clustered.xyz"), Path("search_gfnff.xyz"))
                    print(f"Completed CREST sample and cluster with artificial charge {chrg}.")
                    if count_conf(Path("search_gfnff.xyz")) >= int(self.nclust2):
                        print(f"Reached target number of conformers ({self.nclust2}). Stopping further sampling.")
                        return
                    #Overwrite the charge file with the original charge
                    chrg_file=open(".CHRG","w")
                    chrg_file.write(str(self.charge))
                    chrg_file.close()
                except subprocess.CalledProcessError as e:
                    print(f"WARNING: CREST with artificial charge {chrg} failed (possibly due to topology changes).")
                    print(f"Error: {e}")
                    print("Continuing with structures already collected in search_gfnff.xyz...")
                    # Restore original charge file even on failure
                    chrg_file=open(".CHRG","w")
                    chrg_file.write(str(self.charge))
                    chrg_file.close()
                    continue
                except Exception as e:
                    print(f"WARNING: Unexpected error during artificial charge sampling with {chrg}: {e}")
                    print("Continuing with structures already collected in search_gfnff.xyz...")
                    # Restore original charge file even on failure
                    chrg_file=open(".CHRG","w")
                    chrg_file.write(str(self.charge))
                    chrg_file.close()
                    continue

            sub_total = sum(v for k,v in timings["CREST_initial_sample_crest"].items())
            timings["CREST_initial_sample_crest"]["total"] = sub_total

    def screen(self):
        set_context("CREST_screen")
        with cd(self.folder_path):
            print("----Starting CREST Screening----")
            #1. energy re-rank (at FF level) because of +/-disp/ES structures (output: crest_ensemble.xyz)
            crest_args = [
                self.crest_bin,
                "-mdopt", "search_gfnff.xyz", # optimize each point on a given trajectory or ensemble file with GFNn–xTB
                *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                "-opt", "normal",
                "-ewin", str(self.ewin),
                "-gfnff",
                "--chrg", str(self.charge),
                "-T", str(self.Threads),
                "-mdlen", self.mdlen
            ]
            print("Running energy re-rank:", crest_args)
            run(crest_args, Path("crest_rerank.out"), step="crest_screening")

            #2. sorting the optimized ensemble (output: overwritten crest_ensemble.xyz)
            sort_args = [
                self.crest_bin,
                "-cregen", "crest_ensemble.xyz",
                "-ewin", str(self.ewin),
                "-ethr", "0.05",                 # ensemble sorting: energy threshold between conformer pairs, kcal/mol
                "-rthr", "0.125",                 # ensemble sorting: RMSD threshold
                "-bthr", "0.01",
                "--chrg", str(self.charge)
            ]
            print("Running sorting:", sort_args)
            run(sort_args, Path("crest_rerank.out"), step="crest_sorting")

            #3. GFN2 full opt (output: overwritten crest_ensemble.xyz)
            opt_args = [
                self.crest_bin,
                "-screen", "crest_ensemble.xyz",        ##opt + cregen. Why not using this in the previous?
                *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
                "-ewin", str(self.ewin),
                "-opt", "vtight",
                "-gfn2",
                "--chrg", str(self.charge)
            ]
            print("Running GFN2 full opt:", opt_args)
            run(opt_args, Path("crest_rerank.out"), step="crest_optimization")

            #4. PCA/k-Means clustering of final ensemble (output: crest_clustered.xyz)

            cluster_args = [
                self.crest_bin,
                "-for", "crest_ensemble.xyz",
                #"-cluster", str(self.nclust),
                "-cluster", "normal",
                "-ewin", str(self.ewin),
                "--chrg", str(self.charge)
            ]
            print("Running PCA/k-means clustering:", cluster_args)
            run(cluster_args, Path("crest_rerank.out"), step="crest_clustering")

            sub_total = sum(v for k,v in timings["CREST_screen"].items())
            timings["CREST_screen"]["total"] = sub_total

    def create_ensemble(self, steps_to_run="all"):
        """
        Run ensemble creation in parts. `steps_to_run` can be:
          - "all"        to run ["part01", "part2", "sp"]
          - one of        "part01", "part2", "sp"
          - a list of    ["part01","sp"] etc.
        """
        set_context("CENSO_refine_ensemble")
        print("Raw steps_to_run:", steps_to_run)
        if isinstance(steps_to_run, str):
            steps = [steps_to_run]
        else:
            steps = list(steps_to_run)

        # expand “all” into the full set of steps
        if "all" in steps:
            steps = ["part01", "part2", "sp"]

        print("Normalized steps:", steps)

        with cd(self.folder_path):
            # 1. CENSO part0/1 for sampling solvent
            if "part01" in steps:
                if not (Path("0_PRESCREENING.out").exists() and Path("1_SCREENING.out").exists()):
                    censo_args = [
                        self.censo_bin,
                        "--input", "crest_clustered.xyz",
                        "--inprc", self.censorc_path,
                        "--charge", str(self.charge),
                        "--solvent", self.crest_solvent,
                        "--omp", str(self.O),
                        "--maxcores", str(self.maxcores),
                        "-u", str(self.uhf),
                        "--bhess", "--sm-rrho", "alpb",
                        *(["--gas-phase"] if self.gas_phase else [])
                    ]
                    print("Running CENSO part0 (CREST solvent):", censo_args)
                    run(censo_args, Path("censo.out"), step="censo_part0")

                    censo_args = [
                        self.censo_bin,
                        "--input", "0_PRESCREENING.xyz",
                        "--inprc", f"{self.censorc_path}_part1",
                        "--charge", str(self.charge),
                        "--solvent", self.crest_solvent,
                        "--omp", str(self.O),
                        "--maxcores", str(self.maxcores),
                        "-u", str(self.uhf),
                        "--bhess", "--sm-rrho", "alpb",
                        *(["--gas-phase"] if self.gas_phase else [])
                    ]
                    print("Running CENSO part1 (CENSO solvent):", censo_args)
                    run(censo_args, Path("censo.out"), step="censo_part1")
                else:
                    print("Skipping part01 (already done).")
            else:
                print("Skipping part01.")

            # 2. CENSO part2 across solvents
            if "part2" in steps:
                base = [
                    self.censo_bin,
                    "--input", "../1_SCREENING.xyz",
                    "--inprc", f"{self.censorc_path}_part2",
                    "--charge", str(self.charge),
                    "-u", str(self.uhf),
                    "--omp", str(self.O),
                    "--maxcores", str(self.maxcores),
                    #"--bhess", "--sm-rrho", "alpb"
                ]
                if self.gas_phase:
                    sd = self.folder_path/"gas"
                    sd.mkdir(exist_ok=True, parents=True)
                    with cd(sd):
                        args = base + ["--gas-phase"]
                        print("Running part2 in gas phase.")
                        run(args, Path(f"censo_gas.out"), step=f"censo_part2_gas")
                        cat(Path("2_OPTIMIZATION.xyz"), Path("../CRENSOconf_ensemble_mult.xyz"))
                else:
                    for solv in self.censo_solvents:
                        sd = self.folder_path/solv      # create a subdir for each solvent
                        sd.mkdir(exist_ok=True, parents=True)
                        with cd(sd):
                            args = base + ["--solvent", solv]
                            print(f"Running part2 for {solv}:", args)
                            run(args, Path(f"censo_{solv}.out"), step=f"censo_part2_{solv}")
                            cat(Path("2_OPTIMIZATION.xyz"), Path("../CRENSOconf_ensemble_mult.xyz"))
            else:
                print("Skipping part2.")

            # 3. Single‐point & clustering
            if "sp" in steps:
                # part1/SP
                sp_args = [
                    self.censo_bin,
                    "--input", "CRENSOconf_ensemble_mult.xyz",
                    "--inprc", f"{self.censorc_path}_part1_sp",
                    "--charge", str(self.charge),
                    "-u", str(self.uhf),
                    "--omp", str(self.O),
                    "--maxcores", str(self.maxcores),
                    "--bhess", "--sm-rrho", "alpb",
                    "--solvent",  self.crest_solvent,
                    *(["--gas-phase"] if self.gas_phase else [])
                ]
                print("Running CENSO part1/SP:", sp_args)
                run(sp_args, Path("censo_finalsampling.out"), step="censo_sp")
                shutil.copy2("1_SCREENING.xyz", "CRENSOconf_ensemble_SP.xyz")

                # append gtot into the XYZ
                with open("1_SCREENING.json") as f:
                    results = json.load(f)["results"]
                with open("CRENSOconf_ensemble_SP.xyz") as fin, open("CRENSOconf_ensemble_SP_with_gtot.xyz","w") as fout:
                    while True:
                        num = fin.readline()
                        if not num: break
                        fout.write(num)
                        comment = fin.readline()
                        m = re.search(r"(CONF\d+)", comment)
                        fout.write(f"{results.get(m.group(1),{}).get('gtot', '')}\n" if m else comment)
                        n = int(num.split()[0])
                        for _ in range(n):
                            fout.write(fin.readline())

                # CREST clustering
                cl_args = [
                    self.crest_bin,
                    "CRENSOconf_ensemble_SP_with_gtot.xyz",
                    "--cregen", "CRENSOconf_ensemble_SP_with_gtot.xyz",
                    "--notopo", "-ethr","0.1","-rthr","0.2","-bthr","0.03",
                    "-ewin",str(self.ewin), "--chrg",str(self.charge),
                    "-cluster",str(self.nclust),
                ]
                print("Running clustering:", cl_args)
                run(cl_args, Path("crest_rerank.out"), step="crest_clustering_final")
                os.rename("crest_clustered.xyz","CRENSOconf_final.xyz")
            else:
                print("Skipping SP + clustering.")

            sub_total = sum(v for k,v in timings["CENSO_refine_ensemble"].items())
            timings["CENSO_refine_ensemble"]["total"] = sub_total

    def compute_vibrational_spectrum(self):
        """
        Runs an xTB Hessian calculation on CRENSOconf_final.xyz,
        writes out the vibrational spectrum to 'vibspectrum', and
        returns a dict with thermodynamic and summary energies.
        """
        from pathlib import Path
        import re

        set_context("xTB_compute_vibrational_spectrum")
        with cd(self.folder_path):
            # 1) build and run xTB Hessian
            xtb_args = [
                self.xtb_bin,
                "CRENSOconf_final.xyz",
                "--bhess",
		         *(["--alpb", self.crest_solvent] if not self.gas_phase else []),
		        "--gfn", "2",
                "--json"
            ]
            print("Running xTB Hessian for vibrational spectrum:", xtb_args)
            run(xtb_args, Path("xtb_hessian.out"), step="xtb_hessian")
            sub_total = sum(v for k,v in timings["xTB_compute_vibrational_spectrum"].items())
            timings["xTB_compute_vibrational_spectrum"]["total"] = sub_total
            # 2) ensure vibspectrum was produced
            if not Path("vibspectrum").exists():
                raise FileNotFoundError("xTB did not produce 'vibspectrum'!")
            print("Done. Your vibrational frequencies + IR intensities are in 'vibspectrum'")

            # 3) parse xtb_hessian.out
            out_lines = Path("xtb_hessian.out").read_text().splitlines()
            thermodynamic = {}
            summary       = {}

            # --- parse the THERMODYNAMIC block ---
            for i, L in enumerate(out_lines):
                if "THERMODYNAMIC" in L:
                    # from the next line until the optimized-geometry message
                    for line in out_lines[i+1:]:
                        if line.strip().startswith("optimized geometry written"):
                            break
                        m = re.search(r"::\s*(.+?)\s+(-?\d+\.\d+)\s+Eh", line)
                        if m:
                            raw_key, val = m.group(1), m.group(2)
                            # sanitize key into a pythonic form
                            key = (raw_key.lower()
                                        .replace(" ", "_")
                                        .replace("(", "")
                                        .replace(")", "")
                                        .replace(".", "")
                                        .replace("/", "_"))
                            thermodynamic[key] = float(val)
                    break

            # --- parse the summary box ---
            for i, L in enumerate(out_lines):
                # look for a dashed line followed by "TOTAL ENERGY"
                if re.match(r"^\s*-{5,}", L) and i+1 < len(out_lines) and "TOTAL ENERGY" in out_lines[i+1]:
                    # read until the next dashed line
                    for line in out_lines[i+1:]:
                        if re.match(r"^\s*-{5,}", line):
                            break
                        m = re.search(r"\|\s*([A-Z \-]+?)\s+(-?\d+\.\d+)\s+[A-Za-z\/α]+", line)
                        if m:
                            raw_key, val = m.group(1), m.group(2)
                            key = raw_key.strip().lower().replace(" ", "_").replace("-", "_")
                            summary[key] = float(val)
                    break

            # 4) return everything
            return {
                "thermodynamic": thermodynamic,
                "summary": summary
            }
            


# ==================== Command-line interface ====================

def _save_timings(folder_path: Path, filename: str = "timings.json") -> Path:
    """Add totals and percentages to the global `timings` dict and write it out."""
    total_time = sum(
        v["total"] for v in timings.values()
        if isinstance(v, dict) and "total" in v
    )
    timings["total"] = total_time
    for data in timings.values():
        if isinstance(data, dict) and "total" in data and total_time:
            data["percentages"] = {
                sub: round((t / total_time) * 100, 2)
                for sub, t in data.items()
                if sub != "total" and not isinstance(t, dict)
            }
            data["percentages"]["total"] = round((data["total"] / total_time) * 100, 2)

    timing_file = folder_path / filename
    timing_file.write_text(json.dumps(timings, indent=2))
    print(f"Saved timings to {timing_file}")
    return timing_file


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="crenso",
        description="Generate a conformer ensemble with CREST and CENSO.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Either --smiles or --xyz is required. See the module docstring "
               "for the environment variables controlling program locations.",
    )
    src = p.add_argument_group("input")
    src.add_argument("--smiles", help="SMILES string of the molecule.")
    src.add_argument("--xyz", help="Path to a starting XYZ structure.")
    src.add_argument("--charge", type=int, default=0, help="Molecular charge.")
    src.add_argument("--uhf", type=int, default=0,
                     help="Number of unpaired electrons.")
    src.add_argument("--folder-name",
                     help="Name of the output folder (default: derived from input).")
    src.add_argument("--base-dir", default=None,
                     help="Directory holding the per-molecule folder "
                          "(default: current directory).")

    solv = p.add_argument_group("solvation")
    solv.add_argument("--crest-solvent", default="hexane",
                      help="ALPB solvent for the CREST search.")
    solv.add_argument("--censo-solvents", nargs="+", default=["h2o"],
                      help="CPCM solvent(s) for the CENSO optimization step.")
    solv.add_argument("--gas-phase", action="store_true",
                      help="Run everything in the gas phase.")

    samp = p.add_argument_group("sampling")
    samp.add_argument("--ewin", type=float, default=6,
                      help="Energy window in kcal/mol.")
    samp.add_argument("--nclust", default="500", help="Target number of clusters.")
    samp.add_argument("--mdlen", default="x1.0", help="CREST MD length multiplier.")
    samp.add_argument("--mdlen-art", default="x2.0",
                      help="MD length multiplier for the artificial-bias runs.")
    samp.add_argument("--skip-artificial-dispersion", action="store_true",
                      help="Skip the artificial-dispersion CREST runs.")
    samp.add_argument("--skip-artificial-charge", action="store_true",
                      help="Skip the artificial-charge CREST runs.")
    samp.add_argument("--nci-mode", action="store_true",
                      help="Use CREST --nci mode (for non-covalent complexes).")
    samp.add_argument("--xtb-preopt", action="store_true",
                      help="GFN-FF loose preoptimization before CREST.")
    samp.add_argument("--noreftopo", action="store_true",
                      help="Disable CREST topology check and post-filter with MolBar.")
    samp.add_argument("--gentle", type=int, default=0, choices=[0, 1, 2, 3, 4],
                      help="Increasingly restrained sampling for fragile molecules.")
    samp.add_argument("--skip-conformer-search", action="store_true",
                      help="Reuse an existing CREST ensemble and go straight to CENSO.")

    run_g = p.add_argument_group("execution")
    run_g.add_argument("--censorc", default=None,
                       help="Path to the 'censo2rc' file; the _part1, _part2 and "
                            "_part1_sp variants must sit next to it.")
    run_g.add_argument("--steps", default="all",
                       help="CENSO steps to run: 'all' or any of part01, part2, sp.")
    run_g.add_argument("--tasks", "-P", type=int,
                       default=int(os.environ.get("SLURM_NTASKS", 1)),
                       help="Number of parallel tasks.")
    run_g.add_argument("--threads", "-O", type=int,
                       default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)),
                       help="Threads per task.")
    run_g.add_argument("--vibspectrum", action="store_true",
                       help="Compute the vibrational spectrum after the ensemble.")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not args.smiles and not args.xyz:
        build_parser().error("one of --smiles or --xyz is required")

    config = CrensoConfig(
        crest_solvent=args.crest_solvent,
        censo_solvents=args.censo_solvents,
        ewin=args.ewin,
        mdlen=args.mdlen,
        mdlen_art=args.mdlen_art,
        nclust=args.nclust,
        artificial_disp=[] if args.skip_artificial_dispersion else [0.5, 0.5, 1.5, 2.0],
        artificial_chrg=[] if args.skip_artificial_charge else [1, 2, -1, -2],
        censorc_path=args.censorc,
        uhf=args.uhf,
        P=args.tasks,
        O=args.threads,
        maxcores=args.tasks,
    )

    structure = CRENSOgen(
        smiles=args.smiles,
        configuration=config,
        base_dir=args.base_dir,
        xyz=args.xyz,
        nci_mode=args.nci_mode,
        folder_name=args.folder_name,
        gas_phase=args.gas_phase,
        noreftopo=args.noreftopo,
        xtb_preopt=args.xtb_preopt,
        gentle_sampling=args.gentle,
    )
    structure.charge = args.charge

    steps = "all" if args.steps == "all" else args.steps.split(",")

    if structure.is_monoatomic():
        print("Monoatomic molecule detected: using the single-point fast path.")
        structure.run_monoatomic_sp()
    else:
        if not args.skip_conformer_search:
            structure.initial_sample_crest()
            if not structure.nci_mode:
                structure.screen()
        structure.create_ensemble(steps)

    if args.vibspectrum:
        print(json.dumps(structure.compute_vibrational_spectrum(), indent=2))

    _save_timings(structure.folder_path,
                  "timings_gas.json" if args.gas_phase else "timings_solvated.json")

    final = structure.folder_path / "CRENSOconf_final.xyz"
    if final.exists():
        print(f"\nDone: {final} ({count_conf(final)} conformers)")
    else:
        print(f"\nFinished, but {final} was not produced. Check the logs above.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
