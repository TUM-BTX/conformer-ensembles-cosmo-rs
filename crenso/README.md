# CRENSO

Automated conformer ensemble generation with [CREST](https://github.com/crest-lab/crest)
and [CENSO](https://github.com/grimme-lab/CENSO).

`crenso.py` runs a complete conformer workflow and produces a final, clustered
ensemble (`CRENSOconf_final.xyz`) for a molecule given as a SMILES string or an
XYZ structure:

1. **CREST sampling** — GFN-FF metadynamics, optionally broadened by additional
   runs with scaled dispersion (`-dispscal`) and with artificial charges, which
   drive the search into regions a single run tends to miss.
2. **CREST screening** — GFN-FF re-ranking, ensemble sorting, full GFN2
   optimization, and PCA/k-means clustering.
3. **CENSO refinement** — four sequential calls: prescreening, screening,
   per-solvent optimization, and a final screening of the merged ensemble.
4. **Final clustering** — CREST `--cregen` clustering on the CENSO-ranked
   ensemble.

Optionally an xTB Hessian is run for the vibrational spectrum.

## Requirements

### Python environment

```bash
conda env create -f environment.yml
conda activate crenso
```

This installs Python 3.13, RDKit and **CENSO 2.1.4** — the versions used for the
published results.

> **CENSO must be 2.x.** `pip install censo` gives CENSO 3.x, whose
> configuration format differs from 2.x; the `config/censo2rc*` files here are
> 2.x files and will not work with it. `environment.yml` therefore installs
> CENSO from the `v2.1.4` GitHub tag.

### Quantum-chemistry programs

These are **installed by you**, not by the environment file, so that you can use
your site's own builds:

| Program | Required | Version used | Where to get it |
|---|---|---|---|
| ORCA | yes | 6.0.0 | [orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de/) (proprietary; licence acceptance required) |
| CREST | yes | 3.0.2 | [crest-lab/crest](https://github.com/crest-lab/crest) or `conda install -c conda-forge crest=3.0.2` |
| xTB | yes | 6.7.1 | [grimme-lab/xtb](https://github.com/grimme-lab/xtb) or `conda install -c conda-forge xtb=6.7.1` |
| MolBar | no | — | Only for `--noreftopo`; install in a **separate** environment (it pins an older NumPy) |

ORCA is reached through CENSO, so its location goes in `config/censo2rc*`
(`orcapath`). CREST, xTB and MolBar are called directly and are found on `PATH`,
or via the environment variables below.

`environment.yml` has commented-out lines for `crest` and `xtb` if you would
rather let conda provide them; ORCA always has to be installed separately.

## Setup

`crest`, `censo` and `xtb` are found on `PATH`. To point at specific builds:

```bash
export CRENSO_CREST_BIN=/path/to/crest
export CRENSO_CENSO_BIN=/path/to/censo
export CRENSO_XTB_BIN=/path/to/xtb
export CRENSO_MOLBAR_BIN=/path/to/molbar   # optional
```

**You must edit `config/censo2rc*` before the first run.** The `[paths]` section
of all four files points at the ORCA and xTB executables of the machine the
paper was produced on:

```ini
[paths]
orcapath = /path/to/orca
xtbpath  = /path/to/xtb
orcaversion = 6.0.0
```

All other settings in those files reproduce the published protocol and should be
left alone unless you intend to change the level of theory. See
`config/README.md` for what each file controls.

## Usage

```bash
# From a SMILES string
python crenso.py --smiles "CC(O)CO" --crest-solvent hexane --censo-solvents h2o

# From an existing structure, in the gas phase
python crenso.py --xyz molecule.xyz --gas-phase

# Charged species, several CENSO solvents, 4 tasks x 8 threads
python crenso.py --smiles "CC(=O)[O-]" --charge -1 \
    --censo-solvents h2o hexane aniline --tasks 4 --threads 8
```

`python crenso.py --help` lists every option. Results are written to
`<base-dir>/<folder-name>/`, together with a `timings_*.json` breakdown of the
wall-clock cost of each step.

### Useful options

- `--crest-solvent` / `--censo-solvents` — ALPB solvent for the search, CPCM
  solvent(s) for the refinement. Each CENSO solvent is optimized separately and
  the results are merged before the final screening.
- `--gentle 1..4` — progressively more restrained sampling for molecules that
  fall apart during metadynamics.
- `--noreftopo` — disable CREST's topology check and instead filter conformers
  afterwards with MolBar. Useful when the reference topology check is too strict.
- `--nci-mode` — CREST `--nci` for non-covalent complexes.
- `--skip-conformer-search` — reuse an existing CREST ensemble in the folder and
  go straight to CENSO.
- `--steps part01,part2,sp` — run only some CENSO stages, for restarts.

## Use as a library

```python
from crenso import CRENSOgen, CrensoConfig

config = CrensoConfig(crest_solvent="hexane", censo_solvents=["h2o"],
                      P=4, O=8, censorc_path="config/censo2rc")

mol = CRENSOgen(smiles="CC(O)CO", configuration=config, base_dir="runs")
mol.initial_sample_crest()
mol.screen()
mol.create_ensemble("all")
```

`CrensoConfig.from_dict()` also accepts a nested dictionary, so existing driver
scripts keep working.

## Citation

If you use this code, please cite the accompanying paper (see `CITATION.cff`).
