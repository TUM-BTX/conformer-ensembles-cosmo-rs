# CRENSO

`crenso.py` uses [CREST](https://github.com/crest-lab/crest) and
[CENSO](https://github.com/grimme-lab/CENSO) to generate conformer ensembles.
It takes a SMILES string or an XYZ file and writes the final ensemble to
`CRENSOconf_final.xyz`.

The calculation has four steps:

1. CREST searches for conformers using GFN-FF metadynamics. Optional runs use
   scaled dispersion (`-dispscal`) and artificial charges.
2. The conformers are ranked with GFN-FF, sorted, optimized with GFN2, and
   grouped using PCA/k-means clustering.
3. CENSO runs prescreening, screening, and optimization in each solvent. The
   results are then combined and screened again.
4. CREST `--cregen` clusters the CENSO-ranked ensemble.

## Requirements

### Python environment

Run these commands from the repository root:

```bash
conda env create -f environment.yml
conda activate crenso
```

This installs Python 3.13, RDKit, and CENSO 2.1.4, as used in the paper.
The configuration files in this repository require CENSO 2.x. The environment
file installs version 2.1.4 from GitHub; these settings do not work with CENSO 3.x.

### Other programs

Install ORCA, CREST, and xTB separately. We used:

| Program | Version | Download |
|---|---|---|
| ORCA | 6.0.0 | [ORCA forum](https://orcaforum.kofo.mpg.de/) (licence acceptance required) |
| CREST | 3.0.2 | [CREST on GitHub](https://github.com/crest-lab/crest) |
| xTB | 6.7.1 | [xTB on GitHub](https://github.com/grimme-lab/xtb) |

CREST and xTB can also be installed with conda. Uncomment their lines in
`environment.yml`, or run `conda install -c conda-forge crest=3.0.2 xtb=6.7.1`.

MolBar is only needed for `--noreftopo`. Install it in a separate environment
because it requires an older NumPy version.

## Setup

`crest`, `censo`, and `xtb` should be on your `PATH`. You can also set their
locations explicitly:

```bash
export CRENSO_CREST_BIN=/path/to/crest
export CRENSO_CENSO_BIN=/path/to/censo
export CRENSO_XTB_BIN=/path/to/xtb
export CRENSO_MOLBAR_BIN=/path/to/molbar   # optional
```

Before the first run, set the ORCA and xTB paths in all four
`crenso/config/censo2rc*` files. The paths currently refer to the machine we used
for the paper. Edit the `[paths]` section to match your installation:

```ini
[paths]
orcapath = /path/to/orca
xtbpath  = /path/to/xtb
orcaversion = 6.0.0
```

Keep the other settings to use the same calculation protocol as the paper.
[config/README.md](config/README.md) describes the four files.

## Running a calculation

Run these examples from the `crenso/` directory:

```bash
# From a SMILES string
python crenso.py --smiles "CC(O)CO" --crest-solvent hexane --censo-solvents h2o

# From an existing structure, in the gas phase
python crenso.py --xyz molecule.xyz --gas-phase

# Charged species, several CENSO solvents, 4 tasks x 8 threads
python crenso.py --smiles "CC(=O)[O-]" --charge -1 \
    --censo-solvents h2o hexane aniline --tasks 4 --threads 8
```

`python crenso.py --help` lists all options. Results are saved in
`<base-dir>/<folder-name>/`. The `timings_*.json` file records how long each
step took.

### Useful options

- `--crest-solvent` sets the ALPB solvent for the CREST search.
- `--censo-solvents` sets the CPCM solvents for CENSO. Conformers are optimized
  in each solvent, then combined for the final screening.
- `--gentle 1..4` reduces the sampling intensity if molecules break apart during
  metadynamics. Higher values apply stronger restrictions.
- `--noreftopo` replaces CREST's topology check with a MolBar check after the
  search. This can help if CREST's check rejects too many conformers.
- `--nci-mode` turns on CREST's `--nci` option for non-covalent complexes.
- `--skip-conformer-search` starts CENSO using an existing CREST ensemble in the
  calculation folder.
- `--steps part01,part2,sp` selects which CENSO stages to run, for example when
  restarting a calculation.

## Running from Python

You can also call the workflow from a Python script:

```python
from crenso import CRENSOgen, CrensoConfig

config = CrensoConfig(crest_solvent="hexane", censo_solvents=["h2o"],
                      P=4, O=8, censorc_path="config/censo2rc")

mol = CRENSOgen(smiles="CC(O)CO", configuration=config, base_dir="runs")
mol.initial_sample_crest()
mol.screen()
mol.create_ensemble("all")
```

Settings can also be passed as a nested dictionary using
`CrensoConfig.from_dict()`.

## Citation

If you use this code, please cite the accompanying paper (see
[CITATION.cff](../CITATION.cff)).
