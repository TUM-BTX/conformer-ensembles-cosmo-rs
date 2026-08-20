# conformer-ensembles-cosmo-rs

Code and settings accompanying the paper:

> **The effect of conformer ensembles and their generation method on COSMO-RS
> predictions of phase equilibria**

This repository holds the scripts and configuration files needed to reproduce
the calculations and analyses reported in the paper. It is a companion to the
publication rather than a general-purpose package, but the CRENSO driver in
particular is written to be reusable.

## Contents

| Directory | Contents |
|---|---|
| [`crenso/`](crenso/) | The CRENSO workflow: a driver that runs CREST conformer sampling followed by four CENSO refinement steps, plus the `censo2rc` configuration files defining the level of theory. |

Further analysis code will be added here.

## Quick start

```bash
conda env create -f environment.yml
conda activate crenso

python crenso/crenso.py --smiles "CC(O)CO" --crest-solvent hexane --censo-solvents h2o
```

The environment file installs the Python dependencies. The quantum-chemistry
programs (ORCA, CREST, xTB) are installed separately; see
[`crenso/README.md`](crenso/README.md) for the versions used and how to point
the workflow at them.

## Citation

If you use this code, please cite the accompanying paper. See `CITATION.cff`.

## License

See `LICENSE`.
