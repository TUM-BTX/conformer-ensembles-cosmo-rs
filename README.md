# conformer-ensembles-cosmo-rs

Supporting code and calculation settings for our paper:

> **The effect of conformer ensembles and their generation method on COSMO-RS
> predictions of phase equilibria**

Here we collect the scripts, settings, and molecule lists used in the study.
The CRENSO workflow can also be used to generate and refine conformer ensembles
for other molecules.

## Contents

- [`crenso/`](crenso/) contains the CRENSO script, which runs CREST conformer
  sampling followed by four CENSO refinement steps. The accompanying `censo2rc`
  files specify the level of theory.
- [`molecule_sets/`](molecule_sets/) lists the molecules used in each part of the
  study, with names, SMILES, CAS numbers, and synonyms. It also includes the
  literature references for the experimental LLE data.

## Running CRENSO

Create and activate the Python environment, then run CRENSO with a SMILES string:

```bash
conda env create -f environment.yml
conda activate crenso

python crenso/crenso.py --smiles "CC(O)CO" --crest-solvent hexane --censo-solvents h2o
```

The environment file installs the Python dependencies. You will also need
separate installations of ORCA, CREST, and xTB. See
[`crenso/README.md`](crenso/README.md) for the versions we used and instructions
for setting their paths.

## Citation

If you use this code, please cite the accompanying paper. The citation details
are in [`CITATION.cff`](CITATION.cff).

## License

See [`LICENSE`](LICENSE).
