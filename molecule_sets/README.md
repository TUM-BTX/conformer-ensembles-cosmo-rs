# Molecule sets

This folder contains the molecule lists used in the paper and the references
for the experimental LLE data.

| File | Entries | Used for |
|---|---|---|
| `initial_molecule_database.json` | 2569 | All molecules collected from the databases listed below. |
| `conformer_generation_comparison_molecules.json` | 2434 | Molecules from the initial database used to compare conformer generation methods. |
| `gnnimplicitsolvent_molecules.json` | 127 | Molecules used for the GNNImplicitSolvent calculations. |
| `crenso_solvents_molecules.json` | 48 | Molecules used to compare CREST and CENSO solvent combinations. |
| `lle_nist_references.json` | 259 | References for the experimental LLE data from NIST. |

The GNNImplicitSolvent and CRENSO solvent lists were prepared separately. They
are not subsets of the initial database.

## Format

The four molecule files use keys such as `molecule_1` and `molecule_2`. Each
entry has the following fields:

```json
{
  "name": "(2-Aminophenyl)methanol",
  "SMILES": "C1=CC=C(C(=C1)CO)N",
  "CAS": ["5344-90-1"],
  "synonyms": ["2-Aminobenzyl alcohol", "o-Aminobenzyl alcohol"]
}
```

`CAS` and `synonyms` are lists. They are empty when no CAS number or synonyms
were found. The initial database and conformer comparison files also include
`databases`, which lists the sources for each molecule, and `num_databases`,
which gives the number of sources.

The initial database includes molecules from:

`CompSol_Pure_298K`, `FlexiSol`, `FreeSolv`, `CompSol_Binary_298.15K`,
`Klamt1998`, `LLE_NIST`, `ReSCoSS`, and `SAMPL0`–`SAMPL9`.

## References file

In `lle_nist_references.json`, each citation key has an entry with the
reference details:

```json
{
  "title": "...",
  "author": ["Shen, Zhipeng", "Wang, Qinbo"],
  "journal": "...", "year": "...", "volume": "...",
  "number": "...", "pages": "...",
  "doi": "...", "doi_link": "https://doi.org/..."
}
```
