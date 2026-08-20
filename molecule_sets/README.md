# Molecule sets

The sets of molecules used in the different parts of the study, and the
literature references for the experimental LLE data.

| File | Entries | Used for |
|---|---|---|
| `initial_molecule_database.json` | 2569 | The full starting database, aggregated from the source databases listed below. |
| `conformer_generation_comparison_molecules.json` | 2434 | The molecules carried through the comparison of conformer generation methods. A strict subset of the initial database. |
| `gnnimplicitsolvent_molecules.json` | 127 | The set used for the GNNImplicitSolvent calculations. |
| `crenso_solvents_molecules.json` | 48 | The benchmark set for screening CREST × CENSO solvent combinations. |
| `lle_nist_references.json` | 259 | Bibliographic references for the experimental LLE data taken from NIST. |

The GNNImplicitSolvent and CRENSO-solvent sets were curated for their specific
sub-studies and are not subsets of the initial database.

## Format

The four molecule files are JSON objects keyed by `molecule_1`, `molecule_2`, …
Each entry carries:

```json
{
  "name": "(2-Aminophenyl)methanol",
  "SMILES": "C1=CC=C(C(=C1)CO)N",
  "CAS": ["5344-90-1"],
  "synonyms": ["2-Aminobenzyl alcohol", "o-Aminobenzyl alcohol"]
}
```

`CAS` and `synonyms` are lists and may be empty when no identifier could be
resolved. `initial_molecule_database.json` and
`conformer_generation_comparison_molecules.json` additionally carry
`databases` (the source databases the molecule appears in) and
`num_databases`.

Source databases represented in the initial database:

`CompSol_Pure_298K`, `FlexiSol`, `FreeSolv`, `CompSol_Binary_298.15K`,
`Klamt1998`, `LLE_NIST`, `ReSCoSS`, and `SAMPL0`–`SAMPL9`.

## References file

`lle_nist_references.json` is keyed by a citation key and holds standard
bibliographic fields:

```json
{
  "title": "...",
  "author": ["Shen, Zhipeng", "Wang, Qinbo"],
  "journal": "...", "year": "...", "volume": "...",
  "number": "...", "pages": "...",
  "doi": "...", "doi_link": "https://doi.org/..."
}
```
