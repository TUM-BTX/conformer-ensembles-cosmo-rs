# CENSO configuration files

The CRENSO workflow makes four sequential CENSO calls. Each is driven by one of
the files here, and each file enables exactly one CENSO stage (`run = True`);
every other stage, including the `[nmr]` and `[uvvis]` property blocks, is
switched off and never executes.

| Call | File | Active stage | Method | Solvation | Window (kcal/mol) |
|---|---|---|---|---|---|
| 1. Prescreening | `censo2rc` | `[prescreening]` | B97-D3/def2-SV(P) | not set | 4.0 |
| 2. Screening | `censo2rc_part1` | `[screening]` | r²SCAN-3c | CPCM | 3.5 |
| 3. Optimization | `censo2rc_part2` | `[optimization]` | r²SCAN-3c | CPCM | 2.5 |
| 4. Final screening | `censo2rc_part1_sp` | `[screening]` | r²SCAN-3c | CPCM | 3.5 |

Steps 1, 2 and 4 run in the CREST search solvent. Step 3 is repeated once per
CENSO solvent (or once in the gas phase), and the optimized structures are
merged into the ensemble that step 4 re-ranks.

Notes:

- r²SCAN-3c is a composite method; its `def2-mTZVPP` basis is intrinsic, and the
  files state it explicitly.
- Steps 1, 2 and 4 evaluate the mRRHO contribution from a single-point Hessian
  (`bhess`, ALPB). Step 3 sets `evaluate_rrho = False`, so the per-solvent
  optimizations are ranked on electronic and solvation energy alone.
- `censo2rc_part1` and `censo2rc_part1_sp` describe the *same* level of theory.
  They differ only in the `func` key of the inactive `[prescreening]` block,
  which has no effect on the calculation.

## Before running

The `[paths]` section of all four files ships with placeholders. Replace them
with the ORCA and xTB executables on your machine:

```ini
[paths]
orcapath = /path/to/orca      <- edit
xtbpath  = /path/to/xtb       <- edit
orcaversion = 6.0.0
```

`orcaversion` must match the ORCA you point at, since CENSO adjusts its input
syntax to the version. The published results were produced with **ORCA 6.0.0**
(shared-library build, OpenMPI 4.1.6, AVX2) and xTB 6.7.1.

Everything outside `[paths]` reproduces the protocol used in the paper and is
unchanged from the production run.
