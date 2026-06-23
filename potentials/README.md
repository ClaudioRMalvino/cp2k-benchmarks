# Reference interatomic potentials

External, published ML-potential parameter sets and their provenance/reference
data, used as physics-validation targets for the CP2K NNP optimisation work.
Consolidated here (2026-06-23) from the former top-level `NNPs/`, `potentials/`,
and `nacl/` directories.

| Subdirectory | Source | Contents | Used for |
|---|---|---|---|
| `morawietz-water-2016/` | Morawietz et al. 2016, *"How van der Waals interactions determine the properties of water"* (SI) | 4 bulk-water NNPs (BLYP, BLYP-vdW, RPBE, RPBE-vdW); n2p2/RuNNer `weights.*`/`scaling.data`. **No `input.nn`** — symmetry functions are in the paper Tables S1–S4. | Fig S4 water reference family |
| `RPBE-vdW-2016/` | figshare | A single, *runnable* RPBE-vdW potential **with `input.nn`**. Its `weights`/`scaling` are byte-identical to `morawietz-water-2016/RPBE-vdW/`; this copy just adds the symmetry-function definitions. | Fig S4 production runs; `diag_n128_chebyshev.sh` gate test |
| `oneill-nacl-water-2024/` | O'Neill et al. 2024 (JPCL), NaCl-in-water (cloned; see its `PROVENANCE.md`) | 6 committee models (MP2, RPA, r2SCAN, revPBE0-D3, revPBE-D3, optB88-vdW) + input templates. | NaCl / ion-transport validation |

## Notes
- `oneill-nacl-water-2024/data-sets/` (~251 MB of training/validation data) is
  **git-ignored** — it stays on disk but out of the repo. `PROVENANCE.md` records
  the exact upstream commit to re-fetch it.
- The committee `models/**/*.data` and the canonical `potentials/**/*.data` are
  force-kept in git (see the repo `.gitignore`).
- The H₂O-64 headline deck (`../H2O-64_NNP_MD.inp`) does **not** use any potential
  here — it references `NNP/bulkH2O-jcp2020-cnnp/` (Schran et al. JCP 2020), staged
  into run directories at runtime.
