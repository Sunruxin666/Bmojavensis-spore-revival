# Worklog

## 2026-07-18 | Initial GitHub release candidate

- Initial objective: separate the reproducible Figure 3–7 workflow from the private manuscript working tree; Figure 2 was subsequently added below after author review.
- Inputs: finalized submission scripts, phenotype tables, processed DEG outputs, KEGG lookup tables and secondary-spore workbooks.
- Changes: created portable `data/`, `R/`, `docs/` and `results/` structure; replaced legacy project paths; added preflight, dependency installer, render orchestrator, release checklist and statistical-decision ledger.
- Initial exclusions: Figures 1–2, raw microscopy fields, raw sequencing reads, RStudio state and legacy rendered outputs. The Figure 2 exclusion was superseded by the later inclusion entry below.
- Validation: all R files parsed; preflight passed; the complete Figure 3–7 and audit workflow ran successfully from the portable repository. Key Figure 3/4 audit tables, Figure 5 processed outputs and Figure 7 AUC/correlation tables matched the submission workspace byte-for-byte by SHA-256. No file exceeds 50 MB, and no private absolute path, host name, password, token or API key was detected in tracked candidate files.
- Remaining release decisions: public-data approval, KEGG lookup-table redistribution, citation metadata, license selection and final Figure 3/4 multiplicity policy.

## 2026-07-18 | Figure 2 inclusion

- Corrected scope after author review: Figure 1 remains excluded because no editable experimental-design artwork is available; Figure 2 is required in the repository.
- Added the author-approved flattened Figure 2 as a versioned visual input. The R script crops only microscopy panel a, regenerates panels b and c from biological-replicate germination data, exports their source data and assembles the submission figure.
- Reproducibility boundary: the 32 individual raw microscopy fields are not in the working archive, so the microscopy montage cannot yet be rebuilt or independently reprocessed from raw images.
- Validation: Figure 2 passed R parsing, dependency/input preflight, 600-dpi visual inspection and a complete Figure 2–7 end-to-end rerun.
