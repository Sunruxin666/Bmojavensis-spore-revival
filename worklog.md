# Worklog

## 2026-09-07 | Manuscript, SI and claim-level reproducibility audit

- Reran the complete Figure 2–7 and Supplementary Figure S1–S3 workflow after aligning the manuscript Methods and Results.
- Replaced static Supplementary Tables S1 and S2 with script-generated exports. Table S1 now carries explicit timepoint and contrast metadata; Table S2 contains the complete unique functional-category membership.
- Corrected missing glyphs in Supplementary Figs. S1 and S3 with device-independent mathematical labels.
- Replaced ambiguous Figure 7 correlation labels with `1 g + X-ray` and `SMG + X-ray` and synchronized panel comments with the submitted a–c layout.
- Added `R/validate_manuscript_claims.R`; all 41 numerical and structural manuscript/SI checks passed.
- Added `renv.lock` for the validated R 4.5.2 dependency environment; retained `requirements.R` as a fallback installer.
- Updated the README, input manifest, project status and validation record to match the generated deliverables.

## 2026-09-05 | Submission-language synchronization

- Standardized reader-facing terminology to `X-ray dose` and `gravity condition (1 g versus SMG)`; RCCS now identifies only the method used to generate SMG.
- Regenerated Figures 2–7, Supplementary Figs. S1–S3, source-data exports and statistical audits without changing the underlying data or statistical models.
- Aligned supplementary output numbering with the manuscript: logistic-fit parameters are now Table S3, while RNA-seq count QC is an unnumbered audit export.
- Verified that all 24 logistic models converged (mean R² = 0.981; range, 0.961–0.993).
- Retained the non-fatal local `glmmTMB`/`TMB` binary-version warning in the validation record for reproducibility.
- Renamed the Logistic coefficient `r` from the misleading internal identifier `mu_max` to `growth_rate_coef` in processed phenotype data, figure scripts and audit exports; the fitted maximum slope remains a separate derived field.
- Rebuilt the guide-revised English and Chinese manuscripts with seven embedded main figures, and rebuilt the English and Chinese Supplementary Information with three embedded supplementary figures.

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
