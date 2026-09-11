# Validation record

## 2026-09-07 manuscript/SI alignment validation

The complete workflow was rerun twice after the manuscript Methods, Results and Supplementary Information were aligned. Figures 2–7, Supplementary Figs. S1–S3, Supplementary Tables S1–S3, source-data exports and statistical audits regenerated successfully.

All 16 R scripts parsed successfully, and dependency/input preflight detected all 14 DEG contrasts and six gene-to-KO mapping files.

Supplementary Table S1 is now generated directly from all 14 processed differential-expression contrasts and contains explicit timepoint, gravity-condition, dose, contrast and source-file metadata. Supplementary Table S2 is generated from the complete unique functional-category membership rather than maintained as a static file. Supplementary Table S3 remains generated from the replicate-level logistic fits.

Unicode chemical and optical-density subscripts in Supplementary Figs. S1 and S3 were replaced with device-independent plotmath labels, eliminating missing-glyph boxes in the PDF and raster exports. Figure 7 panel-c treatment labels now distinguish `1 g + X-ray` from `SMG + X-ray` rather than using the ambiguous label `Radiation`.

`R/validate_manuscript_claims.R` now checks 41 prespecified numerical and structural claims. The checks cover inoculum OD600 and CV, germination contrasts and Holm-adjusted p values, growth-rate contrasts, PERMANOVA/PERMDISP results, functional-category correlations, secondary-spore AUC values, Supplementary Table row counts and all 24 logistic fits. All checks passed on 2026-09-07.

An `renv.lock` file records the validated R 4.5.2 environment and its package dependencies. `requirements.R` is retained as an unpinned fallback for systems where restoring the lockfile is impractical.

The remaining local `glmmTMB`/`TMB` binary-version warning was non-fatal and did not affect any audit or claim-validation result.

Validation date: 2026-09-05

## Current submission validation

The complete workflow was rerun on 2026-09-05 after standardizing reader-facing terminology to X-ray dose and gravity condition (1 g versus SMG). RCCS is identified only as the method used to generate SMG. Figures 2–7, Supplementary Figs. S1–S3, source-data exports and statistical audits regenerated successfully without changes to the underlying data or statistical models.

The supplementary workflow now exports the logistic-fit record as `Table_S3_logistic_model_fit_quality.csv`, matching the Supplementary Information. RNA-seq count-level QC is exported separately as `RNAseq_count_QC_audit.csv` and is not assigned a Supplementary Table number.

The Logistic coefficient `r` is now named `growth_rate_coef` throughout the processed input, scripts and regenerated audit tables. This prevents the coefficient from being confused with the fitted maximum slope `(K - A)r/4`, which remains a separate derived field.

The regenerated figures were embedded in `20260905-Manuscript_guide-revised_clean.docx`, `20260905-Manuscript_guide-revised_marked.docx` and `20260905-Manuscript_guide-revised_Chinese-review.docx`. Supplementary figures were embedded in the corresponding guide-revised Supplementary Information files.

## Structural checks

- All 14 R files parsed successfully.
- Dependency and input preflight passed.
- Fourteen DEG contrast tables and six gene-to-KO mapping files were detected.
- No candidate repository file exceeds 50 MB.
- No private absolute path, compute-host name, credential token or API key was detected in candidate tracked files.

## End-to-end execution

`Rscript R/render_all.R` completed successfully and regenerated Figures 2–7, source-data tables, statistical audits and `sessionInfo.txt` under `results/`.

Figure 2 was additionally checked at its 600-dpi export size (4488 × 4251 px). The fixed microscopy panel remained legible, including both 10-µm scale bars, and the regenerated quantitative panels retained the trends in the author-approved manuscript figure.

Figures 3 and 4 were regenerated without the redundant in-figure ANOVA summary captions; the corresponding inferential results remain in the manuscript text and statistical audit outputs. Figure 5 was regenerated after synchronizing `showtext` with the 600-dpi raster device, correcting the previously undersized text. Visual inspection at manuscript scale confirmed readable pathway names, dose labels, legends and p values.

All Figure 3 and Figure 4 jitter layers use explicit fixed seeds, preventing run-to-run drift in the displayed raw-point positions.

The final Chinese review manuscript was rendered to a 30-page A4 PDF. Seven embedded figure objects were detected; the title page, Figure 6, the application paragraphs in the Discussion and the Conclusion were visually inspected with no clipping or overlap. Figure 6 retains only the two data-derived panels a–b and labels `r` as the fitted growth-rate coefficient.

The Supplementary Information workflow regenerated Figs. S1-S3, Table S3 and the separate RNA-seq QC audit. Fig. S2 contains all 48 RNA-seq samples (24 at 3 h and 24 at 10 h). Fig. S3 displays all 24 replicate-level logistic fits; every fit converged, with mean R² = 0.981 and range 0.961–0.993.

Supplementary references resolve only to Figs. S1-S3 and Tables S1-S3; obsolete provisional figure and table numbering has been removed. The Chinese review Supplementary Information was rendered to a five-page PDF with three embedded figure objects; representative figure pages were visually inspected with no clipping or overlap.

The local validation environment emitted a non-fatal warning that `glmmTMB` had been built against TMB 1.9.19 while TMB 1.9.21 was installed. Figure 2 completed successfully; a clean public release environment should install both packages together through `requirements.R` (or a future lockfile) to avoid this local binary-version warning.

## Output identity checks

The following outputs matched the submission workspace byte-for-byte by SHA-256:

- Figure 3 multiplicity audit.
- Figure 4 fit-QC summary.
- Figure 5 annotated significant-DEG table.
- Figure 5 transcriptome–phenotype bridge.
- Figure 5 KEGG pathway scores.
- Figure 5 DEG-burden table.
- Figure 7 replicate-level net-response table.
- Figure 7 group AUC summary.
- Figure 7 correlation table.

Except for the explicitly disclosed fixed microscopy input in Figure 2 panel a, raster figures were regenerated rather than copied. Some PNG byte hashes and font rasterization differ across graphics-device runs, but dimensions are unchanged and the underlying numerical source data are identical. Figures 6a and 6b were pixel-identical to the submission exports in this validation run.
