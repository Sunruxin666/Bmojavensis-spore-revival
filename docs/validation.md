# Validation record

Validation date: 2026-07-18

## Structural checks

- All 13 R files parsed successfully.
- Dependency and input preflight passed.
- Fourteen DEG contrast tables and six gene-to-KO mapping files were detected.
- No candidate repository file exceeds 50 MB.
- No private absolute path, compute-host name, credential token or API key was detected in candidate tracked files.

## End-to-end execution

`Rscript R/render_all.R` completed successfully and regenerated Figures 2–7, source-data tables, statistical audits and `sessionInfo.txt` under `results/`.

Figure 2 was additionally checked at its 600-dpi export size (4488 × 4251 px). The fixed microscopy panel remained legible, including both 10-µm scale bars, and the regenerated quantitative panels retained the trends in the author-approved manuscript figure.

Figures 3 and 4 were regenerated without the redundant in-figure ANOVA summary captions; the corresponding inferential results remain in the manuscript text and statistical audit outputs. Figure 5 was regenerated after synchronizing `showtext` with the 600-dpi raster device, correcting the previously undersized text. Visual inspection at manuscript scale confirmed readable pathway names, dose labels, legends and p values.

All Figure 3 and Figure 4 jitter layers use explicit fixed seeds, preventing run-to-run drift in the displayed raw-point positions.

The revised Figures 3–5 were embedded in `20260718-Manuscript-v19-Fig345-revised-red.docx`. The document was rendered to PDF (26 pages) and every page was visually inspected. No clipping or overlap was detected; Figures 3 and 4 retain their complete axes, and Figure 5 text is legible in the embedded version.

Figure 6 was reduced to the two data-derived panels a–b. The interpretive panel c, its rendering script, assembly dependency, repository documentation and direct manuscript legend text were removed. The revised `20260718-Manuscript-v20-Fig6c-removed-red.docx` was rendered to PDF (26 pages) and every page was visually inspected. The two-panel Figure 6 is complete and legible, with no clipping or overlap.

The Supplementary Information workflow regenerated Figs. S1-S3 and Tables S1/S4. Fig. S2 contains all 48 RNA-seq samples (24 at 3 h and 24 at 10 h). Fig. S3 displays all 24 replicate-level logistic fits; every fit converged, with mean R² = 0.981 and range 0.961–0.993. The five-page `Supplementary_Information_submission.docx` was rendered and every page was visually inspected with no clipping, overlap or missing figure content.

The synchronized main manuscript `20260718-Manuscript-v21-SI-synchronized-red.docx` was rendered to PDF (26 pages) and every page was visually inspected. Supplementary references now resolve only to Figs. S1-S3 and Tables S1-S4; obsolete provisional Fig. S5 text was removed. Red author-confirmation placeholders remain intentionally visible for completion before submission.

In the subsequent v22 clarification, the redundant Results reference to count-level QC/Table S1 was removed while the Methods reference was retained; Results now points only to the complete differential-expression results in Table S2 and the PCA in Fig. S2. The Figure 5b in-panel pathway-score formula was removed because the formula, DEG threshold and bubble-size definition are already stated in the self-contained legend. Figure 5 was regenerated from the revised R script and re-embedded. The 26-page v22 manuscript was rendered and visually inspected with no clipping, overlap or missing content.

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
- Figure 7 exploratory-correlation table.

Except for the explicitly disclosed fixed microscopy input in Figure 2 panel a, raster figures were regenerated rather than copied. Some PNG byte hashes and font rasterization differ across graphics-device runs, but dimensions are unchanged and the underlying numerical source data are identical. Figures 6a and 6b were pixel-identical to the submission exports in this validation run.
