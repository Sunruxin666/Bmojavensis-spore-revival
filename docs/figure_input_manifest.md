# Figure-to-input manifest

All figure scripts and all processed inputs used directly for plotting are
contained in this repository. Raw sequencing reads are not required to rerun
the submitted plots and are not duplicated here.

| Figure | Script | Direct input(s) | Boundary |
|---|---|---|---|
| Methods QC | `R/inoculum_OD600_uniformity_qc.R` | `data/phenotype/cleaned_OD600_data.csv` | Exports the overall 0 h OD600, all eight treatment-group means, and their coefficient of variation. |
| 1 | Not included | Not applicable | Experimental-design artwork will be rebuilt in BioRender. |
| 2 | `R/figure2_submission.R` | `data/phenotype/cleaned_germination_data.csv`; `data/figure2/Figure2A_microscopy_montage_flattened.png` | Panel a is a fixed flattened montage because the 32 raw fields are held separately. Panels b/c are regenerated. |
| 3 | `R/figure3_submission.R` | `data/phenotype/cleaned_germination_data.csv` | Biological-replicate phenotype table included. |
| 4 | `R/figure4_submission.R` | `data/phenotype/cleaned_OD600_data.csv` | OD600 plotting and fitting table included. |
| 5 | `R/figure5_submission.R` | `data/transcriptomics/`; `data/annotations/`; PERMANOVA and growth summaries under `data/phenotype/` | Processed DEG and gene-to-KO inputs used by the plot are included; raw reads are outside figure-level reproduction scope. |
| 6a | `R/figure6a_submission.R` | Figure 5 tables generated under `results/` | Generated automatically after Figure 5. |
| 6b | `R/figure6b_submission.R` | Figure 5 tables under `results/`; `data/phenotype/cleaned_germination_data.csv` | Generated automatically after Figures 3–5 inputs are available. |
| 6 integrated | `R/figure6_integrate.R` | Generated panel-a and panel-b PNG files | Assembly contains only the two data-derived panels; the former interpretive schematic was removed. |
| 7 | `R/figure7_submission.R` | `data/secondary_spores/*.xlsx`; Figure 5-derived tables under `results/` | Baseline and antimicrobial workbooks plus processed transcriptome inputs are included. |
| S1 | `R/supplementary_figures_submission.R` | `data/secondary_spores/control_1.xlsx`, `control_2.xlsx`, `DDAC.xlsx`, `soda_1.xlsx` | Exports relative-OD trajectories and replicate-level normalized AUC values for DDAC and NaHCO3. |
| S2 | `R/supplementary_figures_submission.R` | `data/phenotype/RNAseq_all_counts_3h.csv`, `RNAseq_all_counts_10h.csv` | PCA contains all 48 RNA-seq samples and exports sample coordinates. |
| S3 / Table S3 | `R/supplementary_figures_submission.R` | `data/phenotype/cleaned_OD600_data.csv` | Regenerates all 24 replicate-level logistic fits and exports the observed points, plotted fitted coordinates, model parameters, goodness of fit and convergence status. |
| Table S1 | `R/figure5_submission.R` | All 14 processed DEG tables under `data/transcriptomics/` | Exports every tested gene with explicit timepoint, gravity-condition, dose, contrast and source-file metadata. |
| Table S2 | `R/figure6a_submission.R` | `results/fig5_sig_annotated.csv` plus the scripted annotation-keyword rules | Exports the complete unique membership of the five functional categories. |
| Claim validation | `R/validate_manuscript_claims.R` | Generated audit and source-data tables under `results/` | Fails if 41 key manuscript/SI numerical or structural claims fall outside prespecified tolerances. |

`R/render_all.R` runs the scripts in dependency order and writes all figures,
source-data exports, statistical audits and `sessionInfo.txt` to `results/`.
