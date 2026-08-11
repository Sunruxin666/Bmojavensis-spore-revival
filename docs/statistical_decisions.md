# Statistical and interpretation decisions

## Figure 3

The plotting script retains the original per-dose Welch tests. `figure34_statistical_audit.R` exports Holm and Benjamini–Hochberg sensitivity analyses. A final multiplicity policy must be selected before public release and manuscript acceptance.

## Figure 4

Growth parameters are obtained from four-parameter logistic fits. Fit-boundary and finite-value checks are exported by `figure34_statistical_audit.R`. The plotting script currently retains the original `emmeans` contrasts; multiplicity-adjusted sensitivity results are exported separately.

## Figures 5 and 6

Pathway and functional-module scores are derived from thresholded significant DEGs. The signed score is `sum(log2 fold change) / sqrt(k)`, where `k` is the number of contributing genes. Displayed pathways and representative genes were selected for the manuscript figure. Missing DEG-derived scores must not be interpreted as measured biological zeros without explicit justification.

## Figure 7

`control_1` and `control_2` are treated as measurement sets from one experiment and merged into one combined baseline trajectory. Their six traces are not counted as six independent biological replicates. Baseline AUC is calculated from the combined group-mean trajectory.

For H2O2, gentamicin, DDAC and NaHCO3, AUC is calculated separately for each biological replicate (`n = 3`). Replicate-level net response is the corresponding group-level combined baseline AUC minus the challenged-replicate AUC.

Pearson correlations in panel c are calculated from three dose-level group observations per estimable stratum. They are descriptive and hypothesis-generating; no significance stars or confirmatory mechanistic claims are used.
