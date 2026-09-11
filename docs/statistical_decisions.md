# Statistical and interpretation decisions

## Figure 3

The endpoint-specific factorial models are the primary analyses. Four planned 1 g-SMG comparisons are made across X-ray doses within each endpoint. Welch tests report the SMG-minus-1 g estimate and 95% confidence interval, and Holm adjustment controls the family-wise error rate across the four comparisons. Figure symbols use the Holm-adjusted p values.

## Figure 4

Growth parameters are obtained from four-parameter logistic fits. The coefficient `r` describes the steepness of the fitted OD increase and is named `growth_rate_coef` throughout the workflow; the fitted maximum slope is a separate quantity, `(K - A)r/4`. Fit-boundary and finite-value checks are exported by `figure34_statistical_audit.R`. Four planned model-based 1 g-SMG contrasts are made across X-ray doses within each endpoint. The contrast estimates are accompanied by 95% confidence intervals, and figure symbols use Holm-adjusted p values across the four comparisons.

## Figures 5 and 6

Pathway and functional-module scores are derived from thresholded significant DEGs. The signed score is `sum(log2 fold change) / sqrt(k)`, where `k` is the number of contributing genes. For Figure 6a, duplicate annotations are collapsed to one entry per gene symbol by retaining the contrast with the largest absolute shrunken log2 fold change; the four largest absolute responses are then retained within each prespecified functional category. The selected rows are exported to `Figure6a_representative_gene_selection.csv`. Missing DEG-derived scores must not be interpreted as measured biological zeros without explicit justification.

`transcriptome_permanova_audit.R` reproduces the Figure 5a PERMANOVA from low-count-filtered, Hellinger-transformed count matrices and reports PERMDISP diagnostics separately at 3 and 10 h. The diagnostic uses the eight dose-by-gravity treatment groups and a declared permutation seed.

## Figure 7

`control_1` and `control_2` are treated as measurement sets from one experiment and merged into one combined baseline trajectory. Their six traces are not counted as six independent biological replicates. Baseline AUC is calculated from the combined group-mean trajectory.

For H2O2, gentamicin, DDAC and NaHCO3, AUC is calculated separately for each biological replicate (`n = 3`). Replicate-level net response is the corresponding group-level combined baseline AUC minus the challenged-replicate AUC.

Pearson correlations in panel c are calculated from three dose-level group observations per estimable stratum. They are descriptive and hypothesis-generating; no significance stars or confirmatory mechanistic claims are used.
