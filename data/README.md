# Data manifest

This directory contains the minimum processed inputs needed by the figure scripts. It intentionally excludes raw microscopy fields, raw sequencing reads, RStudio state, intermediate R objects and legacy rendered figures.

## `figure2/`

- `Figure2A_microscopy_montage_flattened.png`: author-approved 32-view microscopy montage used directly as Figure 2a. The script regenerates panels b and c from `phenotype/cleaned_germination_data.csv`. This is a fixed visual input, not a substitute for the unavailable raw microscopy fields.
- `Figure2_flattened_reference.png`: provenance copy of the complete author-approved Figure 2 from which the panel-a input was isolated; it is not used by the rendering script.

## `phenotype/`

- `cleaned_germination_data.csv`: biological-replicate germination measurements.
- `cleaned_OD600_data.csv`: later growth measurements used for curve fitting.
- `Fig5A_permanova_3h.csv` and `Fig5A_permanova_10h.csv`: PERMANOVA summaries.
- `growth_kinetics_group_summary.csv`: group-level phenotype summary used to connect transcriptome and phenotype analyses.

## `annotations/`

- `ko_definition.tsv`: KEGG orthologue descriptions.
- `ko_to_pathway.tsv`: KO-to-pathway mappings.
- `pathway_names.tsv`: pathway identifiers and names.

These KEGG-derived lookup tables should be reviewed for redistribution terms and may instead be regenerated from KEGG at release time.

## `transcriptomics/`

Processed full differential-expression tables for seven contrasts at 3 h and seven contrasts at 10 h, plus six gene-to-KO mapping files. The `count_matrices/` subdirectory contains the six processed gene-count matrices used to reconstruct the all-sample Supplementary Fig. S2 PCA and Supplementary Table S1 count-level QC summary. Raw sequencing reads are not included.

## `secondary_spores/`

Two baseline measurement workbooks and four antimicrobial-challenge workbooks. The two baseline workbooks are combined for baseline estimation and are not treated as six independent biological replicates.
