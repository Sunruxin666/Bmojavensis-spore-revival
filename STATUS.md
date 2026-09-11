# Project status

Last updated: 2026-09-07 Asia/Shanghai

## Objective

Prepare a portable, auditable GitHub repository that reproduces manuscript Figures 2–7 without relying on the legacy local directory structure.

## Current state

| Layer | Status | Validation | Remaining action |
|---|---|---|---|
| Portable paths | Complete | No private absolute or legacy project-root paths detected | Maintain in future edits |
| Minimal processed inputs | Complete | Preflight found 14 DEG tables and 6 KO maps | Author data-release approval |
| Figure and SI scripts | Validated candidate | Figures 2–7, Supplementary Figs. S1–S3 and Supplementary Tables S1–S3 reran successfully; Figure 2 uses a disclosed fixed microscopy input | Tag accepted-manuscript version |
| Statistical and claim audits | Validated candidate | Holm-adjusted planned contrasts, PERMANOVA/PERMDISP and 41 manuscript/SI claims are checked automatically | Freeze accepted-manuscript outputs |
| Documentation | Candidate | README, input manifest, decision ledger and validation record synchronized through 2026-09-07 | Add final citation and license |

## Claims safe to report

- The repository is a pre-release reproduction workflow for Figures 2–7, with Figure 2 panel a treated as a disclosed fixed raster input.
- Figure 7 uses a combined baseline estimate and replicate-level antimicrobial AUCs.
- Figure 7 panel-c correlations are exploratory and based on three dose-level observations per estimable stratum.

## Release blockers

1. Author approval for public redistribution of processed data and KEGG-derived lookup tables.
2. Final software/data license and citation metadata.
3. Confirm whether KEGG-derived lookup tables will be redistributed or regenerated at release time.
