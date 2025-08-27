# SCRIPT Meta-Omics Analysis

**Transitions in the lung microbiota landscape associate with distinct patterns of pneumonia progression**\
[DOI: 10.1101/2024.08.02.24311426](https://doi.org/10.1101/2024.08.02.24311426) [PMID: PMC11326345](https://pmc.ncbi.nlm.nih.gov/articles/PMC11326345/)

*AKA BAS_manuscript*

------------------------------------------------------------------------

## Overview

This repository contains all code, analysis workflows, and supporting documentation for analysis presented in Sumner et al., *Transitions in the lung microbiota landscape associate with distinct patterns of pneumonia progression*, a part of the SCRIPT study, which investigates the lung microbial ecosystem in patients with severe pneumonia using systems biology approaches. The project integrates metagenomics, metatranscriptomics, host transcriptomics, and quantitative PCR to provide a comprehensive perspective on microbial dynamics and host responses.

The analyses here support the findings presented in the corresponding manuscript submitted for peer review. This README provides guidance for researchers seeking to reproduce results, review methods, or use the codebase for related projects.

------------------------------------------------------------------------

## Table of Contents

-   [Project Structure](#project-structure)
-   [Setup & Installation](#setup--installation)
-   [Reproducibility](#reproducibility)
-   [Data Availability](#data-availability)
-   [Analysis Scripts](#analysis-scripts)
-   [Export Utilities](#export-utilities)
-   [Citation](#citation)
-   [Contact](#contact)
-   [License](#license)

------------------------------------------------------------------------

## Project Structure {#project-structure}

The repository is organized as follows:

-   `config.r` — Configure analysis scripts by pointing to correct files for processing.
-   `src/common/` — Shared utility scripts for loading data, feature processing, clustering, and plotting.
-   `00_Amplicon.Rmd` — Amplicon sequencing analysis.
-   `01_Metagenomic_Metatranscriptomic.Rmd` — Metagenomic & metatranscriptomic processing.
-   `02_Virome.Rmd` — Virome processing and analysis.
-   `03_HostTranscriptomics.Rmd` — Host transcriptomic profiling processing.
-   `05_QuantitativePCR.Rmd` — qPCR processing.
-   `06_Metadata.Rmd` — Metadata integration and QC.
-   `07_PERMANOVAs.Rmd` — Multivariate statistical analyses comparing clinical and other variables with omics profiles.
-   `08_Clusters.Rmd`, - Main pneumotype analyses and visualization notebook.
-   `08_MantelTests.Rmd` - Mantel test analyses.
-   `09_Dysbiosis.Rmd` — Dysbiosis scoring visualization and analysis.
-   `10_ordinations.Rmd` — Ordination methods for dimensionality reduction and visualization.
-   `11_DifferentialFeatures_Unfreeze.Rmd` — Differential feature analysis.
-   `13_halla.Rmd`, `14_HallaNetwork.Rmd` — Hierarchical All-against-All Association testing.
-   `15_DemographicsTable.Rmd` — Cohort descriptive statistics and tables.
-   `figures/` — Output figures for the manuscript.
-   Other supporting files and scripts.

------------------------------------------------------------------------

## Setup & Installation

All analyses are carried out in R. Required packages and version constraints are listed below. For full reproducibility, use the provided session info and install package versions as specified.

### Core Dependencies

-   R (≥ 4.0.0)
-   `tidyverse`
-   `ggplot2`, `dplyr`, `readr`, `vegan`, `phyloseq`, `DESeq2`, `halla`, and others as specified in individual scripts.

------------------------------------------------------------------------

## Reproducibility {#reproducibility}

-   All code is version-controlled and scripts are annotated for clarity.
-   Random seeds are set where applicable.
-   Intermediate data files and outputs are included or referenced for each analysis.
-   A list of required input data files and expected formats is provided in each script header.

------------------------------------------------------------------------

## Data Availability {#data-availability}

------------------------------------------------------------------------

## Analysis Scripts {#analysis-scripts}

## Setup

-   `src/common/`: Common functions for various tasks and analyses. Core helper functions for processing metadata and features tables.

    -   `src/common/bracken_processing_functions.r`:
    -   `src/common/clean_feature_names.r`:
    -   `src/common/consensus_cluster_functions.R`:
    -   `src/common/decontam_functions.r`:
    -   `src/common/drop_samples.r`:
    -   `src/common/feature_tables_functions.r`:
    -   `src/common/get_clusters.r`:
    -   `src/common/get_distances.r`:
    -   `src/common/get_viral.r`:
    -   `src/common/graph_bubble_plot.r`:
    -   `src/common/graph_features.r`:
    -   `src/common/graph_violin.r`:
    -   `src/common/humann_processing_functions.r`:
    -   `src/common/krakenuniq_processing_functions.r`:
    -   `src/common/load_antibiotic_table.r`:
    -   `src/common/load_culture_feature_table.r`:
    -   `src/common/load_features.r`:
    -   `src/common/load_lists.r`:
    -   `src/common/load_metadata.r`:
    -   `src/common/mantel_test_functions.R`:
    -   `src/common/metaphlan_processing_functions.r`:
    -   `src/common/nature_theme.r`:
    -   `src/common/ordination_functions.r`:
    -   `src/common/sbatch_header.sh`:

-   `src/other/`: Assorted analyses and fun visualizations that do not fit into other groups

    -   `src/other/all_feature_heatmaps.R`:
    -   `src/other/alpha_diversity_violins.r`:
    -   `src/other/amylase_heatmap.R`:
    -   `src/other/autostereogram_pcoa.r`:
    -   `src/other/get_spp_kos.R`:
    -   `src/other/metadata_completeness_heatmap.R`:
    -   `src/other/n_values.R`: -`src/other/phyla_differentials.R`:
    -   `src/other/pnt_timeseries_over_hos.R`:
    -   `src/other/prevalence_feature_heatmap_by_pneumotype.R`:
    -   `src/other/temporal_alpha_diversity.r`:

-   `src/overview_figures/`: Some figures showcasing overview of data

    -   `src/overview_figures/overview_basic_plots.R`:
    -   `src/overview_figures/overview_heat.r`:

-   `src/revisions/`: Some analyses requested during revisions

    -   `src/revisions/host_subpermanova_pneumotype_heatmap.r`:
    -   `src/revisions/kegg_enrichment_analysis.r`:
    -   `src/revisions/phylogenetic_tree_gtdbtk.R`:
    -   `src/revisions/pneumotype_anovas.R`:
    -   `src/revisions/pneumotype_deseq2.r`:
    -   `src/revisions/tidyheats.r`:

-   `src/sanity_checks/`: Sanity checks for data quality

    -   `src/sanity_checks/metadata_nmi_hclust.r`:
    -   `src/sanity_checks/prevalence_threhold_sanity_check.r`:
    -   `src/sanity_checks/qc_read_counts.r`:
    -   `src/sanity_checks/quast_repetitive_sequence_assembly_comparison.r`:

-   `src/temporal/`: Feature level temporal analyses

    -   `src/temporal/temporal_heatmap.r`:

## Export Utilities {#export-utilities}

Export scripts and one-liners for compiling figures for manuscript submission (e.g., for Overleaf):

**Main Figures:**

``` sh
pdfunite clinical_data/figure1_demographics_overview.pdf dysbiosis/figure2_patchwork_e.pdf pneumotypes/figure4_all_pneumotype_data.pdf host_deseq/figure_4_mdnp_differentials_long.pdf host_deseq/figure_5_cluster_deseq.pdf main_figures.pdf
```

**Supplemental Figures:**

``` sh
pdfunite demographics_sampling.pdf mantel_heat.pdf supplemental_ordination_* bubble_plot_abundance_facet.pdf bubble_plot_abundance_facet_all.pdf dysbiosis_density.pdf dysbiosis_sofa.pdf supplemental_figures.pdf
```

------------------------------------------------------------------------

## Citation {#citation}

If you use this codebase, please cite the associated manuscript:

> Sumner JT, et al. Transitions in the lung microbiota landscape associate with distinct patterns of pneumonia progression. [MedRxiv], [2024], [10.1101/2024.08.02.24311426].

------------------------------------------------------------------------

## Contact {#contact}

For questions, collaboration, or dataset access requests, contact:

-   [jtsumner](https://github.com/jtsumner)

------------------------------------------------------------------------

## License {#license}

This repository is released under the MIT license. See LICENSE file for details.

------------------------------------------------------------------------

## Acknowledgements

We thank the SCRIPT consortium, funding agencies, clinical collaborators, and all contributors to this project.

------------------------------------------------------------------------
