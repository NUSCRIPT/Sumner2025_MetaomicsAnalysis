# to run: system("./run_analysis.sh")

# Primary Data Cleaning Notebooks
#Rscript -e 'library(rmarkdown); rmarkdown::render("00_Amplicon.Rmd","html_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("01_Metagenomic_Metatranscriptomic.Rmd","html_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("02_Virome.Rmd","html_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("03_HostTranscriptomics.Rmd","html_document")'

# Primary Analyses and Visualizations
Rscript -e 'library(rmarkdown); rmarkdown::render("07_PERMANOVAs.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("08_MantelTests.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("10_ordinations.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("06_Metadata.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("09_Dysbiosis.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("08_Clusters.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("15_DemographicsTable.Rmd","html_document")'
Rscript -e 'library(rmarkdown); rmarkdown::render("11_DifferentialFeatures_Unfreeze.Rmd","html_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("19_Overview_UpSets.Rmd","html_document")'
# Rscript -e 'library(rmarkdown); rmarkdown::render("18_HealthySegalDatasetComparison.Rmd","html_document")'

## HAllA
Rscript -e 'library(rmarkdown); rmarkdown::render("13_halla.Rmd","html_document")'
#Rscript -e 'library(rmarkdown); rmarkdown::render("14_HallaNetwork.Rmd","html_document")'
# Rscript -e 'library(rmarkdown); rmarkdown::render("","html_document")'


# Secondary Analyses and Visualizations
Rscript -e 'source("src/sanity_checks/qc_read_counts.r")'
Rscript -e 'source("src/other/metadata_completeness_heatmap.R")'
Rscript -e 'source("src/other/pnt_timeseries_over_hos.R")'
Rscript -e 'source("src/other/prevalence_feature_heatmap_by_pneumotype.R")'
Rscript -e 'source("src/overview_figures/overview_basic_plots.R")'

Rscript -e 'source("src/revisions/pneumotype_anovas.R")'
Rscript -e 'source("src/revisions/pneumotype_deseq2.r")'
Rscript -e 'source("src/other/alpha_diversity_violins.r")'
# Rscript -e 'source("")'
# Rscript -e 'source("")'
# Rscript -e 'source("")'


# Rscript -e 'source("other/amylase_heatmap.R")' # NEED TO CHECK CORRECT!!!!!!!!


# /path/to/my_data/BAS_manuscript/
# /path/to/my_data/BAS_manuscript/