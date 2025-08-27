library(tidyverse)
library(phyloseq)
# source('src/common/load_metadata.r')
# source('src/common/get_distances.r')
source("src/common/feature_tables_functions.r")
library(janitor)

features_amp <- readRDS("objects/amp/AMP_08_PhyloseqGenusTSS.rds") %>% 
  phyloseq::otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="feature") %>% 
  rename_with(.fn = function(x){str_split_i(x, "-", 1)}) %>% 
  as_tibble() %>% 
  # get_most_prevalent(min_prevalence = 0.01) %>% 
  normalize_feature_table(.,normalize=T,transform="log") 
names(features_amp) <- gsub(names(features_amp), pattern = "PT", replacement = "") 

features_amp_asv <- readRDS("objects/amp/AMP_04_PhyloseqTSS.rds") %>% 
  phyloseq::otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="feature") %>% 
  rename_with(.fn = function(x){str_split_i(x, "-", 1)}) %>% 
  as_tibble() %>% 
  normalize_feature_table(.,normalize=T,transform="log") 
names(features_amp_asv) <- gsub(names(features_amp_asv), pattern = "PT", replacement = "") 

features_mgx <- readRDS("objects/mgx/mgx_01_table_MetaphlanLog.rds")
features_mgx_pfam <- readRDS("objects/mgx/mgx_03_table_HumannPFAMLog.rds") 
features_mtx_bkn <- readRDS('objects/mtx/mtx_05_table_BrackenLog.rds')
features_mtx_mpa <- readRDS('objects/mtx/mtx_01_table_MetaphlanLog.rds')
features_mtx_pfam <- readRDS('objects/mtx/mtx_03_table_HumannPFAMLog.rds')


htx <- read_rds("objects/htx/HTX_01A_COMPFeatureTable.rds") %>% 
  dplyr::rename(feature = "Symbol") %>% 
  get_most_prevalent(min_prevalence = 0.1) %>% 
  normalize_feature_table(.,normalize=T,transform="log") %>% 
  get_bal_samples()

# Load Viral
source('src/common/get_viral.r')
features_viral <- get_viruses()

# Load Culture
source("src/common/load_culture_feature_table.r")
culture_feature <- get_culture_feature_table("genus")

# # Load ABX
# source("src/common/load_antibiotic_table.r")
# abx_feature <- get_abx_table() 

mox_feature_list <- list(
  "AMP [Genus]" = features_amp,
  "AMP [ASV]" = features_amp_asv,
  "DNA [Taxonomy]" = features_mgx,
  "DNA [Viral]" = features_viral,
  "RNA [Taxonomy BKN]" =  features_mtx_bkn,
  "RNA [Taxonomy MPA]" =  features_mtx_mpa,
  "DNA [PFAM]" = features_mgx_pfam,
  "RNA [PFAM]" = features_mtx_pfam,
  "RNA [Host Transcriptomics]" = htx,
  "CFU [Culture]" = culture_feature
  # "ABX [Medication]" = abx_feature
)

source('src/common/drop_samples.r')
drop_samples()

get_phylum_taxa <- function(){
  features_phyla <- readRDS("objects/amp/AMP_11_PhyloseqPhylum.rds")
  drop_lst <- get_drop_lst()
  
  # Make names actual phyla for visualization 
  new_names <- tax_table(features_phyla) %>% as.matrix() %>% as.data.frame() %>% pluck('Phylum')
  taxa_names(features_phyla) <- new_names
  features_phyla <- features_phyla %>% 
    phyloseq::otu_table() %>% 
    as.data.frame() %>% 
    rownames_to_column(var="feature") %>% 
    # rename_with(.fn = renaming_helper, .cols = everything(), map_df = map_amp) %>% 
    rename_with(.fn = function(x){str_split_i(x, "-", 1)})
  
  names(features_phyla) <- gsub(names(features_phyla), pattern = "PT", replacement = "") 
  
  features_phyla <- features_phyla %>% 
    select(!any_of(drop_lst))
  
  return(features_phyla)
}




# Add alpha
get_alpha <- function(){
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_amp %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_amp = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_amp_asv %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_amp_asv = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_mgx %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mgx_mpa = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_viral %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mgx_viral = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_mtx_bkn %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mtx_bkn = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_mtx_mpa %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mtx_mpa = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_mgx_pfam %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mgx_pfam = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(features_mtx_pfam %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_mtx_pfam = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(htx %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_htx = 'value'))
  permmeta <<- permmeta %>% left_join(vegan::diversity(t(culture_feature %>% column_to_rownames('feature')), index='shannon') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(shannon_culture = 'value'))
  
}
