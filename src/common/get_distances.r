
library(tidyverse)
library(vegan)
library(ade4)
library(ggpubr)
library(usedist)
source("src/common/nature_theme.r")
source("src/common/feature_tables_functions.r")
library(RColorBrewer)
# map_dt is two column map file with col names: from, to.
fix_dist <- function(dist_obj, map_dt) {
  
  # Clean map file to be only relevant subset
  dist_labels <- labels(dist_obj)
  
  map_dt <- map_dt %>% filter(from %in% dist_labels)
  
  # Reorder dist obj to be same as map file (from) & then rename (to)
  clean_dist_obj <- dist_obj %>%
    dist_subset(., map_dt$from) %>% 
    dist_setNames(., map_dt$to)
  
  return(clean_dist_obj)
}

# Load map obj
map_amp <- readRDS("objects/amp/AMP_14_BetaBrayCurtis.rds") %>% 
  labels(.) %>% 
  unlist() %>% 
  unique() %>% 
  as_tibble_col(column_name="from") %>% 
  separate(from, into = c("to", NA), sep="-", remove=FALSE) %>%
  mutate(to = str_replace(to, "PT", "")) %>%
  filter(str_detect(to, 'BAL')) # remove BRL (tech rep of BAL)


# Load distance obj 
amp <- readRDS("objects/amp/AMP_16_BetaGenusWUnifrac.rds") %>% as.dist() %>% fix_dist(., map_amp)
amp_asv <- readRDS("objects/amp/AMP_14_BetaBrayCurtis.rds") %>% as.dist() %>% fix_dist(., map_amp)

mgx_tax <- readRDS("objects/mgx/mgx_02_distance_Metaphlan.rds") 
mgx_pfam <- readRDS("objects/mgx/mgx_04_distance_HumannPFAM.rds") 

mtx_tax <- readRDS("objects/mtx/mtx_06_distance_Bracken.rds") 
mtx_mpa_tax <- readRDS("objects/mtx/mtx_02_distance_Metaphlan.rds")
mtx_pfam <- readRDS("objects/mtx/mtx_04_distance_HumannPFAM.rds") 

htx <- read_rds("objects/htx/HTX_01A_COMPFeatureTable.rds") %>% 
  dplyr::rename(feature = "Symbol") %>% 
  get_most_prevalent(min_prevalence = 0.1) %>% 
  normalize_feature_table(.,normalize=T,transform="none")
summary_table <- get_sample_types(htx)

htx_dist <- htx %>% 
  get_bal_samples() %>% 
  get_distance(method="jaccard")

source('src/common/get_viral.r')
viral <- get_viruses(distance = T)

# Load Culture
source("src/common/load_culture_feature_table.r")
culture_dists <- get_culture_feature_table("genus") %>% get_distance(method="jaccard") 

# # Load ABX
# source("src/common/load_antibiotic_table.r")
# abx_dists <- get_abx_table() %>% get_distance(method="gower") 

mox_dist_list <- list(
  "AMP [Genus]" = amp,
  "AMP [ASV]" = amp_asv,
  "DNA [Taxonomy]" = mgx_tax,
  "DNA [Viral]" = viral,
  "RNA [Taxonomy BKN]" =  mtx_tax,
  "RNA [Taxonomy MPA]" =  mtx_mpa_tax,
  # "DNA [Pathways]" = mgx_path,
  "DNA [PFAM]" = mgx_pfam,
  "RNA [PFAM]" = mtx_pfam,
  "RNA [Host Transcriptomics]" = htx_dist,
  "CFU [Culture]" = culture_dists
  # "ABX [Medication]" = abx_dists
)


source('src/common/drop_samples.r')
drop_samples()
