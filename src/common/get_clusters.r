# source('src/common/load_metadata.r')
# source('src/common/get_distances.r')
library(tidyverse)
library(janitor)
# Takes map object used in load_distances and is  to 
# change to common sample value
renaming_helper_cluster <- function(to_map, map_df) {
  clean_cluster <- to_map %>% left_join(map_df, by=c('name'='from')) %>% dplyr::rename(sample_name = 'to') %>% dplyr::select(sample_name, cluster_num)
  return(clean_cluster)
}

mox_cluster_list <- list(
  "Amplicon" = readRDS("objects/amp/AMP_20_ConsensusClusters.rds") %>% renaming_helper_cluster(., map_df = map_amp)
)


get_pretty_name <- function(ugly_name){
  pretty_name_clust <- janitor::make_clean_names(paste0('cluster_num', '_', ugly_name))
  return(pretty_name_clust)
}
add_clusters <- function(my_meta=permmeta){
  cluster_names <- names(mox_cluster_list)
  for (i in seq_along(mox_cluster_list)){
    mox_cluster_list[[i]] <- mox_cluster_list[[i]] %>% 
      dplyr::rename_with(., ~ get_pretty_name(ugly_name=cluster_names[[i]]), .cols=starts_with('cluster_num'))
    my_meta <- my_meta %>% left_join(mox_cluster_list[[i]])
  }
  return(my_meta)
}
permmeta <- add_clusters()
