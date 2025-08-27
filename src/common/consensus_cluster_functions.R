library(ConsensusClusterPlus)

# PURPOSE: wrapper for performing consensus clustering on metaomic data reproducibly and consistently
# <sample_distances> = type dist, distance object for clustering
# <plot> = type NULL or string, NULL if you do not want to save plots, "pdf" if you do. must have title if saving
# <title> = type string, required if plotting
# RETURNS consensus cluster report  object
run_consensus_clustering <- function(sample_distances, plot=NULL, file_name="consensus_clusters"){
  consensus_report = ConsensusClusterPlus(sample_distances, #, upper = FALSE, diag = FALSE
                                maxK=10,
                                reps=1000,
                                pItem=.8,
                                pFeature=1,
                                clusterAlg="pam",
                                plot = plot,
                                title = file_name,
                                seed=521)
  return(consensus_report)
}

# PURPOSE: cuts consensus cluster output at specified k clusters 
# <consensus_report> = type consensus cluster report object
# <nc> = type int, number of clusters to cut to
# RETURNS sample cluster assignment as tibble0
get_consensus_clusters <- function(consensus_report, nc){
  clusters <- consensus_report[[nc]]$consensusClass %>% 
    enframe() %>% 
    dplyr::rename(cluster_num = "value",
                  sample_name = "name") %>%
    mutate(cluster_num = factor(cluster_num))
  
  return(clusters)
}

# PURPPOSE: get cladogram from consensus clustering
# <consensus_report> = type consensus cluster report object
# <nc> = type int, number of clusters to cut to
# RETURNS sample cluster tree as hclust bject
get_consensus_tree <- function(consensus_report, nc){
  clust.col <- consensus_report[[nc]]$consensusTree
  
}

# PURPSE: get percentile ranked features associated with each cluster ; threshold 90th percentile
# <feature_table>
# <cluster_table> = output of get_consensus_clusters
get_feature_cluster_ranks <- function(feature_table, cluster_table){
  ranks <- feature_table %>% 
    tidy_features() %>% 
    left_join(cluster_table) %>% 
    group_by(cluster_num, feature) %>% 
    summarise(means = mean(value)) %>% 
    mutate(rank = percent_rank(means)) %>% 
    arrange(cluster_num, desc(rank)) %>% 
    ungroup()
  
  high_features <<- ranks %>% filter(rank > 0.9) %>% pluck("feature") %>% unique()
  
  ranks <- ranks %>% 
    filter(feature %in% high_features) %>%
    mutate(rank = case_when(rank > 0.9 ~ 1, rank <=.9 ~ 0)) %>% # dim reduce
    select(!means) %>%  
    pivot_wider(names_from = cluster_num, values_from = rank) %>% 
    mutate(sum = rowSums(across(where(is.numeric))))  %>% 
    arrange(sum)
  
  return(ranks)
  
}

# # # EXAMPLE USE CASE: 
# dist_mpa_mgx <- bkn %>% get_sample_subset() %>% get_distance("bray")
# clusters <- dist_mpa_mgx %>%run_consensus_clustering(plot = "pdf", file_name ="consensus_cluster_results/mpa_mgx/")
# cluster_table <- get_consensus_clusters(cluster_list[["RNA [BRACKEN]"]],3)
# mox_dist_list[["RNA [BRACKEN]"]] %>%my_ordination_helper(my_meta = cluster_table) %>%  visualize_ordination("cluster_num", "cluster_num",3,.5)
# get_feature_cluster_ranks(get_sample_subset(bkn), cluster_table) %>% View()
# 
# mpa_mgx %>% tidy_features() %>% filter(sample_name %in% (cluster_table %>% filter(cluster_num==2) %>% pluck("sample_name"))) %>% group_by(feature) %>% summarise(mean = mean(value)) %>% arrange(desc(mean)) %>% filter(mean >0)%>% View()