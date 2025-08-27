library(tidyverse)
source("src/common/load_metadata.r")
source("src/common/load_features.r")

# Tidy Data
dt <- mox_feature_list[["AMP [Genus]"]] %>% tidy_features() %>% left_join(permmeta)

# Calc prevalence + determine percentiles
dt <- dt %>% group_by(cluster_num_amplicon, feature) %>% 
    drop_na(cluster_num_amplicon) %>% 
    summarize(prevalence = sum(value>0),
              perc = sum(value>0)/n()) %>% 
    ungroup() %>% 
    group_by(feature) %>% 
    mutate(cv = sd(perc)/mean(perc),
           mm = max(perc) - min(perc)) %>% 
    ungroup() %>% 
    mutate(qmm = pnorm(scale(mm))) %>% 
    arrange(desc(mm))

ggplot(dt, aes(x=qmm))+
    geom_density() +
    labs(x="Max(abundance)-Min(abundance) Percentile") +
    theme_nature()
ggsave("figures/sanity_checks/max_min_abundance_percentile_distribution.pdf",units="in",height = 1.38, width = 2.41)
# Filter
# dt <- dt %>% filter(mm > quantile(mm, .95,)) 
dt <- dt %>% filter(mm > quantile(mm, .6) & mm < quantile(mm, 1)) 

pdf("figures/pneumotypes/prevalence_percentiled_heatmaps.pdf", height=11.33, width=4.33)
dt %>% pivot_wider(id_cols = feature, names_from = cluster_num_amplicon, values_from = perc) %>% 
    column_to_rownames("feature") %>% 
    pheatmap::pheatmap(show_rownames = T, scale = 'row',border_color = "black",number_color = "black", fontsize = 6, cluster_rows = T,
                       display_numbers = (dt %>% 
                                              pivot_wider(id_cols = feature, names_from = cluster_num_amplicon, values_from = prevalence) %>% 
                                              column_to_rownames("feature")
                                          )
    )

dev.off()                       
