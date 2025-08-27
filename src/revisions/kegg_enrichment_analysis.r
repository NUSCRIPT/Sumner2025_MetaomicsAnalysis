library(tidyverse)
library(clusterProfiler)


run_kegg_enrichment <- function(omic, meta, meta_value=NULL, direction = "up"){
    genes <- mox_maaslin_list[[meta]][[omic]]$results %>%
        filter(qval < 0.05, coef > 0) 
    
    if (!is.null(meta_value)){
        genes <- genes %>% filter(value == meta_value)
    }
    
    genes <- genes %>% 
        pluck("feature") %>% 
        str_split_i(., ":|\\.", i = 1)
    print(genes)
    
    universe <- mox_maaslin_list[[meta]][[omic]]$results %>% 
        pluck("feature") %>% 
        sort() %>% 
        unique() %>% 
        str_split_i(., ":|\\.", i = 1) #%>% str_remove(.,"K") %>% as.integer()
    universe
    
    # universe <- mox_feature_list[["DNA [PFAM]"]] %>% 
    #     get_most_prevalent() %>%
    #     pluck("feature") %>% str_split_i(., ":|\\.", i = 1) #%>% str_remove(.,"K") %>% as.integer()
    
    ego_kegg <- enrichKEGG(gene = genes,
                           universe = universe,,
                           pAdjustMethod = "fdr",
                           pvalueCutoff = .05,
                           qvalueCutoff = .05,
                           keyType = "kegg",
                           # minGSSize = 3,
                           # maxGSSize = 500,
                           organism = "ko") %>% 
        data.frame() %>% 
        arrange(FoldEnrichment) %>% 
        mutate(Description = factor(Description, levels=unique(Description))) %>% 
        arrange(Description) %>% 
        mutate(Omic = omic,
               Meta = meta,
               MetaValue = meta_value,
               .before = 1)
    return(ego_kegg)
}

ego_2 <- run_kegg_enrichment("DNA [PFAM]", "Clusters",meta_value = "2")
ego_3 <- run_kegg_enrichment("DNA [PFAM]", "Clusters", meta_value="3")
ego_4 <- run_kegg_enrichment("DNA [PFAM]", "Clusters",meta_value = "4")

clusters <- do.call("rbind", list(ego_2,ego_3,ego_4)) %>% 
    mutate(annotation = gtools::stars.pval(qvalue)) %>% 
    mutate(MetaValue = case_when(MetaValue == "2" ~ "M",
                     MetaValue == "3" ~ "SP",
                     MetaValue == "4" ~ "OL")) %>% 
    rename(Pneumotype = "MetaValue")



ggplot(clusters %>% 
           complete(Pneumotype,Description) %>%
           arrange(Description, category) %>% 
           fill(category) %>% 
           mutate(median_q = median(qvalue)) %>% 
           arrange(median_q), 
       aes(x=Pneumotype, y=Description, fill=FoldEnrichment, label = annotation)) +
    geom_tile(color = "black") +
    geom_text() +
    ggpubr::theme_pubclean() +
    scale_fill_viridis_c(option = "plasma",na.value = '#999999')+
    facet_grid(rows = vars(category),scales = "free_y", 
               space = "free_y",
               labeller = label_wrap_gen(width = 2)) +
    theme(strip.text = element_text(size=5, face="bold"),
          # strip.clip = "off",
          strip.background=element_rect(linewidth=.2, fill="grey",color="black"))

ggsave("figures/kegg_enrichment/enrichment_analysis.pdf", height=8.28, width=6.44, units ="in")
write_tsv(clusters, "tables/maaslin_enriched_clusters.tsv")



ego_dys <- run_kegg_enrichment("DNA [PFAM]", "Dysbiotic",meta_value = NULL) %>% 
    mutate(annotation = gtools::stars.pval(qvalue)) %>% 
    mutate(Meta = "MDNP>90th")

ggplot(ego_dys,
       aes(x=Meta, y=Description, fill=FoldEnrichment, label = annotation)) +
    geom_tile(color = "black") +
    geom_text() +
    ggpubr::theme_pubclean() +
    scale_fill_viridis_c(option = "plasma",na.value = '#999999')+
    facet_grid(rows = vars(category),scales = "free_y", 
               space = "free_y",
               labeller = label_wrap_gen(width = 2)) +
    theme(strip.text = element_text(size=5, face="bold"),
          # strip.clip = "off",
          strip.background=element_rect(linewidth=.2, fill="grey",color="black")) 
ggsave("figures/kegg_enrichment/dysbiosis_ko_enrichment_analysis.pdf", height=8.28, width=6.44, units ="in")
write_tsv(ego_dys, "tables/maaslin_enriched_dysbiotic.tsv")


ego_qpcr <- run_kegg_enrichment("DNA [PFAM]", "QPCR",meta_value = NULL) %>% 
    mutate(annotation = gtools::stars.pval(qvalue)) %>% 
    mutate(Meta = "Bacterial Load (16S copies)")

ggplot(ego_qpcr,
       aes(x=Meta, y=Description, fill=FoldEnrichment, label = annotation)) +
    geom_tile(color = "black") +
    geom_text() +
    ggpubr::theme_pubclean() +
    scale_fill_viridis_c(option = "plasma",na.value = '#999999')+
    facet_grid(rows = vars(category),scales = "free_y", 
               space = "free_y",
               labeller = label_wrap_gen(width = 2)) +
    theme(strip.text = element_text(size=5, face="bold"),
          # strip.clip = "off",
          strip.background=element_rect(linewidth=.2, fill="grey",color="black")) 
ggsave("figures/kegg_enrichment/qpcr_ko_enrichment_analysis.pdf", height=8.28, width=6.44, units ="in")
write_tsv(ego_qpcr, "tables/maaslin_enriched_qpcr.tsv")



ego_culture<- run_kegg_enrichment("DNA [PFAM]", "Culture",meta_value = NULL) %>% 
    mutate(annotation = gtools::stars.pval(qvalue)) %>% 
    mutate(Meta = "Bacterial Culture Positive")

ggplot(ego_culture,
       aes(x=Meta, y=Description, fill=FoldEnrichment, label = annotation)) +
    geom_tile(color = "black") +
    geom_text() +
    ggpubr::theme_pubclean() +
    scale_fill_viridis_c(option = "plasma",na.value = '#999999')+
    facet_grid(rows = vars(category),scales = "free_y", 
               space = "free_y",
               labeller = label_wrap_gen(width = 2)) +
    theme(strip.text = element_text(size=5, face="bold"),
          # strip.clip = "off",
          strip.background=element_rect(linewidth=.2, fill="grey",color="black")) 
ggsave("figures/kegg_enrichment/culture_ko_enrichment_analysis.pdf", height=8.28, width=6.44, units ="in")
write_tsv(ego_culture, "tables/maaslin_enriched_culture.tsv")


# ego_amy <- run_kegg_enrichment("DNA [PFAM]", "Amylase",meta_value = NULL) %>% 
#     mutate(annotation = gtools::stars.pval(qvalue)) %>% 
#     mutate(Meta = "MDNP>90th")
