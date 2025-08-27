source('src/common/load_metadata.r')
pdf("figures/sanity_checks/metadata_completeness.pdf", height = 25, width = 12)

permmeta %>% pivot_longer(2:last_col(),values_transform = as.character) %>% 
    mutate(value = case_when(is.na(value) ~ 0,
                             !is.na(value) ~ 1)) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 
    column_to_rownames('sample_name') %>% 
    pheatmap::pheatmap(.,cluster_rows = F, cluster_cols = T,
                       display_numbers = permmeta %>% column_to_rownames('sample_name'), 
                       fontsize_number = 2,fontsize_row = 5,
                       color = colorRampPalette(c("white", "grey"))(2))
                          
dev.off()
