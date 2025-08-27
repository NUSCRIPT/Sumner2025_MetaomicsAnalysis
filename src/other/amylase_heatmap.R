get_alpha()
# saved but need to change back to amylase i think 2024.09.25
amy_samples <- permmeta %>% 
    drop_na(amylase_bf_log, cluster_num_amplicon) %>% 
    arrange(amylase_bf_log, dysbiosis_score)  %>% 
    # filter(!sample_name %in% flip_lst) %>%
    pluck('sample_name')
amy_annot <-  permmeta %>% 
    filter(sample_name  %in% amy_samples) %>% 
    select(sample_name, dysbiosis_score, amylase_bf_log, cluster_num_amplicon) %>% 
    arrange(factor(sample_name, levels=amy_samples)) %>% 
    column_to_rownames('sample_name')
amy_features <- mox_maaslin_list$Clusters[["AMP [Genus]"]]$results %>%
    filter(coef>0) %>%
    filter(qval<0.05) %>% 
    arrange(desc(coef)) %>% 
    distinct(feature) %>% 
    pluck("feature") #%>% head(15)
amy_ord_features <- mox_feature_list[["AMP [Genus]"]] %>% 
    tidy_features() %>% 
    mutate(value = case_when(value == 0 ~ NA, .default = value)) %>% 
    wide_features() %>% 
    select(any_of(c('feature', amy_samples))) %>% 
    filter(make.names(feature) %in% amy_features) %>% 
    arrange(factor(feature, levels=amy_features)) %>% 
    mutate(feature = fix_asv(feature)) %>% 
    column_to_rownames('feature')
library(RColorBrewer)
# colors <- c(colorRampPalette(c('white','#a2bbd4', '#6a81a5'))(100))
library(ggplotify)
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(0.5, 5.5), c("white", "blueviolet"))
col_fun2 = colorRamp2(c(0,dysbiosis_cutoff, 1), c("white","#a9c3a4", "#7ea679")) 

column_ha = HeatmapAnnotation(Amylase = amy_annot$amylase_bf_log,
                              MDNP = amy_annot$dysbiosis_score,
                              annotation_name_gp = gpar(fontsize = 7),
                              col = list(Amylase = col_fun,
                                         MDNP = col_fun2),
                              border = T,
                              simple_anno_size = unit(.5, 'cm'),
                              annotation_legend_param = list(
                                  labels_gp = gpar(fontsize = 5),
                                  title_gp=gpar(fontface ="bold",fontsize= 7), 
                                  legend_height= unit(1, "cm"),
                                  # direction = "horizontal",
                                  title_position="leftcenter-rot",
                                  legend_width = unit(.25, "cm"),
                                  border='black')
                              )

pdf('figures/complex_amylase_heat.pdf',width = 4.75, height=2.28)
ht_list <- Heatmap(as.matrix(amy_ord_features), #scale(amy_ord_features,center = T,scale = T) ,
        name='taxa',
        show_column_names = FALSE, 
        row_dend_width = unit(.1, "cm"), 
        show_row_dend = F,
        cluster_columns = F,
        show_column_dend=F,
        cluster_rows = F,
        # clustering_distance_rows = "euclidean",#clustering_method_rows = "complete",
        row_names_gp = gpar(fontsize = 7),
        col = colorRamp2(c(log2(1.000001), log2(101)), c("#a2bbd4", "darkblue")),#viridis(100),#colors,
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 5),
                                    title="Abudance",
                                    title_gp=gpar(fontface ="bold",fontsize= 7), 
                                    legend_height= unit(1, "cm"),
                                    title_position="leftcenter-rot",
                                    direction = "vertical",
                                    legend_width = unit(.25, "cm"),
                                    border='black'),

        border_gp = gpar(col = "black", lty = 1,lwd=.7),
        na_col="white",
        top_annotation = column_ha)
draw(ht_list, merge_legend = T,legend_grouping = "original",heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
