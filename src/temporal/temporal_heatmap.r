source("src/common/load_features.r")
require(tidyHeatmap)

gene_palette = circlize::colorRamp2(c(0.0,1e-20, .05, .1, .2, 6.7), c("black", rev(RColorBrewer::brewer.pal(5, "YlGnBu")))) #c("black", viridis::mako(5)))
taxa_palette = circlize::colorRamp2(c(0.0,1e-9, 2, 4,6,8), c("black", rev(RColorBrewer::brewer.pal(5, "YlGnBu")))) #c("black", viridis::mako(5)))

plot_temporal_heatmap <- function(feature_table, feature_palette, min_covariance=T, num_sd=2){
    
    # Make Feature Table into Metadata Table based on Naming convention
    get_timeslots = function(feature_table){
        md <- get_sample_types(feature_table) %>% 
            filter(sample_type!="CONBAL") %>% 
            separate(sample_name,sep = "BAL",into = c("patient", "day"),remove = F,convert = T) %>% 
            drop_na()%>% 
            mutate(baseline = factor(case_when(day == 0 ~ "baseline", day < 10 ~ "early", day>=10 ~ "late"), 
                                     levels=c("baseline", "early", "baseline_only","late")))
        
        longitudinal_pts <- md %>% filter(day!=0) %>% pluck("patient") %>% unique()
        
        md <- md %>% mutate(baseline = factor(case_when(!patient %in% longitudinal_pts ~ "baseline_only", 
                                                        .default = baseline),
                                              levels=c("baseline_only", "baseline", "early","late")),
                            longitudinal = case_when(patient %in% longitudinal_pts ~ "longitudinal", .default = "baseline only"))
        
        # md <- md %>% filter(patient %in% longitudinal_pts)
        return(md)
    }
    
    
    get_heatmap = function(feature_table, feature_palette, min_covariance=min_covariance, num_sd=num_sd, md){
        temporal_heatmap <- feature_table %>% 
            mutate(feature = str_replace_all(feature, "_", " ")) %>% 
            select(feature, any_of(md$sample_name)) %>%
            get_feature_covariance(., md, "day",min_covariance=min_covariance,absolute=T, num_sd=num_sd) %>% 
            
            # drop_zero_sum_features() %>% 
            # get_most_abundant_n(100) %>% 
            tidy_features() %>% 
            left_join(md) %>% 
            left_join(cov_table) %>% 
            mutate(covariance_direction = factor(case_when(covariance > 0 ~ "up",
                                                           covariance < 0 ~ "down"), levels=c("up","down"))
            ) %>% 
            arrange(day) %>%  
            mutate(#day = factor(day, levels = sort(unique(day))),
                sample_name = factor(sample_name, levels = unique(sample_name)),
                patient = factor(patient, levels = unique(patient)),
                feature = factor(feature, levels = covaried_features)) %>% 
            arrange(sample_name) %>%
            rename(`log2(abundance+1)` = "value") %>% 
            # mutate(value =case_when(value ==0 ~ NaN, .default = value)) %>%
            group_by(baseline,covariance_direction) %>% 
            heatmap(.column = sample_name,
                    .row = feature,
                    .value = `log2(abundance+1)`,
                    scale = "none",
                    cluster_rows = T, 
                    cluster_columns = T,
                    na_col = "black",
                    rect_gp = grid::gpar(col = "black", lwd = .5),
                    border = T,
                    # row_km = 4,
                    show_column_names = F,
                    column_title = "Samples",
                    row_title = "Features",
                    column_title_gp = grid::gpar(fontsize = 7),
                    row_title_gp = grid::gpar(fontsize = 7),
                    row_names_gp = grid::gpar(fontsize = 7, fontface="italic"),
                    show_column_dend = F,
                    show_row_dend = F,
                    # palette_value = circlize::colorRamp2(c(0.0,1e-20, .01, .1, 1, 2,6), c("black", rev(RColorBrewer::brewer.pal(6, "YlGnBu")))), #c("black", viridis::mako(5)))
                    palette_value = feature_palette, 
                    
                    palette_grouping = list(
                        # For second grouping (property_group)
                        c("#4C413FFF", "#4C413FFF"),
                        
                        # c("#b58b4c", "#74a6aa"),
                        # For first grouping (vs)
                        c("#4C413FFF", "#4C413FFF", "#4C413FFF","#4C413FFF")
                        
                        
                    ),
                    heatmap_legend_param = list(
                        legend_direction = "horizontal",
                        border ="black",
                        legend_width = unit(3, "cm"),
                        title_gp = grid::gpar(fontface="bold", fontsize = 8)
                    ),
                    height = length(covaried_features)*unit(2.7, "mm"),
                    width =  0.4*(ncol(feature_table)-1)*unit(2.7, "mm") # .6 & 1.5
                    
                    
            ) %>% 
            annotation_tile(day, 
                            # palette = RColorBrewer::brewer.pal(9, "YlOrRd"),
                            palette = 
                                circlize::colorRamp2(
                                    expm1(seq(log1p(0), log1p(max(md$day)), length.out = 9)),
                                    RColorBrewer::brewer.pal(9, "YlOrRd")),
                            
                            border=T, 
                            size = unit(.3,"cm"),
                            annotation_legend_param = list(
                                legend_direction = "horizontal",
                                border ="black",
                                legend_width = unit(3, "cm"),
                                title_gp = grid::gpar(fontface="bold", fontsize = 8)
                            ),
                            annotation_name_gp= grid::gpar(fontface="bold", fontsize = 6)
            ) %>% 
            annotation_tile(covariance,annotation_legend_side="bottom",
                            palette = 
                                circlize::colorRamp2(
                                    c(seq(min(cov_table$covariance), 0, length.out=5)[1:5],
                                      seq(0, max(cov_table$covariance), length.out=5)[2:5]),
                                    rev(RColorBrewer::brewer.pal(9, "Spectral"))
                                ),
                            
                            border=T, 
                            size = unit(.3,"cm"),
                            annotation_legend_param = list(
                                legend_direction = "horizontal",
                                border ="black",
                                legend_width = unit(3, "cm"),
                                title_gp = grid::gpar(fontface="bold", fontsize = 8)
                            ),
                            annotation_name_gp= grid::gpar(fontface="bold", fontsize = 6)) 
        
        
        return(temporal_heatmap)
        
    }
    
    md <- get_timeslots(feature_table)
    temporal_heatmap <- get_heatmap(feature_table, feature_palette, min_covariance, num_sd, md)
    
    temporal_heatmap <- ComplexHeatmap::draw(temporal_heatmap %>% as_ComplexHeatmap,
                                             heatmap_legend_side="bottom",
                                             annotation_legend_side="top",
                                             merge_legend=T)
    
    return(temporal_heatmap)
}


#Genes
pdf("figures/temporal/covariance_features_heatmap.pdf", width =13, height=15.5, onefile = T)

th <- plot_temporal_heatmap(mox_feature_list[["RNA [PFAM]"]], gene_palette, min_covariance = T, num_sd = 2)
# th %>% tidyHeatmap::save_pdf("figures/temporal/covariance_features_rna_pfam_heatmap.pdf", units="in", width=9.23, height=8)

th <- plot_temporal_heatmap(mox_feature_list[["DNA [PFAM]"]], gene_palette, min_covariance = T, num_sd = 3)
# th %>% save_pdf("figures/temporal/covariance_features_rna_pfam_heatmap.pdf", units="in", width=12.30, height=3.80)

# Taxa
th <- plot_temporal_heatmap(mox_feature_list[["AMP [ASV]"]], taxa_palette, min_covariance = T, num_sd = 2)
# th %>% save_pdf("figures/temporal/covariance_features_asv_heatmap.pdf", units="in", width=12.61, height=3.05)

th <- plot_temporal_heatmap(mox_feature_list[["AMP [Genus]"]], taxa_palette, min_covariance = T, num_sd = 2)
th <- plot_temporal_heatmap(mox_feature_list[["DNA [Taxonomy]"]], taxa_palette, min_covariance = T, num_sd = 2)
th <- plot_temporal_heatmap(mox_feature_list[["RNA [Taxonomy BKN]"]], taxa_palette, min_covariance = T, num_sd = 2)
th <- plot_temporal_heatmap(mox_feature_list[["RNA [Taxonomy MPA]"]], taxa_palette, min_covariance = T, num_sd = 2)
th <- plot_temporal_heatmap(mox_feature_list[["RNA [Host Transcriptomics]"]], gene_palette, min_covariance = T, num_sd = 3)
dev.off()

