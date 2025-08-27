library(tidyHeatmap)
source("src/common/load_features.r")

# Get top N feautures
dsa <- lapply(mox_feature_list, get_most_abundant_n,n_features = 20)

# Exclude Viral and only include one RNA taxonomic
dsa <- dsa[!names(dsa) %in% c("RNA [Taxonomy MPA]",  "DNA [Viral]")] 

# Tidy and make into super-column
dsa <- lapply(dsa, tidy_features)
dsa <- dsa %>% bind_rows(.id = "groups") 
dsa <- dsa %>% pivot_wider(names_from = sample_name, values_from = value, values_fill = NaN) %>% pivot_longer(3:last_col())

head(dsa)

# Get intersection of samples in different omics and sort by intersecion info
omics <- dsa %>% 
  filter(!is.nan(value)) %>% 
  distinct(name, groups) %>% 
  arrange(groups) %>% 
  group_by(name) %>% 
  reframe(ome = list(groups)) %>% 
  rowwise() %>% 
  filter(length(ome)>1) %>% 
  mutate(ome = paste(ome, collapse=" ")) %>% 
  ungroup() %>% 
  add_count(ome) %>% 
  filter(n>10) %>% # AT LEAS TEN SAMPLES TO INCLUDE
  arrange(desc(n),ome) %>%
  group_by(ome) %>% 
  mutate(var_temp = ifelse(row_number()==1,1,0)) %>% 
  ungroup() %>% 
  mutate(ome_rank = cumsum(var_temp))

# Make as bool
bool_table <- dsa %>% 
  filter(!is.nan(value)) %>% 
  distinct(groups, name) %>% 
  mutate(value = factor("Sample available", levels=c("Sample available","Sample unavailable"))) %>% 
  pivot_wider(names_from = groups, values_from = value, values_fill = "Sample unavailable") %>% 
  distinct()

# Combine
omics <- omics %>% left_join(bool_table)

# Factors and things for plotting
dsa <- dsa %>%
  mutate(`profile type` = factor(case_when(str_detect(groups, "AMP|Taxonomy") ~ "Taxonomic", 
                                      str_detect(groups, "PFAM|KO") ~ "Microbial Gene-Level",
                                      str_detect(groups, "Host Transcriptomics") ~ "Host Gene-Expression"), levels = c("Taxonomic", "Microbial Gene-Level", "Host Gene-Expression"))
               ) %>% 
  mutate(`'omic type` = factor(case_when(str_detect(groups, "AMP") ~ "Amplicon", 
                                                        str_detect(groups, "DNA") ~ "Metagenomic",
                                                        str_detect(groups, "RNA") ~ "Total RNA"))
         ) %>% 
  mutate(groups = factor(groups, levels = c("AMP [Genus]", 
                                            "AMP [ASV]", 
                                            "DNA [Taxonomy]", 
                                            # "DNA [Viral]", 
                                            "RNA [Taxonomy BKN]", 
                                            # "RNA [Taxonomy MPA]", 
                                            "DNA [PFAM]", 
                                            "RNA [PFAM]", 
                                            "RNA [Host Transcriptomics]"))
         )


annot_unit <- unit(.2, "cm")
bool_colors <- list("Sample available" = "#4C413FFF","Sample unavailable" = "#999999")
bool_colors <- c("#4C413FFF","#999999")
overview_heat <- dsa %>% 
  left_join(omics) %>% 
  drop_na(ome) %>% # drop valid section in other dropped groups
  drop_na() %>% 
  arrange(desc(n), ome, ) %>% 
  mutate(feature = case_when(groups == "AMP [Genus]" ~ str_replace(feature, "ASV[1234567890]+_", ""), .default = feature),
         feature = str_replace_all(feature, "_", " "),
         feature = str_replace_all(feature, "\\$", " ")) %>% 
  mutate(feature = paste(feature, "[", as.integer(groups), "]"),
         ome = factor(ome, levels=unique(ome)),
         ome_rank = factor(ome_rank, levels=unique(ome_rank))) %>% 
  # left_join(bool_table) %>% 
  group_by(groups, ome_rank) %>% 
  heatmap(.column = name,
          .row = feature,
          .value = value,
          scale = "row",
          cluster_rows = T, cluster_columns = T,
          na_col = "grey",
          # rect_gp = grid::gpar(col = NA, lwd = 0),
          border = T,
          show_column_names = F,
          column_title = "Samples",
          row_title = "Multi'omic features",
          # row_names_gp = grid::gpar(fontsize = 4),
          show_column_dend = F,
          show_row_dend = F,
          palette_grouping = list(
            c('#4C413FFF', '#4C413FFF','#5A6F80FF','#D8AF39FF','#5A6F80FF', '#D8AF39FF','#E8C4A2FF')
          )) %>% 
  annotation_tile(`'omic type`, palette = c("#4C413FFF", "#5A6F80FF", "#D8AF39FF"), border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5))  %>% 
  annotation_tile(`profile type`, palette = c("#999999", "#278B9AFF", "darkgreen"),border=T, size = annot_unit, annotation_name_gp= grid::gpar(fontsize = 5))
overview_heat

overview_heat_bool <- overview_heat %>% 
  annotation_tile(`AMP [ASV]`, palette = bool_colors, border=T, size = annot_unit, annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`AMP [Genus]`, palette = bool_colors, border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`DNA [Taxonomy]`, palette = bool_colors, border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`RNA [Taxonomy BKN]`, palette = rev(bool_colors), border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  # annotation_tile(`RNA [Taxonomy MPA]`, palette = rev(bool_colors), border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`DNA [PFAM]`, palette = bool_colors, border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`RNA [PFAM]`, palette = rev(bool_colors), border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5)) %>% 
  annotation_tile(`RNA [Host Transcriptomics]`, palette = bool_colors, border=T, size = annot_unit,annotation_name_gp= grid::gpar(fontsize = 5), show_legend=F)
  
overview_heat
overview_heat_bool
overview_heat %>% save_pdf("figures/overview_omics/multiomic_heatmap_top20_long.pdf", unit = "in", height=8.63, width=10.28)
overview_heat_bool %>% save_pdf("figures/overview_omics/multiomic_heatmap_top20_bool_long.pdf", unit = "in", height=9.91, width=10.89)
