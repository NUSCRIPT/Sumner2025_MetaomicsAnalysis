require(ggpubr)
source("src/common/nature_theme.r")
source("src/common/load_features.r")
dsa <- lapply(mox_feature_list, get_bal_samples)
dsa <- lapply(dsa, get_most_abundant_n,n_features = 20)
dsa <- lapply(dsa, tidy_features)
dsa <- dsa %>% bind_rows(.id = "groups") 
dsa <- dsa %>% pivot_wider(names_from = sample_name, values_from = value, values_fill = NaN) %>% pivot_longer(3:last_col())

dsa_lollipop <- dsa %>% 
  mutate(groups = factor(groups, levels = c("AMP [ASV]", "AMP [Genus]", "DNA [Taxonomy]", "RNA [Taxonomy MPA]", 
                                            "RNA [Taxonomy BKN]", "DNA [PFAM]", "RNA [PFAM]","RNA [Host Transcriptomics]"))
  ) %>% 
  mutate(feature = case_when(groups == "AMP [Genus]" ~ str_replace(feature, "ASV[1234567890]+_", ""), .default = feature),
         feature = str_replace_all(feature, "_", " "),
         feature = str_replace_all(feature, "\\$", " ")) %>% 
  mutate(feature = paste0("[", as.integer(groups), "] ", feature)) %>% 
  group_by(feature,groups) %>% 
  filter(!is.nan(value)) %>% 
  mutate(mean = mean(value),
         prevalence = sum(value > 0)/n()) %>% 
  distinct(groups, feature, mean, prevalence) %>% 
  pivot_longer(3:4) %>% 
  ungroup() %>% 
  mutate(groups = factor(groups, levels=unique(groups)),
         name = factor(name, levels=unique(name)),
  ) %>% 
  
  arrange(groups,feature,name) %>% 
  filter(name=="mean") %>% 
  ungroup() %>%
  arrange(value) %>% 
  mutate(feature = factor(feature, levels=unique(feature))) %>%
  arrange(feature)

taxa_lolli <- dsa_lollipop %>% 
  filter(!groups %in% c("DNA [PFAM]", "RNA [PFAM]","RNA [Host Transcriptomics]")) %>% 
  # mutate(feature = str_sub(feature, end=-7)) %>% 
  ggplot(aes(x=value, y=feature)) + 
  geom_errorbarh(aes(xmin = 0, xmax = value), color = "black", height=0) +
  geom_point(aes(size=value, fill=value), shape=21) +
  scale_fill_viridis_b() +
  ggh4x::facet_wrap2(~groups,drop = T,axes = "all",scales = "free_y", ncol=6) +
  theme_nature() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_text(face="bold", size=9),
        axis.text.y = element_text(face="italic")) +
  labs(x="Log2 Mean Abundance", y = "Features", title = "Top 20 features by mean abundance in taxnomic profiles")
taxa_lolli
ggsave("figures/overview_omics/top20_taxa_lollipop.pdf", taxa_lolli, units="in", width=11.67, height=2.76)

gene_lolli <- dsa_lollipop %>% 
  filter(groups %in% c("DNA [PFAM]", "RNA [PFAM]","RNA [Host Transcriptomics]")) %>% 
  # mutate(feature = str_sub(feature, end=-7)) %>% 
  ggplot(aes(x=value, y=feature)) + 
  geom_errorbarh(aes(xmin = 0, xmax = value), color = "black", height=0) +
  geom_point(aes(size=value, fill=value), shape=21) +
  scale_fill_viridis_b() +
  ggh4x::facet_wrap2(~groups,drop = T,axes = "all",scales = "free_y", ncol=1) +
  theme_nature() +
  theme(plot.title = element_text(hjust = .85),
        strip.text = element_text(face="bold", size=9),
        axis.text.y = element_text(face="italic")) +
  labs(x="Log2 Mean Abundance", y = "Features", title = "Top 20 features by mean abundance in microbial- and host-gene profiles")
gene_lolli
ggsave("figures/overview_omics/top20_gene_lollipop.pdf", gene_lolli, units="in", width=4.36, height=5.99)

#######################################################################################################################3
# Boxplots 


plot_feature <- function(feature_table){
  feature_table %>% 
    tidy_features() %>% 
    mutate(
           feature = str_replace_all(feature, "_", " "),
           feature = str_replace_all(feature, "\\$", " ")) %>% 
    group_by(feature) %>% 
    mutate(average =mean(value)) %>% 
    ungroup() %>% 
    arrange(average,feature) %>% 
    mutate(feature = factor(feature, levels = unique(feature))) %>% 
    arrange(feature) %>% 
    ggplot(aes(x=value, y=feature)) +
    geom_jitter(width = 0,height = .2, alpha=.5, color="#5A6F80FF") +
    geom_boxplot(fill="darkgrey",width=.5, outlier.shape = NA,alpha=.8) +
    theme_nature() +
    rotate_x_text() +
    labs(x="Log2 abundance+1", y="Feature")+
    theme(plot.title = element_text(hjust = 0.5)
          # axis.text.y = element_text(face="italic")
          ) +
    xlim(0,6.7)
}


pad_features <- function(feature_table){
  feature_table <- feature_table %>% 
    mutate(feature = str_replace_all(feature, ":", " "),
           feature = str_replace_all(feature, "-", " "),
           feature = str_pad(feature, width=100, side="left"))
  return(feature_table)
}

dsa <- lapply(mox_feature_list, get_bal_samples)
dsa <- lapply(dsa, get_most_abundant_n,n_features = 50)

dsa <- lapply(dsa, pad_features)
dsa[["AMP [Genus]"]] <- dsa[["AMP [Genus]"]] %>% mutate(feature = str_replace(feature, "ASV[1234567890]+_", ""))
top_plots <- lapply(dsa, plot_feature)

top_plots[["AMP [ASV]"]] + 
  labs(title = "Amplicon\n(ASV-Level)") 

ggsave('figures/overview_omics/jitter_boxplot_amplicon_asv.pdf', units = "in", height=6.56, width=4.45)

top_plots[["AMP [Genus]"]] + 
  labs(title = "Amplicon\n(Genus-Level)") 

ggsave('figures/overview_omics/jitter_boxplot_amplicon_genus.pdf', units = "in", height=6.56, width=4.09)

top_plots[["DNA [Taxonomy]"]] + 
  labs(title = "Metagenomics\n(Taxonomy)") 

ggsave('figures/overview_omics/jitter_boxplot_mgx_taxonomy.pdf', units = "in", height=6.56, width=4.40)

top_plots[["DNA [PFAM]"]] + 
  labs(title = "DNA\n(Microbial Gene-Level)") 

ggsave('figures/overview_omics/jitter_boxplot_mgx_gene.pdf', units = "in", height=6.56, width=4.82)

top_plots[["RNA [Taxonomy MPA]"]] + 
  labs(title = "RNA\n(MetaPhlAn)") 

ggsave('figures/overview_omics/jitter_boxplot_mtx_taxonomy.pdf', units = "in", height=6.56, width=4.69)

top_plots[["RNA [Taxonomy BKN]"]] + 
  labs(title = "RNA\n(Bracken)") 

 ggsave('figures/overview_omics/jitter_boxplot_mtx_bracken.pdf', units = "in", height=6.56, width=4.41)

top_plots[["RNA [PFAM]"]] + 
  labs(title = "RNA\n(Microbial Gene-Level)") 

ggsave('figures/overview_omics/jitter_boxplot_mtx_gene.pdf', units = "in", height=6.56, width=4.71)

top_plots[["RNA [Host Transcriptomics]"]] + 
  labs(title = "RNA\n(Host Gene-Level)") 

ggsave('figures/overview_omics/jitter_boxplot_htx.pdf', units = "in", height=6.56, width=4.09)


