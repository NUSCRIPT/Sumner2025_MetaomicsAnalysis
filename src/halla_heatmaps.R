library(tidyverse)
library(viridisLite)
library(RColorBrewer)
source("src/common/nature_theme.r")

# Prepare to make omic names pretty
mox_list <- c("AMP [ASV]", "DNA [Taxonomy]", "DNA [Viral]",
              "RNA [Taxonomy BKN]", "DNA [PFAM]", "RNA [PFAM]",
              "RNA [Host Transcriptomics]", "CFU [Culture]", "ABX [Medication]", "Metadata")

type_map <- list(A = mox_list,
                 B = janitor::make_clean_names(mox_list)) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

# Read in data and make pretty names
df <-read_tsv("tables/network/halla_network_full.tsv") %>% 
  left_join(type_map, by = c("Source_Type" = "B")) %>% 
  mutate(Source_Type = A) %>% 
  select(!A) %>%   
  left_join(type_map, by = c("Target_Type" = "B")) %>% 
  mutate(Target_Type = A) %>% 
  select(!A) 

colnames(df)


# Get the following counts per comparison:
# (1) Number of clusters
# (2) Number of significant associations (feature-features)

cluster_n <- df %>% 
  distinct(Source_Type, Target_Type, Cluster) %>% 
  count(Source_Type, Target_Type) %>%
  mutate(Source_Type = factor(Source_Type,levels=mox_list), 
         Target_Type=factor(Target_Type, levels=(mox_list))) %>% 
  complete(Source_Type,Target_Type,fill = list(n=NaN)) %>% # Fill triangle w NA
  arrange(Source_Type, Target_Type) %>% 
  filter(as.integer(Source_Type) <= as.integer(Target_Type)) %>% 
  mutate(n = ifelse(Source_Type == Target_Type, NA, n)) %>% 
  mutate(n_feature = ifelse(is.na(n), "N/A", as.character(n))) %>% 
  mutate(n_feature = ifelse(Source_Type == Target_Type, "", n_feature))


pairwise_n <- df %>% 
  count(Source_Type, Target_Type) %>% 
  mutate(Source_Type = factor(Source_Type,levels=mox_list), 
         Target_Type=factor(Target_Type, levels=(mox_list))) %>% 
  complete(Source_Type,Target_Type,fill = list(n=NaN)) %>% # Fill triangle w NA
  arrange(Source_Type, Target_Type) %>% 
  filter(as.integer(Source_Type) <= as.integer(Target_Type)) %>% 
  mutate(n = ifelse(Source_Type == Target_Type, NA, n)) %>% 
  mutate(n_feature = ifelse(is.na(n), "N/A", as.character(n))) %>% 
  mutate(n_feature = ifelse(Source_Type == Target_Type, "", n_feature))


# Heatmap helper function for plotting count tables
counts_heatmap <- function(count_table){
  
  # Get fill colors and text colors
  set.seed(101)
  alpha <- 1.8
  beta <- .1
  colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
  count_table$lblcolor <- ifelse(pbeta(percent_rank(count_table$n), alpha, beta) < .1, "black", "white")
  
  count_heat <- 
    ggplot(count_table, aes(x=Source_Type, y=Target_Type, fill = n, label=n_feature)) +
    geom_tile(color = "black", linewidth = .3) + #.9
    # geom_text(aes(color = lblcolor), size=6/.pt)+
    geom_text(aes(label=n_feature, color=lblcolor), size=8/.pt, nudge_y = 0) +
    guides(fill = guide_colourbar(label = TRUE,
                                  ticks = TRUE,
                                  title = "",
                                  frame.linewidth=.5,
                                  frame.colour="black",
                                  ticks.colour = "black"),
           color = 'none') +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(position="left",expand=c(0,0)) +
    scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100),
                         values=colorvalues,
                         na.value="white",
                         # limits=c(0, 1),
                         name="No. Associations") +
    scale_color_manual(values=c(white="white", black="black")) +
    theme_nature() +
    theme(legend.position = "right",
          legend.text = element_text(size=8),
          axis.text = element_text(size=8),
          axis.title = element_text(size=8,face = "bold"),
          axis.text.x=element_text(angle = 335, hjust = 0),
          panel.border = element_rect(colour = "black", fill=NA,linewidth = .75),
          plot.title = element_text(size=10, hjust = 0.5, margin = margin(0,0,0,0)),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()
    ) +
    xlab("Source") +
    ylab("Target") 
  return(count_heat)
}


# Plot Upper Triangle
pairwise_heat <- counts_heatmap(pairwise_n)+
  ggtitle("No. Feature Associations")

ggsave("figures/network_and_mox/halla_pairwise_features.pdf",pairwise_heat,units="in",width=5.47, height=2.65)

cluster_heat <- counts_heatmap(cluster_n)+
  ggtitle("No. Clusters")

ggsave("figures/network_and_mox/halla_number_clusters.pdf",cluster_heat, units="in",width=5.47, height=2.65)

library(patchwork)
combined_heat<- (cluster_heat + theme(axis.title.x=element_blank(),
                                      axis.text.x=element_blank())) / pairwise_heat

ggsave("figures/network_and_mox/halla_stacked_heats.pdf",combined_heat, units="in",width=5.78, height=4.72)



# Plot Lower Triangle
pairwise_heat <- pairwise_n %>% 
  mutate(tmp = Source_Type, Source_Type = Target_Type, Target_Type = tmp) %>%  # Swap Triangle
  select(!tmp) %>% 
  counts_heatmap()+
  ggtitle("No. Feature Associations") +
  scale_y_discrete(position="right",expand=c(0,0))  +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02,.8))

ggsave("figures/network_and_mox/halla_pairwise_features_lower.pdf",pairwise_heat,units="in",width=5.47, height=2.65)

cluster_heat <- cluster_n %>% 
  mutate(tmp = Source_Type, Source_Type = Target_Type, Target_Type = tmp) %>%  # Swap Triangle
  select(!tmp) %>% 
  counts_heatmap()+
  ggtitle("No. Clusters") +
  scale_y_discrete(position="right",expand=c(0,0)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02,.8))

ggsave("figures/network_and_mox/halla_number_clusters_lower.pdf",cluster_heat, units="in" ,width=5.47, height=2.65)

combined_heat<- cluster_heat + pairwise_heat + plot_layout(ncol=2,axes = "collect_y", axis_titles = "collect") & theme(plot.margin = unit(c(1,5,1,1), "pt"))
combined_heat
ggsave("figures/network_and_mox/halla_stacked_heats_lower.pdf",combined_heat, units="in",width=8.89, height=2.69)



# Plot bothin one square
rbind(cluster_n,
      pairwise_n %>%   
        mutate(tmp = Source_Type, Source_Type = Target_Type, Target_Type = tmp) %>%  # Swap Triangle
        select(!tmp)) %>% 
  counts_heatmap() +
  ggtitle("No. Clusters") +
  geom_text(x = 10, y=0,# Set text's position to the right end of the plot
            hjust = .1,
            label="No. Feature Associations",
            angle=90,
            vjust=3,
            family="Helvetica",
            fontface="bold",
            size = 7/.pt)+
  coord_cartesian(clip = 'off') +   # This keeps the labels from disappearing
  theme(plot.margin = unit(c(1,3,1,1), "pt"),
        panel.border = element_rect(colour = "black", fill=NA,linewidth = .35),
        legend.text = element_text(size = 6, margin = margin(1,1,1,1)),
        legend.margin = margin(1,1,1,10),
        axis.text.x=element_text(angle = 25, hjust = 1)) # This widens the right margin
ggsave("figures/network_and_mox/halla_number_both_clusters_features.pdf",units="in",width=4.06, height=1.95)

ggsave("figures/network_and_mox/halla_number_both_clusters_features.pdf",units="in",width=3.76, height=1.63)
