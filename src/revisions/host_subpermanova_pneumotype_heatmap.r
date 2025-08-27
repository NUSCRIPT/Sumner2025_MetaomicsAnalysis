source("src/common/get_distances.r")
source("src/common/load_metadata.r")

run_permanova <- function(a_dist, a_meta) { 
  # Instantiate vectors/count
  adonis_out <- c()
  adonis_out_df <- c()
  out_cnt <- 1
  
  # Iterate through desired metadata 
  for (i in permanova_metadata) {
    print(paste("Analyzing:",i))
    
    # Get pretty name from list
    pretty_i <- names(which(i == permanova_metadata)) 
    
    # Drop samples with missing metadata
    a_meta_subset <- a_meta %>% tidyr::drop_na(i)
    
    # Make sure more than one metadata component to prevent errors
    if ( length(unique(a_meta_subset[[i]])) > 1 ) {
      
      # Prevent errors, subset to samples with relevant metadata (no NA) and reorder
      common_labels <- intersect(labels(a_dist), a_meta_subset$sample_name)
      a_meta_subset <- dplyr::filter(a_meta_subset, sample_name %in% common_labels)
      a_dist_subset <- dist_subset(a_dist, a_meta_subset$sample_name)
      
      print(dim(a_meta_subset))
      
      # Execute
      form <- as.formula(paste("a_dist_subset", i, sep="~"))
      permanova <- vegan::adonis2(form, data = a_meta_subset, permutations = 4999, parallel = 10)
      print(permanova)
      
      
      # Save relevant components of each PERMANOVA
      adonis_out[[out_cnt]] <- permanova
      adonis_out_df[[out_cnt]] <- slice_head(as_tibble(permanova)) %>% mutate(var = pretty_i)
      
    }
    
    out_cnt <- out_cnt + 1
    
  }
  
  return(adonis_out_df)
}


cluster_subset <- permmeta %>% 
  dplyr::select(sample_name, cluster_num_amplicon) %>% 
  drop_na() %>% 
  mutate(cluster_num_amplicon = factor(cluster_num_amplicon, c(unique(cluster_num_amplicon), "not_cluster"))) %>% 
  pivot_wider(names_from = cluster_num_amplicon,
              values_from = cluster_num_amplicon,
              names_prefix = "is_cluster_",
              values_fill = "not_cluster")
permanova_metadata <- names(cluster_subset)[-1]
cluster_perm<- run_permanova(mox_dist_list[["RNA [Host Transcriptomics]"]],cluster_subset)

sub_adonis <- do.call("rbind", cluster_perm) %>%
  dplyr::rename(P = "Pr(>F)") %>%
  mutate(FDR = p.adjust(P, method = "fdr")) %>% 
  mutate(percent_variance = paste0(as.character(round(R2*100, 1)), "%"),
         FDR_annotation = case_when(FDR < 0.001 ~ "***",
                                    FDR < 0.01 ~  "**",
                                    FDR < 0.05 ~ "*",
                                    TRUE ~  ""))
sub_adonis$cluster <- permanova_metadata


alpha <- 1.9
beta <- 0.1
colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
sub_adonis$lblcolor <- ifelse(qbeta(sub_adonis$R2, alpha, beta) < 0.7, "black", "white")


pn_permanova <- ggplot(sub_adonis, aes(x="RNA [Host Transcriptomics]", y=cluster,label = FDR_annotation, fill = R2)) + 
  geom_tile(color = "black", width=1) + #.9
  geom_text(aes(label=FDR_annotation, color = lblcolor), size=9/.pt, nudge_y = .12) + 
  geom_text(aes(label=percent_variance, color = lblcolor), size=8/.pt, nudge_y = -.12) +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = TRUE,
                                title = "",
                                frame.linewidth=.4,
                                frame.colour="black",
                                ticks.colour = "black")) +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(position="right",expand=c(0,0)) + 
  scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                       values=colorvalues,
                       na.value="white", 
                       limits=c(0, 1), 
                       name="Variance Explained") +        
  scale_color_manual(values=c(white="white", black="black")) +  
  theme_nature() + 
  theme(axis.text.x=element_text(angle=335,hjust=0),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.5),
        axis.line.y.right = element_blank(), 
        axis.line.x.bottom = element_blank()
  ) +
  guides(color='none') +
  labs(x="", y="",fill="Variance Explained")
pn_permanova

ggsave("figures/host_deseq/sub_permanova.pdf", units="in", height=1.47, width=1.31)
saveRDS(pn_permanova, "objects/figure_objects/sub_permanova.rds")
