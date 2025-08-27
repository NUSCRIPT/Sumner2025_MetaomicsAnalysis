require(ade4)
require(vegan)
library(usedist)
require(RColorBrewer)

set.seed(100)
run_mantel <- function(a_dist, b_dist) { 
  # Subset dist matrix 
  common_labels <- intersect(labels(a_dist), labels(b_dist))
  a_subset <- dist_subset(a_dist, common_labels)
  b_subset <- dist_subset(b_dist, common_labels)
  print(length(common_labels))
  mantel_results <- mantel.randtest(a_subset, b_subset, nrepet = 5000)
  return(mantel_results)
}

# run_mantel2 <- function(a_dist, b_dist) { 
#   # Subset dist matrix 
#   common_labels <<- intersect(labels(a_dist), labels(b_dist))
#   a_subset <<- dist_subset(a_dist, common_labels)
#   b_subset <<- dist_subset(b_dist, common_labels)
#   print(length(common_labels))
#   md <- permmeta %>% filter(sample_name %in% common_labels)
#   h1 <- with(md, how(nperm = 999,
#                      plots = Plots(strata = script_id),
#                      within = Within(type = "free")
#   )
#   ) #, nperm = 199, minperm = 200
#   h1 <- with(md, how(nperm = 999,
#                      blocks=initial_bal_sample,
#   )
#   ) #, nperm = 199, minperm = 200
#   
#   my_strata <- case_when(str_detect(common_labels, "BAL00") ~ 0, .default = 1)
#   mantel_results <- vegan::mantel(a_subset, b_subset,
#                                   method="spearman", 
#                                   permutations = 500
#                                   # strata = my_strata
#   )
#   return(mantel_results)
# }

run_mantel_list <- function(dist_list, dist_names) {
  
  # Instantiate df/cnt for holding mantel test results
  # set.seed(100)
  mantel_dfs <- c()
  cnt <- 1
  
  # Iterate through dist_list
  for (seq_a in seq_along(dist_list)) {
    for (seq_b in seq_along(dist_list)) {
      
      # Get names for plotting 
      a_dist_name <- names(dist_list)[seq_a]
      b_dist_name <- names(dist_list)[seq_b]
      print(paste('Comparing: ', a_dist_name, b_dist_name))
      # Only run mantel on (1) different matrixes (2) one time
      if (seq_a < seq_b) {
        
        
        # Execute mantel test between the two dist. matrixes
        mantel_results <- run_mantel(dist_list[[a_dist_name]], dist_list[[b_dist_name]])
        
        # Create df of mantel results and hold in endpoint df for plotting
        mantel_df <- as_data_frame(rbind(c(mantel_results$expvar, 
                                           mantel_results$pvalue, 
                                           mantel_results$obs, 
                                           a_dist_name, 
                                           b_dist_name)
        )
        )
        mantel_dfs[[cnt]] <- mantel_df
        
        # Update element index for next result
        cnt <- cnt + 1      
      } else if (seq_a == seq_b) {
        print(paste('Comparing: ', a_dist_name, b_dist_name))
        # Dont run mantel on self/self comparison
        # Create fake/na df of mantel so plots in triangle 
        mantel_df <- as_data_frame(rbind(c('Std.Obs'= NaN, 
                                           'Expectation' = NaN, 
                                           'Variance'=NaN,
                                           NaN, 
                                           NaN,
                                           a_dist_name, 
                                           b_dist_name)
        )
        )
        mantel_dfs[[cnt]] <- mantel_df 
        # Update element index for next result
        cnt <- cnt + 1      
      }
    }
  }
  
  return(mantel_dfs)
}

get_clean_mantel <- function(all_mantel_results){
  mantel_clean <- do.call("rbind", all_mantel_results) %>%
    dplyr::rename(P = "V4",
                  MantelStat = "V5",
                  a_dist_name = "V6",
                  b_dist_name = "V7") %>%
    filter(P!=NaN) %>% 
    #filter(a_dist_name != b_dist_name) %>% 
    mutate(FDR = p.adjust(P, method = "fdr"),
           MantelStat = as.double(MantelStat))  %>%
    mutate(
      #FDR = round(FDR, 3),
           R2 = MantelStat^2,
           percent_variance = case_when(is.nan(R2) ~ '', .default = paste0(as.character(round(R2*100, 2)), "%")),
           FDR_annotation = case_when(FDR < 0.001 ~ "***",
                                      FDR < 0.01 ~  "**",
                                      FDR < 0.05 ~ "*",
                                      TRUE ~  ""
           )) %>%
    mutate(a_dist_name = factor(a_dist_name, levels = names(mox_dist_list)),
           b_dist_name = factor(b_dist_name, levels = names(mox_dist_list))) 
  return(mantel_clean)
}

#%>% filter(a_dist_name != b_dist_name)
plot_mantel_heat <- function(mantel_clean){
  alpha <- 2.9
  beta <- 0.1
  colorvalues <- pbeta(seq(0, 1, length=101), alpha, beta)
  mantel_clean$lblcolor <- ifelse(qbeta(mantel_clean$R2, alpha, beta) < 0.8, "black", "white")
  
  mantel_heat <- ggplot(mantel_clean, aes(x=b_dist_name, y=a_dist_name, fill = R2, label=FDR_annotation)) +
    geom_tile(color = "black", size=.4) + #.9
    # geom_text(aes(color = lblcolor), size=6/.pt)+
    geom_text(aes(label=FDR_annotation, color = lblcolor), size=7/.pt, nudge_y = .12) + 
    geom_text(aes(label=percent_variance, color = lblcolor), size=6/.pt, nudge_y = -.12) +
    guides(fill = guide_colourbar(label = TRUE,
                                  ticks = TRUE,
                                  title = "",
                                  frame.linewidth=.5,
                                  frame.colour="black",
                                  ticks.colour = "black"),
           color = 'none') +
    scale_x_discrete(expand=c(0,0)) + 
    scale_y_discrete(position="right",expand=c(0,0)) + 
    scale_fill_gradientn(colours=colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
                         values=colorvalues,
                         na.value="white", 
                         limits=c(0, 1), 
                         name="Variance Explained") + 
    scale_color_manual(values=c(white="white", black="black")) +
    theme_nature() + 
    theme(legend.position = "right",
          axis.text.x=element_text(angle=335,hjust=0),
          #panel.border = element_rect(colour = "black", fill=NA,linewidth = .5)
    ) +
    xlab(NULL) +
    ylab(NULL)
  mantel_heat
  return(mantel_heat)
}

# Note, found out you can do spearman, which may be more applicable for noneuclid Bray and rnaVSdna comparisons in future studies
#vegan::mantel(mox_dist_list[["AMP [Genus]"]], mox_dist_list[["AMP [ASV]"]],method = "spearman")
