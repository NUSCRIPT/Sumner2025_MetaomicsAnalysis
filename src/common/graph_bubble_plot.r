library(tidyverse)
library(ghibli)
source('src/common/load_lists.r')
source('src/common/load_metadata.r')
source('src/common/get_clusters.r')

# GOAL: Function to plot relative abundances of feature data 
# For features: capable of plotting all features, top N, or select features of interest
# For samples: capable of plotting all samples or will summarize to metadata category
# Options:
#   my_feature = feature data frame (feature is leftmost column)
#   my_distance = distance matrix for determining order of X axis
#   foi = vector of strings, when want to subset to specific features
#   plot_all=TRUE
#   filter_top_n=NULL # orNULL or int
#   category_level_plot='cluster_num_amplicon'
# Example:
bubble_plot <- function(my_distance, my_features, plot_all=TRUE, filter_top_n=NULL, category_level_plot=NULL, foi){
  
  library(correlation)
  # Tidy features
  my_features <- my_features %>% pivot_longer(2:last_col()) %>% 
    rename(sample_name='name')
  
  # Use PCoA to find meaningful bubbleplot order
  set.seed(200)
  my_pcoa <- cmdscale(my_distance, k = 4, eig = TRUE, add = FALSE)
  my_ord <- as.data.frame(my_pcoa$points)
  
  # Get spearman corr w pcoa components to order features - Chose axis to align to (V1, V2)
  my_corr <<- my_ord %>% rownames_to_column('sample_name') %>% 
    left_join(my_features) %>% 
    group_by(feature) %>% 
    correlation(method = "spearman")
  
  feature_order <- my_corr %>% filter(Parameter1 %in% c('V1', 'V2'), Parameter2=='value') %>% arrange(abs(rho)) %>% pluck('Group') #%>% arrange(abs(rho))
  
  my_features <- my_features %>%  mutate(feature = factor(feature, levels=unique(feature_order)))
  
  # Subset pcoa to means of groups if desired for arranging groups
  if(is.character(category_level_plot)){
    my_ord <- my_ord %>% rownames_to_column('sample_name') %>% 
      left_join(permmeta) %>% 
      drop_na(category_level_plot) %>% 
      group_by(.data[[category_level_plot]]) %>% 
      summarize(V1 = mean(V1),
                V2 = mean(V2)) %>% 
      column_to_rownames(category_level_plot)
  }
  
  bubble_order <- my_ord %>% arrange(V1, V2) %>% row.names()
  
  
  # Subset to Features Of Interest OR plot all
  if(plot_all==TRUE){
    foi <- my_features$feature
  }else{
    my_features <- my_features %>% filter(feature %in% foi)
  }
  
  my_features <- my_features %>% 
    left_join(permmeta)
  
  # Prevalence
  my_features <- my_features %>% 
    group_by(feature) %>% 
    mutate(prevalence = sum(value>0)/n(),
           median_abund = median(value)) %>% 
    filter(prevalence > .01) %>% #FILTER PREVALENCE
    ungroup() #%>% 
  # filter(prevalence > 0.1) #FILTER PREVALENCE
  
  if(is.character(category_level_plot)){
    my_features <- my_features %>% 
      group_by(.data[[category_level_plot]], feature) %>% 
      summarize(
        prevalence = sum(value>0)/n(),#prevalence,
        median_abund = median(value),
        value = median(value)) %>% #median_abund) %>% 
      # unique() %>% # reframe duplicates on prev/abund
      mutate(sample_name = cluster_num_amplicon)
  }
  
  print(head(my_features))
  print(dim(my_features))
  # Order cols&rows
  my_features <- my_features %>% 
    mutate(sample_name = factor(sample_name, levels=bubble_order)) %>% 
    arrange(desc(prevalence), desc(median_abund))
  # my_features <- my_features %>%  mutate(feature = factor(feature, levels=unique(feature_order)))
  
  
  if(is.double(filter_top_n)){
    my_features <- my_features %>% 
      filter(feature %in% unique(my_features$feature)[1:filter_top_n])
  }
  
  my_features <- my_features %>% 
    arrange(match(sample_name, bubble_order), feature) 
  
  # Dont plot below limit of detection
  my_features <- my_features %>% 
    filter(value>0)
  # filter(median_abund >0)
  
  # Dont plot NA sample_names
  my_features <- my_features %>% 
    drop_na(sample_name) 
  print(dput(unique(as.character(my_features$feature))[1:30]))
  min_size <-  min(my_features$value) 
  max_size <-  max(my_features$value)*1.05
  # Bubble plot
  ggp <- ggplot(my_features, aes(x=sample_name, y=feature)) +
    geom_point(aes(size = value, fill=as.factor(cluster_num_amplicon)), color='black', alpha = .8, shape=21) + #'grey38'
    theme_nature() +
    scale_size_continuous(limits = c(min_size, max_size), range = c(1,5))+
    theme(
      panel.grid.major.y = element_line(colour="grey80", linetype="dashed",linewidth = .2),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = .2),
      axis.line.x = element_blank(),
      axis.line.y = element_blank()) +
    rotate_x_text(90) +
    guides(fill=guide_legend(title="Pneumotype"), size=guide_legend(title="Abundance"))
  return(ggp)
}
