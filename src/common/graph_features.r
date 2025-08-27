library(tidyverse)
library(ggplot2)
library(ggridges)
library(colorspace)
source('src/common/clean_feature_names.r')
# Major graphing helper for visualizing differentially expressed features
#   use this for most of fig2 to create density plots with zero counts in flavor 
#   of Lloyd-Price et al 2019 nature
# <feature_col_name> = string of feature columns name
# <group_col_name> = string of groups column name
# <baseline_group> = string of base comparison group (e.g, NPC)
# <recode_dash> = string to change feature's with - to . to match array (maaslin2)
# <make_names> = bool, to use make.names in case feature_array is from maaslin
graph_features <- function(features_omic, 
                           feature_array, 
                           groups_df, 
                           feature_col_name, 
                           group_col_name, 
                           baseline_group, 
                           recode_dash = ".", 
                           scale_val=4, 
                           logspace=T,
                           specify_colors=TRUE,
                           zscalespace = T,
                           clean_feature=T,
                           make_names = F) {
  
  print(feature_array)
  
  my_groups <- groups_df %>% 
    dplyr::rename(meta_group = group_col_name) %>%
    dplyr::select(sample_name, meta_group) %>%
    mutate(meta_group = factor(meta_group))
  
  # meta_group = fct_rev(meta_group))
  ng <- length(levels(my_groups$meta_group))
  nf <- length(feature_array)
  
  
  # Replace feature column name
  my_features <- features_omic %>% 
    dplyr::rename(feature = feature_col_name)
  
  # Rename features with make.names depending on parameters
  # Note: make.names is needed if feature array is supplied from maaslin2
  if (make_names == T){
    my_features <- my_features %>% 
      mutate(feature = make.names(feature))
  }
  
  # Subset features to feature_array- replace dash (now not needed but kept for dependencies)
  my_features <- my_features %>% 
    mutate(feature = str_replace_all(feature, "-|\\[|\\]", recode_dash)) %>% 
    filter(feature %in% feature_array)
  
  if (clean_feature == T){
    my_features <- my_features %>%  mutate(feature = fix_features(feature))
    feature_array = fix_features(feature_array) # formatted feature names
    
  }
  
  # feature_array = fix_features(feature_array) # formatted feature names
  # print(feature_array)
  
  # Arrange features table into tidy format 
  my_features <- my_features %>% 
    pivot_longer(2:last_col()) %>% 
    left_join(my_groups, by=c("name"="sample_name")) %>%
    mutate(feature = factor(feature, levels=feature_array)) %>% # level features
    drop_na(meta_group)
  
  print(feature_array)
  print(my_features)
  # # Drop zeros
  my_features_nz <- my_features %>% 
    filter(value > 0) 
  # og removed zeros, see green noteboook for thoughts (pg 26)
  
  
  # Get median of baseline group for transformation
  grouped_medians <- my_features_nz %>%
    filter(meta_group==baseline_group) %>%
    group_by(feature
             # meta_group
    ) %>%
    summarize(med = median(value),
              sdev = sd(value)) %>%
    # filter(meta_group==baseline_group) %>%
    ungroup()
  
  # Get median of overall group for in case baseline is 0 prevalence
  grouped_medians_all <- my_features_nz %>%
    group_by(feature
             # meta_group
    ) %>%
    summarize(med_all = median(value),
              sdev_all = sd(value)) %>%
    # filter(meta_group==baseline_group) %>%
    ungroup()
  
  grouped_medians <- grouped_medians_all %>% 
    left_join(grouped_medians) %>% 
    mutate(med = case_when(is.na(med) ~ med_all, .default = med),
           sdev = case_when(is.na(sdev) ~ sdev_all, .default = sdev)) %>% 
    select(!c(med_all, sdev_all))
  # grouped_medians <- grouped_medians %>% select(!meta_group)
  
  if (logspace){
    my_features_nz <- my_features_nz %>% 
      left_join(grouped_medians) %>%
      mutate(value = log2(value/med)) %>%  # og log10(value/med), value-med is centered
      ungroup()
  }
  
  
  # get min of values for calc limits and zero boxes
  my_min = min(my_features_nz$value)*1.5 #-(ng*.4+1)#*1.3 #- abs(min(my_features_nz$value))*.5
  my_max = max(my_features_nz$value)*1.5 # Expand 30%
  d = my_max - my_min
  d_bar = d *.0375
  d_bar_tot = d_bar * ng
  # Make zero data for zero rectangles, height weighted by fraction zero at given feature
  zeros <- my_features %>% 
    group_by(feature, meta_group) %>%
    summarize(num_zero = sum(value<=0)/n()) %>%
    mutate(x="zero") %>%
    ungroup() %>%
    mutate(xmin=(my_min-d_bar_tot)+(as.numeric(meta_group)-1)*d_bar,# .3,  #.1
           xmax=xmin+d_bar*.95,
           ymin=as.numeric(meta_group),
           ymax=(num_zero*(ng+scale_val)+ymin)) %>% 
    # PArenthesees here make more uniform / proportional and easier interpreetation
    
    arrange(feature,desc(meta_group))
  
  my_features_nz <- my_features_nz %>% 
    left_join(zeros) %>%
    mutate(num_nz = 1-num_zero) %>%
    arrange(feature,desc(meta_group))
  
  # Prettify feature names
  pretty_features <- gsub("_ec_[0-9].*","",feature_array) %>% 
    str_replace_all(., "_", " ") %>% 
    str_to_title(.) %>% 
    str_wrap(., whitespace_only = T, width=30)
  names(pretty_features) <- feature_array
  
  print(ng)
  print(nf)
  
  # Spacing Estimators
  psy <- (-1/(nf*4))*1.8
  print(zeros, n=Inf)
  print(my_features_nz)
  
  
  #plot
  ggp <- ggplot(my_features_nz , group=meta_group
  ) + 
    geom_vline(xintercept=0, color='grey', linewidth=.3) +
    
    geom_segment(data = zeros,
                 aes(x = xmin,
                     y = meta_group,
                     xend = my_max,
                     yend = meta_group,
                     color = meta_group
                 ),
                 # color='black',
                 # alpha=.8,
                 linewidth=.35) +
    geom_rect(data=zeros,
              mapping=aes(xmin=xmin,
                          xmax=xmax,
                          ymin=ymin,
                          ymax=ymax,
                          fill=meta_group,
                          color=meta_group,
              ),
              linewidth=.3,
              position="identity",
    ) +
    geom_density_ridges(aes(x=value, 
                            y = meta_group, 
                            fill=meta_group, 
                            color=meta_group,
                            group=meta_group,
                            scale=scale_val*num_nz,
                            rel_min_height=0.001 # If too low, can lead to overlap beyond min X
    ),
    stat='density_ridges',
    linewidth=.3,
    # color = 'black'
    ) +
    
    geom_vline(aes(xintercept=unique(max(zeros$xmax))+.05), 
               linewidth=.3
    )+
    labs(color=group_col_name, 
         fill=group_col_name
    )+
    geom_text(aes(x=my_max, y=0, label=feature),
              nudge_x = 0, nudge_y = 0, #-1
              hjust=1, #1 and my_max for rigth sided
              vjust=.75, #0
              check_overlap = T,
              size = 4/.pt,
              parse = T
    )+
    facet_wrap(~feature, 
               ncol=1, 
               shrink = F, 
               strip.position = "right", 
               # labeller = as_labeller(pretty_features),
               drop = F
    ) +
    theme_nature() +
    theme(strip.background=element_blank(), # Witchcraft but math 
          strip.text=element_text(margin=c(0,0,0,0)),
          strip.text.y.right = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y =element_blank(),
          legend.position = 'top',
          legend.box = "horizontal",
          legend.margin      = margin(0,0,0,0),
          legend.key.size    = unit(5, "pt"),
          legend.box.spacing = unit(1, "pt"),
          panel.spacing = unit(0, "npc"),
          panel.spacing.y = unit(#psy,
            -1/((ng+scale_val)*nf),
            "native") # npc makes non-realtive to save size
    ) +
    scale_y_discrete(expand = expansion(mult=c(0,0),  add = c(0,0))) + # if add > 0, very slightly off on zero reaching to 1
    scale_x_continuous(expand=c(0, 0), limits=c(my_min-d_bar_tot,my_max))+
    coord_cartesian(clip = "off", ylim = c(-1, scale_val+ng), expand=F) +
    guides(fill = guide_legend(title.position = "top"))
  
  print(((ng+scale_val)*nf))
  
  # Add colors if special metadata 
  if (specify_colors){
    ggp <- ggp +
      scale_fill_manual(values=mox_color_lists[[group_col_name]]) +
      scale_colour_manual(values=unlist(lapply(mox_color_lists[[group_col_name]], darken,amount = .6))) 
    
  } else {
    ggp <- ggp +
      scale_colour_hue(l=60)
  }
  return(ggp)
}
#Example: graph_features(features_amp, da_amp, groups_df = permmeta_amp, feature_col_name = "Genus", group_col_name = "success", baseline_group = "2")
#Example: graph_features(features_amp, da_amp, groups_df = permmeta_amp, feature_col_name = "Genus", group_col_name = "pt_category", baseline_group = "NPC")


# Additional Plotting functions for more standard and exploratory plotting of 
#   Distributions and correlation plots, especially those in 09_dysbiosis
# <my_data> = data frame containing columns, x_var and y_var in it 
# <x/y_var> = type string; is a column name in my_data; x must be quantitative, y is category
# <x_lab> = type string or NULL; pretty string name for axis labels or defaults to col name
# <y_lab> = type string or NULL; see above
# <ridge_plot> = type bool; determines if plot in style of ridgeline or density=
plot_distribution <- function(my_data, x_var, y_var=NULL, x_lab=NULL, y_lab=NULL, ridge_plot = T){
  if ( is.null(x_lab) ) {
    x_lab = x_var
  }
  
  if ( ridge_plot ){
    if ( is.null(y_lab) ) {
      y_lab = y_var
    }
    ggp <- ggplot(my_data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[y_var]])) +
      geom_density_ridges(scale = 2, size=.4, alpha=.9, jittered_points = TRUE) +
      labs(x=x_lab, y=y_lab) +
      scale_y_discrete(expand = c(.1, 0)) 
  } else {
    
    ggp <- ggplot(my_data, aes(x = .data[[x_var]], fill = .data[[y_var]])) +
      geom_density(size=.5, alpha=.3) +
      labs(x=x_lab, y="Density")
  }
  
  ggp <- ggp + 
    theme_nature() + 
    theme(aspect.ratio = .6)
  return(ggp)
}
# Example: plot_distribution(permmeta, "dysbiosis_score", "pt_category", "Dysbiosis score", "Disease state", ridge_plot=F)


# Helper function to plot correlations in 09 dysbiosis and other helpers
# <my_data> = data frame containing columns, x_var and y_var in it 
# <x_var> = type string; is column name in my_data; x axis data
# <y_var> = type string; is column name in my_data; y axis data
# <x_lab> = type string or NULL; pretty string name for axis labels or defaults to col name
# <y_lab> = type string or NULL; see abova
plot_correlations <- function(my_data, x_var, y_var=NULL, x_lab=NULL, y_lab=NULL){
  if ( is.null(x_lab) ) {
    x_lab = x_var
  }
  
  if ( is.null(y_lab) ) {
    y_lab = y_var
  }
  
  ggp <- ggplot(my_data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[y_var]])) +
    geom_point() +
    geom_smooth(size=.5, method = "lm") +
    theme_nature() +
    labs(x=x_lab, y=y_lab) +
    theme(aspect.ratio = .6)
  
  return(ggp)
}
# Example: plot_correlations(permmeta, "dysbiosis_score", "Quantity.Mean.Log", "Dysbiosis score")
