require(tidyverse)
require(vegan)
# Function designed to (1) perform ordination analysis and (2) visualize PCoAs 
my_ordination_helper <- function(my_distance, my_meta=permmeta){
  set.seed(200)
  my_pcoa <- cmdscale(my_distance, k = 4, eig = TRUE, add = FALSE)

  # cleanup
  my_ord <- as.data.frame(my_pcoa$points)
  names(my_ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates
  
  # Percent explained variation
  my_eig <- eigenvals(my_pcoa)
  eig.percent <- 100*head(my_eig/sum(my_eig))
  print(eig.percent)
  eig.1 <- paste("PCo1 (", as.character(round(eig.percent[1], digits = 1)), "%)", sep = "")
  eig.2 <-paste("PCo2 (", as.character(round(eig.percent[2], digits = 1)), "%)", sep = "")
  
  ## plot PCoA (FIGURE 1A)
  my_meta_ord<- my_ord %>% 
    rownames_to_column("sample_name") %>% 
    left_join(my_meta) %>% mutate(sham = '21') # CHANGED PERMMETA
  
  ordination_results <- list('ord_obj' = my_meta_ord,
                             'label_eig.1' = eig.1,
                             'label_eig.2' = eig.2)
  return(ordination_results)
}
##################################################33333
visualize_ordination <- function(ordination_results, fill_value, shape_value, 
                                 size_value=2, stroke_value=.4, add_background_points=F, my_alpha=1){
  
  # Make arbitrary entry for shape_value if not present in order to have 'no shape' option
  if (is.null(ordination_results[[1]][[shape_value]])) {
    ordination_results[[1]][[shape_value]] <- 'Map'
  }
  
  # Plot ordination
  ggp <- ordination_results[[1]] %>%
    arrange(!is.na(.data[[fill_value]]), .data[[fill_value]]) %>% # arrange plotting order, sort first by not na to put na on bottom layer
    ggplot(.) +
    {if(add_background_points)geom_point(data = ordination_results[[1]][1:4], aes(x=pcoa1, y=pcoa2), color='grey')}+
    geom_point(aes(x = pcoa1, 
                   y = pcoa2, 
                   fill = .data[[fill_value]], 
                   shape = .data[[shape_value]]
    ), 
    alpha = my_alpha,
    size = size_value, 
    stroke = stroke_value) +
    xlab(ordination_results[[2]]) +
    ylab(ordination_results[[3]]) +
    theme_nature() +
    scale_shape_manual(values = c(21,22,23,24,25)) +
    theme(aspect.ratio = 1,
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) 
  
  # map legend to continuous/categorical fill values
  if (is.numeric(ordination_results[[1]][[fill_value]])) {
    ggp <- ggp + guides(fill = guide_colorbar()) + scale_fill_viridis()
  } else {
    ggp <- ggp + guides(fill = guide_legend(override.aes=list(shape=21))) 
    
  }
  
  return(ggp)
  
}
################################################################################33