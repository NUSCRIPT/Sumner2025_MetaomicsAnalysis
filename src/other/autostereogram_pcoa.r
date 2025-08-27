source("src/common/load_features.r")
source("src/common/load_metadata.r")
source("src/common/load_lists.r")
source("src/common/ordination_functions.r")
source("src/common/get_distances.r")
source("src/common/drop_samples.r")
drop_samples()
library(tidyverse)
library(ggforce)
library(viridis)

visualize_ordination_stereo <- function(ordination_results, fill_value, shape_value, 
                                 size_value=2, stroke_value=.4, add_background_points=F, my_alpha=1){
    
    # Make arbitrary entry for shape_value if not present in order to have 'no shape' option
    if (is.null(ordination_results[[1]][[shape_value]])) {
        ordination_results[[1]][[shape_value]] <- 'Map'
    }
    
    # Plot ordination
    ggp <- ordination_results[[1]] %>%
        # drop_na(.data[[fill_value]]) %>%  # DELETE THIS LATER JUST FOR HELPING MAKE PLOTS IN PATCHWORK
        arrange(.data[["pcoa3"]], !is.na(.data[[fill_value]]), .data[[fill_value]]) %>% # arrange plotting order, sort first by not na to put na on bottom layer
        arrange(.data[["pcoa3"]]) %>% # arrange plotting order, sort first by not na to put na on bottom layer
        
        ggplot(.,aes(x = pcoa1, 
                     y = pcoa2, 
                     depth=pcoa3)) +
        {if(add_background_points)geom_point(data = ordination_results[[1]][1:4], aes(x=pcoa1, y=pcoa2), color='grey')}+
        geom_point(aes(
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
    ggp2<<-ggp
    # map legend to continuous/categorical fill values
    if (is.numeric(ordination_results[[1]][[fill_value]])) {
        ggp <- ggp + guides(fill = guide_colorbar()) + scale_fill_viridis()
    } else {
        ggp <- ggp + guides(fill = guide_legend(override.aes=list(shape=21))) 
        
    }
    
    return(ggp)
    
}
################################################################################33
# Amplicon
ordination_results_amp <- my_ordination_helper(mox_dist_list[["RNA [Taxonomy BKN]"]], my_meta = permmeta) 
# ordination_results_amp <- my_ordination_helper(mox_dist_list[['Amplicon ASV']]) 
ordination_results_amp$ord_obj<-ordination_results_amp$ord_obj %>% filter(str_detect(sample_name,"REP|CON",negate = T))

visualize_ordination_stereo(ordination_results_amp, 'pcoa3', 'sham', 
                            5,
                            stroke_value = .8, 
                     add_background_points = F,my_alpha = .9) +
    theme(
        plot.margin = margin(1,1,1.5,1.2, "mm"),
          panel.spacing = unit(1, "lines")
          # legend.position="none"
        ) +
    labs(x="", y="") +
    # scale_fill_manual(values=mox_color_lists[['cluster_num_amplicon']]) +
    facet_stereo(IPD = -63.5) +scale_depth(range=c(-0.3,.3))
    

    
    
    ggplot(ordination_results_amp$ord_obj, aes(pcoa1, pcoa2, color = cluster_num_amplicon, depth=pcoa3)) +
        geom_point(alpha=.8, size=3) + 
        theme_bw() +
        scale_color_brewer(palette='Dark2') +
        guides(colour=guide_legend(override.aes=list(alpha=1))) + 
        facet_stereo(IPD=-63.5) +
        theme(legend.position='bottom') 
    xlim(pc1.r) + ylim(pc2.r)
    facet_stereo(IPD = 63.5, panel.size=100,shrink=T) 
    
    