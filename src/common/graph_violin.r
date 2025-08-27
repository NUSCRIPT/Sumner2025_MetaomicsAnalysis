source('src/common/nature_theme.r')
library(ghibli)
# GOAL: Wrapper for plotting violin plots in consistent manner:
# EXAMPLE: ggplot(my_data, aes(x=my_categorical_column, y=my_quantitative_data, fill=my_categorical_column)) %>% plot_violin(.)
plot_violin <- function(ggp){
    # ggplot(., aes(x=cluster_num_amplicon, y=Quantity.Mean.Log, fill=cluster_num_amplicon)) + #x=outcome
        ggp <- ggp + 
            geom_violin(trim=TRUE, linewidth=.3)+
            geom_boxplot(width=0.1, fill="white", outlier.size = .1,linewidth=.2)+
            stat_summary(fun.data = n_fun, geom = "text", size = 3/.pt,vjust=-.5) + # to add n values
            theme_nature() +
            theme(
                # axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                # panel.grid.major.y = element_line(colour="grey", linetype="dashed"),
                # plot.margin = margin(5, 5, 5, 5),
                legend.position = 'none'
            )
        return(ggp)
}

# GOAL: Add statistical values to basic violin plots in a consistent manner
# EXAMPLE: ggplot(my_data, aes(x=my_categorical_column, y=my_quantitative_data, fill=my_categorical_column)) %>% plot_violin(.) %>%  plot_statistics(.)

plot_statistics <- function(ggp, stat_method='wilcox_test'){
    library(ggpubr)
    ggp <- ggp + 
        geom_pwc(bracket.nudge.y = .15, #.5
                 size = 0.2,
                 label.size = 1,
                 step.increase = .1,#0.12,
                 vjust = .5,#.5
                 tip.length = 0,#0.02
                 label = "p.adj.signif",
                 method=stat_method,
                 p.adjust.method = "fdr",
                 hide.ns = T,
                 symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, Inf), 
                                    symbols = c("***", "**", "*", "ns"))
                 )
    
    return(ggp)
}

n_fun <- function(x){
    return(data.frame(y = max(x), label = paste0("",length(x))))
}
#1.05*max(x)
# # EXAMPLE: 
# # 1. Instantiate plot
# my_plot <- new_metadata %>% 
#     drop_na(Quantity.Mean.Log) %>% 
#     ggplot(., aes(x=cluster_num_amplicon, y=Quantity.Mean.Log, fill=cluster_num_amplicon))  #x=outcome
#
# # 2. Execute in violin and plot stats...make pretty
# plot_violin(my_plot) %>% plot_statistics(.) +  
#     scale_fill_ghibli_d("LaputaMedium", direction = -1) +
#     guides(fill = guide_legend(title = "Pneumotype")) +
#     ylab('Log 16S copies/µL') 
    
