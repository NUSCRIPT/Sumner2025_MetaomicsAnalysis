library(tidyverse)
source("src/common/load_metadata.r")
source("src/common/load_features.r")
source("src/common/load_lists.r")
library(viridis)



                       
draw_pheatmap <- function(my_feature_name="DNA [Taxonomy]"){
    dt <- get_most_prevalent(my_feature_name, .01) 
    
    dt_features <- dt %>% pluck('feature')
    dt_samples <- dt[-1] %>% colnames() %>% sort()
    
    md <- permmeta %>% filter(sample_name %in% dt_samples) %>% arrange(sample_name)
    
    pseudocount <- dt %>% tidy_features() %>% filter(value > 0) %>% summarize(min = sin(min(value))^2) %>% pluck('min')/2
    
    dt <- dt %>% select(feature, md$sample_name) %>% arrange(feature) %>% column_to_rownames('feature')
    dt <- log10((sin(dt)^2+pseudocount)) # Undo AST and replace w Log10
    md_sub <- md  %>% select(sample_name, c("cluster_num_dna_taxonomy","cluster_num_amplicon", "cluster_num_dna_ko","Episode_is_cured", "Episode_etiology", "Episode_category", "Smoking_status","dysbiosis_score", "neutrophils_bodyfluid", "Quantity.Mean.Log", "amylase_bf_log")) %>% column_to_rownames('sample_name')
    
    my_pheat <- pheatmap::pheatmap(dt,
                       annotation_col = md_sub, 
                       show_colnames = F, 
                       show_rownames = T,
                       fontsize_row = 1, 
                       fontsize = 2,
                       angle_col = 90,
                       border_color = "black",
                       color = colorRampPalette(c("white", "black"))(200),
                       cluster_rows = T,
                       cluster_cols = T,
                       treeheight_col = 5,
                       treeheight_row = 0,
                       # clustering_distance_rows = 'correlation',
                       main = "Log10 Transformed Values Prevalence > 1%")
    return(my_pheat)
    }

draw_pheatmap("Amplicon")
draw_pheatmap("Amplicon ASV")
draw_pheatmap("DNA [Taxonomy]")
tmp<-draw_pheatmap("DNA [Viral]")
draw_pheatmap("RNA [Taxonomy]")            
draw_pheatmap("RNA [Taxonomy, mpa]")

draw_pheatmap("DNA [KO]")
draw_pheatmap("RNA [Host Transcriptomics]")




