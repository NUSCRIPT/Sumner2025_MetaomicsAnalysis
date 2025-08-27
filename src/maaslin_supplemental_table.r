library(tidyverse)
mox_maaslin_list <- readRDS('objects/mox/MOX_01_maaslin_results.rds')

full_table_maaslin <- function(mox_maaslin_list){
    df_list <- list()
    cnt <- 0
    categories <- names(mox_maaslin_list)
    # Seq through maaslin differential categories (i = Culture,...)
    for (i in seq_along(mox_maaslin_list)){
        i_category <- categories[i]
        omics <- names(mox_maaslin_list[[i]])
        
        # Seq through omics results within the categories (j = AMP [ASV],...)
        for (j in seq_along(mox_maaslin_list[[i]])){
            cnt <- cnt + 1
            j_omic <- omics[j]
            results <- mox_maaslin_list[[i]][[j]]$results
            
            sig_results <<- results %>% 
                mutate(Model = i_category,
                       DataType = j_omic,
                       .before=1)

            df_list[[cnt]] <- sig_results
        }
    }
    
    df_list <-do.call("rbind", df_list)
    return(df_list)
}

all_maaslin_results <-full_table_maaslin(mox_maaslin_list = mox_maaslin_list)
all_maaslin_results <- all_maaslin_results %>% filter(qval<.1)
write_tsv(all_maaslin_results, "tables/DifferentialAbundanceResults.tsv")
