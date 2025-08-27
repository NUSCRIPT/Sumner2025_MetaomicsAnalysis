require(tidyverse)

# helper function for putting feature tables into tidyformat
tidy_features <- function(my_feature_table){
    my_feature_table <- my_feature_table %>% pivot_longer(2:last_col(), names_to='sample_name')
    return(my_feature_table)
}

# <file_path> = string, path to humann output
# <type> = string, "ko" for kegg. no functionality yet
# <normalize> = BOOL, defaults TRUE - for TSS scaling
# <transform> = string, default is none, other options include "ast" or "log"
# <genes_only> = BOOL, defaults TRUE - removeing unclassified & unmapped AFTER norm
# NOTE: can do log w or without TSS normalizaiton. with TSS, will be mult by 100. Both cases give pseudolog of 1
read_humann <- function(file_path, type="ko", normalize=TRUE, transform="none", genes_only=TRUE){
    
    # Initiate DF via readin
    d <- read_tsv(file_path)
    
    # 1: CANALIZE DATA STRUCTURE / SYNTAX
    d = d %>% dplyr::rename(feature=1) 
    names(d) = gsub(names(d), pattern = "_Abundance-CPM", replacement = "") 
    names(d) = gsub(names(d), pattern = "_Abundance-RELAB", replacement = "") 
    names(d) <- gsub(names(d), pattern = "PT_", replacement = "") 
    
    d <- d %>% tidy_features() 
    
    if (normalize == TRUE){
        d <- d %>% 
            group_by(sample_name) %>% 
            mutate(value = value/sum(value)) %>% 
            ungroup()
        
        if (transform == "ast") {
            d <- d %>% 
                mutate(value = asin(sqrt(value)))
            
        } else if (transform == "log") {
            d <- d %>% 
                mutate(value = log10(value*100+1))
            
        }
    } else {
        if (transform == "log") {
            d <- d %>% 
                mutate(value = log10(value+1))
        }
    }
    
    if (genes_only == TRUE){
        d <- d %>% 
            filter(!feature %in% c("UNGROUPED", "UNMAPPED"))
    }
    
    d <- d %>% pivot_wider(names_from = "sample_name", values_from = "value")
        
    return(d)
}
#a_sample_sums <- function(x) {colSums(Filter(is.numeric, x)) %>% enframe() %>% arrange(., value)}
# Some examples:
# read_humann(FILE_mgx_ko, normalize = T, transform = "ast", genes_only = T)
# read_humann(FILE_mgx_ko, normalize = T, transform = "log", genes_only = T)
# read_humann(FILE_mgx_ko, normalize = F, transform = "log", genes_only = T)
# read_humann(FILE_mgx_ko, normalize = T, transform = "none", genes_only = T)
# read_humann(FILE_mgx_ko, normalize = F, transform = "none", genes_only = T)


# For very basic read-in function if you wwant to use feature_table_functions.r
read_humann_basic <- function(file_path){
  
  # Initiate DF via readin
  d <- read_tsv(file_path)
  
  # 1: CANALIZE DATA STRUCTURE / SYNTAX
  d = d %>% dplyr::rename(feature=1) 
  names(d) = gsub(names(d), pattern = "_Abundance-CPM", replacement = "") 
  names(d) = gsub(names(d), pattern = "_Abundance-RELAB", replacement = "") 
  names(d) <- gsub(names(d), pattern = "PT_", replacement = "") 
  return(d)
}
  



