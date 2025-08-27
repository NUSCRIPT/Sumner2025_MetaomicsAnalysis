# FOR PROCESSING BRACKEN FILES. NOTABLY (1) READING THEM IN 

# Takes a bracken file and processes it to meet feature_table format specifications
# <bracken_file> =  string instance with /path/to/metaphlan_all.tsv
# <sample_suffix> = type string, matching ends of column names to be removed
# <counts> = type BOOL, select proportional data (FALSE , default) or reads counts (TRUE)
# returns 
#   tibble object with [ rows = species, columns = samples ]. 

read_bracken <- function(bracken_file, counts=FALSE, sample_suffix=".bracken.tsv_frac"){
  if (counts){
    quantification_suffix = "_frac"
  } else {
    quantification_suffix = "_num"
  }
  
  feature_table <- read_tsv(bracken_file) %>% 
    rename(feature=1) %>% 
    select(!c(taxonomy_id, taxonomy_lvl)) %>%
    select(!ends_with(quantification_suffix))
  
  names(feature_table) <- gsub(names(feature_table), pattern = sample_suffix, replacement = "") 
  names(feature_table) <- gsub(names(feature_table), pattern = "PT_", replacement = "")
  return(feature_table)
}

