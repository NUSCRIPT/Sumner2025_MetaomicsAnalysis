require(decontam)

#grepl(paste(patterns, collapse="|"), Letter)
# PURPOSE: Use decontam to (1) identify and (2) remove contaminants from feature tables
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <md> = type tibble, samples (columns) by variables (rows). first column is sample_name and matches  feature_tables columns
# <md_column> = type string, corresponds to column name with concentration metadata OR with control metadata (T=negative, F=sample)
# <method> = type string, either "prevalence" or "frequency" to match isContaminant functional
# RETURNS decontaminated feature table. 
# NOTE also puts contaminant dataframe in global env. 
decontaminate_feature_tables <- function(feature_table, md, md_column, method="frequency", my_threshold=0.1){
  feature_mat <- feature_table %>% 
    column_to_rownames("feature") %>% 
    t(.)
  
  if (method == "frequency"){
    concentrations <- md %>% 
      filter(sample_name %in% colnames(feature_table)[-1]) %>%
      mutate(sample_name = factor(sample_name, colnames(feature_table)[-1])) %>% 
      arrange(sample_name) %>% 
      pluck(md_column)
    
    length(concentrations) == length(rownames(feature_mat)) # SANITY CHECK
    
    contaminants <<- isContaminant(seqtab = feature_mat, 
                                   method = method,
                                   conc = concentrations,
                                   threshold = my_threshold)
    
    
  } else if (method == "prevalence"){
    negatives <<- md %>% 
      filter(sample_name %in% colnames(feature_table)[-1]) %>%
      mutate(sample_name = factor(sample_name, colnames(feature_table)[-1])) %>% 
      arrange(sample_name) 
    
    nsn <<- negatives$sample_name
    fsn <<- colnames(feature_table)[-1]
    
    print(paste("SANITY CHECK, METADATA AND FEATURE COLUMNS ORDERED IDENTCALLY:",
                identical(as.character(nsn), fsn)))
    
    negatives <<- negatives %>% pluck(md_column)
    
    print(paste("SANITY CHECK, SAME LENGTH:",
                length(negatives) == length(rownames(feature_mat)))) # SANITY CHECK
    
    contaminants <<- isContaminant(seqtab = feature_mat, 
                                   method = method,
                                   neg = negatives,
                                   threshold = my_threshold)
    
  }
  
  
  
  contaminant_vector <<- contaminants %>%  filter(contaminant==T) %>% rownames()
  
  print(paste(length(contaminant_vector), "contaminants found and removed from dataset"))
  
  feature_table <- feature_table %>% filter(!feature %in% contaminant_vector)
  return(feature_table)
}

