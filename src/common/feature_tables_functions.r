require(tidyverse)

# PURPOSE: helper function for putting feature tables into tidyformat
# <feature_table> =  type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# RETURNS feature table as tidy 3-column format (feature, sample_name, value)
tidy_features <- function(feature_table){
  feature_table <- feature_table %>% 
    pivot_longer(2:last_col(), names_to='sample_name')
  return(feature_table)
}

# PURPOSE: helper function for reversing tidy_features 
# <feature_table> =  type tibble, 3-column tidy-formatted feature table as made by tidy_features
# RETURNS feature table as tidy 3-column format (feature, sample_name, value)
wide_features <-  function(feature_table){
  feature_table <- feature_table %>% 
    pivot_wider(names_from = "sample_name", values_from = "value")
  return(feature_table)
}

# PURPOSE: wrapper to easily and reproducibly filter tables by prevalence at a given abundance
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <min_prevalence> = type double, minimum prevalence in range 0-1 
# <min_abundance> = type double, minimum abundance threshold for an instance to count toward prevalence 
# returns filtered table
# instantiates <top_features>, type character vector, into global env. Is list of top prevalent features
get_most_prevalent <- function(feature_table, min_prevalence=0.1, min_abundance=0){
  # get prev features
  top_features <<- feature_table %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarize(prevalence = sum(value > min_abundance)/n()) %>% 
    arrange(desc(prevalence)) %>%
    filter(prevalence>=min_prevalence) %>%
    pluck('feature')
  
  print(paste(length(top_features), "features passing prevalence filtering"))
  # filter table
  my_features_table <- feature_table %>% 
    filter(feature %in% top_features)
  
  return(my_features_table)
}

# PURPOSE: wrapper to get top N most abundant features by mean in a feature table 
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <min_prevalence> = type double, minimum prevalence in range 0-1 
# returns filtered table
# instantiates <top_features>, type character vector, into global env. Is list of top prevalent features
get_most_abundant_n <- function(feature_table, n_features=20){
  # get prev features
  top_features <<- feature_table %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarize(abund = mean(value)) %>% 
    arrange(desc(abund)) %>%
    slice_max(abund, n = n_features) %>%
    pluck('feature')
  
  print(paste(length(top_features), "features passing abundance filtering"))
  # filter table
  my_features_table <- feature_table %>% 
    filter(feature %in% top_features)
  
  return(my_features_table)
}

# PURPOSE: wrapper function for easily & reproducibly performing TSS normalization with optional AST or LOG10P transformation.
# <d> = type tibble (feature_table), feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <normalize> = BOOL, defaults TRUE - for TSS scaling
# <transform> = string, default is none, other options include "ast" or "log"
# NOTE: can do log w or without TSS normalizaiton. with TSS, will be mult by 100. Both cases give pseudolog of 1
normalize_feature_table <- function(d, normalize=TRUE, transform="none"){
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
        mutate(value = log2(value*100+1))
      
    }
  } else {
    if (transform == "log") {
      d <- d %>% 
        mutate(value = log2(value+1))
      
    } else if (transform == "log100") {
      d <- d %>% 
        mutate(value = log2(value*100+1))
    } else if (transform == "ast_no_tss") {
      d <- d %>% 
        mutate(value = asin(sqrt(value)))
    } 
  }
  
  d <- d %>% pivot_wider(names_from = "sample_name", values_from = "value")
  return(d)
}

# PURPOSE: wrapper to drop unclassified feature
# <d> = type tibble (feature_table), feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <unclassified_terms> = type character vector, array of feature terms that represent unclassified, unamapped values in feature table
# RETURNS feature table without features defined in unclassified_terms vector
# NOTE: default vector includes terms for humann (UNGROUPED, UNMAPPED) and metaphlan (unclassified)
drop_unclassified <- function(d, unclassified_terms = c("UNGROUPED", "UNMAPPED", "UNCLASSIFIED", "UNINTEGRATED")){
  d <- d %>% 
    filter(!feature %in% unclassified_terms)
  
  return(d)
  
}

# PURPOSE: wrapper to drop selected feature, similar to drop_unclassified but for more general selection
# <d> = type tibble (feature_table), feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <feature> = type character vector, array of feature terms that represent unclassified, unamapped values in feature table
# RETURNS feature table without features defined in feature_terms vector
# NOTE: default vector includes terms for humann (UNGROUPED, UNMAPPED) and metaphlan (unclassified)
drop_features <- function(d, feature_terms = c("")){
  d <- d %>% 
    filter(!feature %in% feature_terms)
  
  return(d)
}
# PURPOSE: select sample columns from a feature table using a regex string
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <pattern> = type string, regex expression used to subset sample columns. e.g., "BAL" to get only samples with BAL in them
# <negate> = type BOOL, whether to apply positive (default) or negative selection
# RETURN feature table with selected samples only
get_sample_subset <- function(feature_table, pattern="BAL", negate=FALSE){
  if (!negate){
    feature_table <- feature_table %>% 
      dplyr::select("feature", matches(pattern))
    
  } else if (negate){
    feature_table <- feature_table %>% 
      dplyr::select("feature", !matches(pattern))
    
  }
  
  return(feature_table)
}

# PURPOSE: wrapper for getting distance metric from feature table. works by transforming feature table to appropriate transposition before plugging in.
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# <method> = type string, for chosing distance metric taken by vegan
get_distance <- function(feature_table, method = "bray"){
  sample_distances <- feature_table %>% 
    column_to_rownames("feature") %>%
    as.matrix() %>% 
    t() %>% 
    vegan::vegdist(method=method)
  
  return(sample_distances)
  
}

# PURPOSE: get a dataframe describing sample cateogories based on sample_name
# <feature_table>
# RETURNS tibble with  2 columns (sample_name, sample_type)
# NOTE this is dataset specific and not a totally generalizable function. change if neccesary
get_sample_types <- function(feature_table){
  summary_table <- feature_table %>% 
    colnames() %>% 
    as_tibble() %>% 
    dplyr::rename(sample_name = "value") %>% 
    dplyr::slice(-1) %>% # drop feature
    mutate(sample_type = case_when(str_detect(sample_name, "ZYMO") ~ "Zymo", # MGX
                               str_detect(sample_name, "hURNA") ~ "hURNA", # MTX
                               str_detect(sample_name, "CON_BAL") ~ "CONBAL", # MTX
                               str_detect(sample_name, "REP") ~ "REP", # MTX
                               str_detect(sample_name, "BAL") ~ "BAL", # MTX
                               .default = "Other"))
  return(summary_table)
}

# PURPOSE: subset a feature table to only BAL samples (default) or other subserseded tyepe
# <feature_table>
# <sample_type> = type string, matches string found in sample_type column output in get_sample_types
get_bal_samples <- function(feature_table, my_sample_type = "BAL"){
  summary_table <- get_sample_types(feature_table) 
  
  samples <- summary_table %>% 
    filter(sample_type == my_sample_type) %>% 
    pluck("sample_name")
  
  feature_table <- feature_table %>% 
    dplyr::select(1, all_of(samples))
  
  return(feature_table)
  
}

# PURPOSE: get tibble with column/sample sum from a feature table. useful for checking if any samples have zero features values after filtering
get_sample_sums <- function(feature_table){
  sample_sums <- enframe(colSums(feature_table[,-1]))
  return(sample_sums)
}

# PURPOSE: drop sample where all feature values are zero - helpful after filtering
# returns subset feature table
# NOTE: zero_samples is character vector listing samples with no values, put into global env.
drop_zero_sum_samples <- function(feature_table){
  sample_sums <- feature_table %>% 
    get_sample_sums()
  
  zero_samples <<- sample_sums %>% 
    filter(value == 0) %>% 
    pluck("name") 
  if (length(zero_samples) > 0){
    print(paste(as.character(length(zero_samples)), "zero-sum samples detected"))
    }
  feature_table <- feature_table %>% 
    dplyr::select(!all_of(zero_samples))
  
  return(feature_table)
}

# PURPOSE: print a summary of a given feature table, specifically # samples and # features, in the middle of a pipe series. To be used in middle of pipes for tracking
# <feature_table> =  type tibble, feature (rows) by samples (columns) table with feature column as first row with name "feature"
# RETURNS same feature table
get_intermediate_summary <- function(feature_table){
  print(paste("SUMMARY: Table has", dim(feature_table)[1], "features (rows) and", dim(feature_table)[2]-1, "samples (columns)"))
  return(feature_table)
}


# PURPOSE: drop features that have 0 abundance across all samples (i.e., 0 prevalence) in feature table; 
# RETURNS feature table
# NOTE: This often occurs when you subset samples
drop_zero_sum_features <- function(feature_table){
  zero_features <<- feature_table %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarise(sum = sum(value)) %>% 
    filter(sum == 0) %>% 
    distinct(feature) %>% 
    pluck("feature")
  
  if (length(zero_features) > 0){
    print(paste(as.character(length(zero_features)), "zero-sum features detected"))
  }
  
  feature_table <- feature_table %>% 
    filter(!feature %in% zero_features)
  
  return(feature_table)
  
}

# PURPOSE filter feature table to include only features with N min variance
# <feature_table>
# RETURNS feature table filtered
# NOTE global varaibles set too
get_feature_variance <- function(feature_table, min_variance=0.1){
  var_table <-feature_table %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarise(variance = var(value)) %>% 
    arrange(desc(variance))
  
  varied_features <<- var_table %>% 
    filter(variance>min_variance) %>% 
    pluck("feature")
  
  print(paste(length(varied_features), "features passing variance filtering"))
  
  feature_table <- feature_table %>% 
    filter(feature %in% varied_features)
  
  return(feature_table)
  
}

# PURPOSE filter feature table to include only features that covary with a  metadata value. Useful for plotting heatmaps
# <feature_table>
# <md> = type tibble, contains metadata
# <md_column> = type string, corresponds to column in md
# <min_covariance> = type double, [0-1], minimum level of covariance to include - if TRUE, takes standard-deviation/data-driven based approach
# <absolute> = type bool, include negative values too
# RETURNS feature table filtered
# NOTE global varaibles set too
get_feature_covariance <- function(feature_table, md, md_column, min_covariance=0.1, absolute=FALSE, num_sd=2){
  cov_table <<-feature_table %>% 
    dplyr::select(feature, any_of(md$sample_name)) %>% 
    tidy_features() %>% 
    left_join(md) %>% 
    group_by(feature) %>% 
    summarise(covariance = cov(value, .data[[md_column]])) %>% 
    arrange(desc(covariance))
  
  if (min_covariance==TRUE){
    # Get <num_sd> standard deviations above the mean
    min_covariance <- sd(cov_table$covariance)*num_sd + abs(mean(cov_table$covariance))
    print(paste("Data-driven selection - using min covariance of", min_covariance))
  }
  
  if (absolute){
    covaried_features <<- cov_table %>% 
      filter(abs(covariance)>min_covariance) %>% 
      pluck("feature")
  } else {
    covaried_features <<- cov_table %>% 
      filter(covariance>min_covariance) %>% 
      pluck("feature")
    
  }
  
  print(paste(length(covaried_features), "features passing covariance filtering"))
  
  feature_table <- feature_table %>% 
    filter(feature %in% covaried_features)
  
  return(feature_table)
  
}



# PURPOSE: Get feature tables as binaries for when presence/abscence testing is best suited
# <feature_table>
# RETURNS feature table with 1 for present and 0 for absent
get_presence_absence_table <- function(feature_table){
  feature_table <- feature_table %>% 
    tidy_features() %>% 
    mutate(value = case_when(value > 0 ~ 1, value == 0 ~ 0)) %>% 
    wide_features()
  return(feature_table)
}



# PURPOSE: Fuzzy search a feature table based on string
# <pattern> = type string, search term for finding features
# RETURNS Character vector features matching pattern
search_features <- function(feature_table, pattern){
  features <- feature_table %>% 
    filter(str_detect(feature, pattern)) %>%  # regex(pattern, ignore_case=T)
    pluck("feature") %>% 
    sort()
  
  return(features)
}


normalize_rna_dna = function(rna, dna){
  rna <- rna %>% tidy_features() %>% rename(value_rna = "value")
  dna <- dna %>% tidy_features() %>% rename(value_dna = "value")
  
  included_samples <- union(rna$sample_name, dna$sample_name)
  
  rna <- rna %>% filter(sample_name %in% included_samples)
  dna <- dna %>% filter(sample_name %in% included_samples)
  
  dna_smoother <- min(dna$value_dna[dna$value_dna>0])*0.5 # for replacing zeros
  
  dna <- dna %>% 
    mutate(value_dna = value_dna + dna_smoother)
  
  rna_dna_ratio <- full_join(rna, dna, by = c("sample_name", "feature"))
  
  rna_dna_ratio <<- rna_dna_ratio %>% 
    replace_na(list(value_dna =  dna_smoother,
                    value_rna = 0)) %>% 
    mutate(value = value_rna/value_dna)
  
}


# PURPOSE reproducibly plot 2 features from one feature table
plot_2_features <- function(feature_table, feature_1, feature_2){
  feature_subset <- feature_table %>% 
    tidy_features() %>% 
    filter(feature %in% c(feature_1, feature_2)) %>%
    pivot_wider(names_from = "feature", values_from = "value")
  
  ggplot(feature_subset, aes(x=.data[[feature_1]], y=.data[[feature_2]])) + 
    geom_smooth(color="grey", method="lm") +
    geom_point() +
    ggpubr::theme_pubr()
}


# PURPOSE reproducibly plot 2 features from 2 different features tables 
plot_2_mox <- function(feature_table_1, feature_table_2, feature_1, feature_2,
                       my_alpha=1, my_fill="black", my_color="black", my_method = "lm"){
  feature_subset_1 <- feature_table_1 %>% 
    tidy_features() %>% 
    filter(feature %in% c(feature_1)) %>%
    pivot_wider(names_from = "feature", values_from = "value")
  
  feature_subset_2 <<- feature_table_2 %>% 
    tidy_features() %>% 
    filter(feature %in% c(feature_2)) %>%
    pivot_wider(names_from = "feature", values_from = "value")
  
  # Avoid collisions
  if (feature_1 == feature_2){
    feature_subset_1 <- feature_subset_1 %>% 
      dplyr::rename_with(., ~ sprintf("%s [X]", .x), any_of(feature_1))
    
    feature_1 = sprintf("%s [X]", feature_1)
    
    feature_subset_2 <- feature_subset_2 %>% 
      dplyr::rename_with(., ~ sprintf("%s [Y]", .x), any_of(feature_2))
    
    feature_2 = sprintf("%s [Y]", feature_2)
    
  }
  
  feature_subset <<- inner_join(feature_subset_1, feature_subset_2, by = "sample_name") # Prevents collisions
  
  ggplot(feature_subset, aes(x=.data[[feature_1]], y=.data[[feature_2]])) + 
    geom_smooth(color="grey", 
                method=my_method) +
    geom_point(alpha = my_alpha,
               fill = my_fill,
               color = my_color) +
    ggpubr::theme_pubr()
}
# EXAMPLE plot_2_mox(mox_feature_list[["AMP [Genus]"]], mox_feature_list[["AMP [ASV]"]], "ASV22086_Corynebacterium", "ASV22086_Corynebacterium")


# PURPOSE: wrapper to easily and reproducibly get a prevalence table from a feature table
# <feature_table> = type tibble, feature (rows) by columns (samples) table with feature column as first row with name "feature"
# RETURNS a prevalence summary table with 13-column format 
#   (feature, prevalence, n_not_zero, n, mean_abundance, standard_deviation, median_abundance, median_absolute_deviation, 
#     maximum_abundance, minumum_abundance, minumum_abundance_above_zero, quantile_25, quantile_75) 
# EXAMPLE: get_prevalence_table(mox_feature_list[["AMP [Genus]"]])
get_prevalence_table <- function(feature_table){
  # get prev features
  prevalence_table <- feature_table %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarize(prevalence = sum(value > 0)/n(),
              n_not_zero = sum(value > 0),
              n = n(),
              mean_abundance = mean(value),
              standard_deviation = sd(value),
              median_abundance = median(value),
              median_absolute_deviation = mad(value),
              maximum_abundance = max(value),
              minumum_abundance = min(value),
              minumum_abundance_above_zero = min(if (prevalence == 0) {0} else {value[value != 0]}), # Catches zero prevalence features
              quantile_25 = quantile(value, probs=.25),
              quantile_75 = quantile(value, probs=.75)
              ) %>% 
    ungroup() %>% 
    arrange(desc(prevalence), desc(mean_abundance), desc(median_abundance))

  return(prevalence_table)
}
