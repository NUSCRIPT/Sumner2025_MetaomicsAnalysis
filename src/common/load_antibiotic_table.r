
# PURPOSE: get an antibiotic table from metadata
# NOTE: This contains is a direct copy of lines 10-12 of load_metadata. 
get_abx_table <- function(){
    
    report <- read_csv('~/share_with_Jack/20240528_script_all_events_mapped_undated.csv') %>%  
        mutate(bal_barcode=str_remove_all(bal_barcode, '-')) %>% 
        mutate_at(vars(contains(c('_quantity'))), ~as.numeric(str_remove_all(.,',|>')))  # Fix first two col that are characters
    
    abx_table <- report %>% 
        distinct(bal_barcode, base_medication_name,administered_dose) %>% 
        rename(sample_name = bal_barcode, feature = base_medication_name, value = administered_dose) %>% 
        select(feature, sample_name, value) %>% 
        mutate(value = 1) %>% # drop this if you want to take into account the dosage
        complete(sample_name, feature) %>% 
        # drop_na(feature) %>% # Want to take into account instance of no listed ABX?
        replace_na(replace = list(value = 0,
                                  feature = "NoMed"
        )) %>% 
        arrange(sample_name,feature) %>% 
        group_by(sample_name, feature) %>% 
        reframe(value = sum(value)) %>%  # Assume that multiple values means multiple administrations in the same ~day
        ungroup() %>% 
        wide_features() %>% 
        drop_zero_sum_samples() %>%
        drop_zero_sum_features()
        
    
    return(abx_table)
}

# Some Test cases and Examples
# source("src/common/ordination_functions.r")
# source("src/common/load_metadata.r")
# 
# get_abx_table() %>%
#     drop_zero_sum_samples() %>%
#     drop_zero_sum_features() %>%
#     get_distance(method="jaccard") %>%     my_ordination_helper() %>%
#     visualize_ordination(., "bacteria_positive", "bacteria_positive")
