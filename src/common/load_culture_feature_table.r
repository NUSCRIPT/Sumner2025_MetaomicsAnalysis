# Culture table to feature table

source("src/common/load_lists.r")
source("src/common/feature_tables_functions.r")

# PURPOSE: Make presence absence feature table from culture data <report_organisms> created by load_metadata.r
# RETURNS: feature table at genus or species level.
# <species_or_genus> = string, "species" or "genus" to determine which to return
get_culture_feature_table <- function(species_or_genus = "species"){
    
    # Make list of pretty-names joinable for culture table
    organisms_tbl <- enframe(unlist(organisms), name = "key", value = "value") %>% dplyr::rename(name = "key", new = "value")
    
    # Check that report organisms exists - loaded from load_metadata to global env
    if (!exists("report_organism")){
        get_report_for_culture_table()
    }
    
    # Report organisms to basic TIDY feature table
    feature_culture <- report_organism %>% 
        dplyr::rename(sample_name = "bal_barcode") %>% 
        # Give No Culture a value of 1
        mutate(organism_null = case_when(culture_positive == "Negative" ~ 100, 
                                         culture_positive != "Negative" ~ 0)) %>% 
        select(!c("culture_positive", "fungal_positive", "bacteria_positive")) %>% 
        pivot_longer(2:last_col())
    
    
    # Clean names and convert to proper feature table 
    # Also - SUM value from cases where multiple organisms are in new list 
    #   (e.g. Viridians_streptococcus_#1 & Viridians_streptococcus_#2 both to Streptococcus_species )
    feature_culture <- feature_culture %>% 
        left_join(organisms_tbl) %>% 
        mutate(name = new) %>% 
        select(!new) %>% 
        group_by(sample_name, name) %>% 
        reframe(value = sum(value)) %>% 
        ungroup() %>% 
        arrange(sample_name) %>% 
        wide_features() %>% 
        dplyr::rename(feature = name)  %>%
        get_presence_absence_table()
    
    
    feature_culture_genus <- feature_culture %>% 
        tidy_features() %>% 
        separate(feature, into = c("feature", NA), sep="_",extra = "merge") %>% 
        group_by(sample_name, feature) %>% 
        reframe(value = sum(value)) %>% 
        ungroup() %>% 
        arrange(sample_name) %>% 
        wide_features() %>%
        get_presence_absence_table()
    
    if (species_or_genus == "species"){
        return(feature_culture)
    } else if (species_or_genus == "genus"){
        return(feature_culture_genus)
    } else {
        print("INVALID ARGUMENT SELECT SPECIES OR GENUS")
    }
}

# This is a direct copy of lines 10-33 of load_metadata. This is copied to avoid a dependency crashout. 
# Yes, this should be separate function that metadata calls but this is current simplest solution
# ...woe betide the software developer who reads a biologists code...
get_report_for_culture_table <- function(){
    
    report <- read_csv('~/share_with_Jack/20240528_script_all_events_mapped_undated.csv') %>%  
        mutate(bal_barcode=str_remove_all(bal_barcode, '-')) %>% 
        mutate_at(vars(contains(c('_quantity'))), ~as.numeric(str_remove_all(.,',|>')))  # Fix first two col that are characters
    
    # Subset report to relevant nums
    report_organism <<-  report %>% 
        distinct(bal_barcode, pick(contains("organism"))) %>% 
        pivot_longer(contains('organism'),names_to = c("Var", ".value"),names_sep=11) %>% 
        distinct(bal_barcode, name, quantity) %>% 
        mutate(name = case_when(is.na(name) ~ 'organism_null', .default = name),
               quantity = case_when(is.na(quantity) ~ 0, .default = quantity)) %>% 
        # quantity = as.double(str_remove_all(quantity, ">|,"))) %>% 
        group_by(bal_barcode) %>% 
        mutate(culture_positive = case_when(any(quantity > 0) ~ 'Positive', .default = 'Negative'),
               fungal_positive = case_when(any(str_detect(name, 'Yeast|Candida')) ~ 'Positive', .default = 'Negative'),
               bacteria_positive = case_when(any(str_detect(name, 'Yeast|Candida|organism_null',negate = T)) ~ 'Positive', .default = 'Negative')) %>% 
        ungroup() %>% 
        mutate(name = str_remove_all(name, ',|\\)|\\('),
               name = str_replace_all(name, ' ', '_')) %>% 
        add_count(name) %>% # arrange by prevalence of culture result
        arrange(desc(n)) %>% #distinct(name, n) %>% View()
        select(!n) %>% 
        pivot_wider(names_from = 'name', values_from = 'quantity', values_fill = 0)
}


## TEST EXAMPLES
# get_culture_feature_table("species") %>% get_distance(method="jaccard") %>% my_ordination_helper() %>% visualize_ordination(., "bacteria_positive", "bacteria_positive")
# get_culture_feature_table("genus") %>% get_distance(method="jaccard") %>% my_ordination_helper() %>% visualize_ordination(., "bacteria_positive", "bacteria_positive")
