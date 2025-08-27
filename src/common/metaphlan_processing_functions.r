require(tidyverse)

# FOR PROCESSING METAPHLAN FILES. NOTABLY (1) READING THEM IN AND (2) BASIC PLOTTING

# Takes a metaphlan file (all stratifications) and processes it to species level a species level (s__) tibble
# <mpa_file> =  string instance with /path/to/metaphlan_all.tsv
# returns 
#   tibble object with [ rows = species, columns = samples ]. 
#   Samples with unknown feature name returns as HighestClassified_incertae$GGB123_SGB123 - note the $ 
read_metaphlan <- function(mpa_file){
    
    ds <- read_tsv(mpa_file, skip = 1) %>% rename(feature = "clade_name")
    names(ds) <- gsub(names(ds), pattern = "PT_", replacement = "") 
    
    # dst <- ds %>% filter(str_detect(feature, "t__|UNCLASSIFIED"))
    dss <- ds %>% 
        filter(str_detect(feature, "s__|UNCLASSIFIED")) %>%
        filter(str_detect(feature, "t__", negate=T))
    
    dss <- dss %>% 
        # get most defined by deleting NGBs besides SGBs. 
        mutate(feature = str_replace_all(feature, 
                                         c("\\|[fg]__[FG]GB[0-9]+\\b" = "", 
                                           "\\|[pco]__[PCO]FGB[0-9]+\\b" = ""))) %>%  
        # alt separators for final |s__ so can go by end | rather than first | --- | for known species and $ for SGBs
        mutate(feature = case_when(str_detect(feature, "\\|s__GGB") ~ str_replace(feature, "\\|s__", "_incertae$"), 
                                   .default = str_replace(feature, "\\|s__", "|")
        )) %>% 
        # clips off known features
        mutate(feature = str_replace_all(feature,
                                         c("k__Bacteria\\|" = "",
                                           "k__Eukaryota\\|" = "",#))) %>% 
                                           "[pcofg]__[[:alnum:]_]+\\|\\b" = ""))) %>%
        # clean up final remannt
        mutate(feature = str_replace_all(feature,
                                         c("[pcofg]__" = ""))) 
    return(dss)
}
# example: 
# ds_fi <- "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
# dss <- read_metaphlan(ds_fi)
# head(dss)


# For plotting metaphlan barplots
clean_mpa4bar <- function(mpa){
    mpa_tmp <- mpa %>% 
        pivot_longer(cols = 2:last_col(), names_to="Sample", values_to="Abundance") %>%
        separate(Sample, into="Sample", sep=".metaphlan_profile", extra = "drop")
    
    mpa_tmp_grouped <- mpa_tmp %>%
        group_by(feature) %>% 
        mutate(sum_gen = sum(Abundance)) %>%
        filter(sum_gen > 0) %>%
        select(!sum_gen) %>%
        mutate(meanGeneraAbundancePerSample = mean(Abundance))
    
    top_genera <- mpa_tmp_grouped %>% filter(meanGeneraAbundancePerSample > 1) %>% 
        select(feature) %>% 
        unique() %>% 
        pluck(1)
    top_genera
    
    
    top_abund <- mpa_tmp %>% 
        filter(feature %in% top_genera)
    
    
    low_abund <- mpa_tmp %>% filter(!feature %in% top_genera) %>% mutate(feature = "Z_RareTaxa")
    
    low_abund_melt <- low_abund %>% 
        group_by(Sample) %>% 
        mutate(AbundanceOtherTotal = sum(Abundance)) %>% 
        ungroup() %>% 
        select(!Abundance) %>%
        unique() %>%
        rename(Abundance = AbundanceOtherTotal) %>%
        relocate(Sample)
    
    
    mpa_clean <- rbind(top_abund, low_abund_melt) 
    return(mpa_clean)
}
