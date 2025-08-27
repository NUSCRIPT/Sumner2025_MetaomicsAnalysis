

read_kraken_uniq <- function(directory, kmers_n=1000){
    # function to read in a tsv and add the file name as a column
    customized_read_tsv <- function(file){
        read_tsv(file, skip=2) %>%
            mutate(fileName = basename(file))
    }
    
    
    kraken <- list.files(path = directory,
                         pattern = "*_krakenuniq_report.tsv", full.names = TRUE, recursive = TRUE) %>%  # list all the files
        lapply(customized_read_tsv) %>% # read them all in with our custom function
        reduce(bind_rows) %>% 
        mutate(fileName = str_remove(fileName, "_krakenuniq_report.tsv")) %>% 
        mutate(fileName = str_remove(fileName, "PT_"))
    
    
    spp <- kraken %>% filter(rank == "genus"|taxName =="unclassified")
    spp %>% filter(taxName=="unclassified") %>% summarize(m = mean(`%`)) %>% print()
    
    
    
    feature_table <- spp %>% 
        filter(str_detect(taxName, "Plasmodium|Toxoplasma|Homo", negate = T)) %>% # Known Database Contaminants https://doi.org/10.1371/journal.pcbi.1006277
        filter(kmers>kmers_n) %>% # Filter by Kmers as recommended by manuscript
        rename(feature = taxName,
               sample_name = fileName,
               value = reads) %>% 
        select(feature, sample_name, value) %>% 
        wide_features() %>% 
        tidy_features() %>% 
        replace_na(list(value = 0)) %>% 
        wide_features()
    
    return(feature_table)
}
# 
# bkn <- read_kraken_uniq(DIRECTORY_mtx_kraken_uniq)
# bkn <- bkn %>% 
#     drop_zero_sum_samples() %>% 
#     normalize_feature_table() %>% 
#     filter(feature != "unclassified") %>% 
#     get_most_prevalent(min_prevalence = .001,min_abundance = 0.01) %>% 
#     drop_zero_sum_features() %>% 
#     normalize_feature_table(normalize = F, transform = "log100")
# 
# pheatmap::pheatmap(bkn %>% get_most_abundant_n(20) %>% column_to_rownames("feature"))

