
get_spp_kos <- function(){
    # Get genes that only strongly, positively associate with only a single taxa
    # Remove so taxa and gene content profile aren't directly confounding
    file_ko2taxa <-  "./halla/outputs/dna_taxonomy.dna_ko/all_associations.txt"
    
    ko2taxa <-read_halla_results(file_ko2taxa)
    ko2taxa <- ko2taxa %>% filter(association > .6) 
    ko2taxa_counts <- ko2taxa %>% count(Y_features) 
    
    p.hist <-ggplot(ko2taxa_counts, aes(x=n)) + 
        geom_histogram(color="black", fill="white")
    
    spp_ko <- ko2taxa_counts %>% filter(n==1) %>% pluck('Y_features') 
    
    # # Sanity check
    # ko2taxa %>% filter(Y_features %in% spp_ko) %>% count(X_features) %>% arrange(desc(n)) %>% print(n=Inf)
    # ko2taxa %>% filter(!Y_features %in% spp_ko) %>% count(X_features) %>% arrange(desc(n)) %>% print(n=Inf)
    
    return(spp_ko)
}



