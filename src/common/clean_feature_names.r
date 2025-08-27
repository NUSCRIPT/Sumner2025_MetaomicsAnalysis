library(KEGGREST)

ko_detect <- "K[1234567890]+: |K[1234567890]+\\.\\." # "K[1234567890]+_" -> now handle updated format + maaslin
asv_detect <- "ASV[1234567890]+_"
host_detect <- 'ENSG[1234567890]'
votu_detect <- "vOTU[1234567890]+-"


fix_asv <- function(feature){
  
  feature = str_replace(feature, asv_detect, '')
  feature = str_replace_all(feature, "_", "-")
  feature = str_replace_all(feature, "\\[|\\]", "")
  
  return(feature)
}

# For keeping ASV value around put pretty
fix_asv_2 <- function(feature){
  
  asv = feature
  feature = str_replace(feature, asv_detect, '')
  asv = str_replace(asv, feature, '') %>% str_replace(., "_", "")
  feature = str_replace_all(feature, "_", "-")
  feature = str_replace_all(feature, "\\[|\\]", "")
  feature = paste0(feature, " (", asv, ")")
  
  return(feature)
}

fix_ko <- function(feature){
  
  feature = str_split_i(feature, ":|\\.", i = 1) # "_"
  feature = str_replace(feature, 'k', 'K')
  feature_sanity_check =  as_tibble(list(kegg = feature))
  feature = keggList(feature) %>%
    as_tibble(rownames='kegg') %>%
    separate(value, into=c('gene_name', 'value'), sep='; ') %>%
    mutate(value = str_split_i(value, " \\[EC", i=1), # rm EC num,
           gene_name = case_when(gene_name == kegg ~ value, .default = gene_name),
           value = paste0(gene_name, ' (', kegg, ')'))
  # Sometimes KEGGs in humann3 maps no longer in KO database - keep features as unlistedKEGG (KO)
  feature = feature_sanity_check %>% 
    left_join(feature) %>% 
    mutate(value = case_when(is.na(value) ~ paste0("unlistedKEGG ", '(',kegg,')'), .default = value))
  feature = feature$value
  # feature = str_split_i(feature, " \\[EC", i=1) # rm EC num
  
  return(feature)
}


fix_general <- function(feature){
  
  feature = str_replace_all(feature, "_", " ")
  feature = str_replace_all(feature, "\\[|\\]", "")
  
  return(feature)
}


fix_host_genes <- function(feature){
  
  feature <- paste0(
    str_remove(feature, ".*: "), # gene
    ' (', 
    str_split_i(feature, ":", i = 1), # ensembl
    ')')
  
  return(feature)
}


fix_votu <- function(feature){
  
  feature <- paste0(
    str_remove(feature, ".*-"), # gene
    ' (', 
    str_split_i(feature, "-", i = 1), # ensembl
    ')')
  
  return(feature)
}


fix_features <- function(feature){
  print(paste("INPUT:", feature))
  italicize = TRUE
  
  if ( str_detect(feature[1], asv_detect) ) {
    feature <- fix_asv_2(feature)
    
  } else if ( str_detect(feature[1], ko_detect) ){
    feature <- tryCatch(
      {fix_ko(feature)}, 
      error = function(msg){
        message("KEGG not found error, printing errors...")
        message(msg)
        print(msg)
        # return(feature)
      },finally = feature)
    italicize = FALSE
    
  } else if ( str_detect(feature[1], host_detect) ){
    feature <- fix_host_genes(feature)
    italicize = TRUE
    
  } else if ( str_detect(feature[1], votu_detect) ){
    feature <- fix_votu(feature)
    italicize = TRUE
    
  } else {
    print("general")
    feature <- fix_general(feature)
    
  }
  
  if ( italicize ){
    feature <- paste0("italic(\"", feature, "\")")
    
  } else {
    feature <- paste0("\"", feature, "\"")
    
  }
  
  return(feature)
}


## Test cases
# asv_features <- head(features_amp$feature)
# ko_features <- head(features_mgx_ko$feature)
# tax_features <- head(features_mgx$feature)
# 
# fix_ko(ko_features)
# fix_asv(asv_features) # no italic
# fix_general(tax_features)
# 
# fix_features(ko_features)
# fix_features(asv_features) # yes italc
# fix_features(tax_features)
