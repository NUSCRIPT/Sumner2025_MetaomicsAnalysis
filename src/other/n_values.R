library(tidyverse)
library(gtools)
source("src/common/clean_feature_names.r")
source("src/common/get_distances.r")
source("src/common/load_metadata.r")
source("src/common/load_features.r")
drop_samples()

mox_maaslin_list <- readRDS("objects/mox/MOX_01_maaslin_results.rds")
mox_maaslin_list$Category$`DNA [Taxonomy]`$results %>% 
    arrange(feature) %>% 
    filter(qval < 0.05) %>%
    mutate(stars = stars.pval(qval)) %>% 
    dplyr::select(feature, value, coef,qval,stars)

mox_maaslin_list$Category$`AMP [Genus]`$results %>% 
    arrange(feature) %>% 
    filter(qval < 0.05) %>%
    mutate(stars = stars.pval(qval)) %>% 
    dplyr::select(feature, value, coef,qval,stars)

mox_maaslin_list$Category$`DNA [PFAM]`$results %>% 
    mutate(feature = factor(feature, levels = unique(feature ))) %>% 
    arrange(feature) %>%
    filter(qval < 0.05) %>%
    mutate(stars = stars.pval(qval)) %>% 
    dplyr::select(feature, value, coef,qval,stars) %>% 
    head(40) %>% 
    # group_by(feature) %>% 
    mutate(feature = fix_ko(feature)) %>% 
    head(40)

mox_maaslin_list$Clusters$`AMP [Genus]`$results %>% filter(str_detect(feature, 'Parv')) %>% 
    # mutate(feature = factor(feature, levels = unique(feature ))) %>% 
    arrange(feature) %>%
    filter(qval < 0.05) %>%
    mutate(stars = stars.pval(qval)) %>% 
    dplyr::select(feature, value, coef,qval,stars)

mox_maaslin_list$Amylase$`RNA [Taxonomy BKN]`$results %>% 
    arrange(feature) %>% 
    filter(qval < 0.05) %>%
    mutate(stars = stars.pval(qval)) %>% 
    dplyr::select(feature, value, coef,qval,stars)

mox_maaslin_list$Category$`DNA [PFAM]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 100
mox_maaslin_list$Category$`AMP [Genus]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 6
mox_maaslin_list$Category$`DNA [Taxonomy]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() #2



mox_maaslin_list$Category$`DNA [PFAM]`$results %>%  ggplot(aes(x=coef, y=-log(qval))) + geom_point(size=.5) + facet_wrap(~value) + geom_hline(yintercept=-log(0.05))



mox_maaslin_list$Dysbiotic$`DNA [PFAM]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 929
mox_maaslin_list$Dysbiotic$`AMP [Genus]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 41
mox_maaslin_list$Dysbiotic$`DNA [Taxonomy]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() #6



mox_maaslin_list$Amylase$`DNA [PFAM]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 83
mox_maaslin_list$Amylase$`AMP [Genus]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 16
mox_maaslin_list$Amylase$`DNA [Taxonomy]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() #1

mox_maaslin_list$Amylase$`DNA [PFAM]`$results %>%  ggplot(aes(x=coef, y=-log(qval))) + geom_point(size=.5) + facet_wrap(~value) + geom_hline(yintercept=-log(0.05)) + theme(aspect.ratio = 1)

# Clusters

mox_maaslin_list$Clusters$`DNA [PFAM]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 1743
mox_maaslin_list$Clusters$`AMP [Genus]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() # 63
mox_maaslin_list$Clusters$`DNA [Taxonomy]`$results %>% filter(qval<0.05) %>% distinct(feature) %>% dim() #14

mox_maaslin_list$Clusters$`DNA [PFAM]`$results %>%  ggplot(aes(x=coef, y=-log(qval))) + geom_point(size=.5) + facet_wrap(~value) + geom_hline(yintercept=-log(0.05)) + theme(aspect.ratio = 1)
mox_maaslin_list$Clusters$`AMP [Genus]`$results %>% filter(qval<0.05) %>%  mutate(stars = stars.pval(qval)) %>% View()


mox_maaslin_list$Clusters$`AMP [Genus]`$results %>% filter(qval<0.05) %>% filter(value == 2, coef>0) %>%  mutate(stars = stars.pval(qval))
# Line for firmicutes tradeoffs in pneumotype staph
p1 <- features_amp %>% 
    select(feature, permmeta %>% filter(cluster_num_amplicon==2) %>% pluck('sample_name')) %>%  # get staph custer
    filter(feature %in% (mox_maaslin_list$Clusters$`AMP [Genus]`$results %>% filter(qval<0.05) %>% filter(value == 3, coef>0) %>% pluck('feature'))) %>%  # get sig clusters
    tidy_features() %>% 
    group_by(sample_name) %>% 
    mutate(replete = case_when(feature == 'ASV12894_Staphylococcus' & value < .75 ~ 'Staphylococcus Replete', 
                               feature == 'ASV12894_Staphylococcus' & value > .75 ~ 'Staphylococcus Predominant',
                               .default = NA)) %>% 
    fill(replete) %>%
    ungroup() %>% 
    mutate(feature = fix_asv(feature),
           feature = factor(feature, levels = sort(unique(feature), decreasing=T))) %>% 
    ggplot(aes(x=feature, y = value, group=sample_name, fill=replete)) + geom_line(aes(color=replete), alpha=.8) + geom_point(shape=21)  + scale_colour_manual(values=c('grey30', 'blue'))+ scale_fill_manual(values=c('grey30', 'blue')) + theme_nature() + ggpubr::rotate_x_text(60) + ylab('AST Abundance')+ xlab('Firmicutes taxa') 
ggsave('figures/staph_cluster_tradeoff.pdf', p1,units='in', width=4.31, height=2.23)



features_amp %>% 
    select(feature, permmeta %>% filter(cluster_num_amplicon==2) %>% pluck('sample_name')) %>%  # get staph custer
    tidy_features() %>% 
    group_by(feature) %>% 
    summarise(mean = mean(value), median = median(value), sd = sd(value), prev = sum(value > 0)) %>% arrange(desc(mean))


# halla

read_tsv('halla_June27/halla_network.txt') %>% count(x_type, y_type)

# N values for each omics at end of pipeline
for (i in seq_along(mox_dist_list)){
    permmeta %>% 
        filter(sample_name %in% 
               (mox_dist_list[[i]] %>%  labels() %>% grep('BAL', .,value = T))
               ) %>% 
        pluck('sample_name') %>% 
        unique() %>% 
        length() %>% 
        paste(., names(mox_dist_list)[[i]]) %>% 
        print(.)
}

# BAL with at least 1 good omics profiles (sanity check)
lapply(mox_dist_list[names(mox_dist_list)!="CFU [Culture]"], labels) %>% 
    unlist() %>% 
    enframe() %>% 
    distinct(value)

filter(permmeta, sample_name %in% (lapply(mox_feature_list, colnames) %>% unlist(use.names = F) %>% unique() %>% grep('BAL', .,value = T,invert = F))) %>% dim()

# N values for amp seq
readRDS("objects/AMP_08A_PhyloseqGenusAST.rds") %>% 
    phyloseq::otu_table() %>% 
    as.data.frame() %>% 
    rownames_to_column(var="feature") %>% 
    rename_with(.fn = renaming_helper, .cols = everything(), map_df = map_amp) %>% select(any_of(get_drop_lst())) %>% colnames()


# Pneumotype M prevalence cutibacterium 
pm <- permmeta %>% filter(cluster_num_amplicon ==2) %>% pluck("sample_name")
mox_feature_list[["AMP [Genus]"]] %>% tidy_features() %>% filter(feature=="ASV19300_Cutibacterium") %>% filter(value > 0) %>% filter(sample_name %in% pm)