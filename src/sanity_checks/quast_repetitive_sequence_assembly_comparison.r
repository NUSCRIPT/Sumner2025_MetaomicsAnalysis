
quast <- read_tsv("/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MGX/quast_results/results_2025_05_07_11_02_45/report.tsv")
quast <- quast %>% 
    rename(feature="Assembly") %>% 
    tidy_features() %>% 
    mutate(condition = case_when(str_detect(sample_name, "megahit_paired2") ~ "+Paired +16S +OE",
                                 str_detect(sample_name, "megahit_paired3") ~ "+Paired +16S -OE",
                                 str_detect(sample_name, "megahit_paired4") ~ "+Paired -16S +OE",
                                 str_detect(sample_name, "megahit_paired") ~ "+Paired -16S -OE",
                                 str_detect(sample_name, "megahit") ~ "-Paired -16S -OE"),
           sample_name = str_replace(sample_name, "megahit_", ""),
           sample_name = str_replace(sample_name, "paired[2,3,4]*_", ""),
           sample_name = str_replace(sample_name, "paired_", ""),
           sample_name = str_replace(sample_name, ".contigs", ""),
    )
unique(quast$feature)
# [1] "# contigs (>= 0 bp)"        "# contigs (>= 1000 bp)"     "# contigs (>= 5000 bp)"    
# [4] "# contigs (>= 10000 bp)"    "# contigs (>= 25000 bp)"    "# contigs (>= 50000 bp)"   
# [7] "Total length (>= 0 bp)"     "Total length (>= 1000 bp)"  "Total length (>= 5000 bp)" 
# [10] "Total length (>= 10000 bp)" "Total length (>= 25000 bp)" "Total length (>= 50000 bp)"
# [13] "# contigs"                  "Largest contig"             "Total length"              
# [16] "GC (%)"                     "N50"                        "N75"                       
# [19] "L50"                        "L75"                        "# N's per 100 kbp"         
quast %>% 
    # filter(feature == "Total length (>= 10000 bp)" ) %>% 
    ggplot(aes(x=condition,
               # y=log10(value+1)
               y=value
    )) +
    geom_violin(scale = "width", fill="grey")+
    geom_jitter(width =.1, alpha=.7) +
    geom_boxplot(width =.1,alpha=.9) +
    ggpubr::stat_pwc(ref.group = "-Paired -16S -OE",label.size = 1)+
    theme_pubr(base_size = 5) +
    rotate_x_text(45) +
    ggh4x::facet_wrap2(~feature, scales = "free_y", axes = "all")


quast %>% 
    filter(feature == "Largest contig") %>%
    ggplot(aes(x=condition,
               y=log10(value+1)
               # y=value
    )) +
    geom_violin(scale = "count", fill="grey")+
    geom_jitter(width =.1, alpha=.7) +
    geom_boxplot(width =.1,alpha=.9) +
    ggpubr::stat_pwc(ref.group = "-Paired -16S -OE")+
    theme_pubr() +
    rotate_x_text(45) 


quast %>% 
    filter(feature == "Total length (>= 1000 bp)") %>%
    ggplot(aes(x=condition,
               y=log10(value+1)
               # y=value
    )) +
    geom_violin(scale = "count", fill="grey")+
    geom_jitter(width =.1, alpha=.7) +
    geom_boxplot(width =.1,alpha=.9) +
    ggpubr::stat_pwc(ref.group = "-Paired -16S -OE")+
    theme_pubr() +
    rotate_x_text(45) 

quast %>% 
    filter(feature == "Total length (>= 1000 bp)") %>% 
    filter(value > 0) %>% 
    count(condition) %>% 
    ggplot(aes(x=condition,
               y=n
               # y=value
    )) +
    geom_bar(stat = "identity")+
    theme_pubr() +
    rotate_x_text(45) +
    ylab("No. w/ at least 1 contig > 1kb")


quast %>% pivot_wider(names_from = "feature", values_from = "value") %>% 
    filter(str_detect(sample_name, "ZYMO", 
                      negate = T
    )) %>% 
    ggplot(aes(log10(`# contigs (>= 1000 bp)`+1), log10(`Total length (>= 1000 bp)`+1), 
               color=condition,fill=condition)) +
    geom_point(alpha=.8) +
    geom_smooth(method="lm",se=F)+
    # geom_vline(xintercept = log10(1),linetype = "dashed", color="red") +
    theme_pubr() +
    rotate_x_text(45) +
    theme(legend.position = "left") 
# ggh4x::facet_grid2(~condition,axes="all",remove_labels = "all")

ggExtra::ggMarginal(last_plot(),groupFill = TRUE,)
