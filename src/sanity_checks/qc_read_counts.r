# For cleaning, graphing read QC tables produced from kneaddata 

library(tidyverse)
library(scales)
source('src/common/nature_theme.r')
mgx_qc = "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MGX/kneaddata/mgx_kneaddata_count_table.tsv"
mtx_qc = "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MTX/kneaddata/mtx_kneaddata_count_table.tsv"

# Function that takes kneaddata count table and converts to tidy format
# Also will sum pairs and orphans when data are paired
get_tidy_counts <- function(qc_file){
    
    counts <- read_tsv(qc_file)
    
    colnames(counts)
    # sum pairs and orphans
    tidy_counts <- counts %>% 
        mutate("raw total" = select(., starts_with("raw")) %>% rowSums(),
               "trimmed total" = select(., starts_with("trimmed")) %>% rowSums(),
               "decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA total" = select(., starts_with("decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA")) %>% rowSums(),
               "decontaminated chm13.draft_v1.0_plusY total" = select(., starts_with("decontaminated chm13.draft_v1.0_plusY")) %>% rowSums(),
               "decontaminated human_hg38_refMrna total" = select(., starts_with("decontaminated human_hg38_refMrna")) %>% rowSums(),
               "trimmed total" = select(., starts_with("trimmed")) %>% rowSums(),
               "final total" = select(., starts_with("final")) %>% rowSums(),
               ) %>% 
        pivot_longer(2:last_col())
    
    # separate qc step and read type for easy plotting 
    # note that final and decontam silva are same value when plottng 'total'
    tidy_counts <- tidy_counts %>% 
        separate(name, into = c("qc_step", "read"), sep = ' (?=[^ ]+$)') %>%  # final space
        mutate(qc_step = factor(qc_step, levels = c("raw", "trimmed", "decontaminated chm13.draft_v1.0_plusY", "decontaminated human_hg38_refMrna",
                                                    "decontaminated SILVA_128_LSUParc_SSUParc_ribosomal_RNA", "final")),
               Sample = factor(Sample)) %>% 
        mutate(control = case_when(str_detect(Sample, "ZYMO") ~ "Zymo", # MGX
                                   str_detect(Sample, "hURNA") ~ "hURNA", # MTX
                                   str_detect(Sample, "CON_BAL") ~ "CON_BAL", # MTX
                                   str_detect(Sample, "REP") ~ "REP", # MTX
                                   .default = "BAL")) %>% 
        arrange(Sample, desc(qc_step)) 
    return(tidy_counts)
}

mgx_counts <- get_tidy_counts(mgx_qc)
tidy_counts_total <- mgx_counts %>% filter(read == "total") %>% filter(qc_step != "decontaminated human_hg38_refMrna")


# Plot Read Counts
ggplot(tidy_counts_total , aes(x=value, y=qc_step)) +
    geom_line(
        aes(group=Sample),
        position=position_jitter(height = .3, width=0, seed = 100), 
        alpha=0.25, color="grey") +
    geom_point(aes(color=control),
        position = position_jitter(height = .3, width=0, seed = 100),
        size = 1, alpha = .5) +
    geom_vline(xintercept=20000, linetype="dashed", 
               color = "red", size=.5) +
    annotate(
        "text", 
        x = 20000, 
        y = 0, 
        color='red',
        size=2,
        hjust=-0.15,
        vjust=-.5, 
        label = as.character(expression(2%*%10^{4})), 
        parse=T) +
    labs(
        x="Read count", 
        y = "Processing step",
        title = "MGX Read Quality Control") +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_discrete(
        labels = function(x) 
            str_replace_all(x,"_",replacement = " ") %>% 
            str_wrap(., width = 20)) +
    theme_nature() +
    theme(
        text = element_text(size=10),
        axis.text = element_text(size=8),
        plot.margin = margin(10,10,10,10),
        axis.line.x = element_line(size=.4),
        axis.line.y = element_line(size=.4),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.title=element_blank()
    ) +
    scale_color_manual(values = c("black", "gold"))

ggsave(filename = "figures/qc_reads/mgx_qc_read_counts.pdf", units = "in", width = 4.38, height=3.41, device = cairo_pdf)

mtx_qc = "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MTX/kneaddata/mtx_kneaddata_count_table.tsv"
mtx_counts <- get_tidy_counts(mtx_qc)
tidy_counts_total <- mtx_counts %>% filter(read == "total") %>% 
    filter(str_detect(Sample,"CON_BAL", negate = T))

# Plot MTX Read Counts
ggplot(tidy_counts_total , aes(x=value, y=qc_step)) +
    geom_line(
        aes(group=Sample),
        position=position_jitter(height = .3, width=0, seed = 100), 
        alpha=0.25, color="grey") +
    geom_point(aes(color=control),
               position = position_jitter(height = .3, width=0, seed = 100),
               size = 1, alpha = .5) +
    geom_vline(xintercept=1e5, linetype="dashed", 
               color = "red", size=.5) +
    annotate(
        "text", 
        x = 1e5, 
        y = 0, 
        color='red',
        size=2,
        hjust=-0.15,
        vjust=-.5, 
        label = as.character(expression(1%*%10^{5})), 
        parse=T) +
    labs(
        x="Read count", 
        y = "Processing step",
        title = "MTX Read Quality Control") +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_discrete(
        labels = function(x) 
            str_replace_all(x,"_",replacement = " ") %>% 
            str_wrap(., width = 20)) +
    theme_nature() +
    theme(
        text = element_text(size=10),
        axis.text = element_text(size=8),
        plot.margin = margin(10,10,10,10),
        axis.line.x = element_line(size=.4),
        axis.line.y = element_line(size=.4),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.title=element_blank()
    ) +
    scale_color_manual(values = c("black", "gold", "red", "purple"))

ggsave(filename = "figures/qc_reads/mtx_qc_read_counts.pdf", units = "in", width = 4.38, height=3.41, device = cairo_pdf)

# correlate mgx mtx reads
mgx_final <- mgx_counts %>% filter(qc_step == "final", read == "total")
mtx_final <- mtx_counts %>% filter(qc_step == "final", read == "total")

mgx_mtx <- full_join(mgx_final,mtx_final, by = join_by(Sample, qc_step, read))




ggplot(mgx_mtx, (aes(x=value.x, y = value.y))) +
    geom_point() +
    theme_nature() +
    theme(
        text = element_text(size=10),
        axis.text = element_text(size=8),
        plot.margin = margin(10,10,10,10),
        axis.line.x = element_line(size=.4),
        axis.line.y = element_line(size=.4),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.title=element_blank()
    ) +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x="MGX final reads", y="MTX final reads")

ggsave(filename = "figures/qc_reads/mgx_mtx_reads_scatter.pdf", units = "in", width = 2.52, height=2.05, device = cairo_pdf)

source("src/common/load_metadata.r")
mgx_final <- mgx_final %>% 
    mutate(Sample = str_remove(Sample, "PT_"))
mgx_final_md<- mgx_final %>% 
    left_join(permmeta, by=c("Sample" = "sample_name"))

mgx_final_md %>% 
    # filter(is.na(Episode_etiology)) %>% 
    filter(!Sample %in% get_drop_lst()) %>% 
    print(n=100)

ggplot(mgx_final_md %>% filter(!Sample %in% get_drop_lst()), aes(x=value, y=Episode_etiology)) + 
    geom_point(aes(color=control),
               position = position_jitter(height = .3, width=0, seed = 100),
               size = 1, alpha = .5) +
    geom_vline(xintercept=20000, linetype="dashed", 
               color = "red", size=.5) +
    annotate(
        "text", 
        x = 20000, 
        y = 0, 
        color='red',
        size=2,
        hjust=-0.15,
        vjust=-.5, 
        label = as.character(expression(2%*%10^{4})), 
        parse=T) +
    labs(
        x="Read count", 
        y = "Episode etiology",
        title = "MGX Read Count by Etiology") +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_discrete(
        labels = function(x) 
            str_replace_all(x,"_",replacement = " ") %>% 
            str_wrap(., width = 20)) +
    theme_nature() +
    theme(
        text = element_text(size=10),
        axis.text = element_text(size=8),
        plot.margin = margin(10,10,10,10),
        axis.line.x = element_line(size=.4),
        axis.line.y = element_line(size=.4),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.title=element_blank()
    ) +
    stat_summary(fun.data = median_hilow, shape = 18, geom = "pointrange", color="grey") +
    scale_color_manual(values = c("black", "gold"))

ggsave(filename = "figures/qc_reads/mgx_final_read_count_etiology.pdf", units = "in",width = 4.38, height=3.41, device = cairo_pdf)


# 
# 
# mgx_final <- mgx_final %>% 
#     mutate(Sample = str_remove(Sample, "PT_"))
# mgx_final_md<- mgx_final %>% 
#     left_join(permmeta, by=c("Sample" = "sample_name"))
# 
# mgx_final_md %>% 
#     # filter(is.na(Episode_etiology)) %>% 
#     filter(!Sample %in% get_drop_lst()) %>% 
#     print(n=100)
# 
# ggplot(mgx_final_md %>% filter(!Sample %in% get_drop_lst()), aes(x=value, y=Episode_etiology)) + 
#     geom_point(aes(color=control),
#                position = position_jitter(height = .3, width=0, seed = 100),
#                size = 1, alpha = .5) +
#     geom_vline(xintercept=20000, linetype="dashed", 
#                color = "red", size=.5) +
#     annotate(
#         "text", 
#         x = 20000, 
#         y = 0, 
#         color='red',
#         size=2,
#         hjust=-0.15,
#         vjust=-.5, 
#         label = as.character(expression(2%*%10^{4})), 
#         parse=T) +
#     labs(
#         x="Read count", 
#         y = "Episode etiology",
#         title = "MGX Read Count by Etiology") +
#     scale_x_log10(
#         breaks = scales::trans_breaks("log10", function(x) 10^x),
#         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     scale_y_discrete(
#         labels = function(x) 
#             str_replace_all(x,"_",replacement = " ") %>% 
#             str_wrap(., width = 20)) +
#     theme_nature() +
#     theme(
#         text = element_text(size=10),
#         axis.text = element_text(size=8),
#         plot.margin = margin(10,10,10,10),
#         axis.line.x = element_line(size=.4),
#         axis.line.y = element_line(size=.4),
#         plot.title = element_text(size=10, hjust = 0.5),
#         axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
#         axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
#         legend.title=element_blank()
#     ) +
#     stat_summary(fun.data = median_hilow, shape = 18, geom = "pointrange", color="grey") +
#     scale_color_manual(values = c("black", "gold"))
# 




mtx_final <- mtx_final %>% 
    mutate(Sample = str_remove(Sample, "PT_"))
mtx_final_md<- mtx_final %>% 
    left_join(permmeta, by=c("Sample" = "sample_name"))

mtx_final_md %>% 
    filter(is.na(Episode_etiology)) %>%
    filter(!Sample %in% get_drop_lst()) %>% 
    print(n=100)

ggplot(mtx_final_md %>% filter(!Sample %in% get_drop_lst()), aes(x=value, y=Episode_etiology)) + 
    geom_point(aes(color=control),
               position = position_jitter(height = .3, width=0, seed = 100),
               size = 1, alpha = .5) +
    geom_vline(xintercept=1e5, linetype="dashed", 
               color = "red", size=.5) +
    annotate(
        "text", 
        x = 1e6, 
        y = 0, 
        color='red',
        size=2,
        hjust=-0.15,
        vjust=-.5, 
        label = as.character(expression(1%*%10^{6})), 
        parse=T) +
    labs(
        x="Read count", 
        y = "Episode etiology",
        title = "MTX Read Count by Etiology") +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_discrete(
        labels = function(x) 
            str_replace_all(x,"_",replacement = " ") %>% 
            str_wrap(., width = 20)) +
    theme_nature() +
    theme(
        text = element_text(size=10),
        axis.text = element_text(size=8),
        plot.margin = margin(10,10,10,10),
        axis.line.x = element_line(size=.4),
        axis.line.y = element_line(size=.4),
        plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.title=element_blank()
    ) +
    stat_summary(fun.data = median_hilow, shape = 18, geom = "pointrange", color="red") +
    scale_color_manual(values = c("black", "gold", "red", "purple"))

ggsave(filename = "figures/qc_reads/mtx_final_read_count_etiology.pdf", units = "in", width = 4.38, height=3.41, device = cairo_pdf)




