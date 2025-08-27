#Preparation ----
source('src/common/load_metadata.r')
source("src/common/load_features.r")

library(tidyverse)

set.seed(100)

read_htx <- function(file_path){
    print(file_path)
    raw_htx <-  read_tsv(file_path,
                         skip = 1) %>%
        rename_with(~ str_replace(., "-BAL-", "BAL"), everything()) %>%
        dplyr::select(!Gene_Symbol) %>%
        column_to_rownames("Symbol")
    return(raw_htx)
}

make_map <- function(file_path){
    print(file_path)
    raw_map <-  read_tsv(file_path,
                         skip = 0) %>%
        filter(Symbol=="Symbol") %>% 
        dplyr::select(3:last_col()) %>%
        pivot_longer(1:last_col()) %>%
        mutate(value = str_replace(value, "-BAL-", "BAL")) %>% 
        return(raw_map)
}

## Read Files ----
htx_base_path <- "/path/to/my_data/htx/raw_counts/"

htx_files <- list.files(path = htx_base_path, full.names=TRUE,recursive=FALSE, pattern = ".txt") %>%
  str_subset(., ".xlsx", negate=TRUE) # Remove metadata file

htx_data_list<- lapply(htx_files, read_htx)

htx_raw <- do.call("cbind", htx_data_list) %>% 
  subset(., select=which(!duplicated(rev(names(.)))))  # remove the double-sequenced samples


htx_map_list <- lapply(htx_files, make_map)
htx_map <- do.call("rbind", htx_map_list)

head(htx_raw)
############################

# Get and remove protein coding genes
get_protein_genes <- function(htx_raw){
    library("biomaRt")
    yourListofGenes<- rownames(htx_raw)
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    myResult <- getBM(values = yourListofGenes, attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"), mart=mart)
    protein_genes <- myResult %>% filter(gene_biotype == 'protein_coding') 
    # For converting later
    ensembl2symbol <<- protein_genes %>% 
      dplyr::rename(feature = "ensembl_gene_id") %>% 
      mutate(hgnc_symbol = paste0(feature, ": ", hgnc_symbol)) %>% 
      distinct(feature, hgnc_symbol)
    # Get protein genes
    protein_genes <- protein_genes %>% pluck('ensembl_gene_id') 
    return(protein_genes)
}

protein_genes <- get_protein_genes(htx_raw) 
detach("package:biomaRt", unload = TRUE)
htx_raw <- htx_raw %>% 
  rownames_to_column('feature') %>% 
  filter(feature %in% protein_genes) %>% 
  left_join(ensembl2symbol) %>% 
  mutate(feature = hgnc_symbol) %>% 
  dplyr::select(!hgnc_symbol) %>% 
  column_to_rownames('feature')

dim(htx_raw)
head(htx_raw)
################################
#amplicon samples that will be used to subset HTX data prior to pre-filtering
amp_samples <- permmeta %>% 
  drop_na(dysbiosis_score) %>% 
  pluck("sample_name")

# Subset to only AMP samples & determine prevalence @ N reads
htx_clean <- htx_raw %>%
    select(any_of(amp_samples)) %>% 
    rownames_to_column("Symbol") %>%
    select(contains("Symbol")|contains("BAL")) %>%
    select(!contains(c("CON", "hURNA", "REDACTEDBALID1", "REDACTEDBALID2"))) %>% # drop samples with few genes reads(<800) / few reads(<1e6)
    rowwise() %>%
    mutate(
        prev = sum(c_across(2:last_col()) > 10), # at least ten reads, common for bulk RNAseq
        .before = 2
    ) %>%
    ungroup()

# number of samples for to represent 10% of sampled population for prevalence filtering
n_prevalence_10 <- (length(colnames(htx_clean))-2)*.1

htx_clean <- htx_clean %>% 
    filter(prev >= n_prevalence_10) %>%
    arrange(desc(prev)) %>%
    select(!c(prev))

# Drop to only relevant BAL samples
htx_clean <- htx_clean %>% 
    select(contains("Symbol")|contains("BAL")) %>%
    select(!contains("CON")) %>%
    column_to_rownames("Symbol")

# Remove technical duplicates 
htx_clean <- htx_clean %>% select(!contains(".1")) # NEW

# Generate metadata for DESeq2 based on permmeta
sample_names <- colnames(htx_clean)

exp_matrix <- as_tibble(sample_names) %>% 
    filter(str_detect(value, "\\.",negate=TRUE)) %>%
    mutate(condition = case_when(str_ends(value, "00") ~ "Baseline",
                                 T ~ "BAL")) %>%
    dplyr::rename(sample_name = "value") %>% 
    left_join(permmeta, by=c('sample_name' = "sample_name")) %>% 
    mutate(cluster_num_amplicon = factor(cluster_num_amplicon, levels=c("1","2","3","4"))) #as.factor(cluster_num_amplicon))

# Add Batch Info
exp_matrix <- exp_matrix %>% 
  left_join(htx_map %>% filter(str_detect(value, "\\.",negate = T)), by = c("sample_name" = "value"),multiple = "first") %>% 
  separate(name, into = c(NA, "batch", NA), sep = "_") 

# Subset to only samples with Pneumotype data & 
# Prepare for DESeq2
exp_matrix_pneumo <- exp_matrix %>% 
  drop_na(cluster_num_amplicon) %>%
  mutate(sample_name_dup = sample_name) %>% 
  column_to_rownames("sample_name_dup")

htx_deseq <- htx_clean %>% 
  select(exp_matrix_pneumo$sample_name)
htx_deseq %>% rownames_to_column("feature") %>% get_sample_sums() %>% arrange(value)
library(DESeq2)
#### NEW CODE

#DESeq ----
## Pneumotypes ----
# Instantiate DESeq Object
dds <- DESeqDataSetFromMatrix(countData=htx_deseq, 
                              colData=exp_matrix_pneumo,
                              design= ~ batch + cluster_num_amplicon)

# dds <- DESeq(dds, fitType = "local")
dds <- DESeq(dds,
             betaPrior =T, sfType = "poscount",
             test="Wald"
             )

table(colData(dds)$cluster_num_amplicon)
# 1  2  3  4 
# 47 25 15 19 

resultsNames(dds)
# Should've just made this a for-loop...
pairwise_res<- rbind(
  results(dds,name =  "cluster_num_amplicon1",lfcThreshold = 0, altHypothesis = "greaterAbs") %>% 
    as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "1"),
  
  results(dds,name = "cluster_num_amplicon2",altHypothesis = "greaterAbs") %>% 
    as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "2"),
  
  results(dds,name =  "cluster_num_amplicon3",altHypothesis = "greaterAbs") %>% 
    as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "3"),
  
  results(dds,name = "cluster_num_amplicon4",altHypothesis = "greaterAbs") %>% 
    as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "4")
  
)

# Set betaprior to FALSE  and uncomment this if you really want to do pairwise. 
# Logically, I think the comparison to the mean / not cluster makes more sense but is influenced by the distribution of the clusters [ie more common clusters will be more powerful which is better distilled in pairwise...tradeoffs]
# Also Log2FC are way higher in this one but i dnot think its meaningfully different - actually plotting the factors are outlier skewed
# # With or without lfcShrink shows some difference, but largely same message
# pairwise_res <- rbind(
#   results(dds, contrast=c("cluster_num_amplicon", "4", "1")) %>% 
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "4", "1"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "4vs1"),
#   
#   results(dds, contrast=c("cluster_num_amplicon", "3", "1")) %>% 
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "3", "1"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "3vs1"),
#   
#   results(dds, contrast=c("cluster_num_amplicon", "2", "1")) %>% 
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "2", "1"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "2vs1"),
#   
#   results(dds, contrast=c("cluster_num_amplicon", "4", "2")) %>%
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "4", "2"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "4vs2"),
# 
#   results(dds, contrast=c("cluster_num_amplicon", "3", "2")) %>% 
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "3", "2"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "3vs2"),
#   
#   results(dds, contrast=c("cluster_num_amplicon", "4", "3")) %>% 
#     # lfcShrink(dds, contrast=c("cluster_num_amplicon", "4", "3"), res=., type = 'ashr') %>%
#     as_tibble(rownames='ENSEMBL') %>% arrange(pvalue) %>% mutate(comparison = "4vs3")
#   )


# Subset to significant q values and fold change
res_sig <- pairwise_res %>% 
  filter(padj < 0.05, abs(log2FoldChange) > .5) %>% 
  # Be able to get top N in each group wise comparison
  arrange(padj,desc(abs(log2FoldChange))) %>% 
  group_by(comparison) %>% 
  mutate(idx = row_number()) %>% 
  ungroup()

pairwise_res <- pairwise_res %>% 
  mutate(significance = 
           case_when(padj < 0.05 & abs(log2FoldChange) > .5 ~ "|log2FC| > 0.5 & FDR P < 0.05", 
                     .default= "not significant")
         ) %>% 
  arrange(significance, (padj))

volcano <- ggplot(pairwise_res, 
       aes(x=log2FoldChange, y=-log10(padj), fill=significance)) +
  geom_vline(xintercept = 0.5,linetype="dashed", color="grey") +
  geom_vline(xintercept = -0.5,linetype="dashed", color="grey") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey") +
  geom_point(color="black", shape=21, alpha=.9, size=1,stroke=.25) +
  scale_fill_manual(values = c("red", "black")) +
  theme_nature() +
  ggh4x::facet_wrap2(~comparison,axes = "all", scale="free_x") +
  xlim(-1.5, 1.5)
volcano

ggsave("figures/host_deseq/volcano_plot.pdf", units="in", width=2.60, height=1.82)

volcano_1_4 <- ggplot(pairwise_res %>% filter(comparison %in% c(1,4)),
                  aes(x=log2FoldChange, y=-log10(padj), fill=significance)) +
  geom_vline(xintercept = 0.5,linetype="dashed", color="grey") +
  geom_vline(xintercept = -0.5,linetype="dashed", color="grey") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey") +
  geom_point(color="black", shape=21, alpha=.8, size=2,stroke=.25) +
  scale_fill_manual(values = c("red", "black")) +
  theme_nature() +
  ggh4x::facet_wrap2(~comparison,axes = "all", ncol = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title  = element_text(margin = margin(b=1)),
        legend.text = element_text(margin = margin(l=1)),
        legend.title.position = "top",
        legend.position = "top",
        legend.key.spacing.y = unit(1, "pt"),
        legend.margin = margin(b=5))
volcano_1_4

ggsave("figures/host_deseq/volcano_plot_1_4.pdf", units="in", width=1.35, height=2.75)
write_tsv(pairwise_res, file = "tables/deseq_pneumotypes_results.tsv")

ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)

###Visual ----
library(pheatmap)
library(viridis)
## ALL SAMPLES
mat <- assay(ntd) %>% 
  as_tibble(rownames='ENSEMBL') %>% 
  filter(ENSEMBL %in% res_sig$ENSEMBL) %>% 
  column_to_rownames('ENSEMBL')

df<- permmeta %>% 
  filter(sample_name %in% colnames(mat)) %>% 
  select(sample_name,dysbiotic, cluster_num_amplicon) %>% 
  mutate(dysbiotic = paste0('dysbiotic', dysbiotic)) %>% 
  column_to_rownames('sample_name') %>% 
  arrange(cluster_num_amplicon)

mat <- mat %>% select(rownames(df))

pheatmap::pheatmap(mat, cluster_rows=T, show_rownames=T,scale='none',
         cluster_cols=T, annotation_col = df, color = viridis(100),
         border_color = "black",fontsize_row = 2,
         show_colnames = F,fontsize = 5,
         filename = "figures/host_deseq/pheatmap_significant.pdf")

# Visualize differentially expressed genes across samples 
gene_palette = circlize::colorRamp2(c(1e-30,1e-20, .05, .1, .2, 6.7), c("black", rev(RColorBrewer::brewer.pal(5, "YlGnBu")))) #c("black", viridis::mako(5)))

mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
  filter(feature %in% res_sig$ENSEMBL) %>% 
  tidy_features() %>% 
  left_join(permmeta) %>% 
  drop_na(cluster_num_amplicon) %>% 
  group_by(cluster_num_amplicon) %>% 
  
  tidyHeatmap::heatmap(
    .row=feature,
    .column = sample_name,
    palette_value=gene_palette,
    .value = value,
    scale = "row"
  )

source("src/common/graph_features.r")
source("src/common/load_lists.r")
top_de_cluster <- graph_features(mox_feature_list[["RNA [Host Transcriptomics]"]],
               feature_array = res_sig %>% filter(idx < 10) %>% distinct(ENSEMBL) %>% pluck("ENSEMBL"),
               baseline_group = "1",
               groups_df = permmeta %>% drop_na(dysbiotic),
               group_col_name = "cluster_num_amplicon",feature_col_name = "feature",
               # feature_col_name = "feature",
               recode_dash = "-",clean_feature = T) +
  xlab("Log2 Expression")

top_de_cluster
ggsave("figures/host_deseq/top10_DE_host_genes_pneumotypes.pdf", units="in", width=2.61, height=4.14)

# Save for later
pairwise_res_pneumotype <- pairwise_res



##Dysbiotic ----
## SUMMARY 
# DYSBIOTIC
dds <- DESeqDataSetFromMatrix(countData=htx_deseq, 
                              colData=exp_matrix_pneumo,
                              design= ~  batch + dysbiosis_score)
dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds,name = "dysbiosis_score")

summary(res)

pairwise_res <- res %>% 
  as_tibble(rownames='ENSEMBL') %>% 
  mutate(comparison = "MDNP Score")


res_sig <- pairwise_res %>% 
  arrange(pvalue, log2FoldChange) %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 1.5)

pairwise_res <- pairwise_res %>% 
  mutate(significance =
           case_when(padj < 0.05 & abs(log2FoldChange) > 1.5 ~ "|log2FC| > 1.5 &\nFDR P < 0.05", 
                     .default= "not significant")
         ) %>%
  arrange(desc(significance), padj)

head(pairwise_res, 20)

###Visual ----

mdnp_volcano <- ggplot(pairwise_res %>% arrange(significance, padj),
       aes(x=log2FoldChange, y=-log(padj), fill=significance)) +
  geom_vline(xintercept = 1.5,linetype="dashed", color="grey") +
  geom_vline(xintercept = -1.5,linetype="dashed", color="grey") +
  geom_hline(yintercept = -log(0.05), linetype="dashed", color="grey") +
  geom_point(color="black", shape=21, alpha=.9, size=1,stroke=.25) +
  scale_fill_manual(values = c("red", "black")) +
  theme_nature() 

mdnp_volcano

ggsave("figures/host_deseq/volcano_plot_mdnp.pdf",mdnp_volcano,  units="in", width=2.56, height=1.04)
write_tsv(pairwise_res, file = "tables/deseq_mdnp_results.tsv")

# Plot top dysbiotic-associated genes @ Sample-level resolution
gene_palette = circlize::colorRamp2(c(0.0,1e-20, .05, .1, .2, 6.7), c("black", rev(RColorBrewer::brewer.pal(5, "YlGnBu")))) #c("black", viridis::mako(5)))

dysbiotic_heatmap <- mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
  filter(feature %in% head(res_sig$ENSEMBL,20)) %>% 
  tidy_features() %>% 
  left_join(permmeta) %>% 
  mutate(feature = str_remove(feature, ".*: "),
         feature = tolower(feature)) %>% 
  drop_na(cluster_num_amplicon) %>% 
  arrange(dysbiotic,desc(dysbiosis_score)) %>% 
  group_by(dysbiotic) %>% 
  arrange(dysbiotic, dysbiosis_score) %>% 
  dplyr::rename(MDNP = "dysbiosis_score") %>% 
  
  tidyHeatmap::heatmap(
    .row=feature,
    .column = sample_name,
    .value = value,
        scale = "none",
        cluster_rows = T, 
        cluster_columns = T,
        na_col = "black",
        rect_gp = grid::gpar(col = "black", lwd = .5),
        border = T,
        show_column_names = F,
        column_title = "Samples",
        row_title = "Features",
        column_title_gp = grid::gpar(fontsize = 7),
        row_title_gp = grid::gpar(fontsize = 7),
        row_names_gp = grid::gpar(fontsize = 4, fontface="italic"),
        show_column_dend = F,
        show_row_dend = F,
        palette_value = gene_palette,
        
        palette_grouping = list(
          # For second grouping (property_group)
          c("#4C413FFF", "#4C413FFF"),
          
          # c("#b58b4c", "#74a6aa"),
          # For first grouping (vs)
          c("#4C413FFF", "#4C413FFF", "#4C413FFF","#4C413FFF")
          
          
        ),
        heatmap_legend_param = list(
          title="Expr.",
          border ="black",
          labels_gp = grid::gpar(fontsize = 4),
          legend_width = unit(3, "cm"),
          title_gp = grid::gpar(fontface="bold", fontsize = 7)
        )

        ) %>% 
  tidyHeatmap::annotation_tile(MDNP, 
                  palette = 
                    circlize::colorRamp2(
                      seq(0,1, by=1/8),
                      RColorBrewer::brewer.pal(9, "Greys")),
                  border=T, 
                  size = unit(.3,"cm"),
                  annotation_legend_param = list(
                    title="MDNP",
                    border ="black",
                    legend_width = unit(3, "cm"),
                    labels_gp = grid::gpar(fontsize = 4),
                    title_gp = grid::gpar(fontface="bold", fontsize = 7)
                  ),
                  annotation_name_gp= grid::gpar(fontface="bold", fontsize = 6))
dysbiotic_heatmap
tidyHeatmap::save_pdf(dysbiotic_heatmap,filename="figures/host_deseq/de_host_genes_dysbiosis_score_heatmap.pdf", units = "in", 
                      height=1.73, width=7.34)

# Plot top features by dysbiotic vs not dybiotic
# source("src/common/graph_features.r")
# source("src/common/load_lists.r")
zlr_dysbiotic <- graph_features(mox_feature_list[["RNA [Host Transcriptomics]"]],
               feature_array = head(res_sig$ENSEMBL,15),
               baseline_group = "FALSE",
               groups_df = permmeta %>% drop_na(dysbiotic),
               group_col_name = "dysbiotic",feature_col_name = "feature",
               # feature_col_name = "feature",
               recode_dash = "-",clean_feature = F)
zlr_dysbiotic 
ggsave("figures/host_deseq/top10_DE_host_genes_dysbiosis_score.pdf", units="in", width=2.17, height=3.27)

# REGRESSIONS
# Plot subset of samples for correlation

p_ninj <- permmeta %>% 
  left_join(mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
              tidy_features() %>% 
              filter(feature =="ENSG00000131669: NINJ1")) %>% 
  ggplot(aes(dysbiosis_score, value)) + 
  geom_point(fill="grey", color="black", shape=21,alpha = .8, size=1) + 
  geom_smooth(color="black",method="lm",linewidth = .3) +
  theme_nature() +
  ylab("NINJ1\n(Log2 Expression)") +
  xlab("MDNP Score")
p_ninj

ggsave("figures/host_deseq/feature_plot_mdnp_ninj.pdf", units="in", width=1.59, height = 1.36)

  
p_pfkb3 <- permmeta %>% 
  left_join(mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
              tidy_features() %>% 
              filter(feature =="ENSG00000170525: PFKFB3")) %>% 
  ggplot(aes(dysbiosis_score, value)) + 
  geom_point(fill="grey", color="black", shape=21,alpha = .8, size=1) + 
  geom_smooth(color="black",method="lm",linewidth = .3) +
  theme_nature() +
  ylab("PFKFB3\n(Log2 Expression)") +
  xlab("MDNP Score")
p_pfkb3

ggsave("figures/host_deseq/feature_plot_mdnp_pfkb3.pdf", units="in", width=1.59, height = 1.36)


p_il1b <- permmeta %>% 
  left_join(mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
              tidy_features() %>% 
              filter(feature =="ENSG00000125538: IL1B")) %>% 
  ggplot(aes(dysbiosis_score, value)) + 
  geom_point(fill="grey", color="black", shape=21,alpha = .8, size=1) + 
  geom_smooth(color="black",method="lm",linewidth = .3) +
  theme_nature() +
  ylab("IL1B\n(Log2 Expression)") +
  xlab("MDNP Score")
p_il1b
ggsave("figures/host_deseq/feature_plot_mdnp_il1b.pdf", units="in", width=1.59, height = 1.36)

p_nfkb2 <- permmeta %>% 
  left_join(mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
              tidy_features() %>% 
              filter(feature =="ENSG00000077150: NFKB2")) %>% 
  ggplot(aes(dysbiosis_score, value)) + 
  geom_point(fill="grey", color="black", shape=21,alpha = .8, size=1) + 
  geom_smooth(color="black",method="lm",linewidth = .3) +
  theme_nature() +
  ylab("NFKB2\n(Log2 Expression)") +
  xlab("MDNP Score")
p_nfkb2
ggsave("figures/host_deseq/feature_plot_mdnp_nfkb2.pdf", units="in", width=1.59, height = 1.36)


###########
#Enrichment Analysis ----
## MDNP ----
library(clusterProfiler)
library(org.Hs.eg.db)

gene_dyb <- pairwise_res %>% 
  filter(significance != "not significant") %>% 
  filter(log2FoldChange > 0) %>% 
  mutate(ENSEMBL= str_remove(ENSEMBL, ":.*")) %>% 
  pluck("ENSEMBL")
universe <- mox_feature_list[["RNA [Host Transcriptomics]"]] %>% 
  mutate(feature= str_remove(feature, ":.*")) %>% 
  pluck("feature")

set.seed(100)

ego <- enrichGO(gene          = gene_dyb,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = "ENSEMBL",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write_tsv(data.frame(ego), "tables/deseq_enriched_mdnp.tsv")

dotplot(ego, showCategory=15) + 
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno")

data.frame(ego) %>% names()

enriched_dyb <- ggplot(data.frame(ego) %>% 
                         filter(Count >= 5, qvalue<.05) %>% 
                         slice_max(FoldEnrichment,n = 15) %>%
                         arrange(FoldEnrichment) %>%
                         mutate(Description = str_wrap(Description,width = 40)) %>% 
                         mutate(Description = factor(Description, levels=Description)), 
       aes(x=FoldEnrichment, y=Description, size=Count,fill=-log10(qvalue))) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25),
                       limits = c(-log10(1), max(-log10(data.frame(ego)$qvalue)*1.05))) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.y = element_text(lineheight = .75)) +
  labs(subtitle = "MDNP score enriched biological processes") + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.6),
        panel.grid.major.y = element_line(color="grey"),
        panel.grid.minor.y = element_line(color="grey"),
        panel.grid.major.x = element_line(color="grey"),
        strip.text = element_text(face="bold",hjust = 0),
        strip.clip = "off",
        axis.line.x = element_blank(),
        axis.line.y = element_blank())


enriched_dyb

ggsave("figures/host_deseq/enriched_bp_mdnp_score.pdf", width=4.28, height=2.38)



enriched_dyb_wide <- ggplot(data.frame(ego) %>% 
                              filter(Count >= 5, qvalue<.05) %>% 
                              slice_max(FoldEnrichment,n = 15) %>% 
                              arrange(FoldEnrichment) %>%
                              mutate(Description = factor(Description, levels=Description)), 
                            aes(y=FoldEnrichment, x=Description, size=Count,fill=-log10(qvalue))) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25),
                       limits = c(-log10(1), max(-log10(data.frame(ego)$qvalue)*1.05))) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 335, hjust = 0))+
  theme(legend.position = "left", 
        legend.position.inside = c(0.6,.4),
        legend.box = "vertical",
        legend.title = element_text(margin = margin(b=1)),
        legend.text = element_text(margin = margin(l=1)),
        legend.spacing.y = unit(3, "pt"),
        axis.text.y = element_text(lineheight = .6,size = 5),
        plot.subtitle = element_text(face="bold",hjust=1,margin=margin(2,0,0,0)))  + 
  scale_x_discrete(position = "bottom") + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.6),
        panel.grid.major.y = element_line(color="grey"),
        panel.grid.minor.y = element_line(color="grey"),
        panel.grid.major.x = element_line(color="grey"),
        strip.text = element_text(face="bold",hjust = 0),
        # axis.text.y = element_text(lineheight = .75),
        strip.clip = "off",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right")

enriched_dyb_wide
ggsave("figures/host_deseq/enriched_bp_mdnp_wide.pdf", width=3.18, height=1.85)

## SL Pneumotype ----
###
# Not for pneumotype
#### 
# Cluster 4 Positive
gene4 <- pairwise_res_pneumotype %>% 
  filter(significance != "not significant") %>% 
  filter(log2FoldChange > 0) %>% 
  filter(comparison=="4") %>% 
  mutate(ENSEMBL= str_remove(ENSEMBL, ":.*")) %>% 
  pluck("ENSEMBL")

ego4 <- enrichGO(gene         = gene4,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = "ENSEMBL",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write_tsv(data.frame(ego4), "tables/deseq_enriched_bp_oral_like.tsv")

enriched_4 <- ggplot(data.frame(ego4) %>% 
                       filter(Count >= 5, qvalue<.05) %>% 
                       slice_max(FoldEnrichment,n = 15) %>% 
                       arrange(FoldEnrichment) %>%
                       mutate(Description = str_wrap(Description,width = 40)) %>% 
                       mutate(Description = factor(Description, levels=Description)), 
                     aes(x=FoldEnrichment, y=Description, size=Count,fill=-log10(qvalue))) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25)) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.y = element_text(lineheight = .75))+
  labs(subtitle= "PnOL enriched biological processes") + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.6),
        panel.grid.major.y = element_line(color="grey"),
        panel.grid.minor.y = element_line(color="grey"),
        panel.grid.major.x = element_line(color="grey"),
        strip.text = element_text(face="bold",hjust = 0),
        strip.clip = "off",
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

enriched_4
ggsave("figures/host_deseq/enriched_bp_oral_like.pdf", width=4.19, height=2.11)


# Cluster 4 negative is too small, only 2 de

## OL Pneumotype ----

# Cluster 4 Positive
gene1 <- pairwise_res_pneumotype %>% 
  filter(significance != "not significant") %>% 
  filter(log2FoldChange < 0) %>% 
  filter(comparison=="1") %>% 
  mutate(ENSEMBL= str_remove(ENSEMBL, ":.*")) %>% 
  pluck("ENSEMBL")

ego1 <- enrichGO(gene         = gene1,
                 universe      = universe,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 keyType       = "ENSEMBL",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

write_tsv(data.frame(ego1), "tables/deseq_enriched_bp_skin_like.tsv")

enriched_1 <- ggplot(data.frame(ego1) %>% 
                       filter(Count >= 5, qvalue<.05) %>% 
                       # slice_min(qvalue,n = 15) %>% 
                       slice_max(FoldEnrichment,n = 15) %>% 
                       arrange(FoldEnrichment) %>%
                       mutate(Description = str_wrap(Description,width = 40)) %>% 
                       mutate(Description = factor(Description, levels=Description)), 
                     aes(x=FoldEnrichment, y=Description, size=Count,fill=-log10(qvalue))) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25)) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.y = element_text(lineheight = .75))+
  labs(subtitle= "PnSL negatively enriched biological processes")

enriched_1
ggsave("figures/host_deseq/enriched_bp_skin_like.pdf", width=4.28, height=2.38)

# Plot both on same scale
faceted <- rbind(data.frame(ego1) %>% 
                       filter(Count >= 5, p.adjust<.05) %>% 
                       slice_max(FoldEnrichment,n = 15) %>%
                       arrange(FoldEnrichment) %>%
                   mutate(Pneumotype = " Underexpressed in Skin-like"),
                 data.frame(ego4) %>% 
                   filter(Count >= 5, p.adjust<.05) %>% 
                   slice_max(FoldEnrichment,n = 15) %>%
                   arrange(FoldEnrichment) %>%
                   mutate(Pneumotype = "Overexpressed in Oral-like")) %>% 
  arrange(Pneumotype, FoldEnrichment) %>% 
  mutate(Description = factor(Description, levels=unique(Description))) %>% 
  arrange(Description)
  


enriched_1_4 <- ggplot(faceted, 
       aes(x=FoldEnrichment, y=Description, size=Count,fill=-log10(qvalue))
       ) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25),
                       limits = c(-log10(1), max(-log10(faceted$qvalue)*1.05))) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.y = element_text(lineheight = .75))+
  ggh4x::facet_wrap2(~Pneumotype,ncol = 1,scale = "free_y",axes="all")+
  theme(legend.position = "left", 
        legend.position.inside = c(0.6,.4),
        legend.box = "vertical",
        legend.title = element_text(margin = margin(b=1)),
        legend.text = element_text(margin = margin(l=1)),
        legend.spacing.y = unit(3, "pt"),
        axis.text.y = element_text(lineheight = .6,size = 5),
        plot.subtitle = element_text(face="bold",hjust=1,margin=margin(2,0,0,0)))  + 
  scale_y_discrete(position = "right") + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.6),
        panel.grid.major.y = element_line(color="grey"),
        panel.grid.minor.y = element_line(color="grey"),
        panel.grid.major.x = element_line(color="grey"),
        strip.text = element_text(face="bold",hjust = 0),
        # axis.text.y = element_text(lineheight = .75),
        strip.clip = "off",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right")

enriched_1_4
ggsave("figures/host_deseq/enriched_bp_both_pneumotypes.pdf", width=6.3, height=2.45)

enriched_1_4_wide <- ggplot(faceted, 
                            aes(x=Description, y=FoldEnrichment,  size=Count,fill=-log10(qvalue))
) +
  geom_point(shape=21) +
  theme_nature() + 
  scale_fill_viridis_c(option = "inferno",
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black",
                                              frame.linewidth = 0.25),
                       limits = c(-log10(1), max(-log10(faceted$qvalue)*1.05))) +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 335, hjust = 0))+
  ggh4x::facet_wrap2(~Pneumotype,ncol = 2,scale = "free_x",axes="all")+
  theme(legend.position = "left", 
        legend.position.inside = c(0.6,.4),
        legend.box = "vertical",
        legend.title = element_text(margin = margin(b=1)),
        legend.text = element_text(margin = margin(l=1)),
        legend.spacing.y = unit(3, "pt"),
        axis.text.y = element_text(lineheight = .6,size = 5),
        plot.subtitle = element_text(face="bold",hjust=1,margin=margin(2,0,0,0)))  + 
  scale_x_discrete(position = "bottom") + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=.6),
        panel.grid.major.y = element_line(color="grey"),
        panel.grid.minor.y = element_line(color="grey"),
        panel.grid.major.x = element_line(color="grey"),
        strip.text = element_text(face="bold",hjust = 0),
        # axis.text.y = element_text(lineheight = .75),
        strip.clip = "off",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right")

enriched_1_4_wide
ggsave("figures/host_deseq/enriched_bp_both_pneumotypes_wide.pdf", width=5.59, height=2.10)

venn <- ggVennDiagram::ggVennDiagram(list("OL\n(Up)" = gene4,"SL\n(Down)" = gene1),
                      label="count",
                      label_alpha = 0,
                      label_size =  5/.pt,
                      set_size = 5/.pt,
                      edge_size = .6,
                      label_color = c("OL" = darken("#9fb3b6"),"SL" = darken("#e18256"),"Intersect" = "black"),
                      set_color = c("OL" = darken("#9fb3b6"),"SL" = darken("#e18256"))
) +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = TRUE,
                                title = "",
                                frame.linewidth=.4,
                                frame.colour="black",
                                barwidth = .6,
                                barheight = 2,
                                ticks.colour = "black")) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme(legend.position = "left") + 
  scale_x_continuous(expand = expansion(mult = .4))
venn
ggsave("figures/host_deseq/venn_diagram.pdf", units="in", width=2.60, height = 1.55)

#Patchwork ----
## MDNP ----

# Convert to ggplot
dysbiotic_heatmap_heat <- dysbiotic_heatmap
dysbiotic_heatmap <- 
  ComplexHeatmap::draw(dysbiotic_heatmap_heat %>% tidyHeatmap::as_ComplexHeatmap(),
                     heatmap_legend_side="right",show_heatmap_legend = FALSE,show_annotation_legend=FALSE,
                     annotation_legend_side="right",
                     merge_legend=F) %>% 
  grid::grid.grabExpr(.) %>% 
  as_ggplot()

library(patchwork)

design = '
AABBCC
AABBCC
AADDEE
AADDEE
FFFFFF
FFFFFF
FFFFFF
FFFFFF
GGGGGG
GGGGGG
'

free(mdnp_volcano + theme(  legend.position = "top",
                       legend.title.position = "top",
                       legend.box = "vertical",
                       legend.direction = "vertical",
                       # legend.position = "inside",
                          legend.position.inside = c(0,.8),
                          legend.title = element_text(margin = margin(b=1)),
                          legend.text = element_text(margin = margin(l=1))
                     )) +
  p_ninj + 
  p_pfkb3 + 
  p_il1b + 
  p_nfkb2 + 
  free(dysbiotic_heatmap + theme(plot.margin=unit(c(1,1,1,1), "pt")),
       type = "panel", side="l")+ 
  enriched_dyb_wide + theme(
                       legend.box = "horizontal",
                       legend.title = element_text(margin = margin(b=1)),
                       legend.text = element_text(margin = margin(l=1)),
                       legend.spacing.y = unit(3, "pt"),
                       axis.text.y = element_text(lineheight = .75),
                       plot.subtitle = element_text(face="bold",hjust=1,margin=margin(2,0,0,0))
  ) +
  plot_layout(design = design,guides="keep",axis_titles = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.margin = margin(3,3,3,3),
        plot.tag = element_text(size = 8,hjust=1)
        )

ggsave(filename="figures/host_deseq/figure_4_mdnp_differentials_long.pdf", 
       width =4.55,  height=5.45, units="in")
       # width= 10.30, height=3.48, units="in") #9.69

#####

dysbiotic_heatmap <- 
  ComplexHeatmap::draw(dysbiotic_heatmap_heat %>% tidyHeatmap::as_ComplexHeatmap(),
                       heatmap_legend_side="right",
                       annotation_legend_side="top",
                       merge_legend=T) %>% 
  grid::grid.grabExpr(.) %>% 
  as_ggplot()
design = '
ABBBB
ABBBB
CBBBB
CDEFG
CDEFG
'

free(mdnp_volcano + theme(legend.position = "inside",
                          legend.position.inside = c(0,.8),
                          legend.title = element_text(margin = margin(b=1)),
                          legend.text = element_text(margin = margin(l=1))
), 
type="space", side="l"
) +
  free(dysbiotic_heatmap + theme(plot.margin=unit(c(1,1,1,1), "pt")
  ),
  type="panel",side = "tbl"
  ) + 
  enriched_dyb + theme(legend.position = "inside", 
                       legend.position.inside = c(0.6,.4),
                       legend.box = "vertical",
                       legend.title = element_text(margin = margin(b=1)),
                       legend.text = element_text(margin = margin(l=1)),
                       legend.spacing.y = unit(3, "pt"),
                       axis.text.y = element_text(lineheight = .75),
                       plot.subtitle = element_text(face="bold",hjust=1,margin=margin(2,0,0,0))
  ) +
  p_ninj + 
  p_pfkb3 + 
  p_il1b + 
  p_nfkb2 + 
  plot_layout(design = design,guides="collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.margin = margin(0,8,0,0),
        plot.tag = element_text(size = 8,hjust=1)
  )
ggsave(filename="figures/host_deseq/figure_4_mdnp_differentials.pdf", 
       width =8.67,  height=3.11, units="in")

## Title 3 ----

# source("src/revisions/host_subpermanova_pneumotype_heatmap.r")

pn_permanova <- readRDS("objects/figure_objects/sub_permanova.rds")
pn_permanova <- pn_permanova +  
  scale_x_discrete(expand = c(0,0), 
                   labels= c("RNA [Host Transcriptomics]"="Gene Expression\nPERMANOVA")) +
  scale_y_discrete(expand = c(0,0), 
                   labels= c("is_cluster_1"="SL",
                             "is_cluster_2"="M",
                             "is_cluster_3"="SP",
                             "is_cluster_4"="OL"
                             ),
                   limits=rev) +
  
  theme(axis.text.x = element_text(angle=0,hjust=.5,size=7),
        axis.text.y = element_text(angle=0,hjust=.5,size=7),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        legend.text = element_text(margin = margin(l=1))
  )


design = '
#ABBDDDEE
CCCCDDDEE
'

free(pn_permanova,type = "label", side="l") +

  free(volcano_1_4 + theme(legend.position = "top",
                           legend.box = "vertical"),
       type = "label", side="b") + 
  ((guide_area() + venn + plot_layout(guides = "collect"))) +
  top_de_cluster + #theme(legend.position = "bottom")+
  enriched_1_4 +
  plot_layout(design = design,
              guides="collect"
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.margin = margin(2,5,2,2),
        plot.tag = element_text(size = 8,face = "bold")
  )
ggsave("figures/host_deseq/figure_5_cluster_deseq.pdf", units="in",width=8.22, height=3.97)

