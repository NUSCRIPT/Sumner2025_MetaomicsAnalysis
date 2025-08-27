source("src/common/nature_theme.r")
source("src/common/metaphlan_processing_functions.r")
require(tidyverse, vegan, tidyHeatmap)

source("src/common/nature_theme.r")
source('src/common/load_metadata.r')
source('src/common/load_lists.r')

# FOR FINDING JUSTIFICATION FOR PREVALENCE THRESHOLDS
ds_fi <- "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
dss <- read_metaphlan(ds_fi)

# GET FEATURE x PREVALENCE AT N% ABUNDANCE
prev_distribution <- dss %>% 
    tidy_features() %>% 
    group_by(feature) %>% 
    summarise(n = n(),
              prev0 = sum(value>0),
              prev0.001 = sum(value>1e-3),
              prev0.01 = sum(value>1e-2),
              prev0.1 = sum(value>1e-1),
              prev1 = sum(value>1),
              prev10 = sum(value>10),
              prev50 = sum(value>50),
              prev99 = sum(value>99),
              prev100 = sum(value>100)) %>% 
    pivot_longer(3:last_col())

# number of features found in at least 1% abundance in 1 sample
prev_features <- prev_distribution %>% filter(name == "prev1", value > 1) %>% distinct(feature) %>% pluck('feature')

# HEATMAP
my_heat <- prev_distribution %>% 
    mutate(value = log10(value+1)) %>%
    mutate(name = factor(name, levels=unique(name)),
           value = value/n*100) %>% 
    # arrange(desc(value)) %>% 
    heatmap(
        .column = name,
        .row = feature,
        .value = value,
        scale = "none",
        # cluster_rows=FALSE,
        cluster_columns = FALSE
    )
my_heat
save_pdf(my_heat, filename = "figures/sanity_checks/heatmap_prevalence_abundance_counts.pdf",units="in", width=3.78, height=9.02)

# LINE CHART
prev_distribution %>% 
    mutate(value = value+1e-1,
           name = str_replace(name, "prev", ""),
           name = as.double(name)) %>% 
    ggplot(aes(x=name+1e-10, y=value, group=feature)) + 
    geom_point(alpha=.2) +
    geom_line(alpha=.7) +
    scale_y_log10() +
    theme_nature() +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x="Abundance Threshold", y = "Feature Prevalence", title = "Species Prevalence at Varied Abundance Thresholds")

ggsave("figures/sanity_checks/scatter_abundance_by_prevalence_in_mgx.pdf",units="in", width=3.78, height=2.25)


prev_distribution %>% 
    group_by(name) %>% 
    summarize(total = sum(value > 0)) %>% 
    mutate(name = str_replace(name, "prev", ""),
           name = as.double(name)) %>% 
    
        ggplot(aes(x=name+1e-4, y=total)) + 
        geom_point(alpha=.2) +
        geom_line(alpha=.7) +
        geom_text(aes(label=total), hjust=-.5, vjust = 0, size=3) +
        geom_vline(xintercept = 5e-4, linetype="dashed", color = "red") +
        scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(x="Abundance Threshold", 
             y = "Total Feature with Prevalence > 0", 
             title = "No. Species with Prevalence > 0 at Varied Abundance Thresholds") +
        ylim(0,800)+s
        theme_nature() 
ggsave("figures/sanity_checks/scatter_total_species_prevalence.pdf",units="in", width=3.97, height=2.3)

