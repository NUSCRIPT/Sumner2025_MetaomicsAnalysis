library(tidyverse)
source("src/common/load_metadata.r")
source("src/common/load_features.r")
source("src/common/load_lists.r")
library(viridis)

ggplot(permmeta %>% drop_na(Episode_category, cluster_num_amplicon), 
       aes(x=log10(BAL_day_after_hos_admission+1), y=dysbiosis_score, group=cluster_num_amplicon,color=cluster_num_amplicon,fill=cluster_num_amplicon)) +
    geom_smooth(alpha=.5) +
    geom_point(shape=21) + 
    theme_nature()+
    # facet_wrap(~cluster_num_amplicon, nrow = 4)+
    scale_fill_manual(values=mox_color_lists[['cluster_num_amplicon']]) +theme(axis.text.x        = element_text(margin=margin(1, 0, 0, 0)))+
    scale_color_manual(values=mox_color_lists[['cluster_num_amplicon']]) +theme(axis.text.x        = element_text(margin=margin(1, 0, 0, 0))) +
    ylab("MDNP score")
ggsave("figures/dysbiosis/dysbiosis_score_temporal_hospital_stay_pneumotypes.pdf", units="in",width=5.30, height=2.41)
ggplot(permmeta, aes(x=baseline_or_bal, y=dysbiosis_score)) +
    geom_boxplot() +
    geom_jitter() +
    theme_nature() +
    ylab("MDNP score")

ggsave("figures/dysbiosis/dysbiosis_score_bal_baseline_boxplot.pdf", units="in",width=2.20, height=2.15)
