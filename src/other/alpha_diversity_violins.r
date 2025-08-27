source("src/common/get_distances.r")
source("src/common/load_features.r")
source("src/common/load_metadata.r")
source("src/common/load_lists.r")
source("src/common/graph_violin.r")
source("src/common/drop_samples.r")
get_alpha()
drop_samples()

# SHANNON DIVERSITY
shannon_amp_cat <-permmeta %>% 
    drop_na(shannon_amp) %>%
    ggplot(., aes(x=Episode_category, y=shannon_amp, fill=Episode_category)) %>%   #x=outcome
    plot_violin(.) %>% plot_statistics(.) +
    scale_fill_manual(values=mox_color_lists[['Episode_category']])  +
    guides(fill = guide_legend(title = "Pneumonia category")) +
    # guides(fill = guide_legend(title = "Pneumonia category")) +
    ylab('Shannon diversity') +
    xlab('Pneumonia category')

ggsave(filename = "figures/alpha_diversity/shannon_amp_Episode_category.pdf", shannon_amp_cat, units = "in", width=1.41, height=1.15)

shannon_amp_out <-permmeta %>% 
    drop_na(shannon_amp, Episode_is_cured) %>%
    ggplot(., aes(x=Episode_is_cured, y=shannon_amp, fill=Episode_is_cured)) %>%   #x=outcome
    plot_violin(.) %>% plot_statistics(.) +
    scale_fill_manual(values=mox_color_lists[['Episode_is_cured']])  +
    guides(fill = guide_legend(title = "Pneumonia category")) +
    # guides(fill = guide_legend(title = "Pneumonia category")) +
    ylab('Shannon diversity') +
    xlab('Pneumonia outcome')

ggsave(filename = "figures/alpha_diversity/shannon_amp_Episode_is_cured.pdf", shannon_amp_out, units = "in", width=1.41, height=1.15)


permmeta <- permmeta %>% left_join(
    vegan::diversity(t(mox_feature_list[['AMP [Genus]']] %>% column_to_rownames('feature')), index='invsimpson') %>% as_tibble(rownames = 'sample_name') %>% dplyr::rename(invsimpson_amp = 'value')
)

###################
# INVERSE SIMPSON #
###################

isimpson_amp_cat <-permmeta %>% 
    drop_na(invsimpson_amp) %>%
    ggplot(., aes(x=Episode_category, y=invsimpson_amp, fill=Episode_category)) %>%   #x=outcome
    plot_violin(.) %>% plot_statistics(.) +
    scale_fill_manual(values=mox_color_lists[['Episode_category']])  +
    guides(fill = guide_legend(title = "Pneumonia category")) +
    # guides(fill = guide_legend(title = "Pneumonia category")) +
    ylab('Inv. Simpson') +
    xlab('Pneumonia category')

ggsave(filename = "figures/alpha_diversity/invsimpson_amp_Episode_category.pdf", isimpson_amp_cat, units = "in", width=1.41, height=1.15)


isimpson_amp_out <-permmeta %>% 
    drop_na(invsimpson_amp, Episode_is_cured) %>%
    ggplot(., aes(x=Episode_is_cured, y=invsimpson_amp, fill=Episode_is_cured)) %>%   #x=outcome
    plot_violin(.) %>% plot_statistics(.) +
    scale_fill_manual(values=mox_color_lists[['Episode_is_cured']])  +
    guides(fill = guide_legend(title = "Pneumonia category")) +
    # guides(fill = guide_legend(title = "Pneumonia category")) +
    ylab('Inv. Simpson') +
    xlab('Pneumonia outcome')

ggsave(filename = "figures/alpha_diversity/invsimpson_amp_Episode_is_cured.pdf", isimpson_amp_cat, units = "in", width=1.41, height=1.15)
