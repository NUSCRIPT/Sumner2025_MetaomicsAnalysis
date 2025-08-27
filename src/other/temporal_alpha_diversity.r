source("src/common/load_metadata.r")
source("src/common/load_features.r")
source("src/common/load_lists.r")

get_alpha()
permmeta <- permmeta %>%
    mutate(is_longitudinal = 
               case_when(script_id %in% (permmeta %>% 
                                             filter(initial_bal_sample == 0) %>% 
                                             distinct(script_id) %>% 
                                             pluck("script_id")
                                         ) ~ "longitudinal", 
                         .default="baseline_only")
           )

permmeta %>% drop_na(shannon_amp) %>% filter(is_longitudinal=="longitudinal") %>% 
    ggplot(aes(x=sqrt(days_since_bal), y=shannon_amp, group=script_id)) +
    geom_line(alpha=.5,linewidth=.25) +
    
    geom_point(shape =21, fill="grey",size=1,stroke=.5,alpha=.9) + 
    theme_nature() +
    labs(x = "Square root of days since first BAL", y = "Shannon diversity")

ggsave("figures/alpha_diversity/timesries_shannon_amp_raw.pdf", units="in", height= 1.04, width =1.83)



permmeta %>% drop_na(shannon_amp) %>% filter(is_longitudinal=="longitudinal") %>% 
    arrange(script_id, days_since_bal) %>% 
    # select(sample_name, script_id, shannon_amp, days_since_bal) %>%
    group_by(script_id) %>% 
    mutate(shannon_amp = shannon_amp - first(shannon_amp)) %>%
    # mutate(shannon_amp = scale(shannon_amp)) %>% 
    mutate(shannon_direction = case_when(last(shannon_amp) > first(shannon_amp) ~ "diversity increases",
                                         last(shannon_amp) < first(shannon_amp) ~ "diversity decreases"),
           shannon_direction = factor(shannon_direction, levels = c("diversity increases", "diversity decreases"))) %>% 
    ungroup() %>% 
    arrange(shannon_direction) %>% 
    drop_na(shannon_direction) %>%
    mutate(script_id = factor(script_id, levels = unique(script_id))) %>% 
    
    ggplot(aes(x=sqrt(days_since_bal), y=shannon_amp, group=script_id, fill=shannon_direction)) +
        geom_line(color="black") +
        geom_point(shape =21,size=2,stroke=1,alpha=.9) + 
        theme_nature() +
        ggh4x::facet_wrap2(~Episode_is_cured, axes="all",scales = "fixed") +
        scale_fill_brewer(palette="Set2") +
        labs(x = "Square root of days since first BAL", y = "Shannon diversity\n(H'i-H'1)")

ggsave("figures/alpha_diversity/timesries_shannon_amp_episode_outocome.pdf", units="in", height= 1.29, width =4.23)


permmeta %>% drop_na(shannon_amp) %>% filter(is_longitudinal=="longitudinal") %>% 
    arrange(script_id, days_since_bal) %>% 
    # select(sample_name, script_id, shannon_amp, days_since_bal) %>%
    group_by(script_id) %>% 
    mutate(shannon_amp = shannon_amp - first(shannon_amp)) %>%
    # mutate(shannon_amp = scale(shannon_amp)) %>% 
    mutate(shannon_direction = case_when(last(shannon_amp) > first(shannon_amp) ~ "diversity increases",
                                         last(shannon_amp) < first(shannon_amp) ~ "diversity decreases"),
           shannon_direction = factor(shannon_direction, levels = c("diversity increases", "diversity decreases"))) %>% 
    ungroup() %>% 
    arrange(shannon_direction) %>% 
    drop_na(shannon_direction) %>%
    mutate(script_id = factor(script_id, levels = unique(script_id))) %>% 
    
    ggplot(aes(x=sqrt(days_since_bal), y=shannon_amp, group=script_id, fill=shannon_direction)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_line(color="black") +
    geom_point(shape =21,size=1,stroke=1,alpha=.9) + 
    theme_nature() +
    ggh4x::facet_wrap2(~script_id, axes="all",scales = "fixed") +
    scale_fill_brewer(palette="Set2") +
    labs(x = "Square root of days since first BAL", y = "Shannon diversity\n(H'i-H'1)")

ggsave("figures/alpha_diversity/timesries_shannon_amp_individuals.pdf", units="in", width= 5.29, height =4.52)

permmeta %>% drop_na(shannon_amp) %>% filter(is_longitudinal=="longitudinal") %>% 
    arrange(script_id, days_since_bal) %>% 
    # select(sample_name, script_id, shannon_amp, days_since_bal) %>%
    group_by(script_id) %>% 
    mutate(shannon_amp = shannon_amp - first(shannon_amp)) %>%
    # mutate(shannon_amp = scale(shannon_amp)) %>% 
    mutate(shannon_direction = case_when(last(shannon_amp) > first(shannon_amp) ~ "diversity increases",
                                         last(shannon_amp) < first(shannon_amp) ~ "diversity decreases"),
           shannon_direction = factor(shannon_direction, levels = c("diversity increases", "diversity decreases"))) %>% 
    ungroup() %>% 
    arrange(shannon_direction) %>% 
    drop_na(shannon_direction) %>%
    mutate(script_id = factor(script_id, levels = unique(script_id))) %>% 
    
    ggplot(aes(x=sqrt(days_since_bal), y=shannon_amp, group=script_id, fill=cluster_num_amplicon)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_line(color="black") +
    geom_point(shape =21,size=1.5,stroke=.5) + 
    theme_nature() +
    ggh4x::facet_wrap2(~script_id, axes="all",scales = "fixed") +
    scale_fill_manual(values = mox_color_lists[["cluster_num_amplicon"]]) +

    labs(x = "Square root of days since first BAL", y = "Shannon diversity\n(H'i-H'1)")
ggsave("figures/alpha_diversity/timesries_shannon_amp_individuals_pneumotype.pdf", units="in", width= 5.29, height =4.52)
    
