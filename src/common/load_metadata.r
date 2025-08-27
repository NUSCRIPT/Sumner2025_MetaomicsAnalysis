library(tidyverse)
library(janitor)

get_new_metadata <- function(){
  
  # Load metadata base
  meta <- read_csv('~/share_with_Jack/CarpeDiem_dataset.csv')
  map <- read_csv('~/share_with_Jack/201022_script_anonymized_ids_master.csv') %>% mutate(script_id = as.double(script_id), anonymized_id = as.double(anonymized_id))
  dict <- read_csv('~/share_with_Jack/CarpeDiem_data_dictionary.csv')
  report <- read_csv('~/share_with_Jack/20240528_script_all_events_mapped_undated.csv') %>%  
    mutate(bal_barcode=str_remove_all(bal_barcode, '-')) %>% 
    mutate_at(vars(contains(c('_quantity'))), ~as.numeric(str_remove_all(.,',|>')))  # Fix first two col that are characters
  
  # Subset report to relevant nums
  report_organism <<-  report %>% 
    distinct(bal_barcode, pick(contains("organism"))) %>% 
    pivot_longer(contains('organism'),names_to = c("Var", ".value"),names_sep=11) %>% 
    distinct(bal_barcode, name, quantity) %>% 
    mutate(name = case_when(is.na(name) ~ 'organism_null', .default = name),
           quantity = case_when(is.na(quantity) ~ 0, .default = quantity)) %>% 
    # quantity = as.double(str_remove_all(quantity, ">|,"))) %>% 
    group_by(bal_barcode) %>% 
    mutate(culture_positive = case_when(any(quantity > 0) ~ 'Positive', .default = 'Negative'),
           fungal_positive = case_when(any(str_detect(name, 'Yeast|Candida')) ~ 'Positive', .default = 'Negative'),
           bacteria_positive = case_when(any(str_detect(name, 'Yeast|Candida|organism_null',negate = T)) ~ 'Positive', .default = 'Negative')) %>% 
    ungroup() %>% 
    mutate(name = str_remove_all(name, ',|\\)|\\('),
           name = str_replace_all(name, ' ', '_')) %>% 
    add_count(name) %>% # arrange by prevalence of culture result
    arrange(desc(n)) %>% #distinct(name, n) %>% View()
    select(!n) %>% 
    pivot_wider(names_from = 'name', values_from = 'quantity', values_fill = 0)
  report_organism_minimal <- report_organism %>% select(1:4)
  # Note: the organisms are quite separated (e.g., multiple "Viridans streptococci" and Staph groups)
  
  report_bal <- report %>% 
    distinct(pt_study_id, bal_barcode, hospital_los_days,BAL_day_after_hos_admission, 
             culture_flag, `amylase bf`, antibiotics_flag, total_episodes, pna_category, 
             assess_stday, episode_rank, icu_rank, ICU_Day) %>% 
    dplyr::rename(amylase_bf = `amylase bf`,
                  sample_name = bal_barcode) %>% 
    mutate(amylase_bf = as.numeric(str_replace(amylase_bf, "<10", "5")), # half of presumptive LoQ
           amylase_bf_log = log10(amylase_bf), 
           pna_category = case_when(str_detect(pna_category, 'HAP') ~ 'HAP',
                                    str_detect(pna_category, 'VAP') ~ 'VAP',
                                    str_detect(pna_category, 'CAP') ~ 'CAP',
                                    str_detect(pna_category, 'Non-pneumonia control') ~ 'NPC',
                                    .default = pna_category),
           .after = culture_flag)
  # Note: has one entry per bal barcode (n=375)
  
  
  # Join and parse new md by
  my_new_meta <- meta %>% 
    filter(has_bal == T) %>%
    fill( Episode_category, .direction = 'down') %>% # Expand categories, must early in case don't have initial bal after join
    fill( Episode_etiology, .direction = 'down' ) %>% 
    fill( Episode_is_cured, .direction = 'down' ) %>% 
    fill( Episode_duration, .direction = 'down' ) %>% 
    left_join(map, by=c('Patient_id' = 'anonymized_id'), unmatched = 'drop', relationship='many-to-one') %>%
    relocate(script_id, has_bal) %>% 
    inner_join(report_bal, by=c('script_id' = 'pt_study_id', 
                                'ICU_stay' = 'icu_rank', 
                                'ICU_day' = 'ICU_Day'), 
               relationship = 'one-to-one') %>% 
    # pluck('sample_name')
    # filter(.data = report_bal, !sample_name %in% my_new_meta)
    relocate(sample_name)  %>% 
    arrange(sample_name) %>% 
    group_by(script_id) %>%
    relocate(contains('episode'), .after=has_bal) %>% 
    mutate(Episode_category = case_when(Episode_category == 'Non-PNA-ctrl' ~ 'NPC', .default=Episode_category), 
           Episode_etiology = case_when(Episode_category == 'NPC' ~ 'NPC', .default=Episode_etiology),
           Episode_is_cured = case_when(Episode_category == 'NPC' ~ 'NPC', .default=Episode_is_cured)) %>% 
    fill( episode_rank, .direction = 'down' ) %>% # Expand categories
    fill( pna_category, .direction = 'down' ) %>% 
    ungroup() %>%  
    mutate(Episode_category = case_when(is.na(Episode_category) ~ pna_category,
                                        .default = Episode_category)) %>% 
    left_join(
      # Subset to dereplicate means
      readRDS('src/common/QPCR_02_CleamqPCRTable.rds') %>% 
        mutate(bal.id = str_replace(bal.id, "-BAL-", "BAL"),
               Quantity.Mean = case_when((Quantity.Mean < 10 & Quantity.Mean > 0) ~ 5, # LoQ adjustment
                                         T ~ Quantity.Mean )
        ) %>% 
        dplyr::rename(sample_name ="bal.id") %>%
        select(sample_name, Quantity.Mean, Ct.Mean) %>% 
        unique() %>% 
        drop_na(Quantity.Mean) # tmp fix bug from qpcr sample n=1sample redacted, fixed now but code here
    ) %>% 
    mutate(Quantity.Mean.Log = log10(Quantity.Mean),
           days_since_bal = as.double(str_split_i(sample_name, 'BAL', 2)),
           initial_bal_sample = case_when(days_since_bal == 0 ~ 1, days_since_bal != 0 ~ 0),
           baseline_or_bal = case_when(str_ends(sample_name, 'BAL00') ~ "Baseline",.default = 'BAL'),
    ) %>% 
    group_by(script_id) %>% 
    mutate(bal_order = row_number()) %>% 
    ungroup()
  
  
  get_dysbiosis_score <- function(my_distance, npc_names){
    dysbiosis_score <- my_distance %>%
      as.matrix() %>%
      as_tibble(rownames = "sample_name") %>%
      pivot_longer(2:last_col()) %>%
      filter(name %in% npc_names) %>%
      group_by(sample_name) %>%
      summarize(dysbiosis_score = mean(value)) %>%
      filter(str_detect(sample_name, "BRL", negate = TRUE)) # Remove the one tech dup
    
    # Threshold at 90th quartile of NPC 
    dysbiosis_cutoff <<- dysbiosis_score %>% 
      filter(sample_name %in% npc_names) %>% 
      summarize(dysbiosis_cutoff = as.double(quantile(dysbiosis_score, probs=0.9, na.rm=TRUE))) %>%
      pluck("dysbiosis_cutoff")
    
    print(dysbiosis_cutoff)
    
    # Score thresholds
    dysbiosis_score <- dysbiosis_score %>% 
      mutate(dysbiotic = case_when(dysbiosis_score >= dysbiosis_cutoff ~ T,
                                   dysbiosis_score < dysbiosis_cutoff ~ F), 
             .before = 3)
    return(dysbiosis_score)
  }
  source('src/common/drop_samples.r')
  drop_list <- get_drop_lst()
  
  # vector list of only NPC that stay NPC to calc dysbiosis
  npc_names <- my_new_meta %>% 
    filter(Episode_category == "NPC") %>%     
    filter(!sample_name %in% drop_list) %>% 
    select(sample_name) %>% 
    pluck(1)
  source("src/common/get_distances.r")
  dysbiosis_dt <- get_dysbiosis_score(amp, npc_names)
  
  my_new_meta <- my_new_meta %>% 
    left_join(dysbiosis_dt) %>% 
    left_join(report_organism_minimal, by = c('sample_name' = 'bal_barcode'))
  
  return(my_new_meta)
}

# Execute
permmeta <- get_new_metadata() 
colnames(permmeta)
# Set factors/levels
permmeta <- permmeta %>% 
  mutate(Episode_category = factor(Episode_category, levels=c("HAP", "VAP","CAP", "NPC")),
         Binary_outcome = factor(Binary_outcome),
         Intubation_flag = factor(Intubation_flag),
         Hemodialysis_flag = factor(Hemodialysis_flag),
         Norepinephrine_flag = factor(Norepinephrine_flag),
         Episode_etiology = case_when(
           Episode_etiology == 'Indeterminate' ~ 'Culture-negative', # N=3 - checked w dmbi
           .default = Episode_etiology), 
         antibiotics_flag = factor(antibiotics_flag)
  )

source('src/common/get_clusters.r')

source('src/common/drop_samples.r')

neutrophil_report <- read_csv("~/20240628_script_all_events_mapped_neutrophils_vectors.csv") %>% 
  dplyr::rename(amylase_bf = `amylase bf`,
                sample_name = bal_barcode) %>% 
  mutate(amylase_bf = as.numeric(str_replace(amylase_bf, "<10", "5"))) %>%  
  distinct(BAL_day_after_hos_admission, sample_name, amylase_bf, neutrophils_bodyfluid)

permmeta <- permmeta %>% left_join(neutrophil_report, by = c("BAL_day_after_hos_admission","sample_name","amylase_bf"))
drop_samples()
