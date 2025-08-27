# categorical_metadata <- c( 
#   #"pt_study_id", # Individual id
#   "Episode_category", "Episode_etiology", # Diagnosis 
#   "initial_bal_sample", "bal_method", # BAL methodology
#   "success", "cause_failure", "clin_impression_d78", # Clinical hallmarks, disease severity
#   "prior_noninvase_ventilation_rc_flag", "initial_antibiotic",'icu_readmission_flg',
#   "bacteria_positive", 'dysbiotic') # indep
# permmeta %>% select(!where(is.numeric)) %>% colnames() %>% dput()
categorical_metadata <- c("Episode_category", "Episode_etiology", 
  "Episode_is_cured",
  "Discharge_disposition", "Binary_outcome", "Global_cause_failure", 
  "Patient_category", "COVID_status", 
  "Smoking_status", "Intubation_flag", "Hemodialysis_flag", "Norepinephrine_flag", 
  "pna_category", "baseline_or_bal", "dysbiotic", "culture_positive", 
  "fungal_positive", "bacteria_positive", "cluster_num_amplicon"
  # "cluster_num_dna_taxonomy", "cluster_num_dna_ko"
  )

quantitative_metadata <- c("icu_rank", "icu_day", "p_f_ratio_points", "pf_ratio", "platelet_points", 
                           "platelet", "bilirubin_points", "bilirubin", "htn_points", "map", 
                           "htn_med_name", "eye_opening_score", "best_motor_response_score", 
                           "best_verbal_response_score", "gcs_points", "verbal_score_estimated_ind", 
                           "renal_points", "urine_output", "creatinine", "hd_or_crrt_flag", 
                           "hd_flag", "crrt_flag", "sofa", "intub_flag", "copd_emphyasema", 
                           "copd_bronchitis", "bronchitis_other", "asthma", "copd_other", 
                           "bronchiectasis", "other_chronic_pulmonary_disease", "copd", 
                           "fungal_positive", "bacteria_positive", "accession_numbers_id", 
                           "rbc_body_fluid", "wbc_body_fluid", "neutrophils_body_fluid", 
                           "absolute_neutrophils", "macrophage_bf", "monocyte_bf", "lymph_bf", 
                           "absolute_lymphocytes", "eosinophils_body_fluid", "absolute_eosinophils", "other_cells_body_fluid", "amylase_bf", "white_blood_cell_count", 
                           "c_reactive_protein", "ldh", "creatine_kinase_total", "procalcitonin", 
                           "ferritin", "troponin_i", "dysbiosis_score","Quantity.Mean")


permanova_metadata <- list( 
  # 'A: Individuals' = "script_id", # Individual id
  'D: Pneumonia Type' = "Episode_category",  # Diagnosis 
  'D: Pneumonia Etiology' = "Episode_etiology",
  'D: Pneumonia Resolution' = "Episode_is_cured", 
  'D: Binned Outcome' = "Binary_outcome",
  # 'BM: Baseline/Follow Up BAL' = "initial_bal_sample", # BAL methodology
  # 'BM: BAL Method' = "bal_method", 
  
  'CH: Discharge Disposition' = 'Discharge_disposition', # Clinical hallmarks, disease severity
  'CH: Amylase Level (Log10, U/L)' = "amylase_bf_log",
  'CH: FiO2' = 'FiO2',
  'CH: PEEP' = 'PEEP',
  'CH: Plateau_Pressure' = 'Plateau_Pressure',
  'CH: Lung_Compliance' = 'Lung_Compliance',
  # 'CH: SOFA Score' = "SOFA_score", 
  # 'CH: Bilirubin' = 'Bilirubin',
  # 'CH: Albumin' = 'Albumin',
  # 'CH: Lactic_acid' = 'Lactic_acid',
  # 'CH: Bicarbonate' = 'Bicarbonate',
  
  'PH: Antibiotic Use' = "antibiotics_flag",
  'PH: Hospital Length of Stay' = 'hospital_los_days',
  'PH: Days Since Admission' = 'BAL_day_after_hos_admission',
  'PH: Smoking Status' = 'Smoking_status',
  'PH: Race' = 'Race',
  'PH: Age' = 'Age',
  # 'PH: BMI' = 'BMI',
  # 'PH: Multiple ICU Visits' = 'icu_readmission_flg', # weak correlation
  # 'PH: Prior NIV' = "prior_noninvase_ventilation_rc_flag",  # PT history
  
  'CB: Respiratory Culture Results' = "bacteria_positive", # CB
  "CB: Fungal Culture Results" = 'fungal_positive',
  'CB: Neutrophil Level (%)' = "neutrophils_bodyfluid", 
  'CB: Hemoglobin' = 'Hemoglobin',
  'CB: Platelets' = 'Platelets',
  # 'CB: Macrophage Level (%)' = "macrophage_bf",
  # 'CB: Lymphoid Cell Level (%)' = 'Lymphocytes', 
  # 'CB: White Blood Cell Level (%)'= "WBC_count", # redundant
  
  'IB: 16S rRNA Gene Copies (Log10 µL)' = "Quantity.Mean.Log", # independent biomarkers
  'IB: MDNP Score' = 'dysbiosis_score',
  'IB: Pneumotype' = 'cluster_num_amplicon'
  # 'IB: Diversity' = 'invsimpson_amp', # redundant
) 

distance_of_interest <- c("AMP [Genus]", "AMP [ASV]", "DNA [Taxonomy]", "DNA [Viral]", 
                          "RNA [Taxonomy BKN]", "DNA [PFAM]", "RNA [PFAM]", 
                          "RNA [Host Transcriptomics]", "CFU [Culture]")


mox_color_lists <- list(
  'dysbiotic' = c('FALSE' = '#a9c3a4', #86C2DAFF 
                  'TRUE' = '#7ea679'), #06425AFF
  # 'pt_category' =c('NPC' = "#999999",
  #                  'CAP' = '#278B9AFF',
  #                  'HAP' = '#E75B64FF',
  #                  'VAP' = '#D8AF39FF'), #rgb(230, 133, 59, maxColorValue = 255)),
  'Episode_category' =c('NPC' = "#999999",
                   'CAP' = '#a2bbd4',
                   'HAP' = '#f0ba7d',
                   'VAP' = '#a9c3a4'), 
  'pt_category_dysbiotic' =c('NPC' = "#999999",
                             # 'NPCFALSE' = "#999999",
                             # 'NPCTRUE' = '#4C413FFF',
                             'CAPFALSE' = '#a2bbd4',
                             'CAPTRUE' = '#6a81a5',
                             'HAPFALSE' = '#f0ba7d',
                             'HAPTRUE' = '#e18256',
                             'VAPFALSE' = '#a9c3a4',
                             'VAPTRUE' = '#7ea679'
  ),
  # 'pt_category_dysbiotic' =c('NPC' = "#999999",
  #                            # 'NPCFALSE' = "#999999",
  #                            # 'NPCTRUE' = '#4C413FFF',
  #                            'CAPFALSE' = '#278B9AFF',
  #                            'CAPTRUE' = '#5A6F80FF',
  #                            'HAPFALSE' = '#E75B64FF',
  #                            'HAPTRUE' = '#DE7862FF',
  #                            'VAPFALSE' = '#D8AF39FF',
  #                            'VAPTRUE' = '#E8C4A2FF'
  #                            ),
  # 'pt_category_dysbiotic' =c('NPC' = "#999999",
  #                            # 'NPCFALSE' = "#999999",
  #                            # 'NPCTRUE' = '#4C413FFF',
  #                            'CAPFALSE' = '#278B9AFF',
  #                            'HAPFALSE' = '#E75B64FF',
  #                            'VAPFALSE' = '#D8AF39FF',
  #                            'CAPTRUE' = '#5A6F80FF',
  #                            'HAPTRUE' = '#DE7862FF',
  #                            'VAPTRUE' = '#E8C4A2FF'
  # ),
  'bacteria_positive' = c('Negative' = '#f0ba7d', #B50A2AFF
                          'Positive' = "#e18256"), #D98594FF
  'Episode_etiology' = c("NPC" = '#999999',
                         "Bacterial" = '#f0ba7d',
                         "Bacterial/viral" = '#f5e8d0',
                         "Viral" = '#e18256',
                         "Culture-negative" = '#7ea679',
                         "Indeterminate" = '#64897a'
                         ),
  'Episode_is_cured' = c('Cured' = '#f0ba7d',     # rgb(159, 159, 163, maxColorValue = 255), 
                'Indeterminate' = '#2d441e',              #rgb(255, 128, 0, maxColorValue = 255), 
                'Not cured' = '#a9c3a4',
                'NPC' = '#999999'),
  # 'cluster_num_amplicon' =  c('1' = "#F0D77BFF",
  #                             '2' = "#B4DAE5FF",
  #                             '3' = "#AE93BEFF",
  #                             '4' = "#5C5992FF",
  #                             '5' = "#403369FF"),
  'cluster_num_amplicon' =  c('1' = "#e18256",
                              '2' = "#f0ba7d",
                              '3' = "#64897a",
                              '4' = "#9fb3b6",
                              '5' = "#403369FF"),
  'success_bacteria_positive' =c(
                             '0FALSE' = "#999999",
                             '0TRUE' = '#4C413FFF',
                             '1FALSE' = '#278B9AFF',
                             '1TRUE' = '#5A6F80FF',
                             '2FALSE' = '#E75B64FF',
                             '2TRUE' = '#DE7862FF'
                             # 'VAPFALSE' = '#D8AF39FF',
                             # 'VAPTRUE' = '#E8C4A2FF'
  ),
  "antibiotics_flag" = c('0' = '#547b80', '1' = '#67ada9')
    
)


'#9fb3b6'

"#e18256"
'#f0ba7d'

'#7ea679'
'#a9c3a4'


'#6a81a5'
'#a2bbd4'


'#64897a'
'#9fb3b6'

'#547b80'
'#67ada9'

'#2d441e'
'#3A7345'

# library(hues)
# pal <- c('#64897a','#9fb3b6',"#e18256",'#f0ba7d','#7ea679','#a9c3a4','#2d441e','#3A7345','#64897a','#547b80','#67ada9','#6a81a5','#a2bbd4')
# swatch(pal)

# For renaming bugs in culture 
organisms <- list(
  "Achromobacter_species" = "Achromobacter_species",
  "Achromobacter_xylosoxidans" = "Achromobacter_xylosoxidans",
  "Acinetobacter_baumannii" = "Acinetobacter_baumannii",
  "Acinetobacter_ursingii" = "Acinetobacter_ursingii",
  "Beta_Hemolytic_Streptococci_Group_C" = "Streptococcus_species",
  "Beta_Hemolytic_Streptococci_Group_G" = "Streptococcus_species",
  "Burkholderia_cepacia_complex" = "Burkholderia_cepacia",
  "Burkholderia_cepacia_complex_presumptive" = "Burkholderia_cepacia",
  "Candida_albicans" = "Candida_albicans",
  "Candida_glabrata" = "Candida_glabrata",
  "Candida_tropicalis" = "Candida_tropicalis",
  "Citrobacter_freundii_group_ESBL_Positive_Note:_This_organism_produces_the_KPC_carbapenemase._Consultation_with_Infectious_Disease_Service_is_recommended." = "Citrobacter_freundii",
  "Citrobacter_koseri" = "Citrobacter_koseri",
  "Corynebacterium_species" = "Corynebacterium_species",
  "Elizabethkingia_meningoseptica" = "Elizabethkingia_meningoseptica",
  "Enterobacter_aerogenes" = "Enterobacter_aerogenes",
  "Enterobacter_aerogenes_ESBL_Positive" = "Enterobacter_aerogenes",
  "Enterobacter_cloacae_complex" = "Enterobacter_cloacae_complex",
  "Enterococcus_faecalis" = "Enterococcus_faecalis",
  "Enterococcus_faecium_Vancomycin_Resistant" = "Enterococcus_faecium",
  "Enterococcus_species" = "Enterococcus_species",
  "Escherichia_coli" = "Escherichia_coli",
  "Escherichia_coli_ESBL_Positive" = "Escherichia_coli",
  "Haemophilus_influenzae_Beta_Lactamase_Negative" = "Haemophilus_influenzae",
  "Hafnia_alvei" = "Hafnia_alvei",
  "Klebsiella_oxytoca" = "Klebsiella_oxytoca",
  "Klebsiella_pneumoniae" = "Klebsiella_pneumoniae",
  "Klebsiella_pneumoniae_ESBL_Positive" = "Klebsiella_pneumoniae",
  "Lactobacillus_species" = "Lactobacillus_species",
  "Methicillin-Resistant_Staphylococcus_aureus" = "Staphylococcus_aureus",
  "Morganella_morganii" = "Morganella_morganii",
  "Neisseria_meningitidis" = "Neisseria_meningitidis",
  "Neisseria_species" = "Neisseria_species",
  "Pantoea_species" = "Pantoea_species",
  "Proteus_mirabilis" = "Proteus_mirabilis",
  "Providencia_stuartii" = "Providencia_stuartii",
  "Pseudomonas_aeruginosa" = "Pseudomonas_aeruginosa",
  "Pseudomonas_aeruginosa_#2_Encapsulated_Strain" = "Pseudomonas_aeruginosa",
  "Pseudomonas_nitroreducens" = "Pseudomonas_nitroreducens",
  "Raoultella_ornithinolytica_ESBL_Positive" = "Raoultella_ornithinolytica",
  "Serratia_marcescens" = "Serratia_marcescens",
  "Staphylococcus_aureus" = "Staphylococcus_aureus",
  "Staphylococcus_aureus_#2" = "Staphylococcus_aureus",
  "Staphylococcus_coagulase_negative" = "Staphylococcus_coagulasenegative",
  "Stenotrophomonas_maltophilia" = "Stenotrophomonas_maltophilia",
  "Stomatococcus_species" = "Stomatococcus_species",
  "Streptococcus_agalactiae_Group_B" = "Streptococcus_agalactiae",
  "Streptococcus_pneumoniae" = "Streptococcus_pneumoniae",
  "Streptococcus_pseudopneumoniae" = "Streptococcus_pseudopneumoniae",
  "Viridans_streptococcus" = "Streptococcus_species",
  "Viridans_streptococcus_#2" = "Streptococcus_species",
  "Viridans_streptococcus_#3" = "Streptococcus_species",
  "Yeast_Not_Cryptococcus_Species" = "Yeast_Not_Cryptococcus_Species",
  "Yeast_Not_Cryptococcus_Species_#2" = "Yeast_Not_Cryptococcus_Species",
  "organism_null" = "organism_null"
)



