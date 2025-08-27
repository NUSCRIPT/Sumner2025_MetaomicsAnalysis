library(tidyverse)

source('src/common/load_metadata.r')

# Snippet from rstatix that helped me decide ordering of data structure.
# Homogeneity of proportions between groups
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# H0: the proportion of smokers (rows) is similar in the four groups (cols)
# Ha:  this proportion is different in at least one of the populations (cols).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [JACK'S]
# H0: the proportion of clinical group X is similar in the four pneumotypes
# Ha: this proportion is different in at least one of the populations. 

# Step 1: Transform metadata to pneumotype list of count matrix (Clinical Data = Rows; Pneumotype = Columns)
mdCategorical <- permmeta %>% 
  drop_na(cluster_num_amplicon) %>% 
  select(!where(is.numeric)) %>% # get categorical only
  filter(baseline_or_bal == 'Baseline')
myCatVars <- mdCategorical %>% select(!c("cluster_num_amplicon", "sample_name", 'has_bal',
                                         'COVID_status', 'Intubation_flag', "Patient_id/ICU_stay/ICU_day",
                                         'pna_category', 'Patient_category', 
                                         'baseline_or_bal')) %>% 
  colnames()
# LIST STRUCTURE 
# list(aCatVar = list(mdMat = c('tests'), mdFischer = c('testa')), ... = ..., )
results <- list()

for (i in myCatVars){
  print(is.character(i))
  mdMati <- mdCategorical %>% 
    select(cluster_num_amplicon, i) %>% 
    group_by(.data[[i]]) %>% 
    count(cluster_num_amplicon) %>%
    ungroup() %>%  
    drop_na(.data[[i]]) %>%
    pivot_wider(names_from = 2, values_from = 3,values_fill = 0) %>% 
    column_to_rownames(var = i)
  
  results[[i]] = list(name = i, mdMat = mdMati, mdFischer = NULL)
  
  
}

# Step 2: Run fisher test w/ simulated p values on list
library(rstatix)
set.seed(100)
for (i in myCatVars){
  print(i)    
  results[[i]]$mdFischer = fisher_test(results[[i]]$mdMat, results[[i]]$mdMat, simulate.p.value = T, detailed = T,B = 5000)
}


# Step 3: Compile, FDR correct, and plot
tidyList <- function(object){object$mdFischer %>% mutate(name = object$name, .before = n)}
tidyRes <- do.call('rbind', lapply(results, tidyList)) %>% 
  mutate(q =  p.adjust(p),
         q.symbol = gtools::stars.pval(q)) %>% 
  arrange(q) 
source('src/common/nature_theme.r')
ggplot(tidyRes, aes(x='Pneumotype', y = name, fill = -log10(q), label=q.symbol)) +
  geom_tile(color='black') +
  geom_text(color = 'white')+
  theme_nature() +
  ggtitle("Baseline\npneumotype \nassociations \n(Fischer's test)")
ggsave("figures/pneumotypes/pneumotype_anovas_baseline.pdf", units = 'in', width=1.82, height =3.46)

