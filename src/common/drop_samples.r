# -----------------------------------------------------------------------------
# Script Name: drop_samples.r
# Author: Jack T. Sumner
# Contact: jacksumner2026@u.northwestern.edu
#
# Description:
# Purpose of this script is to drop samples that were either removed from 
# the IRB or for which their is ambiguity for the tubes true identity.
#
# Acknowledgments:
# Developed by Jack Sumner. Please cite appropriately if used in published work.
#
# Usage Restrictions:
# This script is provided "as is" without warranty of any kind. It may be used
# and modified for academic and non-commercial purposes only. Redistribution
# or commercial use requires written permission from the author.
#
# If you use or adapt this code, please acknowledge the original author.
# -----------------------------------------------------------------------------

# Note that these did actually have BAL ids in them. These samples are not integral to reported analyses.
# Sample IDs redacted for security reasons. Script left in to show functionality in other parts of analysis.
get_drop_lst <- function(){
  # Tubes for which there is ambiguity
  flip_lst <- c("INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs")
  
  # Tubes that were removed
  rm_lst <- c("INTERNALBALIDs", "INTERNALBALIDs")
  
  # HTX & MGX technical replicates - added 7/16/25
  htx_replicates <- c("INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", 
                      "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs")
  
  # Transplant samples
  transplant <- c("INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs", "INTERNALBALIDs")
  
  
  drop_lst <- c(flip_lst, rm_lst, htx_replicates, transplant)
  
  return(drop_lst)
  
}

drop_samples <- function(){
  
  drop_lst <- get_drop_lst()
  
  if (exists('permmeta')) {
    print('DROPPING SAMPLES FROM PERMMETA')
    permmeta <<- permmeta %>% filter(!sample_name %in% drop_lst)
  }
  
  if (exists('mox_feature_list')) {
    print('DROPPING SAMPLES FROM FEATURE')
    mox_feature_list <<- lapply(mox_feature_list, select, !any_of(drop_lst))
    
  }
  
  if (exists('mox_dist_list')) {
    library(usedist)
    dist_subset_negative <- function(my_dist, my_lst){
      new_dist <- dist_subset(my_dist, labels(my_dist)[!labels(my_dist) %in% my_lst])
      return(new_dist)
    }
    
    print('DROPPING SAMPLES FROM DISTANCE')
    mox_dist_list <<- lapply(mox_dist_list, dist_subset_negative, my_lst=drop_lst)
  }
  
}

