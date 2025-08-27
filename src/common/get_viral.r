
get_viruses <- function(distance=FALSE){
   
    
    features_viral <- readRDS("objects/mgx/mgx_07_table_vOTUsLog.rds")
    
    if (distance==T){
        
        distance_viral <- readRDS("objects/mgx/mgx_08_distance_vOTUsASTJaccard.rds")
        return(distance_viral)
    } else {
        
        return(features_viral)
        
    }

}
