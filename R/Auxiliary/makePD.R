library(tidyverse)
library(ape)

makePD <- function(tree, pam, grid) { 
  
  gridNames <- grid$grids
  
  srGrid <- numeric(length = nrow(grid))
  pdGrid <- numeric(length = nrow(grid))
  
  for(i in 1:nrow(grid)) { 
    svMisc::progress(value = i, max.value = nrow(grid), progress.bar = TRUE) 
    
    sppSel <- pam %>% 
      filter(grids == gridNames[i])
    
    treeSel <- drop.tip(tree, setdiff(tree$tip.label, sppSel$species))
    
    srGrid[i] <- Ntip(treeSel)
    pdGrid[i] <- sum(treeSel$edge.length)
    
  } 
  
  resPD <- data.frame(grids = gridNames, SR = srGrid, PD = pdGrid) 
  
  res <- full_join(resPD, grid, by = "grids")
  
  return(res)
  
}
