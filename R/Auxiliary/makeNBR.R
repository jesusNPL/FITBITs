library(ape)

##### Function that counts the number of branches in a tree #####
getNumberBranches <- function(phy, phyPC, phyRDM, clade) { 
  
  ## Inits
  vecParents <- numeric(length = 3)
  names(vecParents) <- c("All_spp", "PC", "RDM")
  vecChilds <- numeric(length = 3)
  names(vecChilds) <- c("All_spp", "PC", "RDM")
  
  ## All species
  tabPhy <- data.frame(phy$edge)
  names(tabPhy) <- c("parent", "child")
  ## Phylogenetic completeness
  tabPhyPC <- data.frame(phyPC$edge)
  names(tabPhyPC) <- c("parent", "child")
  ## Random loss
  tabPhyRDM <- data.frame(phyRDM$edge)
  names(tabPhyRDM) <- c("parent", "child")
  
  ## Fill vectors parents
  vecParents[1] <- length(unique(tabPhy$parent))
  vecParents[2] <- length(unique(tabPhyPC$parent))
  vecParents[3] <- length(unique(tabPhyRDM$parent)) 
  
  ## Fill vectors childs
  vecChilds[1] <- length(unique(tabPhy$child))
  vecChilds[2] <- length(unique(tabPhyPC$child))
  vecChilds[3] <- length(unique(tabPhyRDM$child))
  ## Final table
  res <- data.frame(rbind(vecParents, vecChilds)) 
  res$Clade <- clade
  
  return(res)
}

##### Function to select single branches from PC #####

make_PC_unique <- function(treePC) { 
  
  tips <- treePC$tip.label # species list
  genera <- (sapply(strsplit(tips, "_"), function(x) x[1])) # genera list
  
  genera <- data.frame(sort(genera)) 
  names(genera) <- "genus"
  
  genera_count <- genera %>% # count genera
    count(genus) 
  
  genera_count <- genera_count[order(genera_count$n), ] # order counts
  
  genera_unique <- genera_count %>% # filter unique branches
    filter(n == 1)
  
  genera_double <- genera_count %>% # filter double branches
    filter(n == 2)
  
  ii <- sapply(genera_unique$genus, function(x, y) grep(x, y)[1], y = tips)
  treeUnique <- drop.tip(treePC, setdiff(treePC$tip.label, tips[ii]))
  
  iii <- sapply(genera_double$genus, function(x, y) grep(x, y)[1], y = tips)
  treeDouble <- drop.tip(treePC, setdiff(treePC$tip.label, tips[iii]))
  
  res <- list(branchesUnique = treeUnique$tip.label, 
              treeUnique = treeUnique, 
              branchesDouble = treeDouble$tip.label, 
              treeDouble = treeDouble)
  
  return(res)
}

##### Make unique branches #####
#unique_birds <- make_PC_unique(treePC = phy_birds_PC)
