
library(ape)
library(PDcalc)

## Function to correct tree

check_and_fix_ultrametric <- function(phy){
  
  if (!is.ultrametric(phy)){
    
    vv <- vcv.phylo(phy)
    dx <- diag(vv)
    mxx <- max(dx) - dx
    for (i in 1:length(mxx)){
      phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
    }
    if (!is.ultrametric(phy)){
      stop("Ultrametric fix failed\n")
    }	
  }
  
  return(phy)
}
##### Load data #####

### Seed plants
tree_seed_otb <- read.tree("Data/Phylos/big_seed_plant_trees/ALLOTB.tre")
tree_seed_mb <- read.tree("Data/Phylos/big_seed_plant_trees/ALLMB.tre")

### solve polytomies
tree_seed_otb_check <- check_and_fix_ultrametric(tree_seed_otb)
tree_seed_mb_check <- check_and_fix_ultrametric(tree_seed_mb)


tree_seed_otb_check <- multi2di(tree_seed_otb)
is.ultrametric(tree_seed_otb_check)
is.binary(tree_seed_otb_check)

tree <- bifurcatr(tree_seed_mb, runs = 1)
plot(tree, show.tip.label = F)

#write.nexus(tree, file = "Data/Phylos/big_seed_plant_trees/ALLOTB_binary.nex")
write.nexus(tree, file = "Data/Phylos/big_seed_plant_trees/ALLMB_binary.nex")
