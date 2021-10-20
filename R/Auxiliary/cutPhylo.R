
# http://blog.phytools.org/2014/11/pruning-trees-to-one-member-per-genus.html

#setwd("/Users/jesusnpl/Dropbox/Collaborations/FITBITs")

##### lumberjack function #####
## This function leaves one representative of each lineage at a specific time, 
## for example, representatives at 10 million years.

### Arguments
## tree: a phylogenetic tree
## ages: a vector of ages

lumberjack <- function(tree, ages, plot = FALSE) { 
  
  if ( ! ("phytools" %in% installed.packages())) 
  {install.packages("phytools", dependencies = TRUE)}
  
  library(phytools)
  
  Phylos <- list()
  
  for(i in 1:length(ages)) {
    print(paste0("Cutting phylogeny at MY ", ages[i]))
    ## get all node heights
    H <- phytools::nodeHeights(tree)
    ## time from the root
    t <- max(H)-ages[i]
    ## identify all edges that cross mybp
    h1 <- which(H[, 1] < t)
    h2 <- which(H[, 2] > t)
    ii <- intersect(h1, h2)
    ## all daughter nodes of those edges
    nodes <- tree$edge[ii, 2]
    ## internal phytools function
    getDescendants <- phytools:::getDescendants
    ## find all descendants from each edge
    tips <- lapply(nodes, getDescendants, tree = tree)
    ## find all tips to keep
    tips2 <- tree$tip.label[sapply(tips, function(x, y) x[x <= Ntip(y)][1], y = tree)]
    ## drop all the others
    treeCut <- drop.tip(tree, setdiff(tree$tip.label, tips2))
    ## Plot the cut tree
    if(plot == TRUE) {
      plot(treeCut, show.tip.label = FALSE)
      axisPhylo()
      Sys.sleep(1)
    } 
    ## Store cut tree
    Phylos[[i]] <- treeCut
  }
  return(Phylos)
}


#library(ape)
#treeMAM <- read.nexus("Data/PHYLACINE/Phylogenies/Complete_phylogeny.nex")

#age <- seq(from = 1, to = 100, by = 1)

#MAM_cut <- lumberjack(tree = treeMAM[[3]], ages = age, plot = FALSE)
#MAM_cut[[50]]

##### lumberPD #####

## Function that estimate phylogenetic diversity and number of species on a specific age
## Tree = phylogenetic tree
## ages =  numeric vector of ages
lumberPD <- function(tree, ages, saveTree = FALSE) {
  
  PDobs <- sum(tree$edge.length)
  SR <- Ntip(tree)
  age0 <- 0
  
  orig <- data.frame(PDobs, SR, age0)
  names(orig) <- c("PD", "SR", "My")

  cuts <- lumberjack(tree, ages, plot = FALSE)

  PDcuts <- numeric(length = length(ages))
  SRcuts <- numeric(length = length(ages))
  
  for(i in 1:length(ages)) {
    PDcuts[i] <- sum(cuts[[i]]$edge.length)
    SRcuts[i] <- Ntip(cuts[[i]])
    
  }
  
  chops <- data.frame(PDcuts, SRcuts, ages)
  names(chops) <- c("PD", "SR", "My")
  chops <- rbind(orig, chops)
  #reps <- rands
  #randtrees <- list()
  #randPD <- numeric(length = reps)
  
  #for(i in 1:reps) {
    #temp <- tree
    #randtips <- sample(temp$tip.label, length(temp$tip.label))
    #temp$tip.label <- randtips
    #randPD[i] <- sum(temp$edge.length)
    #randtrees[[i]] <- temp 
  #}
  #class(randtrees) <- "multiPhylo"
  #PDmean <- mean(randPD)
  
  
  
  # Store trees
  if(saveTree == TRUE) { 
    Phylos <- cuts
    results <- list(chops, Phylos)
    return(results)
  } else {
    return(chops)
  }
  
}

check_phylo <- function(tree, brlens_tol = 1e-9){
  #phy <- ape::drop.tip(tree, setdiff(tree$tip.label, sp_keep))
  phy <- ape::multi2di(tree, random = FALSE)
  phy <- phytools::force.ultrametric(phy)
  phy$edge.length[phy$edge.length < brlens_tol-9] <- brlens_tol
  ladderize(phy)
}

##### lumberPD_random #####
## Function that estimates PD by removing N species randomly
## Tree = phylogenetic tree
## nKeep = number of species to keep
## nSim = number of times the PD estimation is repeated

lumberPD_random <- function(tree, nKeep = 4000, nSims = 100) { 
  spp <- tree$tip.label
  
  PDrand <- numeric(length = nSims)
  SRrand <- numeric(length = nSims)
  #PDsing <- numeric(length = nSims)
  
  for(i in 1:nSims) { 
    svMisc::progress(i, max.value = nSims)
    
    sppKeep <- sample(spp, nKeep)
    #sppKeep <- spp[spp %in% sppDrop]
    tree2 <- drop.tip(tree, setdiff(tree$tip.label, sppKeep))
    #tree3 <- phytools::drop.tip.singleton(tree = tree, tip = )
    PDrand[i] <- sum(tree2$edge.length)
    SRrand[i] <- Ntip(tree2)
    #PDsing[i] <- sum(tree3$edge.length)
  }
  
  rands <- data.frame(cbind(PDrand, SRrand))
  means <- mean(PDrand) 
  medians <- median(PDrand)
  quants <- quantile(PDrand, probs = c(0.05, 0.95))
  
  res <- list(means, medians, quants, rands)
  names(res) <- c("PD_rand_mean", "PD_rand_median", 
                  "PD_rand_quantiles", "PD_rands")
  return(res)

}

#tree_mam <- read.nexus("Data/Phylos/Mammal_tree/mcc_mammal_tree.nex")

#x <- lumberPD_random(tree = tree_mam, nKeep = 100, nSims = 100)

##### lumberPD_clade #####
## Function that estimate PD by removing nodes of specific species 

lumberPD_clade <- function(tree, sppDrop) {
  tree$node.label <- paste0("node", seq(1:tree$Nnode)) 
  ## Selecting a node of interest (e.g. MRCA)
  mrca_node <- getMRCA(tree, tip = sppDrop)
  ## The node label for the mrca_node
  nodo <- tree$node.label[mrca_node-Ntip(tree)]
  ## Extract clade
  tree2 <- extract.clade(tree, nodo)
  ## Species to keep
  spp2keep <- setdiff(tree$tip.label, tree2$tip.label)
  ## Drop species from the clade
  tree3 <- drop.tip(tree, setdiff(tree$tip.label, spp2keep))
  # Calculations
  PDclade <- sum(tree3$edge.length)
  SRclade <- Ntip(tree3)
  PDtot <- sum(tree$edge.length)
  SRtot <- Ntip(tree)
  dPD <- PDtot-PDclade
  dSR <- SRtot-SRclade

  res <- round(c(PDclade, SRclade, PDtot, SRtot, dPD, dSR))
  names(res) <- c("PD_clade", "SR_clade", "PDtotal", "SR_total", "deltaPD", "deltaSR")
  return(res)
}

## Usage
#dropSpp <- c("Dasypus_yepesi", "Dasypus_hybridus", "Dasypus_pilosus", 
             #"Dasypus_novemcinctus", "Dasypus_sabanicola")

#tt <- lumberPD_clade(tree = tree_mam, sppDrop = dropSpp)
#tt
