
setwd("Dropbox/Collaborations/FITBITs/")

library(ape)

source("R/cutPhylo.R")

##### Load adata #####

### Amphibians
tree_amph <- read.nexus("Data/Phylos/Amphibian_tree/mcc_amphibian_tree.nex")
tree_amph <- drop.tip(tree_amph, "Homo_sapiens")

### Squamata
tree_squa <- read.nexus("Data/Phylos/Squamate_tree/mcc_squamate_tree.nex")

### Birds
tree_bird_hack <- read.nexus("Data/Phylos/Bird_hackett_tree/mcc_birds_hackett_tree.nex")
#tree_bird_mayr <- read.nexus("Data/Phylos/Bird_mayr_hackett_tree/mcc_birds_mayr_hackett_tree.nex")

### Mammals
tree_mam <- read.nexus("Data/Phylos/Mammal_tree/mcc_mammal_tree.nex")

##### Run lumberjack #####
### Inits
age <- seq(from = 0.1, to = 100, by = 0.1)

phy_lst <- list(tree_amph, tree_squa, 
                tree_bird_hack, tree_mam)

clades <- c("amphibians", "squamata", 
            "birds", "mammals")

times <- length(clades)

treeChops <- list()

### Run
for(i in 1:times) {
  
  print(paste0("Starting lumberjack for ", clades[i], " clade..."))
  
  chopTMP <- lumberPD(tree = phy_lst[[i]], ages = age)

  treeChops[[i]] <- chopTMP
}

names(treeChops) <- clades

save(treeChops, file = "output/Vertebrates/tree_cuts_VerTer_01_My.RData")

### Plots
plot(treeChops$amphibians[, 1], col = "black", pch = 16, frame.plot = 0, 
     ylab = "PD", xlab = "MY")

lines(treeChops$amphibians[, 2], col = "red", lwd = 2, ylab = "PD", xlab = "MY")

##### Cut Phylogeny #####

### Inits
age <- seq(from = 0.1, to = 100, by = 0.1)

phy_lst <- list(tree_amph, tree_squa, 
                tree_bird_hack, tree_mam)

clades <- c("amphibians", "squamata", 
            "birds_hack", "mammals")

times <- length(clades)

Chops <- list()

### Run
for(i in 1:times) {
  
  print(paste0("Starting lumberjack for clade ", clades[i]))
  
  chopTMP <- lumberPD(tree = phy_lst[[i]], ages = age, saveTree = TRUE)
  
  Chops[[i]] <- chopTMP[[2]]
}
Chops[[4]]

save(Chops, clades, age, phy_lst, 
     file = "output/Vertebrates/tree_chops_VerTer_01_MY.RData")

## add extra space to right margin of plot within frame
par(mar = c(5, 4, 4, 6) + 0.1)

## Plot first set of data and draw its axis
plot(treeChops$birds[, 3], treeChops$birds[, 1], pch = 16, axes = FALSE, 
     ylim = c(0, max(treeChops$birds[, 1])), 
     xlab = "", ylab = "", type = "b", col = "black", main = NULL)
axis(2, ylim = c(0, treeChops$birds[, 1]), col = "black", las = 1)  ## las=1 makes horizontal labels
mtext("Phylogenetic diversity", side = 2, line = 2.5)
box()

## Allow a second plot on the same graph
par(new = TRUE)

## Plot the second plot and put axis scale on right
plot(treeChops$birds[, 3], treeChops$birds[, 2], pch = 15,  xlab = "", ylab = "", 
     ylim = c(0, max(treeChops$birds[, 2])), 
     axes = FALSE, type = "b", col = "red")
## a little farther out (line=4) to make room for labels
mtext("Species richness", side = 4, col = "red", line = 4) 
axis(4, ylim = c(0, max(treeChops$birds[, 2])), col = "red", col.axis = "red", las = 1)

## Draw the time axis
axis(1, pretty(range(treeChops$birds[, 3]), 10))
mtext("Time (Mya)", side = 1, col = "black", line = 2.5)  

## Add Legend

library(tidyverse)

scl = 10

treeChops[[1]] %>% 
  ggplot(aes(x = My)) + 
  geom_line(aes(y = PD, colour = "PD"), size = 2.5, color = "black", alpha = 0.7) + 
  geom_line(aes(y = SR*scl, colour = "SR"), size = 2.5, color = "red", alpha = 0.7) +
  scale_y_continuous(name = "Phylogentic diversity", 
                     sec.axis = sec_axis(~./scl, name = "Number of species")) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.text.y.right = element_text(colour = "red"),
        axis.ticks.y.right = element_line(colour = "red"),
        axis.title.y.right = element_text(colour = "red"), 
        axis.line.x = element_line(size = 1), 
        axis.line.y = element_line(size = 1), 
        axis.line.y.right = element_line(size = 1, color = "red")
  ) + 
  #viridis::scale_color_viridis(discrete = TRUE) +
  labs(x = "Time (Mya)")

