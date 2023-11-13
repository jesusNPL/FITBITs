
library(tidyverse)
library(patchwork)
library(cowplot)
library(scales)

### Load data

load("output/RCP/My01/Seed_plants/cps_PD_SEEDS.RData")

### Inits
clades <- c("Seed", "woody")

cpsPD <- cps_PD[[1]]

##### Plot chops #####
scl = 10

## Seed plants
lim_seed <- c(cpsPD[1, 2], cpsPD[1, 3], cpsPD[1, 4])

plot_seed <- treeChops$Seed %>% 
  mutate(mya = (My+1)*1000000) %>% 
  ggplot(aes(x = (mya))) + 
  geom_line(aes(y = PD, colour = "PD"), size = 2.5, color = "black", alpha = 0.7) + 
  geom_line(aes(y = SR*scl, colour = "SR"), size = 2.5, color = "red", alpha = 0.7) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(name = "Phylogentic diversity", 
                     breaks = seq(0, max(treeChops$Seed$PD)+10000, 
                                  by = round(max(treeChops$Seed$PD)/7)), 
                     sec.axis = sec_axis(~./scl, name = "Number of species", 
                                         breaks = seq(0, max(treeChops$Seed$SR)*scl, 
                                                      by = round((max(treeChops$Seed$PD)/7)/10)))) +
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
        axis.line.y.right = element_line(size = 1, color = "red"), 
        plot.title = element_text(size = 20)
  ) + 
  #viridis::scale_color_viridis(discrete = TRUE) +
  labs(x = "Time (My)") + 
  ggtitle("Seed plants - OTB") + 
  geom_vline(xintercept = (lim_seed[1]*1000000), # species cut
             color = "black", linetype = "dashed", size = 1, alpha = 0.7) + 
  geom_vline(xintercept = (lim_seed[2]*1000000), 
             color = "black", linetype = "dotted", size = 0.8, alpha = 0.7) + 
  geom_vline(xintercept = (lim_seed[3]*1000000), 
             color = "black", linetype = "dotted", size = 0.8, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$Seed[6, ][2]*scl), 
             color = "red", linetype = "dashed", size = 1, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$Seed[1, ][2]*scl), 
             color = "red", linetype = "F1", size = 1, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$Seed[4, ][2]*scl), 
             color = "red", linetype = "dotted", size = 0.8, alpha = 0.7)

plot_seed
head(treeChops$Seed, 10)

### Woody plants
lim_wood <- c(cpsPD[2, 2], cpsPD[2, 3], cpsPD[2, 4])

plot_wood <- treeChops$woody %>% 
  mutate(mya = (My+1)*1000000) %>% 
  ggplot(aes(x = (mya))) + 
  geom_line(aes(y = PD, colour = "PD"), size = 2.5, color = "black", alpha = 0.7) + 
  geom_line(aes(y = SR*scl, colour = "SR"), size = 2.5, color = "red", alpha = 0.7) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(name = "Phylogentic diversity", 
                     breaks = seq(0, max(treeChops$woody$PD), 
                                  by = round(max(treeChops$woody$PD)/7)), 
                     sec.axis = sec_axis(~./scl, name = "Number of species", 
                                                    breaks = seq(0, max(treeChops$woody$SR)*scl, 
                                                                 by = round((max(treeChops$woody$PD)/7)/10)))) +
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
        axis.line.y.right = element_line(size = 1, color = "red"), 
        plot.title = element_text(size = 20)
  ) + 
  #viridis::scale_color_viridis(discrete = TRUE) +
  labs(x = "Time (My)") + 
  ggtitle("Woody plants") + 
  geom_vline(xintercept = (lim_wood[1]*1000000), # species cut
             color = "black", linetype = "dashed", size = 1, alpha = 0.7) + 
  geom_vline(xintercept = (lim_wood[2]*1000000), 
             color = "black", linetype = "dotted", size = 0.8, alpha = 0.7) + 
  geom_vline(xintercept = (lim_wood[3]*1000000), 
             color = "black", linetype = "dotted", size = 0.8, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$woody[29, ][2]*scl), 
             color = "red", linetype = "dashed", size = 1, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$woody[1, ][2]*scl), 
             color = "red", linetype = "F1", size = 1, alpha = 0.7) + 
  geom_hline(yintercept = as.numeric(treeChops$woody[23, ][2]*scl), 
             color = "red", linetype = "dotted", size = 0.8, alpha = 0.7)

plot_wood
head(treeChops$woody, 10)

##### Plot plants#####

save(plot_seed, plot_wood, 
     file = "output/figures/plots_phylocuts_SR_SEED_LOGS.RData")

pdf("output/figures/fig_xx_PD_SR_cuts_SeedPlants_clades_horizontal.pdf", 
    width = 12, height = 4)

plot_grid(plot_seed, plot_wood,
          nrow = 1, ncol = 2, #labels = c("A)", "B)", "C)", "D)"), 
          label_size = 17, hjust = -25, vjust = 2)

dev.off()

