
library(mcp)
library(tidybayes)
library(tidyverse)

##### Load and prepare data #####
load("output/Vertebrates/tree_cuts_VerTer_01_My.RData")

clade <- "birds"


treeChops_Birds_log <- treeChops[[3]] %>% 
    filter(My <= 30) %>% 
    mutate(SR_log = log(SR), PD_log = log(PD)) 

##### Inits #####
iter = 25000
chains = 4

## Model PD

mod_PD <- list(
  PD ~ 1 + My,  # intercept + slope
  1 ~ 0 + My,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
) 

# Null model
mod_PD_null <- list(PD ~ 1)

##### Run Bayesian Regressions with Change Points #####
## Run model
fit_pd <- mcp(data = treeChops_Birds_log, 
              model = mod_PD, 
              chains = chains, 
              iter = iter) 

pps_PD <- plot(fit_pd)
pps_PD

fit_pd$loo <- loo(fit_pd) 

## Run null model
fit_pd_null <- mcp(mod_PD_null, 
                   data = treeChops_Birds_log, 
                   par_x = "My")  # fit null model here 

fit_pd_null$loo <- loo(fit_pd_null)

### Save models 
# PD
save(fit_pd, fit_pd_null, pps_PD, 
     file = "output/RCP/BRCP_PD_Birds.RData") 

##### Model comparison #####
print(loo::loo_compare(fit_pd_null$loo, fit_pd$loo))
  
##### Extract changing points #####
cps1_PD <- summary(fit_pd)[1, ]
cps2_PD <- summary(fit_pd)[2, ]
cps3_PD <- summary(fit_pd)[3, ] 

### Plot Birds over 30 MY
pBirds <- pps_PD + labs(x = "Time (My)", y = "Total PD loss  - Birds") + 
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), 
        legend.position = "none", 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15), 
        axis.line.x = element_line(size = 1), 
        axis.line.y = element_line(size = 1), 
        plot.title = element_text(size = 20)
  ) + 
  geom_vline(xintercept = as.numeric(cps1_PD[2]), color = "darkgray", linetype = "dashed", size = 1) + 
  geom_vline(xintercept = as.numeric(cps1_PD[3]), color = "darkgray", linetype = "dashed", size = 1) +
  geom_vline(xintercept = as.numeric(cps1_PD[4]), color = "darkgray", linetype = "dashed", size = 1) 

pBirds  

pdf(file = "output/RCP/BRCP_birds.pdf", 
    height = 7, width = 10)

pBirds 

dev.off()
