
library(mcp)
library(tidybayes)
library(tidyverse)

##### Load and prepare data #####
load("output/Vertebrates/tree_cuts_VerTer_01_My.RData")

clades <- c("amphibians", "squamata", 
            "birds", "mammals")

treeChops_log <- list() 

for(i in 1:length(clades)) {
  treeChops_log[[i]] <- treeChops[[i]] %>% 
    #filter(My <= 30) %>% 
    mutate(SR_log = log(SR), PD_log = log(PD)) 
}

##### Inits #####
iter = 25000
chains = 4

## Model SR
mod_SR <- list(
  SR ~ 1 + My,  # intercept + slope
  1 ~ 0 + My,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
)

mod_SR_log <- list(
  SR_log ~ 1 + My,  # intercept + slope
  1 ~ 0 + My,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
)
# Null model
mod_SR_null <- list(SR ~ 1)

## Model PD

mod_PD <- list(
  PD ~ 1 + My,  # intercept + slope
  1 ~ 0 + My,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
) 

mod_PD_log <- list(
  PD_log ~ 1 + My,  # intercept + slope
  1 ~ 0 + My,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
)
# Null model
mod_PD_null <- list(PD ~ 1)

##### Run Bayesian Regressions with Change Points #####
## SR models
models_SR <- list()
models_SR_null <- list()

## PD models
models_PD <- list()
models_PD_null <- list()

### LOGS

## SR models
models_SR_LOG <- list()

## PD models
models_PD_LOG <- list()

## Plots
pps_SR <- list()
pps_PD <- list()

pps_SR_LOG <- list()
pps_PD_LOG <- list()

for(j in 1:length(clades)) {
  print(paste0("Starting Bayesian Regressions with CP for ", clades[j], " clade"))
  
  ### Models
  fit_sr <- mcp(data = treeChops_log[[j]], 
                     model = mod_SR, 
                     chains = chains, 
                     iter = iter) 
  pps_SR[[j]] <- plot(fit_sr)
  fit_sr$loo <- loo(fit_sr) 
  models_SR[[j]] <- fit_sr
  
  fit_pd <- mcp(data = treeChops_log[[j]], 
                model = mod_PD, 
                chains = chains, 
                iter = iter) 
  pps_PD[[j]] <- plot(fit_pd)
  fit_pd$loo <- loo(fit_pd) 
  models_PD[[j]] <- fit_pd
  
  ### LOG 
  fit_sr_log <- mcp(data = treeChops_log[[j]], 
                model = mod_SR_log, 
                chains = chains, 
                iter = iter) 
  pps_SR_LOG[[j]] <- plot(fit_sr_log)
  fit_sr_log$loo <- loo(fit_sr_log) 
  models_SR_LOG[[j]] <- fit_sr_log
  
  fit_pd_log <- mcp(data = treeChops_log[[j]], 
                model = mod_PD_log, 
                chains = chains, 
                iter = iter) 
  pps_PD_LOG[[j]] <- plot(fit_pd_log)
  fit_pd_log$loo <- loo(fit_pd_log) 
  models_PD_LOG[[j]] <- fit_pd_log
  
  ### Null models
  fit_sr_null <- mcp(mod_SR_null, 
                          data = treeChops_log[[j]], 
                          par_x = "My")  # null model here
  fit_sr_null$loo <- loo(fit_sr_null)
  models_SR_null[[j]] <- fit_sr_null 
  
  fit_pd_null <- mcp(mod_PD_null, 
                     data = treeChops_log[[j]], 
                     par_x = "My")  # fit null model here
  fit_pd_null$loo <- loo(fit_pd_null)
  models_PD_null[[j]] <- fit_pd_null

}

### Save models 
# SR
save(models_SR, 
     file = "output/RCP/My01/BRCP_SR_30my.RData") 

save(models_SR_LOG,  
     file = "output/RCP/My01/BRCP_SR_LOG_30my.RData")

save(models_SR_null, 
     file = "output/RCP/My01/BRCP_SR_NULL_30my.RData")
# PD
save(models_PD, 
     file = "output/RCP/My01/BRCP_PD_30my.RData") 

save(models_PD_LOG,  
     file = "output/RCP/My01/BRCP_PD_LOG_30my.RData") 

save(models_PD_null,  
     file = "output/RCP/My01/BRCP_PD_NULL_30my.RData") 
### Save plots
# SR
save(pps_SR,  
     file = "output/RCP/My01/plots_BRCP_SR_30my.RData")
save(pps_SR_LOG, 
     file = "output/RCP/My01/plots_BRCP_SR_LOG_30my.RData")
# PD
save(pps_PD, clades, 
     file = "output/RCP/My01/plots_BRCP_PD_30my.RData")
save(pps_PD_LOG,  
     file = "output/RCP/My01/plots_BRCP_PD_LOG_30my.RData")

##### Model comparison #####

for(k in 1:length(clades)) {
  
  print(loo::loo_compare(models_SR[[k]]$loo, models_SR_null[[k]]$loo))
  print(loo::loo_compare(models_PD[[k]]$loo, models_PD_null[[k]]$loo))
  
}

##### Extract changing point #####
## SR
cps1_SR <- list()
cps2_SR <- list()
cps3_SR <- list()

## PD
cps1_PD <- list()
cps2_PD <- list()
cps3_PD <- list()

for(i in 1:length(clades)) {
  ## SR
  cps1_SR[[i]] <- summary(models_SR[[i]])[1, ]
  cps2_SR[[i]] <- summary(models_SR[[i]])[2, ]
  cps3_SR[[i]] <- summary(models_SR[[i]])[3, ] 
  ## PD
  cps1_PD[[i]] <- summary(models_PD[[i]])[1, ]
  cps2_PD[[i]] <- summary(models_PD[[i]])[2, ]
  cps3_PD[[i]] <- summary(models_PD[[i]])[3, ] 
}

## SR
cps1SR <- do.call(rbind, cps1_SR)
cps1SR$Clade <- clades
cps2SR <- do.call(rbind, cps2_SR)
cps2SR$Clade <- clades
cps3SR <- do.call(rbind, cps3_SR)
cps3SR$Clade <- clades

cps_SR <- list(cps1SR, cps2SR, cps3SR)
save(cps_SR, file = "output/RCP/My05/cps_SR_all.RData")

## PD
cps1PD <- do.call(rbind, cps1_PD)
cps1PD$Clade <- clades
cps2PD <- do.call(rbind, cps2_PD)
cps2PD$Clade <- clades
cps3PD <- do.call(rbind, cps3_PD)
cps3PD$Clade <- clades

cps_PD <- list(cps1PD, cps2PD, cps3PD)
save(cps_PD, file = "output/RCP/My05/cps_PD_all.RData")

pps[[6]] + labs(x = "Time (My)", y = "Total PD loss  - Plants") + 
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
  geom_vline(xintercept = 5, color = "darkgray", linetype = "dashed", size = 1)
  
predict(fit)
