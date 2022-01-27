setwd("Z:/3-Personal Research Folders/Jesus/FITBITs")

# load packages
library(dplyr)
#library(ggmap)
library(rgdal)
library(raster)
library(sf)

birds <- st_read(dsn = "Ranges_birds/BOTW.gdb", layer = "All_Species")

head(birds)


spp <- unique(birds$binomial)

spp2 <- gsub(" ", "_", spp)

nSPP <- length(spp)

for(i in 1:nSPP) {
  print(spp[i])
  
  sel_spp <- birds %>% 
    filter(binomial == spp[i]) 
  
  st_write(sel_spp, paste0("Ranges_birds/Ranges_ind_birds/", spp2[i], ".shp"), append = FALSE)
  
}


