# First, install and load "CommEcol" package.

#install.packages("CommEcol")
library(CommEcol)

# Now copy and paste the following function into R.

radiusextract <- function(lon, lat, ras, radius = 0.5) {
  coor <- coordinates(ras)
  lonpar <- coor[, 1][which(abs(coor[, 1] - lon) %in% min(abs(coor[, 1] - lon)))]
  latpar <- coor[, 2][which(abs(coor[, 2] - lat) %in% min(abs(coor[, 2] - lat)))]
  data <- data.frame(x = coor[, 1], y = coor[, 2], val = values(ras))
  return(select.window(xf = lonpar[1], yf=latpar[1], radius = radius, xydata = data)) 
}
# Now we can use the function to search for the temperature values 
# of not only the raster cell that the coordinates fall, 
#but also the raster cells within a distance of 0.167 degree (approximately 10 arcminutes). 

## Arguments: 
# "lon" is the longitude value of a place, 
# "lat" is the latitude value, 
# "gridshp" is the raster with the values you want (for example, temperature) 
# "radius" is the distance radius which will be used to search for values.

getData('worldclim', var = 'bio', res = 10)

bio1 <- raster("wc10/bio1.bil")
res <- radiusextract(lon = -116.867, lat = 36.462, ras = bio1, radius = 0.167)
res
