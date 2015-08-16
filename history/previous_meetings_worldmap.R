## plot all GEOSTATS using plotGoogleMaps:

library(plotGoogleMaps)
library(sp)
library(RColorBrewer)
library(rgdal)

## download the table:
geostat <- read.csv("GEOSTAT_previous_meetings.csv")
str(geostat)
coordinates(geostat) <- ~ Longitude + Latitude
proj4string(geostat) <- CRS("+proj=longlat +datum=WGS84")

m1 = plotGoogleMaps(geostat, iconMarker=iconlabels(geostat$Name, colPalette='white'), zoom=12, mapTypeId='TERRAIN', control.width=0, add=T, zIndex = 2)
geostat$url <- paste("<a href='",geostat$Photo_URL,"'>", geostat$Name," </a>",sep='' )
geostat.s <- geostat[,c("Name", "url")]
m2 = plotGoogleMaps(geostat.s, iconMarker=geostat$Icon_URL, zoom=12,   mapTypeId='TERRAIN', clickable=TRUE, zIndex=1, control.width= 0, previousMap=m1, filename='GEOSTAT_history.htm')