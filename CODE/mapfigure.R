#Map creation for Beetle Manuscript with Nicola Rammell
#Created 2 March 2020
#Using Hakai Institute's guide to mapping in R, with coastal focus
#Code modified from: https://hecate.hakai.org/rguide/mapping-in-r.html#site-maps

####PACKAGE LOADING####
library(raster)
library(maps) 
library(mapdata)
library(mapproj)
library(maptools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(ggsn)
library(tidyverse)
library(here)
library(sf)
library(ggspatial)

####MAP CODE####

###Create a BC coast map
#load the world map shapefile database
#this database has a lower resolution (which is fine for large scale map)
m <- map_data("world", c("usa", "Canada"))

#make a basic BC coast map with Vancouver labelled
ggplot() + 
  geom_polygon(data = d, aes(x = long, y = lat, group = group),
                        fill = 'grey80') + 
  geom_point(data = BC.df, aes(x = -123.141759, y = 49.251556), size = 2) +
  geom_text(aes(x = -123.141759, y = 49.0, label = "Vancouver, BC")) +
  theme_classic() +
  coord_map("conic", lat0 = 18, xlim=c(227, 240), ylim=c(46,55)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank())
ggsave("FIGURES/BCcoastmap.png", dpi = 1000)


###Create a smaller inset with detailed coastlines of Denny and Campbell Island
#load shapefiles, which were downloaded from the Hakai Institute
BC.shp <- readOGR(here("data","05_shapefiles", "COAST_TEST2.shp"))

#chose the lat/long extent you want to show
kunsoot <- extent(-128.2, -127.5, 52.09, 52.25)

#crop your shapefile polygons to the extent defined
#Heads up, this takes a WHILE to run! ~30 mins
BC.shp2 <- crop(BC.shp, kunsoot)

### project and fortify (i.e. turn into a dataframe)
BC.df <- fortify(BC.shp2)

###Creating a map of Campbell and Denny Island with Bella Bella and Kunsoot
ggplot() + 
  theme_bw() +
  geom_polygon(data = BC.df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.1, fill = 'grey95')+
  coord_cartesian(xlim = c(-128.18, -128), ylim = c(52.11, 52.2)) +
  geom_point(data = BC.df, aes(x = -128.143314, y = 52.163193), size = 2) +
  geom_text(aes(x = -128.143314, y = 52.158, label = "Bella Bella", size = 3)) +
  geom_point(data = BC.df, aes(x = -128.008487, y = 52.149128), size = 2) +
  geom_text(aes(x = -128.025, y = 52.143, label = "Kunsoot River", size = 3)) +
  scale_y_continuous(breaks = c(52.12, 52.15, 52.18)) +
  scale_x_continuous(breaks = c(-128.15, -128.09, -128.03)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10)) +
  ggsn::scalebar(BC.df, dist = 1, transform = TRUE, model = "WGS84",
                 dist_unit = "km", location = "bottomright", st.bottom = FALSE,
                 st.dist = 0.03, border.size = 0.5, 
                 anchor = c(x = -128.01, y = 52.110))
ggsave("FIGURES/dennyislandmap.png", dpi = 1000)

#The BC coastline map and the Denny Island inset map were overlaid on one another
#after they were exported from R to create manuscript figure 1
