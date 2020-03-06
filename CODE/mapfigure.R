#Map creation for Beetle Manuscript with Nicola Rammell
#Created 2 March 2020
#Using Hakai Institute's guide to mapping in R, with coastal focus
#Code modified from: https://hecate.hakai.org/rguide/mapping-in-r.html#site-maps

####PACKAGE LOADING####

#Original code pkgs
#library(ggplot2)
#library(sf)
#library(rnaturalearth)
#library(rnaturalearthdata)
#library(rgeos)
#library(ggspatial)

#Hakai new code libraries
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

####MAP CODE####
#load shapefile date
BC.shp <- readOGR(here("data","05_shapefiles", "COAST_TEST2.shp"))

# chose the lat/long extent you want to show
kunsoot <- extent(-128.3, -127.94, 51.5, 52.5)

### crop your shapefile polygons to the extent defined
#Heads up, this takes a WHILE to run! ~10-15 mins
BC.shp2 <- crop(BC.shp, kunsoot)

### project and fortify (i.e. turn into a dataframe)
BC.df <- fortify(BC.shp2)

# (IF DESIRED) Load .csv file with your specific study site lat/longs
# this file is a dataframe with 4 columns: site_name, otterOcc(Y or N), lat, long  
# EXPTsites <- read.csv("/Users/jennb/Dropbox/Simple_BC_map/EXPTsites.csv", header = T)

ggplot() + 
  theme_bw() +
  geom_polygon(data = BC.df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.1, fill = 'grey95') +
  coord_cartesian(xlim = c(-128.18, -128), ylim = c(52.11, 52.2)) +
  #geom_point(data=EXPTsites, aes(x=long, y=lat, shape=otter), size=4, 
  #colour="blue", stroke=1.5) +  
  #add this to plot site locations
  scale_y_continuous(breaks = c(52.12, 52.15, 52.18)) +
  scale_x_continuous(breaks = c(-128.15, -128.09, -128.03)) +
  north(data = BC.df, scale = 0.05, symbol = 3, 
        anchor = c(x = -127.999, y = 52.218)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10)) +
  ggsn::scalebar(BC.df, dist = 1, transform = TRUE, model = 'WGS84', 
                 dist_unit = "m") +
  annotation_scale(location = "bl", width_hint = 0.5) +
 
ggplot(data = BC.df, aes(x = long, y = lat) + 
  geom_sf() +
  theme_bw() +
  geom_polygon(data = BC.df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.1, fill = 'grey95') +
  coord_cartesian(xlim = c(-128.18, -128), ylim = c(52.11, 52.2)) +
  #geom_point(data=EXPTsites, aes(x=long, y=lat, shape=otter), size=4, 
  #colour="blue", stroke=1.5) +  
  #add this to plot site locations
  scale_y_continuous(breaks = c(52.12, 52.15, 52.18)) +
  scale_x_continuous(breaks = c(-128.15, -128.09, -128.03)) +
  north(data = BC.df, scale = 0.05, symbol = 3, 
        anchor = c(x = -127.999, y = 52.218)) +
  theme(panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA),
        axis.title.y= element_blank(), axis.title.x = element_blank(),
        axis.text.y= element_text(size=10), axis.text.x = element_text(size=10)) +
  ggsn::scalebar(BC.df, dist = 1, transform = TRUE, model = 'WGS84', 
                 dist_unit = "m")

