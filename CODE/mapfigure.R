#Map creation for Beetle Manuscript with Nicola Rammell
#Created 2 March 2020

####PACKAGE LOADING####
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)

####MAP CODE####
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-135.0, -126.0), ylim = c(50.0, 54.0), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)
