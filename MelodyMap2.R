## This script create a map of the data. 

library(ggplot2)
library(dplyr)
library(sf) 
library(rnaturalearth) 
library(patchwork)

#### Data ####
data = read.csv('MelodicEvoSeq.csv')

song_frequency = data %>%
  filter(!is.na(NAME_1)) %>%
  group_by(NAME_1) %>% 
  summarise(n_songs = n())

orange_colours = c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", 
                   "#FD8D3C", "#F16913", "#D94801", "#A63603",
                   "#7F2704")

#### United Kingdom + Ireland ####
uk_sf <- ne_states(country = c("united kingdom", "ireland"), 
                   returnclass = "sf")

# If the administrator is Ireland, the geonunit is Ireland
uk_sf$geonunit[uk_sf$admin == "Ireland"] = "Ireland"

uk_sf = left_join(uk_sf, song_frequency, by = c("geonunit" = "NAME_1"))

# We have no songs from Wales, 
# but indicate count as 0 rather than missing
uk_sf$n_songs[uk_sf$geonunit == "Wales"] = 0 

hist(uk_sf$n_songs)
uk_sf$nsongs_discrete = cut(uk_sf$n_songs, 5)

uk_ireland = ggplot() + 
  geom_sf(data = uk_sf, aes(fill = n_songs), color = NA) +
  theme_minimal() +
  scale_fill_gradientn(colors = orange_colours)

#### Japan ####
japan_sf <- ne_states(country = "japan", 
                   returnclass = "sf")

japan_sf = left_join(japan_sf, song_frequency, by = c("name" = "NAME_1"))

# Code areas with no songs as 0 rather than missing
japan_sf$n_songs[is.na(japan_sf$n_songs)] = 0

japan = ggplot() + 
  geom_sf(data = japan_sf, aes(fill = n_songs), color = NA) +
  theme_minimal() +
  xlim(c(127, 150)) + 
  ylim(c(30, 46)) + 
  scale_fill_gradientn(colors = orange_colours)

#### USA ####
usa_sf <- ne_states(country = "United States of America", returnclass = "sf")

usa_sf = left_join(usa_sf, song_frequency, by = c("name" = "NAME_1"))

# Code areas with no songs as 0 rather than missing
usa_sf$n_songs[is.na(usa_sf$n_songs)] = 0

usa = ggplot() + 
  geom_sf(data = usa_sf, aes(fill = n_songs), color = NA) +
  theme_minimal() +
  coord_sf(crs = st_crs(2163), xlim = c(-2500000, 2500000), 
           ylim = c(-2300000, 730000)) + 
  scale_fill_gradientn(colors = orange_colours)

ggsave(filename = 'figures/ukireland_map.jpeg', plot = uk_ireland)
ggsave(filename = 'figures/japan_map.jpeg', plot = japan)
ggsave(filename = 'figures/usa_map.jpeg', plot = usa)
