### EXPLORING ANALYSIS OF CTMM MODELS ###

##################-
#### 1. Setup #####
##################-

###  LIBRARIES ###

library(data.table)
library(tidyverse)
library(ctmm)
library(sf)
library(ggplot2)
#for parallel processing
library(parallel)
library(foreach)
library(doParallel)

### DIRECTORIES ###
ctmm_path = "./data/ctmm_fits/"
data_filter_path = "./data/tracks_filtered/"
telem_path = "./data/telem_obj/"
akdes_path = "./data/akdes/"
lake_polygon_path = "./data/lake_coords/"

# >>> 1.  Load muddyfoot datasets and ctmm models ####

#muddyfoot tracks
muddyfoot_sub <- readRDS(paste0(data_filter_path, 'muddyfoot_sub.rds'))

#isolate pike
pike_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Northern Pike') %>% 
  as.data.frame()

#isolate perch
perch_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Perch') %>% 
  as.data.frame()

#check treatment order
perch_muddyfoot %>% 
  select(individual_ID, Treatment) %>% 
  distinct()
#[1:15] - Control
#[16:30] - Mix (i.e Exposed)

#isolate roach
roach_muddyfoot <- muddyfoot_sub %>% 
  dplyr::filter(Species == 'Roach') %>% 
  as.data.frame()

roach_muddyfoot %>% 
  select(individual_ID, Treatment) %>% 
  distinct()
#[1:15] - Control
#[16:29] - Exposed

#telemetry objects
pike_muddyfoot_tel = readRDS(paste0(telem_path, 'pike_muddyfoot_tel.rds'))
perch_muddyfoot_tel = readRDS(paste0(telem_path, 'perch_muddyfoot_tel.rds'))
roach_muddyfoot_tel = readRDS(paste0(telem_path, 'roach_muddyfoot_tel.rds'))

#ctmm lists
pike_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_pike_fits/muddyfoot_pike_OUF_models.rds"))
perch_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_perch_fits/muddyfoot_perch_OUF_models.rds"))
roach_muddyfoot_ctmm_fits = readRDS(paste0(ctmm_path, "muddyfoot_roach_fits/muddyfoot_roach_OUF_models.rds")) 

# >>> 2. Seperate ctmm fits into species by treatment groups ####

### Pike ###
mf_pike_control_fits <- pike_muddyfoot_ctmm_fits[1:3]
mf_pike_mix_fits <- pike_muddyfoot_ctmm_fits[4:6]
# Combining 'control' and 'mix' back into a 'total' list
mf_pike_total_fits <- list(mf_pike_control_fits = mf_pike_control_fits, mf_pike_mix_fits = mf_pike_mix_fits)

### Perch ###
mf_perch_control_fits <- perch_muddyfoot_ctmm_fits[1:3]
mf_perch_mix_fits <- perch_muddyfoot_ctmm_fits[4:6]
# Combining 'control' and 'mix' back into a 'total' list
mf_perch_total_fits <- list(mf_perch_control_fits = mf_perch_control_fits, mf_perch_mix_fits = mf_perch_mix_fits)

### Roach ###
mf_roach_control_fits <- roach_muddyfoot_ctmm_fits[1:3]
mf_roach_mix_fits <- roach_muddyfoot_ctmm_fits[4:6]
# Combining 'control' and 'mix' back into a 'total' list
mf_roach_total_fits <- list(mf_roach_control_fits = mf_roach_control_fits, mf_roach_mix_fits = mf_roach_mix_fits)


############################### #
#### 2. AKDE meta-analysis ######
############################### #

#Load in muddyfoot bounding box
muddyfoot_polygon = sf::st_read(paste0(lake_polygon_path, "lake_muddyfoot_polygon.gpkg"))
plot(muddyfoot_polygon)

#Calculate AKDEs

#Setting boundry with 'sp' function
#Need telemetry object
pike_akdes <- akde(pike_muddyfoot_tel, pike_muddyfoot_ctmm_fits, weights=T, sp = muddyfoot_polygon)
#This function calculates individual and population-level autocorrelated kernel density home range estimates (AKDEc)
#sp argument forces hard boundries based on a spatial object

#saveRDS(pike_akdes, file= paste0(akdes_path, "pike_muddyfoot_akdes.rds"))
#pike_akdes = readRDS(paste0(akdes_path, "pike_muddyfoot_akdes.rds"))

#perch
#may need to run again on subset of perch with irregular sampling intervals
perch_akdes <- akde(perch_muddyfoot_tel, perch_muddyfoot_ctmm_fits, weights=F, sp = muddyfoot_polygon)
saveRDS(perch_akdes, file= paste0(akdes_path, "perch_muddyfoot_akdes.rds"))

#roach
roach_akdes <- akde(roach_muddyfoot_tel, roach_muddyfoot_ctmm_fits, weights=F, sp = muddyfoot_polygon)
saveRDS(roach_akdes, file= paste0(akdes_path, "roach_muddyfoot_akdes.rds"))


#Below is code to calculate it all on a consistent grid if necessary
# calculate AKDES on a consistent grid
# estimate aKDE home ranges for each individual 
# first estimate aKDE for one individual whos HR is representative enough to use as a reference for the grid size for all others
# names(pike_muddyfoot_tel)
# pike_muddyfoot_tel[[1]]
# 
# akde_ref <- akde(pike_muddyfoot_tel[[1]], muddyfoot_pike_fits[[1]], weights = FALSE)
# 
# #Estimate AKDE's of all individuals
# pike_akdes <- list()
# for(i in 1:length(pike_muddyfoot_tel))
# {
#   print(names(pike_muddyfoot_tel)[i])
#   pike_akdes[[i]] <- akde(pike_muddyfoot_tel[[i]],muddyfoot_pike_fits[[i]],weights=FALSE,
#                               grid = list(dr = akde_ref$dr, 
#                                           align.to.origin = T))
# }

#Effective sample sizes
# Loop through each element and print the value of DOF.H
for (i in names(pike_akdes)) {
  # Construct the full name to access DOF.H
  value <- pike_akdes[[i]]$DOF.H
  
  # Print the element name and the value of DOF.H
  print(paste("DOF.H for", i, ":", value))
}

#Separating the treatment groups
pike_akde_control <- pike_akdes[1:3]
pike_akde_mix <- pike_akdes[4:6]

# Combining 'control' and 'mix' back into a 'total' list
pike_akdes_total <- list(pike_akde_control = pike_akde_control, pike_akde_mix = pike_akde_mix)

#Mean home-range for each treatment
meta(pike_akdes_total,col='black',sort=F, verbose = T, level.UD = 0.95)

COL <- color(pike_akdes, by='individual')

# color to be spatially distinct
COL_control <- COL[1:3]

# color to be spatially distinct
COL_mix <- COL[4:6]

#Calculate extent
EXT <- extent(pike_akdes_total,level=0.95)

# plot AKDEs
#control fish
plot(pike_akde_control,
     col.DF = COL_control,
     col.level=COL_control,
     col.grid=NA,
     level=NA,
     main="control AKDE", 
     level.UD = 0.95,
     xlim=EXT$x,
     ylim=EXT$y)

#exposed fish
plot(pike_akde_mix,
     col.DF=COL_mix,
     col.level=COL_mix,
     col.grid=NA,
     level=NA,
     main="Mix AKDE", 
     level.UD = 0.95,
     xlim=EXT$x,
     ylim=EXT$y)




#####################################-
#### 2. Meta-analysis - diffusion ####
#####################################-

#by individuals
meta(mf_pike_total_fits,variable="diffusion",level=0.95,level.UD=0.95,
     method="MLE",IC="AICc",boot=FALSE, error=0.01,debias=TRUE,
     verbose=T,units=TRUE,plot=TRUE,sort=FALSE,mean=T,col="black")

#mean here represents that mean diffution of all fish (or the population diffusion)

#by treatment groups
meta(mf_pike_total_fits,variable="diffusion",level=0.95,level.UD=0.95,
     method="MLE",IC="AICc",boot=FALSE, error=0.01,debias=TRUE,
     verbose=T,units=TRUE,plot=TRUE,sort=FALSE,mean=T,col="black")

#diffusion refers to the rate at which an animal location spreads out over time due to random movement
#It is used as a measure to esimate the uncertainty of animal location as time progresses
#Biologically, a higher diffusion rate indicates that the animal tends to move more unpredictably and cover more ground 


meta(mf_pike_total_fits,variable="speed",level=0.95,level.UD=0.95,
     method="MLE",IC="AICc",boot=FALSE, error=0.01,debias=TRUE,
     verbose=T,units=TRUE,plot=TRUE,sort=FALSE,mean=T,col="black")

meta(mf_pike_total_fits,variable="area",level=0.95,level.UD=0.95,
     method="MLE",IC="AICc",boot=FALSE, error=0.01,debias=TRUE,
     verbose=T,units=TRUE,plot=TRUE,sort=FALSE,mean=T,col="black")


### Plot for each treatment ###
#Separating the treatment groups
# Separating into 'control' and 'mix'
pike_control_tel <- pike_muddyfoot_tel[1:3]
pike_mix_tel <- pike_muddyfoot_tel[4:6]

# Combining 'control' and 'mix' back into a 'total' list
pike_total_tel <- list(pike_control_tel = pike_control_tel, pike_mix_tel = pike_mix_tel)

# Initialize a layout for the plots (optional)
par(mfrow = c(1, 2))  # This will arrange plots in a 2x3 grid

# Plot for each individual
for (i in 1:length(pike_mud_fits)) {
  plot(pike_total_tel[[i]], pike_mud_fits[[i]], error = FALSE, main = paste("Plot for treatment", i))
}
par(mfrow = c(1, 1))

#######################-
#### 3, Encounters rates ####
#######################-

help("proximity")
#Pairwise separation distances
DISTS <- distances(pike_muddyfoot_tel[c("F59880","F59884")],
                   pike_mud_fits[c("F59880","F59884")])


#Visualise the separation distances
plot(DISTS$est ~ DISTS$timestamp,
     type = "l",
     col = "#5e548e")
#not sure what unit the y-axis is in. Im guessing metres


#Empirical encounters per day
DISTS$encounter <- ifelse(DISTS$est <= 1, 1, 0)
#According to this code, an encounter occurs when the two individuals are less the 1 metre from each other

#Visualise the results
par(mfrow = c(1,1))
plot(DISTS$encounter ~ DISTS$timestamp)
cdplot(as.factor(DISTS$encounter) ~ DISTS$timestamp)
#conditional density plot

#Empirical Encounter rate (n/day
# Taking difference so that it only counts it as one encounter if it enters into that radius and stays there
unique_encounter <- diff(DISTS$encounter)
unique_encounter <- length(which(unique_encounter == 1))
t <- "day" %#% (DISTS$t[nrow(DISTS)] - DISTS$t[1])
unique_encounter_per_day <- unique_encounter/t
unique_encounter_per_day
#These individuals encountered each other close to 10 times a day. 





##########################-
#### 5. Encounter distribution ####
#########################-

#Plot the data and HR estimates
plot(pike_muddyfoot_tel[c("F59880", "F59884")],
     UD=pike_akdes[c("F59880", "F59884")], error = FALSE,
     col = NA,
     col.DF=c("#f4a261", "#2a9d8f"),
     col.grid = NA)


#Estimate the home range overlap
(over <- ctmm::overlap(pike_akdes[c("F59880", "F59884")]))
#This function calculates a measure of similarity between distributions known as the Bhattacharyya coefficient
#Here we have the overlsop of their autocorrelated kernel density estimates
#DOF referes to the effective sample size
#Values close to 1 mean that their distribution was very similar - perhaps not surprising in a small lake


#Estimate the CDE
CDE <- encounter(pike_akdes[c("F59880", "F59884")])
summary(CDE)

#Visualise the CDE
ctmm::plot(UD = CDE,
     col.DF="red",
     pike_muddyfoot_tel[c("F59880", "F59884")],
     col= c("#e76f51", "#264653"),
     error = FALSE,
     col.grid = NA)
#This code doesn't work

########################### -
#### 6. Population-level AKDE's
############################# -

#Put shapefile in with so argument

PKDE_control <- pkde(pike_control_tel, pike_akde_control, sp = muddyfoot_polygon)
#calculates individual and population-level autocorrelated kernel density home range estimates

#Background map
library(raster)
library(move2)


muddyfoot_pike_mv <- mt_as_move2(pike_muddyfoot, 
                            coords = c("Long","Lat"),
                            crs = "WGS84",
                            time_column = "timestamp",
                            track_id_column = "individual_ID",
                            na.fail = F) # allows or not empty coordinates


#Get outline of lake and bound akde by this
fish1 <- filter_track_data(muddyfoot_pike_mv, .track_id = "59880")
muddyfoot_map <- mapview::mapView(mt_track_lines(fish1)$geometry)
muddyfoot_bbox <-  mapedit::drawFeatures(map = muddyfoot_map)


muddyfoot_background <- 
  pike_muddyfoot %>% 
  ungroup() %>% 
  # transform to EPSG 3857 to match basemap
  st_transform(crs = st_crs(3857)) %>% 
  ggplot() + 
  basemap_gglayer(ext = muddyfoot_bbox,
                  map_service = "mapbox", 
                  map_type = "satellite",
                  map_token = "pk.eyJ1IjoibXNtaTAwMDUiLCJhIjoiY20wMXN5a2dzMHgxNzJyb3FhZ2NwbHo4YyJ9.lzWVygEPwdT4vfLEnlhuEw") +
  geom_path(data = pike_muddyfoot, 
            aes(geometry = geometry,
                x = after_stat(x), 
                y = after_stat(y),), 
            alpha = 0, 
            stat = "sf_coordinates") + 
  scale_fill_identity() + coord_sf() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_light()

# Save the plot
ggsave(plot=muddyfoot_background, "muddyfoot_background.tiff", device = "tiff")

# Create a StackedRaster object from the saved plot
raster <- raster("muddyfoot_background.tiff") # OR stack("my_ggplot.tiff") for colored images

# Get the GeoSpatial Components
lat_long <- ggplot_build(muddyfoot_background)$layout$panel_params[[1]][c("x.range","y.range")] 

# Supply GeoSpatial  data to the StackedRaster 
raster::extent(raster) <- c(lat_long$x.range,lat_long$y.range)
projection(raster) <- CRS("+proj=longlat +datum=WGS84")


Plot PKDE
## Note: can export these as a raster - use in models for roach and perch in RSFs ## as Probability mass function
```{r}
plot(pike_control,PKDE_control,main="Control PKDE", error = FALSE, R = raster)
```

TO DO: tomorrow:

# Making Raster of lake with habitat vs. no habitat # Could also do depth?
```{r}
library(terra)
library(sf)

muddyfoot_bbox <-  mapedit::drawFeatures(map = muddyfoot_map)
st_centroid(muddyfoot_bbox)
polyProj <- st_transform(muddyfoot_bbox, crs="EPSG:32633")

r <- rast(ext(polyProj), res=2, crs="EPSG:32633", vals=0)
r_latlong <- terra::project(r, "EPSG:4326")

#habitat
#Habitat1 = 63.771120, 20.048381
#Habitat2 = 63.771092, 20.048189
habitats <- data.frame(type=c("habitat1","habitat2"),
                       x=c(20.048381, 20.048189),
                       y=c(63.771120, 63.771092))

#r_latlong[cellFromXY(r_latlong, st_coordinates(habitats$geometry))] <- 1
r_latlong <- rasterize(as.matrix(habitats[,c("x","y")]), r_latlong, value=1, background=0)
r_latlong <- mask(r_latlong, muddyfoot_bbox)

plot(r_latlong)

#Need to convert it from terra to raster
r_latlong_r <- raster::raster(r_latlong)
plot(r_latlong_r)
```

RSF (not working - points outside of raster)
```{r}
#' The integrator = "Riemann" option is still in testing, we use it here because it is much faster
control_rsf_riemann <- list()

for(i in 1:length(pike_control))
{
  control_rsf_riemann[[i]] <- rsf.fit(pike_control[[i]], control_AKDES[[i]], R = list(habitat=r_latlong_r), integrator = "Riemann", trace = 3)
}

summary(control_rsf_riemann[[1]])

control_rsf <- mean(control_rsf_riemann)
summary(control_rsf)
```

- Figure out tracks on a map 
```{r}
library(basemaps)

fish1 <- filter_track_data(pike_muddyfoot_mv, .track_id = "H170-1802-59880")
#Creating map of the lake
muddyfoot_map <- mapview::mapView(mt_track_lines(fish1)$geometry)

muddyfoot_bbox <-  mapedit::drawFeatures(map = muddyfoot_map)
# plot static map
plot_tracks <- Pike_muddyfoot %>% ungroup() %>% 
  # transform to EPSG 3857 to match basemap
  st_transform(crs = st_crs(3857)) %>% 
  ggplot() +
  basemap_gglayer(ext = muddyfoot_bbox,
                  map_service = "esri", 
                  map_type = "world_imagery",
                  map_res = 1, dpi = 1200) +
  geom_path(data = Pike_muddyfoot, aes(geometry = geometry,
                                       x = after_stat(x), y = after_stat(y),
                                       color = individual_id, group = Treatment), alpha = 0.5, stat = "sf_coordinates") + 
  facet_wrap(~ individual_id) +
  scale_fill_identity() + coord_sf() +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_light() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", color = "#63605F"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(linewidth = 2, fill = '#395B5F'))
```

- Figure out how to code up the distance matrix
- calculate daily distance travelled & maximum displacement
- Recursion
- Run akde's in list specifying same grid (use code from maina - recalculate overlap)
- Bound AKDE's by the lake outline


### Perch swim distance and speed ###
```{r}
perch_muddyfoot <- muddyfoot_sub %>%
  filter(Species == 'Perch')

perch_muddyfoot_mv <- mt_as_move2(perch_muddyfoot, time_column = 'timestamp', track_id_column = 'individual_id')
```

Trajectory analysis
```{r}
# Annotate speed, azimuth and turning angle to the trajectory.
perch_muddyfoot_mv <- perch_muddyfoot_mv %>% mutate(azimuth = mt_azimuth(perch_muddyfoot_mv), 
                                                    speed = mt_speed(perch_muddyfoot_mv), 
                                                    turnangle = mt_turnangle(perch_muddyfoot_mv),
                                                    distance = mt_distance(perch_muddyfoot_mv))
head(perch_muddyfoot_mv)

#Change units
perch_muddyfoot_mv$speed <- units::set_units(perch_muddyfoot_mv$speed, cm/s)
perch_muddyfoot_mv$distance <- units::set_units(perch_muddyfoot_mv$distance, cm)
str(perch_muddyfoot_mv$speed)

perch_muddyfoot_mv <- as.data.frame(perch_muddyfoot_mv)

#Max speeds for perch are 0.977 m/s. Need to filter out observations greater than this
perch_muddyfoot_mv <- perch_muddyfoot_mv %>% dplyr::filter(speed <= units::set_units(97.7, "cm/s"))

#summarise average daily speed

#Get the date
perch_muddyfoot_mv$date <- strftime(perch_muddyfoot_mv$timestamp, format="%Y/%m/%d")

#Summarise per day
perch_muddyfoot_daily <- perch_muddyfoot_mv %>%
  dplyr::group_by(individual_id, date) %>%
  dplyr::mutate(
    avg_daily_speed = mean(speed),
    total_distance = sum(distance)) %>% 
  dplyr::distinct(individual_id, date, .keep_all = TRUE) %>%
  ungroup()

str(perch_muddyfoot_daily)
```

Plot
```{r}
ggplot(perch_muddyfoot_daily, aes(x = factor(Treatment), y = avg_daily_speed, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Mean daily speed",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
head(perch_muddyfoot_daily)
perch_muddyfoot_daily$Date <- as.factor(perch_muddyfoot_daily$date)

perch_speed <- lmer(avg_daily_speed ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = perch_muddyfoot_daily)

performance::check_model(perch_speed)

summary(perch_speed)
# 2. Obtain estimated means
perch_speed_means <- estimate_means(perch_speed)


# 3. Obtain estimated contrasts
perch_speed_contrast <- estimate_contrasts(perch_speed, contrast = "Treatment", p_adjust = "tukey")

perch_speed_means
perch_speed_contrast
```

Plot
```{r}
perch_muddyfoot_daily$total_distance <- units::set_units(perch_muddyfoot_daily$total_distance, m)

ggplot(perch_muddyfoot_daily, aes(x = factor(Treatment), y = total_distance, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Total distance",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
perch_dist <- lmer(total_distance ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = perch_muddyfoot_daily)

# 2. Obtain estimated means
perch_dist_means <- estimate_means(perch_dist)

# 3. Obtain estimated contrasts
perch_dist_contrast <- estimate_contrasts(perch_dist, contrast = "Treatment", p_adjust = "tukey")

perch_dist_means
perch_dist_contrast
```

### Perch swim distance and speed ###
```{r}
roach_muddyfoot <- muddyfoot_sub %>%
  filter(Species == 'Roach')

roach_muddyfoot_mv <- mt_as_move2(roach_muddyfoot, time_column = 'timestamp', track_id_column = 'individual_id')
```

Trajectory analysis
```{r}
# Annotate speed, azimuth and turning angle to the trajectory.
roach_muddyfoot_mv <- roach_muddyfoot_mv %>% mutate(azimuth = mt_azimuth(roach_muddyfoot_mv), 
                                                    speed = mt_speed(roach_muddyfoot_mv), 
                                                    turnangle = mt_turnangle(roach_muddyfoot_mv),
                                                    distance = mt_distance(roach_muddyfoot_mv))

#Change units
roach_muddyfoot_mv$speed <- units::set_units(roach_muddyfoot_mv$speed, cm/s)
roach_muddyfoot_mv$distance <- units::set_units(roach_muddyfoot_mv$distance, cm)
str(perch_muddyfoot_mv$speed)

roach_muddyfoot_mv <- as.data.frame(roach_muddyfoot_mv)

#Max speeds for perch are 0.977 m/s. Need to filter out observations greater than this
roach_muddyfoot_mv <- roach_muddyfoot_mv %>% dplyr::filter(speed <= units::set_units(84.1, "cm/s"))

#summarise average daily speed

#Get the date
roach_muddyfoot_mv$date <- strftime(roach_muddyfoot_mv$timestamp, format="%Y/%m/%d")

#Summarise per day
roach_muddyfoot_daily <- roach_muddyfoot_mv %>%
  dplyr::group_by(individual_id, date) %>%
  dplyr::mutate(
    avg_daily_speed = mean(speed),
    total_distance = sum(distance)) %>% 
  dplyr::distinct(individual_id, date, .keep_all = TRUE) %>%
  ungroup()
```

Plot
```{r}
ggplot(roach_muddyfoot_daily, aes(x = factor(Treatment), y = avg_daily_speed, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Mean daily speed",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
head(perch_muddyfoot_daily)
roach_muddyfoot_daily$Date <- as.factor(roach_muddyfoot_daily$date)

roach_speed <- lmer(avg_daily_speed ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = roach_muddyfoot_daily)

performance::check_model(perch_speed)

summary(roach_speed)
# 2. Obtain estimated means
roach_speed_means <- estimate_means(roach_speed)


# 3. Obtain estimated contrasts
roach_speed_contrast <- estimate_contrasts(roach_speed, contrast = "Treatment", p_adjust = "tukey")

roach_speed_means
roach_speed_contrast
```

Plot
```{r}
roach_muddyfoot_daily$total_distance <- units::set_units(roach_muddyfoot_daily$total_distance, m)

ggplot(roach_muddyfoot_daily, aes(x = factor(Treatment), y = total_distance, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Total distance",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
roach_dist <- lmer(total_distance ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = roach_muddyfoot_daily)

summary(roach_dist)
# 2. Obtain estimated means
roach_dist_means <- estimate_means(roach_dist)


# 3. Obtain estimated contrasts
roach_dist_contrast <- estimate_contrasts(roach_dist, contrast = "Treatment", p_adjust = "tukey")

roach_dist_means
roach_dist_contrast
```

### Perch swim distance and speed ###
```{r}
pike_muddyfoot <- muddyfoot_sub %>%
  filter(Species == 'Northern Pike')
pike_muddyfoot_mv <- mt_as_move2(pike_muddyfoot, time_column = 'timestamp', track_id_column = 'individual_id')
```

Trajectory analysis
```{r}
# Annotate speed, azimuth and turning angle to the trajectory.
pike_muddyfoot_mv <- pike_muddyfoot_mv %>% mutate(azimuth = mt_azimuth(pike_muddyfoot_mv), 
                                                  speed = mt_speed(pike_muddyfoot_mv), 
                                                  turnangle = mt_turnangle(pike_muddyfoot_mv),
                                                  distance = mt_distance(pike_muddyfoot_mv))

#Change units
pike_muddyfoot_mv$speed <- units::set_units(pike_muddyfoot_mv$speed, cm/s)
pike_muddyfoot_mv$distance <- units::set_units(pike_muddyfoot_mv$distance, m)


#Max speeds for perch are 0.977 m/s. Need to filter out observations greater than this
pike_muddyfoot_mv <- pike_muddyfoot_mv %>% dplyr::filter(speed <= units::set_units(82.3, "cm/s"))


pike_muddyfoot_mv <- as.data.frame(pike_muddyfoot_mv)

#summarise average daily speed

#Get the date
pike_muddyfoot_mv$date <- strftime(pike_muddyfoot_mv$timestamp, format="%Y/%m/%d")

#Summarise per day
pike_muddyfoot_daily <- pike_muddyfoot_mv %>%
  dplyr::group_by(individual_id, date) %>%
  dplyr::mutate(
    avg_daily_speed = mean(speed),
    total_distance = sum(distance)) %>% 
  dplyr::distinct(individual_id, date, .keep_all = TRUE) %>%
  ungroup()
```

Plot
```{r}
ggplot(pike_muddyfoot_daily, aes(x = factor(Treatment), y = avg_daily_speed, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Mean daily speed",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
pike_muddyfoot_daily$Date <- as.factor(pike_muddyfoot_daily$date)

pike_speed <- lmer(avg_daily_speed ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = pike_muddyfoot_daily)

performance::check_model(perch_speed)

summary(pike_speed)
# 2. Obtain estimated means
pike_speed_means <- estimate_means(pike_speed)


# 3. Obtain estimated contrasts
pike_speed_contrast <- estimate_contrasts(pike_speed, contrast = "Treatment", p_adjust = "tukey")

pike_speed_means
pike_speed_contrast
```
```{r}
pike_muddyfoot_daily %>% 
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(mean = mean(avg_daily_speed, na.rm = T),
                   sd = sd(avg_daily_speed, na.rm = T),
                   max = max(avg_daily_speed, na.rm = T),
                   min = min(avg_daily_speed, na.rm = T),
                   num_unique_ID = length(unique(individual_id)),
                   num_data_pts = sum(!is.na(avg_daily_speed)),
                   mean_datapts_ID = num_data_pts/num_unique_ID) %>%
  ungroup()
```

Plot
```{r}
ggplot(pike_muddyfoot_daily, aes(x = factor(Treatment), y = total_distance, fill = factor(Treatment))) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.75,
    # move to the right
    justification = -0.1,
    # remove the slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Total distance",
  ) +
  coord_flip()+
  theme(legend.position="none")
```

linear mixed model
```{r}
perch_dist <- lmer(total_distance ~ Treatment + scale(Std_length) + (1|Date) + (1|individual_id), data = pike_muddyfoot_daily)

summary(perch_dist)
# 2. Obtain estimated means
pike_dist_means <- estimate_means(perch_dist)


# 3. Obtain estimated contrasts
pike_dist_contrast <- estimate_contrasts(perch_dist, contrast = "Treatment", p_adjust = "tukey")

pike_dist_means
pike_dist_contrast
```













```