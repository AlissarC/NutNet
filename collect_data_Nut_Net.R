# create leaf trait dataset for analyses

#########################################
# packages necessary
#########################################
library(tidyr)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(multcompView)
library(ggplot2)

getwd()
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")

#################### load data containing leaf traits (Jen Firn) ####################################################################
#####################################################################################################################################

leaf_traits = read.csv("NutNet-foliar-traits-7JAN2017.csv")
names(leaf_traits)
nrow(leaf_traits) ## 2764
length (unique(leaf_traits$site_code)) # 27 sites
length (unique(leaf_traits$Taxon)) # 244 species 
length (unique(leaf_traits$trt)) # 10 trt
## each site contains one, two or three blocks, each block contains different plots corresponding to different trt.
## in each trt we have leaf data for different species, the number and type of species vary between plots 

##################### load data containing SLA version 2 ##################### 
leaf_v2 =  read.csv("foliar_cover_updated_2.csv")
names(leaf_v2)
nrow(leaf_v2) ## 2664
length (unique(leaf_v2$site_code)) # 27 sites
length (unique(leaf_v2$Taxon)) # 243 species 
length (unique(leaf_v2$trt)) # 10 trt

leaf_new_sla = leaf_v2[, c(1:5, 13)]  # selection of site_code, plot, year, Taxon, block, SLA_v2 
leaf = left_join(leaf_traits, leaf_new_sla) ### join by Taxon, site_code, year, block, plot 
names(leaf)

leaf[leaf == 'NULL'] <- NA ## replace cells with NULL by NA

leaf[, 8:30] <- sapply(leaf[, 8:30], as.numeric) ## set columns with measurements to numeric 
sapply(leaf, class) ## verify the data types for each column

leaf$Ntrt = 0 ## creation of a column with Ntrt= 0 
leaf$Ntrt[leaf$trt == 'N' | leaf$trt == 'NP' | leaf$trt == 'NK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1 ## set Ntrt=1 if any point receives N
### repeat the same for the other treatments 
leaf$Ptrt = 0
leaf$Ptrt[leaf$trt == 'P' | leaf$trt == 'NP' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ktrt = 0
leaf$Ktrt[leaf$trt == 'K' | leaf$trt == 'NK' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1

leaf$Nfix = 'no' ## creation of a column with Nfix as no
leaf$Nfix[leaf$Family == 'Fabaceae'] = 'yes' ## when the families are Fabaceae, set them as N fixers, the other families as non fixers

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf, "leaf.csv")
#######################################################################################################################################

#################### load data containing biomass and cover percentage ####################################################################
#####################################################################################################################################
#### This file contains information on species richness, percentage cover, biomass for each plot, PAR, etc....: not per species but per plot

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
core =  read.csv("comb-by-plot-09-April-2018.csv")
names(core) 
nrow(core) ## 14939
length (unique(core$site_code)) # 100 sites


core[core == 'NULL'] <- NA ## replace NULL by NA

core[, 26:45] <- sapply(core[, 26:45], as.numeric) # set columns with measurments as numeric
sapply(core, class)

core$par_rat = core$Ground_PAR / core$Ambient_PAR ### calculate the ratio of ground PAR over Ambiant PAR 

#### selection of the columns: site_code, year, block, plot, trt, from sum_INT_cover to the end 
core_data = core[, c(2, 25, 15, 16, 22, 26:46)]  

#### selection of the columns: site_code, block, plot, trt, site_name, lat,long annd other columns containing information
### about each site
core_info = core[, c(2, 15, 16, 22, 1, 3:14, 17:21, 23, 24)]


############## merge leaf with core data for each site/year/block/plot/trt  ##########################
leaf_core_data = join(leaf, core_data, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')

leaf_core = join(leaf_core_data, core_info, by = c("site_code", "block", "plot", "trt"), type = 'left', match = 'first')

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf_core, "leaf_core.csv")
##### here for each species in the same plot: same data of core for the whole plot, no data of biomass or PAR per species 
##########################################################################################################################


########################################### load and add climate data to leaf_core #########################################################
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
spei = read.csv("CRU-annual_2018-07-06.csv")
names(spei)   ## this file contains for each site the annual data from 1903 to 2016 of precip, PET and SPEI-6,SPEI-12, SPEI-24
table(spei$year)
spei$p_pet = spei$precip / spei$PET  ### Calculate aridity Index 

sapply(spei, class)
####### Average climate data between 1903 and 1970, and between 1970 and 2016  
spei$year <- as.numeric(spei$year)

### Spei, precip, pet for each site and year (not averaged fore 1901-2015). 
## precip_mean ans pet_mean in the file SPEI are for averages for the period 1975-2016
leaf_core_spei = join(leaf_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')
names(leaf_core_spei)
#leaf_core_spei$TimePeriod = NULL

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf_core_spei, "leaf_core_spei.csv")

# check
nrow(leaf) # 2788
nrow(leaf_core) ## 2788
nrow(leaf_core_spei) ## 2788

######################################## load growing season climate data: temperature, par, vpd, z ###################
## these data are for each longitude and latitude, I don't know for which years they are.... ########################### 
# read in gs climate (>0°C)
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
tmp_globe = read.csv("cru_tmp_climExtract_growingseason_globe.csv")
par_globe = read.csv("cru_par_climExtract_growingseason_globe.csv")
vpd_globe = read.csv("cru_vpd_climExtract_growingseason_globe.csv")
z_globe =  read.csv("z_globe.csv")

leaf_core_spei$tmp = NA
leaf_core_spei$par = NA
leaf_core_spei$vpd = NA
leaf_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(leaf_core_spei)){
  
  currentLat = leaf_core_spei$latitude[i]
  currentLon = leaf_core_spei$longitude[i]
  
  clim_comb = tmp_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  tmp = subset(clim_comb, lat==best_lat & lon==best_lon)$tmp
  
  clim_comb = par_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  par = subset(clim_comb, lat==best_lat & lon==best_lon)$par
  
  clim_comb = vpd_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  vpd = subset(clim_comb, lat==best_lat & lon==best_lon)$vpd
  
  clim_comb = z_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  z = subset(clim_comb, lat==best_lat & lon==best_lon)$z
  
  climate_df = rbind(climate_df, c(tmp, par, vpd, z, currentLat, currentLon, best_lat, best_lon, i))
  
  
}

plot(climate_df[,5] ~ climate_df[,7])
plot(climate_df[,6] ~ climate_df[,8])

leaf_core_spei$tmp = climate_df[, 1]
leaf_core_spei$par = climate_df[, 2]
leaf_core_spei$vpd = climate_df[, 3]
leaf_core_spei$z = climate_df[, 4]
names(leaf_core_spei)

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf_core_spei, "leaf_core_spei.csv")

########################################################
# add species information
########################################################
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
info = read.csv("full-species-info-25-January-2019.csv")
nrow(info) ## 5897
length(unique(info$site_code)) # 108 sites
### this file contains information about species wether they are grass or Forb, C3 or C4, perennial or annual 

code_taxon_data = paste(toupper(leaf_core_spei$site_code), toupper(leaf_core_spei$Taxon), sep = ' ')

write.csv(code_taxon_data, "code_taxon_data.csv")

info_caps = mutate_all(info, .funs = toupper) ### transform all lower cases to upper cases 
names(info_caps)
info_caps$code_taxon = paste(info_caps$site_code, info_caps$standard_taxon, sep = ' ') #bind site_code to santadard_taxon with space between them


n_info = NULL

for (i in 1:length(code_taxon_data)){
  
  ancillary_data = subset(info_caps, code_taxon == code_taxon_data[i])[, c(6:10)]
  
  ancillary_data$n = i
  
  n_info = rbind(n_info, ancillary_data)
  
}

leaf_core_spei_info = cbind(leaf_core_spei, n_info)
leaf_core_spei_info


#############################
## add in new growing season climate data (March 11, 2022)
#############################
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
gs_climate <- read.csv("cru_gs_climate.csv")
gs_climate

leaf_core_spei_info_gsclimate <- left_join(leaf_core_spei_info, gs_climate, by = c("site_code"))

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf_core_spei_info_gsclimate, "leaf_core_spei_info_gsclimate.csv")
nrow(leaf_core_spei_info_gsclimate) #2788


###################################### add in species composition data ###############################################
######################################################################################################################

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
cover = read.csv("full-cover-31-August-2020.csv")
names(cover)
nrow(cover) ## 239683
#cover = subset(cover, live==1)
#nrow(cover) ## 212131
#min(cover$live)

cover_select = select(cover, site_code, plot, year, Taxon, max_cover)

core_leaf_spei_cover = left_join(leaf_core_spei_info_gsclimate, cover_select)
nrow(core_leaf_spei_cover) # 2788


####################################### load biomass data ##############################

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
full_biomass<-read.csv("full-biomass-nutrients-06-December-2017.csv")
names(full_biomass)
length(unique(full_biomass$site_code)) ## 27
full_biomass$trt
########### for each site, block and plot, we have in this file data per category
#### forb, graminoid etc... not per species 

full_biomass[full_biomass == 'NULL'] <- NA # replace empty values by NA
full_biomass$pct_N <- as.numeric(full_biomass$pct_N)
sapply(full_biomass, class)

full_biomass$Ntrt = 0
full_biomass$Ntrt[full_biomass$trt == 'N' | full_biomass$trt == 'NP' | full_biomass$trt == 'NK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when N is added: Ntrt=1

full_biomass$Ptrt = 0
full_biomass$Ptrt[full_biomass$trt == 'P' | full_biomass$trt == 'NP' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when P is added: Ptrt=1

full_biomass$Ktrt = 0
full_biomass$Ktrt[full_biomass$trt == 'K' | full_biomass$trt == 'NK' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1 # when K is added: Ktrt=1

full_biomass$Ntrt_fac = as.factor(full_biomass$Ntrt)
full_biomass$Ptrt_fac = as.factor(full_biomass$Ptrt)
full_biomass$Ktrt_fac = as.factor(full_biomass$Ktrt)


full_biomass$pct_C = as.numeric(as.character(full_biomass$pct_C))
full_biomass$pct_N = as.numeric(as.character(full_biomass$pct_N))
full_biomass$pct_P = as.numeric(as.character(full_biomass$pct_P))
full_biomass$pct_K = as.numeric(as.character(full_biomass$pct_K))
full_biomass$pct_Mg = as.numeric(as.character(full_biomass$pct_Mg))
full_biomass$pct_Ca = as.numeric(as.character(full_biomass$pct_Ca))


## AGN = above ground N!, mass is leaf mass
full_biomass$AGN <- full_biomass$mass*(full_biomass$pct_N*0.01) # calculation of aboveground N quatity

full_biomass$Nfix = 'no'
full_biomass$Nfix[full_biomass$category == 'LEGUME'] = 'yes' # set legumes as N fixers 

nrow(full_biomass) # 2102

biomass_core_leaf_spei = join(core_leaf_spei_cover, full_biomass, by = c("site_code", "year"), type = 'left', match = 'first')
names(biomass_core_leaf_spei)
nrow(biomass_core_leaf_spei) # 2788


### add growing season length data#######################
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
gs_information <- read.csv("Weather_site_inventory_20200921.csv")
head(gs_information)
colnames(gs_information)

gs_length <- gs_information[,c(2, 20)] ## selection of the colum site_code and gs_len

biomass_core_leaf_spei <- join(biomass_core_leaf_spei, gs_length, by = c("site_code"), type = 'left', match = 'first')

biomass_core_leaf_spei$site_code[is.na(biomass_core_leaf_spei$gs_len)]
biomass_core_leaf_spei$gs_len[biomass_core_leaf_spei$site_code == 'gilb.za'] <- 6
biomass_core_leaf_spei$gs_len[biomass_core_leaf_spei$site_code == 'summ.za'] <- 8

biomass_core_leaf_spei$gs_frac <- biomass_core_leaf_spei$gs_len / 12

###### add soil moisture data#######################
soil <- read.csv("soil_data.csv")
biomass_core_leaf_spei <- join(biomass_core_leaf_spei, soil, by = c("site_code"), type = 'left', match = 'first')

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(biomass_core_leaf_spei, "biomass_core_leaf_spei.csv")
biomass_core_leaf_spei$max_cover

#################################### get PET and aridity index : mean annual variables between 1970 and 2000  ###############
#######################################################################################################################
library(sf)
library(raster)
library(tmap)
library(rgdal)
library(mapview)

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/input")
l = list.files(pattern = "v3_yr") #reading world maps containing the string mask in the crus_rasters (maps from jeff work)
l
################# loop to extract all the pixels (raster) from the world maps of PET data ############################# 
######## https://cgiarcsi.community/data/global-aridity-and-pet-database/################

l = list.files(pattern = "v3_yr")
l = c(l, list.files(pattern = "proj_elev", full.names = T)) # bring in new, projected elev data

r = list()
for(i in l){
  r1 = raster(i)
  r[[i]] = r1
  print(r1)
}

r = raster::stack(l)
crs(r) = 4326


records =  biomass_core_leaf_spei %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = crs(r))

records = st_transform(records, crs = st_crs(r))
st_crs(records)==st_crs(r)


df1 = cbind(records,
            raster::extract(r, st_coordinates(records), method = 'simple'))
df1

df1$Long = st_coordinates(df1)[,1]
df1$Lat = st_coordinates(df1)[,2]
st_geometry(df1) = NULL
df1$crs = 4326
names(df1)

colnames(df1)[colnames(df1) == "awi_pm_sr_yr"] <- "aridity"
df1$aridity= df1$aridity * 0.0001
hist(df1$aridity)
plot(df1$aridity, df1$SPEI_12)

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(df1, "df1.csv")
df1$aridity
