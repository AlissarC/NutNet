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
library(ggstatsplot)

library(tidyverse)
library(dplyr)
library(ggstatsplot)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(nlme)
library(marginaleffects)
library(piecewiseSEM)
library(rstantools)
library(multcomp)
library(treemapify)
library(relaimpo)
library(r2glmm)
library(patchwork)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(MuMIn)
library(boot)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(lavaan)
library(patchwork)


################################################# Processing  ##########################################
########################################################################################################
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
traits <- read.csv("df1.csv")
names(traits)
nrow(traits)
traits$aridity

traits$Ntrt_fac <- as.factor(traits$Ntrt)
traits$Ptrt_fac <- as.factor(traits$Ptrt)
traits$Ktrt_fac <- as.factor(traits$Ktrt)
traits$block_fac <- as.factor(traits$block)
traits$trt_fac <- as.factor(traits$trt)

## add in photosynthetic pathway information
levels(as.factor(traits$Family)) # check Amaranthaceae, Asteraceae, Boraginaceae, Caryophyllaceae, Cyperaceae, Euphorbiaceae,
# Polygonaceae, Poaceae, Scrophulariaceae
# only one!
# C4
#traits$photosynthetic_pathway ='NULL'
traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'
                              & traits$Family == 'Cyperaceae' & traits$Taxon == 'FIMBRISTYLIS DICHOTOMA'] <- 'C4'

traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'] <- 'C3'

table(traits$photosynthetic_pathway)

### calculate SLA #############
traits$sla_m2_g = traits$SLA_v2 * (1/1000000)
hist(traits$sla_m2_g)

### calculate LMA #############
traits$lma = 1/traits$sla_m2_g
hist(traits$lma) # some extremely high values

### calculate narea #############
traits$narea = (traits$leaf_pct_N / 100) * (traits$lma)
hist(traits$narea) # some extremely high values

### calculate nmass #############
traits$nmass = (traits$leaf_pct_N)
hist(traits$nmass)

### calculate N: P ratio in leaves #############
traits$leaf_pct_P = traits$leaf_ppm_P/10000 ## transform into percentage
traits$leaf_N_P = traits$leaf_pct_N/traits$leaf_pct_P

### calculate C: N ratio in leaves #############
traits$leaf_N_C = traits$leaf_pct_C/traits$leaf_pct_N

### calculate lai per plot #############
traits$lai = -log(traits$Ground_PAR / traits$Ambient_PAR) / 0.86 # from: http://manuals.decagon.com/Manuals/10242_Accupar%20LP80_Web.pdf page 41
hist(traits$lai)

### calculate par per leaf area to assume par absorbed is reduced in dense canopies: calculation per plot
traits$par_per_leaf_area = traits$par2_gs * ((1 - exp(-0.5 * traits$lai)) / traits$lai) # from Dong et al. (2007) eqn 2
hist(traits$par_per_leaf_area)

### calculate biomass for each species depending on the percentage max-cover for each species ###############
traits$spp_live_mass = traits$live_mass * (traits$max_cover / 100)
hist(traits$spp_live_mass)


## calculate big delta C13 from small delta 
traits$delta = ((-0.008 - traits$leaf_C13_delta_PDB * 0.001) / (1 + traits$leaf_C13_delta_PDB * 0.001)) * 1000
hist(traits$delta)

##### calculate chi for C3 plants #################### 
traits$chi[traits$photosynthetic_pathway == 'C3'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C3'] * 0.001 - 0.0044) / (0.027 - 0.0044)
hist(traits$chi)

##### calculate chi for C4 plants #################### 
traits$chi[traits$photosynthetic_pathway == 'C4'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C4'] * 0.001 - 0.0044) / ((-0.0057 + 0.03*0.4) - 0.0044)


hist(traits$chi)

#### create a new column for fence #################
traits$fence <- 'no'
traits$fence[traits$trt == 'Fence' | traits$trt == 'NPK+Fence'] <- 'yes'

################### transform some variables into log to get normal distribution #####################################
traits$logpar2_gs <- log(traits$par2_gs)

traits$logpar_per_leaf_area <- log(traits$par_per_leaf_area)
hist(traits$logpar_per_leaf_area)

traits$loglma <- log(traits$lma)
hist(traits$loglma)

traits$lognarea <- log(traits$narea)
hist(traits$lognarea)

traits$lognmass <- log(traits$nmass)
hist(traits$lognmass)

names(traits)
traits$log_spp_live_mass <- log(traits$spp_live_mass)
hist(traits$log_spp_live_mass)

############################### calculate Beta from chi ################################################################################
#### load functions ####
############################
### functions to calculate vcmax and jmax 
setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper")
source("optimal_vcmax_R/calc_optimal_vcmax.R")
sourceDirectory("optimal_vcmax_R/functions", modifiedOnly = FALSE)

traits$nstar <- calc_nstar(traits$tmp2_gs, traits$z)

traits$gammastar <- calc_gammastar_pa(traits$tmp2_gs, traits$z)

traits$km <- calc_km_pa(traits$tmp2_gs, traits$z)

traits$vpd_adj <- calc_vpd(traits$tmp2_gs, traits$z, traits$vpd2_gs)

traits$patm <- calc_patm(traits$z)

traits$ca <- 400 * 1e-6 * traits$patm

Kp_25 <- 16 # Pa
Ea_Kp <- 36300 # J mol^-1, for the Pa parameter ## Boyd et al 2015
traits$tmpK <- traits$tmp2_gs + 273.15
traits$kp <- Kp_25 * exp((Ea_Kp * (traits$tmpK - 298.15))/(298.15 * 8.3145 * traits$tmpK))

traits$beta_num <- 1.6 * traits$nstar * traits$vpd_adj * 1000 * ((traits$chi) ^ 2)
traits$beta_denom <- ((1 - traits$chi)^2) * (traits$km + traits$gammastar)

#### Beta denom C4 ################# 
traits$beta_denom[traits$photosynthetic_pathway == 'C4'] <- ((1 - traits$chi[traits$photosynthetic_pathway == 'C4'])^2) * 
   (traits$kp[traits$photosynthetic_pathway == 'C4'] + traits$gammastar[traits$photosynthetic_pathway == 'C4'])


################ calculate beta #########
traits$beta <- traits$beta_num / traits$beta_denom
hist(traits$beta)

traits$logbeta = log(traits$beta)
hist(traits$logbeta)

#################### create new column pft based on c3/c4, fixers or non fixers ##############
traits$pft[traits$photosynthetic_pathway == 'C4'] <- 'c4'
traits$pft[traits$photosynthetic_pathway == 'C3' & traits$Nfix == 'no'] <- 'c3_noNfix'
traits$pft[traits$photosynthetic_pathway == 'C3' & traits$Nfix == 'yes'] <- 'c3_yesNfix'

nrow(traits)# 2788
length(unique(traits$Taxon)) # 244 species 
length(unique(traits$site_code)) # 27 sites

#################### exclude C points where chi is higher than 0.3 and lower than 0.98 ######### 
traits_chi_sub <- subset(traits, chi > 0.1 & chi < 0.95)
nrow(traits_chi_sub) ### 2106
table(traits_chi_sub$pft) ### 162 C4 remains 
hist(traits_chi_sub$chi) # quasi normal distribution 
hist(traits_chi_sub$aridity) # quasi normal distribution 

nrow(traits_chi_sub) # 2106
min(traits_chi_sub$chi) # 0.11 
max(traits_chi_sub$chi) # 0.99

length(unique(traits_chi_sub$Taxon)) # 207 species 
length(unique(traits_chi_sub$site_code)) # 26 sites

################# exclusion of fence plots and subset outliers using MAD ###########################
unique(traits_chi_sub$trt)
unique(traits_chi_sub$fence)

## exclusion of fence treatments 
traits_nofence_chi_sub = subset(traits_chi_sub, fence == "no")
traits_nofence_chi_sub$fence
nrow(traits_nofence_chi_sub) ## 1752
length(unique(traits_nofence_chi_sub$Taxon)) # 208 species 
length(unique(traits_chi_sub$site_code)) # 26 sites

min(traits_chi_sub$aridity) # 0.14
max(traits_chi_sub$aridity) # 2.32

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(traits, "traits.csv")
write.csv(traits_chi_sub, "traits_chi_sub.csv")
write.csv(traits_nofence_chi_sub, "traits_nofence_chi_sub.csv")


names(traits_nofence_chi_sub)
################################## Narea predictions from theoretical vcmax ###############################################
####################################################################################################

### calculate vcmax, jmax, and vpmax with known chi for C3 plants
leaf_chi_subset_c3 <- subset(traits_nofence_chi_sub, photosynthetic_pathway == 'C3')
table(leaf_chi_subset_c3$photosynthetic_pathway) # 1607
gas_exchange_pred_c3 <- calc_optimal_vcmax(pathway = "C3",
                                           tg_c = leaf_chi_subset_c3$tmp, 
                                           paro = leaf_chi_subset_c3$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c3$vpd, 
                                           z = leaf_chi_subset_c3$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c3$chi,
                                           lma = leaf_chi_subset_c3$lma)

### calculate vcmax, jmax, and vpmax with known chi for C4 plants
leaf_chi_subset_c4 <- subset(traits_nofence_chi_sub, photosynthetic_pathway == 'C4')
table(leaf_chi_subset_c4$photosynthetic_pathway) # 145
gas_exchange_pred_c4 <- calc_optimal_vcmax(pathway = "C4",
                                           tg_c = leaf_chi_subset_c4$tmp, 
                                           paro = leaf_chi_subset_c4$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c4$vpd, 
                                           z = leaf_chi_subset_c4$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c4$chi,
                                           lma = leaf_chi_subset_c4$lma)

## add C3 model results to leaf dataset
npred_c3 <- bind_cols(leaf_chi_subset_c3, gas_exchange_pred_c3[ ,39:51])
npred_c3$model_lma <- gas_exchange_pred_c3$lma

## add C4 model results to leaf dataset
npred_c4 <- bind_cols(leaf_chi_subset_c4, gas_exchange_pred_c4[ ,39:51])
npred_c4$model_lma <- gas_exchange_pred_c4$lma

# combine C3 and C4 subsets
leaf_chi_subset_npred <- bind_rows(npred_c3, npred_c4)

leaf_chi_subset_npred$lognphoto <- log(leaf_chi_subset_npred$nphoto)
leaf_chi_subset_npred$lognstructure <- log(leaf_chi_subset_npred$nstructure)

names(leaf_chi_subset_npred)
hist(leaf_chi_subset_npred$lognphoto)
hist(leaf_chi_subset_npred$lognstructure)

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
write.csv(leaf_chi_subset_npred, "leaf_chi_subset_npred.csv")
