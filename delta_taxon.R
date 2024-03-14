#########################################
# packages necessary
#########################################
install.packages("semPlot")
install.packages("piecewiseSEM")

library(tidyr)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(multcompView)
library(ggstatsplot)
library(tidyverse)
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
library(semPlot)
library(lme4)
library(piecewiseSEM)
library(emmeans)
library(car)
library(tidyverse)
library(MuMIn)
library(multcomp)
library(multcompView)
library(nlme)
library(semPlot)
library(cowplot)


####################### script to explore and analyse delta N mass and delta above ground biomass between plots
#### that received N and plots that did not : N-control, NK-K, NP-P, NPK-PK ####################################

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/data_output")
traits_delta<- read.csv("leaf_chi_subset_npred.csv")
nrow(traits_delta) # 1752

#################################### calculate Delta Nmass per Taxon  ######################################
##### if the taxa change between treatments, deltas are calculated only for the Taxon which is kept across treatments

### subset plots which did not receive N (they can be control, K, PK, P) ################## 
traits_lowN <- subset(traits_delta, Ntrt_fac == '0')
nrow(traits_lowN) # 894
table(traits_lowN$trt)


### subset plots which receive N (they can be N, NP, NPK, NK, NPK)
traits_highN <- subset(traits_delta, Ntrt_fac == '1')
nrow(traits_highN ) # 858


#################### combine low N and high N subsets##########################
Delta <- left_join(traits_lowN, traits_highN, 
                   by = c('site_code', 
                          'block_fac','Ptrt_fac', 'Ktrt_fac', 'Taxon'))
names(Delta)
nrow(Delta) ## 894

################ calculate delta: NPK-NP, N-control, NK-K, NP-P ################
Delta$delta_nmass <- ((Delta$nmass.y - Delta$nmass.x) / Delta$nmass.x) * 100
Delta$delta_narea <- ((Delta$narea.y - Delta$narea.x) / Delta$narea.x) * 100
Delta$delta_spp_live_mass <- ((Delta$spp_live_mass.y - Delta$spp_live_mass.x) / Delta$spp_live_mass.x) * 100

#### selection  delta dataframe  ####### 

#### selection of high_N columns, but we will use only site characteristics and climate + taxon characteristics, don't use the other numeric columns concerning plots since they are for high N trt only ############
names(Delta)
delta_subset<- subset(Delta, select = c(1: 186, 368:370))
colnames(delta_subset) <- sub("\\.x$", "", colnames(delta_subset)) ## remove .x  
names(delta_subset)
nrow(delta_subset) # 894

write.csv(delta_subset, "delta_subset.csv")
names(delta_subset)
hist(delta_subset$delta_nmass) #ok
hist(delta_subset$delta_narea) # horrible
hist(delta_subset$delta_spp_live_mass) # horrible

#############  remove instances where any delta values are 2.5 times higher than the MAD ##########################################
#################################################################################################################################
nrow(delta_subset) # 894
delta_nmass_mad <- mad(delta_subset$delta_nmass, na.rm = T)
delta_nmass_mad   # 33.68
delta_nmass_median <- median(delta_subset$delta_nmass, na.rm = T)
delta_nmass_median ### 21.62
delta_nmass_mean<- mean(delta_subset$delta_nmass, na.rm = T)
delta_nmass_mean ## 27.92

deltanmass_mad_data <- subset(delta_subset,
                              delta_nmass < 3 * delta_nmass_mad &
                                delta_nmass > 3 * -delta_nmass_mad)

nrow(deltanmass_mad_data) # 524 
hist(deltanmass_mad_data$delta_nmass) # normal dist
names(deltanmass_mad_data)

######################## delta nmass analysis #########################################################
#### Hyp 1: Cold temperature and drought increase deltanmass, high PAR increases deltanmass #######################
#########################################################################################################################
names(deltanmass_mad_data)
deltanmass_lmer_aridity <- lmer(delta_nmass ~tmp + aridity + par2_gs + Nfix +  
                                  photosynthetic_pathway + Ptrt_fac*Ktrt_fac + 
                                  (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),data = deltanmass_mad_data)
Anova(deltanmass_lmer_aridity, type = 'II')
summary(deltanmass_lmer_aridity)

plot(resid(deltanmass_lmer_aridity) ~ fitted(deltanmass_lmer_aridity))

vif(deltanmass_lmer_aridity) # no colinearity
AIC(deltanmass_lmer_aridity) ## 4983.48

plot(deltanmass_lmer_aridity)
qqnorm(residuals(deltanmass_lmer_aridity))
qqline(residuals(deltanmass_lmer_aridity))

densityPlot(residuals(deltanmass_lmer_aridity))
shapiro.test(residuals(deltanmass_lmer_aridity)) 
outlierTest(deltanmass_lmer_aridity)

residuals <- resid(deltanmass_lmer_aridity)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

plot(fitted(deltanmass_lmer_aridity), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : none

r.squaredGLMM(deltanmass_lmer_aridity) ## 0.280, not terrible 

delta_nmass_model_aridity <- data.frame(Var = c('Tg', 'aridity', 'PAR', 'N fixer', 
                                        'C3/C4', 'Soil P', 'Soil K', 'Soil P x Soil K'))

delta_nmass_model_aridity$df <- as.matrix(Anova(deltanmass_lmer_aridity))[1:8, 2]
delta_nmass_model_aridity$Slope <- c(
  summary(emtrends(deltanmass_lmer_aridity, ~tmp, var = "tmp"))[1, 2],
  summary(emtrends(deltanmass_lmer_aridity, ~aridity, var = "aridity"))[1, 2],
  summary(emtrends(deltanmass_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 2],
  NA, NA, NA, NA, NA)
delta_nmass_model_aridity$SE <- c(
  summary(emtrends(deltanmass_lmer_aridity, ~tmp, var = "tmp"))[1, 3],
  summary(emtrends(deltanmass_lmer_aridity, ~aridity, var = "aridity"))[1, 3],
  summary(emtrends(deltanmass_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 3],
  NA, NA, NA, NA, NA)
delta_nmass_model_aridity$p <- as.matrix(Anova(deltanmass_lmer_aridity))[1:8, 3]
delta_nmass_model_aridity$RelImp <- as.matrix(calc.relip.mm(deltanmass_lmer_aridity)$lmg)[1:8]
delta_nmass_model_aridity$RelImp <- delta_nmass_model_aridity$RelImp * 100
delta_nmass_model_aridity

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/stat_output")
write.csv(delta_nmass_model_aridity, "delta_nmass_model_aridity.csv")

###################### # percentage increase of delta nmass in plots receiving P compared to plots not receiving P (NP-P)
(summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[2,2] - summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[1,2])/
  summary(emmeans(deltanmass_lmer_aridity, ~Ptrt_fac))[1,2]
# - 0.4415 : plots receiving P respond less to N addition than plots not receiving P 
## this effect is significant in the model 

# percentage increase of Nmass in N fixers compared to non-N fixers
(summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[2,2] - summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[1,2])/
  summary(emmeans(deltanmass_lmer_aridity, ~Nfix))[1,2]
# - 0.9136: non fixers increase in deltaNmass is 91.36% greater than  fixers 
## highly significatif in the model

####################### linear regression figures #######################
######### aridity
aridity.regline <- data.frame(emmeans(deltanmass_lmer_aridity,"aridity",
                                      at = list(aridity =seq(min(deltanmass_mad_data$aridity, na.rm = T), max(deltanmass_mad_data$aridity, na.rm = T), 0.01)),
                                      type = "response")) 

(aridity_plot <- ggplot(data = deltanmass_mad_data, aes(x = aridity, y = delta_nmass)) + 
    
    geom_point(aes(color = aridity, size = aridity))+
    scale_size_continuous(range = c(3, 1)) + # Set the range for the size of points
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    
    geom_smooth(data = aridity.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = aridity.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "pink", alpha = 0.3)+ 
    scale_x_continuous(limits = c(0, 2.2)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab(expression (italic('MI'))))

aridity_plot

ggsave("plots_new/aridity_plot.jpeg", plot = aridity_plot, 
       width = 30, height = 20, units = "cm")
######################################################################################################
################ temperature regression ############

tmp.regline <- data.frame(emmeans(deltanmass_lmer_aridity,"tmp",
                                      at = list(tmp =seq(min(deltanmass_mad_data$tmp, na.rm = T), max(deltanmass_mad_data$tmp, na.rm = T), 0.1)),
                                      type = "response")) 

(tmp_plot <- ggplot(data = deltanmass_mad_data, aes(x = tmp, y = delta_nmass)) + 
    
    geom_point(aes(color = tmp, size = tmp))+
    scale_size_continuous(range = c(3, 1)) + # Set the range for the size of points
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    
    geom_smooth(data = tmp.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = tmp.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "pink", alpha = 0.3)+ 
    
    scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['mass'] * ' (%)')) +
    xlab(expression (italic('Tg') ['°C'])))
    
tmp_plot

ggsave("plots_new/tmp_plot.jpeg", plot = tmp_plot, 
       width = 30, height = 20, units = "cm")


########## merge the two plots #############
rel_widths <- c(2, 2)
merged_plot_deltanmass<- plot_grid(aridity_plot, tmp_plot, ncol = 2, 
                                   align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                                   label_x = c(0.9, 0.9), rel_widths = rel_widths)
merged_plot_deltanmass


ggsave("plots_new/merged_plot_deltanmass.jpeg", plot = merged_plot_deltanmass, 
       width = 60, height = 30, units = "cm")


######## Structural equation model for deltanmass and delta ABG corrected by the max cover ##################
#############################################################################################################

delta_subset$PKgroup[delta_subset$Ptrt_fac == '0' & delta_subset$Ktrt_fac == '0'] <- '-P, -K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '1' & delta_subset$Ktrt_fac == '0'] <- '+P, -K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '0' & delta_subset$Ktrt_fac == '1'] <- '-P, +K+µ'
delta_subset$PKgroup[delta_subset$Ptrt_fac == '1' & delta_subset$Ktrt_fac == '1'] <- '+P, +K+µ'

delta_subset$Nfix_level = 0 ## creation of a column with Ntrt= 0 
delta_subset$Nfix_level[delta_subset$Nfix == 'yes'] = 1 
delta_subset$Nfix_level[delta_subset$Nfix == 'no'] = 0 

delta_subset$ps_level = 0 ## creation of a column with Ntrt= 0 
delta_subset$ps_level[delta_subset$photosynthetic_pathway == 'C3'] = 1 
delta_subset$ps_level[delta_subset$photosynthetic_pathway == 'C4'] = 0 

delta_sem_Mad <- subset(delta_subset,
                       delta_nmass < 3 * delta_nmass_mad &
                         delta_nmass > 3 * -delta_nmass_mad)

nrow(delta_sem_Mad) # 524
names(delta_sem_Mad)
subset_SEM = delta_sem_Mad [ , c(9,37,89,97,103,117,118,133,187,188,189,190,191,192)]
names(subset_SEM)
nrow(subset_SEM)

subset_SEM_ready <- na.omit(subset_SEM)
nrow(subset_SEM_ready) # 384
names(subset_SEM_ready)

############### Structural equation model #######################

########################### using psem ####################### 
delta_sem <- psem (
  
  deltaABG = lme(delta_spp_live_mass~ aridity + tmp + par2_gs + Nfix_level + ps_level +Ptrt_fac*Ktrt_fac, 
                 random = ~ 1|Taxon,
                 data = subset_SEM_ready),
  
  deltanmass = lme(delta_nmass~ delta_spp_live_mass + aridity + tmp + par2_gs + Nfix_level + ps_level +Ptrt_fac*Ktrt_fac,
                   random = ~ 1|Taxon,
                   data = subset_SEM_ready),
  int = lm(aridity~tmp , data=subset_SEM_ready)
  
)

summary(delta_sem)
plot(delta_sem)


################################# delta narea analysis ######################################

#############  remove instances where any delta values are 3 times higher than the MAD ##########################################
#################################################################################################################################
nrow(delta_subset) # 894
delta_narea_mad <- mad(delta_subset$delta_narea, na.rm = T)
delta_narea_mad   # 38.71
delta_narea_median <- median(delta_subset$delta_narea, na.rm = T)
delta_narea_median ### 11.42
delta_narea_mean<- mean(delta_subset$delta_narea, na.rm = T)
delta_narea_mean ## 37.07

deltanarea_mad_data <- subset(delta_subset,
                              delta_narea < 3 * delta_narea_mad &
                                delta_narea > 3 * -delta_narea_mad)

nrow(deltanarea_mad_data) # 443 
hist(deltanarea_mad_data$delta_narea) # normal dist

#### Hyp 1: Cold temperature and drought increase deltanmass, N fixation decreases deltanmass #######################
#########################################################################################################################
######### When using aridity index, tmp and aridity have significant effects but not par ##############
deltanarea_lmer_aridity <- lmer(delta_narea ~tmp + aridity + par2_gs + Nfix + 
                                  photosynthetic_pathway + Ptrt_fac*Ktrt_fac + 
                                  (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),
                                data = deltanarea_mad_data)
Anova(deltanarea_lmer_aridity, type = 'II') ## tmp, aridity, p, nfix, photo pathway effects
summary(deltanarea_lmer_aridity)


plot(resid(deltanarea_lmer_aridity) ~ fitted(deltanarea_lmer_aridity))

vif(deltanarea_lmer_aridity) # no colinearity
AIC(deltanarea_lmer_aridity) ## 4464.286

plot(deltanmass_lmer_aridity)
qqnorm(residuals(deltanmass_lmer_aridity))
qqline(residuals(deltanmass_lmer_aridity))

densityPlot(residuals(deltanarea_lmer_aridity))
shapiro.test(residuals(deltanarea_lmer_aridity)) 
outlierTest(deltanarea_lmer_aridity)

residuals <- resid(deltanarea_lmer_aridity)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

plot(fitted(deltanarea_lmer_aridity), residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity : none

r.squaredGLMM(deltanarea_lmer_aridity) ## 0.20, not terrible 

delta_narea_model_aridity <- data.frame(Var = c('Tg', 'aridity', 'PAR', 'N fixer', 
                                                'C3/C4', 'Soil P', 'Soil K', 'Soil P x Soil K'))

delta_narea_model_aridity$df <- as.matrix(Anova(deltanarea_lmer_aridity))[1:8, 2]
delta_narea_model_aridity$Slope <- c(
  summary(emtrends(deltanarea_lmer_aridity, ~tmp, var = "tmp"))[1, 2],
  summary(emtrends(deltanarea_lmer_aridity, ~aridity, var = "aridity"))[1, 2],
  summary(emtrends(deltanarea_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 2],
  NA, NA, NA, NA, NA)
delta_narea_model_aridity$SE <- c(
  summary(emtrends(deltanarea_lmer_aridity, ~tmp, var = "tmp"))[1, 3],
  summary(emtrends(deltanarea_lmer_aridity, ~aridity, var = "aridity"))[1, 3],
  summary(emtrends(deltanarea_lmer_aridity, ~par2_gs, var = "par2_gs"))[1, 3],
  NA, NA, NA, NA, NA)
delta_narea_model_aridity$p <- as.matrix(Anova(deltanarea_lmer_aridity))[1:8, 3]
delta_narea_model_aridity$RelImp <- as.matrix(calc.relip.mm(deltanarea_lmer_aridity)$lmg)[1:8]
delta_narea_model_aridity$RelImp <- delta_narea_model_aridity$RelImp * 100
delta_narea_model_aridity

setwd ("C:/Users/18632/Documents/Alissar/NutNet/NutNet_Paper/stat_output")
write.csv(delta_narea_model_aridity, "delta_narea_model_aridity.csv")

###################### # percentage increase of delta narea in plots receiving P compared to plots not receiving P (NP-P)
(summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[2,2] - summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[1,2])/
  summary(emmeans(deltanarea_lmer_aridity, ~Ptrt_fac))[1,2]
# - 0.607 : plots receiving P respond less to N addition than plots not receiving P 
## this effect is significant in the model 

# percentage increase of Nmass in N fixers compared to non-N fixers
(summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[2,2] - summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[1,2])/
  summary(emmeans(deltanarea_lmer_aridity, ~Nfix))[1,2]
# - 1.124: non fixers increase in deltaNarea is 112.4 % greater than  fixers 
## highly significant in the model


####################### linear regression figures #######################
######### aridity
aridity.regline <- data.frame(emmeans(deltanarea_lmer_aridity,"aridity",
                                      at = list(aridity =seq(min(deltanarea_mad_data$aridity, na.rm = T), max(deltanarea_mad_data$aridity, na.rm = T), 0.01)),
                                      type = "response")) 

(aridity_plot_naraea <- ggplot(data = deltanarea_mad_data, aes(x = aridity, y = delta_narea)) + 
    
    geom_point(aes(color = aridity, size = aridity))+
    scale_size_continuous(range = c(3, 1)) + # Set the range for the size of points
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    
    geom_smooth(data = aridity.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = aridity.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "pink", alpha = 0.3)+ 
    scale_x_continuous(limits = c(0, 2.2)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab(expression (italic('MI'))))

aridity_plot_naraea

ggsave("plots/aridity_plot_naraea.jpeg", plot = aridity_plot_naraea, 
       width = 30, height = 20, units = "cm")
######################################################################################################
################ temperature regression ############

tmp.regline <- data.frame(emmeans(deltanarea_lmer_aridity,"tmp",
                                  at = list(tmp =seq(min(deltanarea_mad_data$tmp, na.rm = T), max(deltanarea_mad_data$tmp, na.rm = T), 0.1)),
                                  type = "response")) 

(tmp_plot_narea <- ggplot(data = deltanarea_mad_data, aes(x = tmp, y = delta_narea)) + 
    
    geom_point(aes(color = tmp, size = tmp))+
    scale_size_continuous(range = c(3, 1)) + # Set the range for the size of points
    scale_color_gradient(low = "darkblue", high = "lightblue") +
    
    geom_smooth(data = tmp.regline, aes(y = emmean), 
                col = 'black', lwd = 2, alpha = 0.8) +
    
    geom_ribbon(data = tmp.regline, aes(y=emmean, ymin = lower.CL, ymax = upper.CL), fill = "pink", alpha = 0.3)+ 
    scale_x_continuous(limits = c(3, 21)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 40, colour = 'black'),
          axis.title.x = element_text(size = 40, colour = 'black'),
          axis.text.x = element_text(size = 30, colour = 'black'),
          axis.text.y = element_text(size = 30, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    ylab(expression('∆'* italic('N')['area'] * ' (%)')) +
    xlab(expression (italic('Tg') ['°C'])))

tmp_plot_narea

ggsave("plots/tmp_plot_narea.jpeg", plot = tmp_plot_narea, 
       width = 30, height = 20, units = "cm")
##########################################################################################################
merged_plot_deltanarea<- plot_grid(aridity_plot_naraea, tmp_plot_narea, ncol = 2, 
                                   align = "vh", labels = c("(a)", "(b)"), label_size = 30, 
                                   label_x = c(0.9, 0.9), rel_widths = rel_widths)
merged_plot_deltanarea


ggsave("plots_new/merged_plot_deltanarea.jpeg", plot = merged_plot_deltanarea, 
       width = 60, height = 30, units = "cm")


