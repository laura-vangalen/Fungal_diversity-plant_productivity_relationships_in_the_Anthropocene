# Code used to conduct analyses and create Figure 2 and Supplementary Figure S1 for "Fungal diversity–plant productivity relationships in the Anthropocene"

library(ggplot2) # v3.5.1
library(ggpubr) # v0.6.0
library(ggeffects) #ggpredict v2.2.0
library(terra) # v1.7-78
library(tidyr) # v1.3.1
library(ggh4x) # facet_grid2 v0.2.8
library(sjPlot) #plot_model v2.8.17
library(effects) # predictorEffects v4.2-2
library(mgcv) # gams v1.9-1
library(gridExtra) # v2.3

######## Methods of defining the fungal groups from Mikryukov et al (2023) ######################################################
# "We used the FungalTraits database v1.3 (24) to categorize fungal species into functional groups based on their ecological 
# and morphological characteristics. These ecological groups are as follows: (i) AM fungi (3.7% of total OTUs), including all 
# Glomeromycota but excluding all Endogonomycetes due to the lack of data separating AM species from free-living species; 
# (ii) EcM fungi (15.4%); (iii) non-EcM Agaricomycetes (9.3%), primarily saprotrophic macrofungi; (iv) molds (5.8%), including 
# Mortierellales, Mucorales, Umbelopsidales, Aspergillaceae and Trichocomaceae of Eurotiales, and Trichoderma of Hypocreales; 
# (v) putative pathogens (11.4%), consisting of plant, animal, and fungal pathogens as primary or secondary lifestyles; (vi) 
# soil borne OHPs (7.5%), excluding Mortierellales; (vii) yeasts (1.3%), excluding dimorphic yeasts; and (viii) other unicellular, 
# nonyeast fungi (9.6%), including Chytridiomycota, Aphelida, Rozellomycota, and other early-diverging fungal lineages."
##################################################################################################################################

# For our analysis, we use (1) for AM fungi, (ii) for EcM fungi, and (v) for pathogens


#################################################################################################################################################################
# extract random values from raster layers to use as data points
#################################################################################################################################################################
# raster layers are from the sources listed in the manuscript methods. All rasters were reprojected to equirectangular projection before use.
AM=rast("./Mikryukov_2023_richness_layers/Alpha_S_AM_GSMc.tif")
ECM=rast("./Mikryukov_2023_richness_layers/Alpha_S_EcM_GSMc.tif")
Path=rast("./Mikryukov_2023_richness_layers/Alpha_S_Path_GSMc.tif")
ph=rast("./Climate_layers/SG_Soil_pH_H2O_005cm_repro_Mikryukov.tif")
soc=rast("./Climate_layers/SG_SOC_Content_005cm_repro_Mikryukov.tif")
temp=rast("./Climate_layers/CHELSA_BIO_Annual_Mean_Temperature_repro_Mikryukov.tif")
precip=rast("./Climate_layers/CHELSA_BIO_Annual_Precipitation_repro_Mikryukov.tif")
ele=rast("./Climate_layers/EarthEnvTopoMed_Elevation_repro_Mikryukov.tif")
man=rast("./Human_layers/ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation_repro_Mikryukov.tif")
urb=rast("./Human_layers/ConsensusLandCoverClass_Urban_Builtup_repro_Mikryukov.tif")
r1=rast("./Human_layers/CSP_Global_Human_Modification_repro_Mikryukov.tif")
r2=rast("./Human_layers/GHS_Population_Density_repro_Mikryukov.tif")
r3=rast("./Human_layers/GRIP4_DistanceToAllRoads_repro_Mikryukov.tif")

stack=c(AM,ECM,Path,ph,soc,temp,precip,ele,man,urb,r1,r2,r3)
set.seed(435)
sample1=spatSample(stack,size=10000,na.rm=T,xy=T) # 10,000 random grid cells
saveRDS(sample1,"raster_sample1.rds")

# now add the logging one. Doing separately because this doesn't differentiate between non-forest land areas and sea areas. So using the other layers to choose the random points, then will know that all points are terrestrial and anything in the logging layer with 128 will be non-forest terrestrial
r1=rast("./Human_layers/FML_v3-2_with-colorbar_repro_Mikryukov.tif")
names(r1)<-"FML_v3.2"

sample1=readRDS("./raster_sample1.rds")
sample1_new=terra::extract(r1,sample1[c("x","y")])
sample1_all=cbind(sample1,sample1_new[c(2:ncol(sample1_new))])
saveRDS(sample1_all,"./raster_sample1.rds")

# 11 – Naturally regenerating forest without any signs of management, including primary forests;
# 20 – Naturally regenerating forest with signs of management, e.g., logging, clear cuts etc;
# 31 – Planted forests;
# 32 – Plantation forests (rotation time up to 15 years);
# 40 – Oil palm plantations;
# 53 – Agroforestry;
# 128 - non-forest area



#################################################################################################################################################################
# correlations between human variables, and tidy data
#################################################################################################################################################################
sample1=readRDS("raster_sample1.rds")
cor(sample1[c("ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation","ConsensusLandCoverClass_Urban_Builtup","CSP_Global_Human_Modification","GRIP4_DistanceToAllRoads","GHS_Population_Density")])
# correlations are fine - highest is 0.74 (ConsensusLandCoverClass_Urban_Builtup and CSP global human modification)

sample1$EarthEnvTopoMed_Elevation=ifelse(sample1$EarthEnvTopoMed_Elevation<0,0,sample1$EarthEnvTopoMed_Elevation) # make negative elevations zero
sample1$FML_v3.2=as.factor(sample1$FML_v3.2) # change to factor
sample1=sample1[sample1$GHS_Population_Density <1500,] # drop population outlier


#################################################################################################################################################################
# Models for richness ~ global human modification index (Figure 2)
#################################################################################################################################################################

# AM
mod1=lm(RES_LOG_arbuscular_mycorrhizal_Richness~CSP_Global_Human_Modification+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1) 
plot(predictorEffects(mod1,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=3,cols=3,main="")

# ECM
mod2=lm(RES_LOG_ectomycorrhizal_Richness~CSP_Global_Human_Modification+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod2)
par(mfrow=c(2,2))
plot(mod2) 
plot(predictorEffects(mod2,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=3,cols=3,main="") 

# pathogens
mod3=lm(RES_LOG_pathogen_Richness~CSP_Global_Human_Modification+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod3)
par(mfrow=c(2,2))
plot(mod3)
plot(predictorEffects(mod3,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=3,cols=3,main="")


### Extract partial residuals ##############################################################################################################################################################

# AM
pred_AM=predict(mod1, newdata = data.frame(CSP_Global_Human_Modification = sample1$CSP_Global_Human_Modification,
                                        SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                        CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(CSP_Global_Human_Modification = sample1$CSP_Global_Human_Modification,
                                        SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                        CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(CSP_Global_Human_Modification = sample1$CSP_Global_Human_Modification,
                                        SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                        CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$CSP_Global_Human_Modification,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long$fungi_type = factor(parital_effects_long$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order

fungi.labs <- c("AM fungi","EcM fungi","Pathogenic fungi") # set label names
names(fungi.labs) <- c('pred_AM','pred_ECM','pred_path')

png("figure_partial_residuals.png",width=7.5,height=3.2,units="in",res=600)
ggplot(parital_effects_long, aes(x=human_orig, y=residuals)) + 
  theme_classic()+
  facet_grid2(~fungi_type,labeller = labeller(fungi_type = fungi.labs),scales="free",independent="y") +
  geom_point(alpha=0.05,colour="blue")+
  geom_smooth(method=lm,colour="blue4")+
  stat_cor(aes(label = ..r.label..),hjust=1,label.x.npc="right",label.y.npc="bottom",vjust=0) +
  theme(strip.background = element_blank(),
  strip.text.x = element_text(hjust=0,size=12))+
  labs(title="",y="Fungal richness (partial residuals)",x="Global human modification index")
dev.off()



#################################################################################################################################################################
# Models for richness ~ other human impact factors (Figure S1)
#################################################################################################################################################################

# AM
mod1=lm(RES_LOG_arbuscular_mycorrhizal_Richness~ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation+FML_v3.2+ConsensusLandCoverClass_Urban_Builtup+GHS_Population_Density+GRIP4_DistanceToAllRoads+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1) 
plot(predictorEffects(mod1,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=4,cols=3,main="")

# ECM
mod2=lm(RES_LOG_ectomycorrhizal_Richness~ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation+FML_v3.2+ConsensusLandCoverClass_Urban_Builtup+GHS_Population_Density+GRIP4_DistanceToAllRoads+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod2)
par(mfrow=c(2,2))
plot(mod2) 
plot(predictorEffects(mod2,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=3,cols=3,main="") 

# pathogens
mod3=lm(RES_LOG_pathogen_Richness~ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation+FML_v3.2+ConsensusLandCoverClass_Urban_Builtup+GHS_Population_Density+GRIP4_DistanceToAllRoads+
          SG_Soil_pH_H2O_005cm+SG_SOC_Content_005cm+CHELSA_BIO_Annual_Mean_Temperature+CHELSA_BIO_Annual_Precipitation+EarthEnvTopoMed_Elevation,data=sample1)
summary(mod3)
par(mfrow=c(2,2))
plot(mod3)
plot(predictorEffects(mod3,partial.residuals=TRUE), partial.residual=list(pch=".", col="#FF00FF80"),rows=3,cols=3,main="")


### Extract partial residuals ##############################################################################################################################################################

######### for cultivated #########
reference_level <- levels(sample1$FML_v3.2)[1] # work out the reference level of the categorical variable

# AM
pred_AM=predict(mod1, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = sample1$ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation,
                                           FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                           SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                           CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = sample1$ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation,
                                            FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                            SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                            CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = sample1$ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation,
                                             FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                             SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                             CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long$fungi_type = factor(parital_effects_long$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order


######### for forest categories #########
reference_level <- levels(sample1$FML_v3.2)[1] # work out the reference level of the categorical variable

# AM
pred_AM=predict(mod1, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                           FML_v3.2=sample1$FML_v3.2,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                           SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                           CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                            FML_v3.2=sample1$FML_v3.2,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                            SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                            CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                             FML_v3.2=sample1$FML_v3.2,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                             SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                             CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$FML_v3.2,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long_tree=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long_tree$fungi_type = factor(parital_effects_long_tree$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order

######### for urban builtup #########
reference_level <- levels(sample1$FML_v3.2)[1] # work out the reference level of the categorical variable

# AM
pred_AM=predict(mod1, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                           FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=sample1$ConsensusLandCoverClass_Urban_Builtup,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                           SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                           CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                            FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=sample1$ConsensusLandCoverClass_Urban_Builtup,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                            SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                            CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                             FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=sample1$ConsensusLandCoverClass_Urban_Builtup,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=0,
                                             SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                             CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$ConsensusLandCoverClass_Urban_Builtup,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long_urban=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long_urban$fungi_type = factor(parital_effects_long_urban$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order

######### for pop density #########
reference_level <- levels(sample1$FML_v3.2)[1] # work out the reference level of the categorical variable

# AM
pred_AM=predict(mod1, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                           FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=sample1$GHS_Population_Density,GRIP4_DistanceToAllRoads=0,
                                           SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                           CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                            FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=sample1$GHS_Population_Density,GRIP4_DistanceToAllRoads=0,
                                            SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                            CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                             FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=sample1$GHS_Population_Density,GRIP4_DistanceToAllRoads=0,
                                             SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                             CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$GHS_Population_Density,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long_pop=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long_pop$fungi_type = factor(parital_effects_long_pop$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order

######### for distance to roads #########
reference_level <- levels(sample1$FML_v3.2)[1] # work out the reference level of the categorical variable

# AM
pred_AM=predict(mod1, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                           FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=sample1$GRIP4_DistanceToAllRoads,
                                           SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                           CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# ECM
pred_ECM=predict(mod2, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                            FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=sample1$GRIP4_DistanceToAllRoads,
                                            SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                            CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant
# Pathogens
pred_path=predict(mod3, newdata = data.frame(ConsensusLandCoverClass_Cultivated_and_Managed_Vegetation = 0,
                                             FML_v3.2=reference_level,ConsensusLandCoverClass_Urban_Builtup=0,GHS_Population_Density=0,GRIP4_DistanceToAllRoads=sample1$GRIP4_DistanceToAllRoads,
                                             SG_Soil_pH_H2O_005cm = 0, SG_SOC_Content_005cm = 0, CHELSA_BIO_Annual_Mean_Temperature = 0,
                                             CHELSA_BIO_Annual_Precipitation = 0,EarthEnvTopoMed_Elevation = 0),se.fit=T) # calculate effect of human variable keeping the others constant

# make df of residuals added to the fit for plotting the partial residuals
partial_effects=data.frame(human_orig=sample1$GRIP4_DistanceToAllRoads,pred_AM=pred_AM$fit+resid(mod1),pred_ECM=pred_ECM$fit+resid(mod2),
                           pred_path=pred_path$fit+resid(mod3))
# ggplot format
parital_effects_long_roads=data.frame(pivot_longer(partial_effects,cols=c(2:ncol(partial_effects)),names_to="fungi_type",values_to="residuals"))
parital_effects_long_roads$fungi_type = factor(parital_effects_long_roads$fungi_type, levels=c('pred_AM','pred_ECM','pred_path')) # set order



## plots ######################################################################################################

fungi.labs <- c("AM fungi","EcM fungi","Pathogenic fungi") # set label names
names(fungi.labs) <- c('pred_AM','pred_ECM','pred_path')


p1=ggplot(parital_effects_long, aes(x=human_orig, y=residuals)) + 
  theme_classic()+
  facet_grid2(~fungi_type,labeller = labeller(fungi_type = fungi.labs),scales="free",independent="y") +
  geom_point(alpha=0.05,colour="blue")+
  geom_smooth(method=lm,colour="blue4")+
  stat_cor(aes(label = ..r.label..),hjust=1,label.x.npc="right",label.y.npc="bottom",vjust=0) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust=0,size=12))+
  labs(title="",y="",x="Cultivated and managed vegetation (%)")
p1

p2=ggplot(parital_effects_long_tree, aes(x=human_orig, y=residuals)) + 
  theme_classic()+
  facet_grid2(~fungi_type,scales="free",independent="y") +
  geom_boxplot(alpha=0.1,colour="blue")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  scale_x_discrete(labels = c('1','2','3',"4","5","6","7"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  labs(title="",y="",x="Forest management classes")
p2

p3=ggplot(parital_effects_long_urban, aes(x=human_orig, y=residuals)) + 
  theme_classic()+
  facet_grid2(~fungi_type,scales="free",independent="y") +
  geom_point(alpha=0.05,colour="blue")+
  geom_smooth(method=lm,colour="blue4")+
  stat_cor(aes(label = ..r.label..),hjust=1,label.x.npc="right",label.y.npc="bottom",vjust=0) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  labs(title="",y="Fungal richness (partial residuals)",x="Urban builtup (%)")
p3

p4=ggplot(parital_effects_long_pop, aes(x=human_orig, y=residuals)) + 
  theme_classic()+
  facet_grid2(~fungi_type,scales="free",independent="y") +
  geom_point(alpha=0.05,colour="blue")+
  geom_smooth(method=lm,colour="blue4")+
  stat_cor(aes(label = ..r.label..),hjust=1,label.x.npc="right",label.y.npc="bottom",vjust=0) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  labs(title="",y="",x="Population density (per 250m2)")
p4

p5=ggplot(parital_effects_long_roads, aes(x=human_orig/1000, y=residuals)) + # in km (original units m)
  theme_classic()+
  facet_grid2(~fungi_type,scales="free",independent="y") +
  geom_point(alpha=0.05,colour="blue")+
  geom_smooth(method=lm,colour="blue4")+
  stat_cor(aes(label = ..r.label..),hjust=1,label.x.npc="right",label.y.npc="bottom",vjust=0) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  labs(title="",y="",x="Distance to nearest road (km)")
p5

png("figure_partial_residuals_ag_forestry.png",width=7.5,height=13,units="in",res=600)
grid.arrange(p1,p2,p3,p4,p5, nrow = 5)
dev.off()


