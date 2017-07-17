# R script to
# - interpolate/subset NPP and ALD data, which has already been computed/cleaned by another script (1-sla.R)
# - read in soil moisture/pH/temp data, which has already been computed/cleaned by another script (FFcore_process.R)
# - read in 2014 and 2015 SLA data
# - exploratory plots with SLA, ALD, NPP, soil temp/moisture/pH

# Carolyn Anderson November 2015

# Load necessary packages
library(ggplot2)
library(reshape2)
library(stringr)

#theme_set(theme_gray(base_size = 30))
#theme_set(theme_gray())

# Need to setwd to Github folder to run 1-sla.R
#setwd("C:/Users/ande352/Documents/GitHub/cpcrw/sla")

# Set wd to get SLA file
#setwd("//pnl/projects/Alaska_Carbon/") #if using Shared Drive
setwd("~/cpcrw-sla/")

# -----------------------------------------------------------------------------
# 2014 SLA data
# -----------------------------------------------------------------------------
sla.2014 <- read.csv("CPCRW_Data_2014/CPCRW_SLA/26Aug2014_SLA.csv")[,-c(1,8)]
sla.2014$Transect <- 5
sla.2014$Year <- 2014
colnames(sla.2014)[c(1:4)] <- c("Tree", "Position.E.to.W_m", "species2", "Leaf")
sla.2014$Position.E.to.W_m[sla.2014$Position.E.to.W_m == 75] <- 72
sla.2014$Species <- ifelse(sla.2014$species2 == "Alder", "ALSP", ifelse(sla.2014$species2 == "Spruce", "PIMA", "BEPA"))

# SLA column, adapted from Ben's github script
# Projected leaf area (PLA) to hemisurface leaf area (HSLA) 
# conversion based on Bond-Lamberty et al. (2003) values
pla_to_hsla <- data.frame(Species=c("ALSP", "BEPA", "PIGL", "PIMA"),
                          p2h=c(1.0, 1.0, 1.55, 1.55))

sla.2014 <- merge(sla.2014, pla_to_hsla)
sla.2014$SLA <- with(sla.2014, p2h * SurfaceArea_cm2 / WeightDry_g)

#write.csv(sla.2014, "CPCRW_Data_2014/CPCRW_SLA/CPCRW_2014_SLA_processed.csv")


## histogram
# alder
ggplot(subset(sla.2014, species2 %in% "Spruce")) +
  aes(x=SLA) +
  geom_histogram() +
  xlab(expression(spruce~SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------
# 2015 SLA data
# -----------------------------------------------------------------------------
SLA_decoder <- read.csv("CPCRW_Data_2014/CPCRW_SLA/CPCRW_2015_SLA_decoder.csv")
sla.2015 <- read.csv("CPCRW_Data_2014/CPCRW_SLA/CPCRW_2015_SLA.csv")[,-c(8,9)]
sla.2015.total <- merge(SLA_decoder, sla.2015, by="Number") #Merge SLA with decoder

sla.2015.total$Year <- 2015
sla.2015.total$species2 <- ifelse(sla.2015.total$Species=="ALSP", "Alder", "Spruce")
sla.2015.total$Transect <- as.numeric(gsub('T',"", sla.2015.total$Transect)); # take out 'T' and change to numbers
colnames(sla.2015.total)[c(1,3,5)] <- c("Tree","Position.E.to.W_m","Leaf")

# SLA column, adapted from Ben's github script
# Projected leaf area (PLA) to hemisurface leaf area (HSLA) 
# conversion based on Bond-Lamberty et al. (2003) values
pla_to_hsla <- data.frame(Species=c("ALSP", "PIGL", "PIMA"),
                          p2h=c(1.0, 1.55, 1.55))
sla.2015.total <- merge(sla.2015.total, pla_to_hsla)

sla.2015.total$SLA <- with(sla.2015.total, p2h * SurfaceArea_cm2 / WeightDry_g)

# Merge with Transect positions
#transect <- read.csv("CPCRW_Data_2014/CPCRW_2014_Transect_Positions.csv") # Transect positions
#sla.2015.total <- merge(sla.2015.total, transect, by="Transect")

#write.csv(sla.2015.total, "CPCRW_Data_2014/CPCRW_SLA/CPCRW_2015_SLA_processed.csv")

# -----------------------------------------------------------------------------
# TOPOGRAPHIC AND SITE-WIDE DATA (ALL 6 TRANSECTS)
# -----------------------------------------------------------------------------
#slope <- read.csv("CPCRW_Data_2014/CPCRW_2014_Slope_Linear_Interp.csv")
slope <- read.csv("CPCRW_Data_2014/CPCRW_2014_slope_avg.for.SLA.csv") #averaged slope
tree <- read.csv("CPCRW_Data_2014/CPCRW_VegSurvey/CPCRW_2014_Tree_Linear_Interp.csv")
npp <- read.csv("CPCRW_Data_2014/CPCRW_NPP/CPCRW_2014_NPP_Linear_Interp.csv")
#ald <- read.csv("CPCRW_Data_2014/CPCRW_ALD/CPCRW_2014_Fall_ALD_Linear_Interp.csv") #orig interpolated ALD
ald <- read.csv("CPCRW_Data_2014/CPCRW_ALD/CPCRW_2014_ALD_avg.for.SLA.csv") #averaged ALD
colnames(ald)[3] <- "ald_cm"

# Merge files WITHOUT SLA
all.merge <- merge(slope, tree, by=c("Position.E.to.W_m", "Transect"))
all.merge <- merge(all.merge, npp, by=c("Position.E.to.W_m", "Transect"))
all.merge <- merge(all.merge, ald, by=c("Position.E.to.W_m", "Transect"))
#write.csv(all.merge, "CPCRW_Data_2014/CPCRW_SLA/23March2016_all.merge.csv")



# Merge files - all transects (2015 SLA)
sla_all.2015 <- merge(sla.2015.total, slope, by=c("Position.E.to.W_m", "Transect"))
sla_all.2015 <- merge(sla_all.2015, tree, by=c("Position.E.to.W_m", "Transect"))
sla_all.2015 <- merge(sla_all.2015, ald, by=c("Position.E.to.W_m", "Transect"))
sla_all.2015 <- merge(sla_all.2015, npp, by=c("Position.E.to.W_m", "Transect"))

# Sort by Species, then Tree, then Leaf
sla_all.2015 <- sla_all.2015[order(sla_all.2015[, 3], sla_all.2015[, 4], sla_all.2015[, 5]), ]

#write.csv(sla_all.2015, "CPCRW_Data_2014/CPCRW_SLA/21Jan2016_SLA_all.processed.csv")

# -----------------------------------------------------------------------------
# SOIL DATA (3 TRANSECTS)
# -----------------------------------------------------------------------------
# Getting interpolated data for soil properties (depths: 1.75cm, 6cm, 12cm)
##FTIR
ftir_1.75cm <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/FTIR/Data/Processed/CPCRW_2014_FFcores_FTIR_hum_Linear_Interp_1.75cm.csv")
ftir_1.75cm$Depth_cm <- 1.75

ftir_6cm <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/FTIR/Data/Processed/CPCRW_2014_FFcores_FTIR_hum_Linear_Interp_6cm.csv")
ftir_6cm$Depth_cm <- 6

tmp1 <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/FTIR/Data/Processed/CPCRW_2014_FFcores_FTIR_hum_Linear_Interp_12cmT5.csv")
tmp2 <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/FTIR/Data/Processed/CPCRW_2014_FFcores_FTIR_hum_Linear_Interp_12cmT6T7.csv")
ftir_12cm <- rbind(tmp1, tmp2)
ftir_12cm$Depth_cm <- 12

##Total C, N, CN
totalc_1.75cm <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/TotalCNS/CPCRW_2014_totalCNS_Linear_Interp_1.75cm.csv")
totalc_1.75cm$Depth_cm <- 1.75

totalc_6cm <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/TotalCNS/CPCRW_2014_totalCNS_Linear_Interp_6cm.csv")
totalc_6cm$Depth_cm <- 6

tmp3 <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/TotalCNS/CPCRW_2014_totalCNS_Linear_Interp_12cmT5.csv")
tmp4 <- read.csv("CPCRW_Data_2014/CPCRW_Carbon/TotalCNS/CPCRW_2014_totalCNS_Linear_Interp_12cmT6T7.csv")
totalc_12cm <- rbind(tmp3, tmp4)
totalc_12cm$Depth_cm <- 12

##Soil mositure, pH, BD
moisture_1.75cm <- read.csv("CPCRW_Data_2014/CPCRW_Coring/Processed/CPCRW_2014_FFcores_Linear_Interp_1.75cm.csv")
moisture_1.75cm$Depth_cm <- 1.75

moisture_6cm <- read.csv("CPCRW_Data_2014/CPCRW_Coring/Processed/CPCRW_2014_FFcores_Linear_Interp_6cm.csv")
moisture_6cm$Depth_cm <- 6

tmp5 <- read.csv("CPCRW_Data_2014/CPCRW_Coring/Processed/CPCRW_2014_FFcores_Linear_Interp_12cmT5.csv")
tmp6 <- read.csv("CPCRW_Data_2014/CPCRW_Coring/Processed/CPCRW_2014_FFcores_Linear_Interp_12cmT6T7.csv")
moisture_12cm <- rbind(tmp5, tmp6)
moisture_12cm$Depth_cm <- 12


##Soil temperature
temp <- read.csv("CPCRW_Data_2014/CPCRW_Temperature/CPCRW_2014_FF_temp_Linear_Interp.csv")
temp12 <- read.csv("CPCRW_Data_2014/CPCRW_Temperature/CPCRW_2014_FF_temp_Linear_Interp_12cm.csv")
temp.all <- merge(temp, temp12)
soil_temp.melt <- melt(temp.all, id.vars=c("Position.E.to.W_m","Transect"))
tmp <- str_split_fixed(soil_temp.melt$variable, "_", 2)
soil_temp.melt2 <- cbind(soil_temp.melt, tmp)
colnames(soil_temp.melt2)[c(4,6)] <- c("Temperature","Depth_cm")
soil_temp.melt2$Depth_cm <- gsub(".{2}$", "", soil_temp.melt2$Depth_cm)
soil_temp.melt2 <- soil_temp.melt2[,c(1,2,4,6)]
soil_temp.melt2$Depth_cm[soil_temp.melt2$Depth_cm==1] <- 1.75 #For temperature, call 1cm --> 1.75cm



## rbind similar data frames
ftir.total <- rbind(ftir_1.75cm, ftir_6cm)
ftir.total <- rbind(ftir.total, ftir_12cm)

moisture.total <- rbind(moisture_1.75cm, moisture_6cm)
moisture.total <- rbind(moisture.total, moisture_12cm)

totalc.total <- rbind(totalc_1.75cm, totalc_6cm)
totalc.total <- rbind(totalc.total, totalc_12cm)



## Merge for total
soil.total <- merge(ftir.total, moisture.total, by=c("Position.E.to.W_m","Transect","Depth_cm"))
soil.total <- merge(soil.total, totalc.total, by=c("Position.E.to.W_m","Transect","Depth_cm"))
soil.total <- merge(soil.total, soil_temp.melt2, by=c("Position.E.to.W_m","Transect","Depth_cm"))

# Merge with moss (no Depth_cm)
moss <- read.csv("CPCRW_Data_2014/CPCRW_VegSurvey/CPCRW_2014_FF_moss_Linear_Interp.csv")[,c(1:3)]
soil.total <- merge(soil.total, moss, by=c("Position.E.to.W_m","Transect"))
# Get only relevant positions (0, 15, 30, 45, 60, 72m)
soil.total.sla <- subset(soil.total, Position.E.to.W_m  %in% c(0, 15, 30, 45, 60, 72))
#write.csv(soil.total.sla, "CPCRW_Data_2014/CPCRW_SLA/25March2016_all.merge.soil.csv")


# Merge with sla
sla_all_soil.2015 <- merge(sla_all.2015, soil.total.sla, by=c("Position.E.to.W_m", "Transect"))
# Sort by Species, then Tree, then Leaf
sla_all_soil.2015 <- sla_all_soil.2015[order(sla_all_soil.2015[, 3], sla_all_soil.2015[, 4], sla_all_soil.2015[, 5]), ]
#write.csv(sla_all_soil.2015, "CPCRW_Data_2014/CPCRW_SLA/25March2016_SLA_subset_soil.processed.csv")

# -----------------------------------------------------------------------------
## Compare 2014 and 2015 SLA data
# -----------------------------------------------------------------------------
# rbind the 2014 and 2015 SLA data
sla.1415 <- rbind(sla.2014, sla.2015.total)

ggplot(subset(sla.1415, species2 %in% "Spruce")) +
  aes(x=Position.E.to.W_m, y=SLA, color=factor(Year)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Transect)

ggplot(subset(sla.1415, Transect %in% 5 & species2 %in% "Spruce")) +
  aes(x=WeightDry_g, y=SurfaceArea_cm2, color=factor(Year)) +
  geom_point() +
  geom_smooth()
