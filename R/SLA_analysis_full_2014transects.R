# R script to
# - read in SLA and transect/soil data
#   - SLA data processed in SLA_data.R
#   - soil data processed in FFcore_process.R
# - analyze SLA data

# Carolyn Anderson November 2015

# Load necessary packages
library(tidyverse)
library(ape)
library(nlme)
library(MuMIn)
library(VGAM)
library(grid)
library(gridExtra)
library(coin)
library(segmented)

# Set wd to get SLA files
#setwd("//pnl/projects/Alaska_Carbon/") #if using Shared Drive
setwd("~/cpcrw-sla/")

# -----------------------------------------------------------------------------
# AVERAGE BASAL AREA AND TREE DENSITY VALUES
# -----------------------------------------------------------------------------
# Run 1-tree_survey.R script from:
# https://github.com/bpbond/cpcrw/blob/master/tree_survey/1-tree_survey.R

# Average Basal area (m2/ha) and Density (/ha), and corresponding standard deviation
mean(d1$`Basal area (m2/ha)`)
sd(d1$`Basal area (m2/ha)`)

mean(d1$`Density (/ha)`)
sd(d1$`Density (/ha)`)

# -----------------------------------------------------------------------------
# PROCESSING DATA (FIGURES/MODELS)
# -----------------------------------------------------------------------------
# Load processed SLA data
sla_all <- read.csv("CPCRW_Data_2014/CPCRW_SLA/21Jan2016_SLA_all.processed.csv")
sla_all_soil <- read.csv("CPCRW_Data_2014/CPCRW_SLA/25March2016_SLA_subset_soil.processed.csv")

# Load processed data WITHOUT SLA (no gaps)
sla_all.nosla <- read.csv("CPCRW_Data_2014/CPCRW_SLA/23March2016_all.merge.csv")
sla_all_soil.nosla <- read.csv("CPCRW_Data_2014/CPCRW_SLA/25March2016_all.merge.soil.csv") #This has depths 1.75cm, 6cm, 12cm
sla_all_soil.nosla2 <- merge(sla_all_soil.nosla, sla_all.nosla, by=c("Position.E.to.W_m","Transect")) #This has depths 1.75cm, 6cm, 12cm

# Merge sla_all.nosla with transect positions (all directions)
tran.pos <- read.csv("CPCRW_Data_2014/CPCRW_2014_Transect_Positions.csv")
sla_all.nosla2 <- merge(sla_all.nosla, tran.pos)

# Add E-W cyclical soil sampling for transects 5-7
core.pos <- read.csv("CPCRW_Data_2014/CPCRW_Coring/19Sept2014_CoreDataSept2014.csv")[,c(2,4)]
colnames(core.pos)[2] <- "Position.E.to.W_m"
core.pos <- merge(core.pos, tran.pos)

# Make ALD classes (based on 2 or 3 categories)
#3 classes
sla_all$ald_class <- ifelse(sla_all$ald_cm < 100, "Low ALD", ifelse(sla_all$ald_cm > 150, "ALD > 150cm", "Mid ALD"))
sla_all$ald_class <- as.factor(sla_all$ald_class)
sla_all_soil$ald_class <- ifelse(sla_all_soil$ald_cm < 100, "Low ALD", ifelse(sla_all_soil$ald_cm > 150, "ALD > 150cm", "Mid ALD"))
sla_all_soil$ald_class <- as.factor(sla_all_soil$ald_class)

#2 classes (140 based on threshold seen in regressions)
sla_all$ald_class2 <- ifelse(sla_all$ald_cm < 140, "Shallow", "Deep")
sla_all$ald_class2 <- as.factor(sla_all$ald_class2)

sla_all_soil$ald_class2 <- ifelse(sla_all_soil$ald_cm < 140, "Shallow", "Deep")
sla_all_soil$ald_class2 <- as.factor(sla_all_soil$ald_class2)

sla_all.nosla$ald_class2 <- ifelse(sla_all.nosla$ald_cm < 140, "Shallow", "Deep")
sla_all.nosla$ald_class2 <- as.factor(sla_all.nosla$ald_class2)

sla_all_soil.nosla2$ald_class2 <- ifelse(sla_all_soil.nosla2$ald_cm < 140, "Shallow", "Deep")
sla_all_soil.nosla2$ald_class2 <- as.factor(sla_all_soil.nosla2$ald_class2)

# -----------------------------------------------------------------------------
# AVERAGE AND MEDIAN ALD VALUES
# -----------------------------------------------------------------------------
range(sla_all.nosla$ald_cm) #min and max ALD values
mean(sla_all.nosla$ald_cm)
median(sla_all.nosla$ald_cm)

# Histogram of ALD values
ggplot(sla_all.nosla) +
  aes(x=ald_cm) +
  geom_histogram()

# -----------------------------------------------------------------------------
# SAMPLING LAYOUT
# -----------------------------------------------------------------------------
# Soil data
sla_all.nosla2$soildata <- -999
sla_all.nosla2$soildata <- ifelse(sla_all.nosla2$Transect < 8, "Soil data", "No soil data")
sla_all.nosla2$Transect <- as.factor(sla_all.nosla2$Transect)

ggplot(subset(sla_all.nosla2)) +
  aes(x=Position.S.to.N_m, y=Position.E.to.W_m) +
#  geom_vline(data=sla_all.nosla2, aes(xintercept=Position.S.to.N_m, linetype="Transect"), show.legend=FALSE) + 
  geom_point(aes(color="value3"), size=8) + #aes(color=soildata)
  geom_point(data=core.pos, aes(x=Position.S.to.N_m, y=Position.E.to.W_m), size=4) +
  xlab("Meters (south to north)") + ylab("Meters (east to west)") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal(xlim=c(-5, 75), ylim=c(-5, 75)) + 
  scale_color_manual(name="", values=c("value3"="gray"))

# -----------------------------------------------------------------------------
## Change some variables to factors for data analysis;
## Create data frames with replicate averages
# -----------------------------------------------------------------------------
#Make Transect, Tree, and Leaf into factors
sla_all_soil$Transect <- as.factor(sla_all_soil$Transect)
sla_all_soil$Tree <- as.factor(sla_all_soil$Tree)
sla_all_soil$Leaf <- as.factor(sla_all_soil$Leaf)

sla_all$Transect <- as.factor(sla_all$Transect)
sla_all$Tree <- as.factor(sla_all$Tree)
sla_all$Leaf <- as.factor(sla_all$Leaf)
#sla_all$Position.E.to.W_m <- as.factor(sla_all$Position.E.to.W_m)

#For alder and spruce SLA
sla_alder <- subset(sla_all_soil, species2 %in% "Alder" & Depth_cm %in% 6)
sla_spruce <- subset(sla_all_soil, species2 %in% "Spruce" & Depth_cm %in% 6)

# Reducing to averages for each tree
#sla_all
colnames(sla_all)
sla_all.red <- sla_all[,c(2:5,8:11,13:22)]
colnames(sla_all.red)
sla_all.ag <- do.call(data.frame, aggregate(cbind(WeightDry_g, SurfaceArea_cm2, SLA) ~  
                                              Position.E.to.W_m + Transect + Species + species2 + Tree + 
                                              SlopePercent + area.cm.2 + Individuals + Spruce.area.cm.2 + spruce.area.frac +
                                              ald_cm + ald_class2 + npp_gC + ald_class + ald_class2,
                                            sla_all.red, function(x) c(mean=mean(x), sd = sd(x))))

sla_all.ag.alder <- subset(sla_all.ag, species2 %in% "Alder")
sla_all.ag.spruce <- subset(sla_all.ag, species2 %in% "Spruce")

#sla_all, only 6cm
colnames(sla_all_soil)
sla_all_soil.red <- subset(sla_all_soil, Depth_cm %in% 6)
sla_all_soil.red <- sla_all_soil.red[,c(2:5,8:11,13:37)]
colnames(sla_all_soil.red)
sla_soil.ag <- do.call(data.frame, aggregate(cbind(WeightDry_g, SurfaceArea_cm2, SLA) ~  
                                               Position.E.to.W_m + Transect + Species + species2 + Tree + Depth_cm + SlopePercent + area.cm.2 + 
                                               Individuals + Spruce.area.cm.2 + spruce.area.frac + ald_cm + ald_class2 + npp_gC + 
                                               Hum_29 + Hum_16 + Hum_15 + Ow.mean + Ow.dry.mean + Ov.mean + pH.mean + BD.mean +
                                               N_per + C_per + S_per + CN + Temperature + MossDepth,
                                             sla_all_soil.red, function(x) c(mean=mean(x), sd = sd(x))))

#sla_all (1.75, 6, 12cm)
colnames(sla_all_soil)
#sla_all_soil.red2 <- subset(sla_all_soil, Depth_cm %in% 6)
sla_all_soil.red2 <- sla_all_soil[,c(2:5,8:11,13:37)]
colnames(sla_all_soil.red2)
sla_soil.ag2 <- do.call(data.frame, aggregate(cbind(WeightDry_g, SurfaceArea_cm2, SLA) ~  
                                               Position.E.to.W_m + Transect + Species + species2 + Tree + Depth_cm + SlopePercent + area.cm.2 + 
                                               Individuals + Spruce.area.cm.2 + spruce.area.frac + ald_cm + ald_class2 + npp_gC + 
                                               Hum_29 + Hum_16 + Hum_15 + Ow.mean + Ow.dry.mean + Ov.mean + pH.mean + BD.mean + 
                                               N_per + C_per + S_per + CN + Temperature + MossDepth,
                                             sla_all_soil.red2, function(x) c(mean=mean(x), sd = sd(x))))


sla_soil.ag.alder <- subset(sla_soil.ag, species2 %in% "Alder")
sla_soil.ag.spruce <- subset(sla_soil.ag, species2 %in% "Spruce")

# -----------------------------------------------------------------------------
## FIGURES
# -----------------------------------------------------------------------------
## Histograms of SLA
# alder
ggplot(sla_alder) +
  aes(x=SLA) +
  geom_histogram() +
  xlab(expression(alder~SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# spruce
ggplot(sla_spruce) +
  aes(x=SLA) +
  geom_histogram() +
  xlab(expression(spruce~SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# SLA across field site (with position reversed to match figure of field site)
ggplot(sla_all) +
  aes(x=Position.E.to.W_m, y=SLA, color=species2, fill=species2) +
  geom_point(position="jitter") +
  xlab("Landscape position (m)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  scale_x_reverse( lim=c(80,0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title=element_blank(),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## Figure: Slope vs. ALD
ggplot(subset(sla_all.nosla)) +
  aes(x=SlopePercent, y=ald_cm) +
  geom_point() +
#  geom_smooth(method="lm") +
  geom_vline(xintercept = 22.9, linetype="solid", color="gray", size=1) +
  xlab("Slope (%)") + ylab("ALD (cm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.nosla, ald_cm < 151)) + #for plot inset
  aes(x=SlopePercent, y=ald_cm) +
  geom_point(size=5) +
#  geom_smooth(method="lm", size=3) +
  xlab("Slope (%)") + ylab("ALD (cm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 35),
#        axis.title = element_text(size = 35),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="red", size=6))

# -----------------------------------------------------------------------------
## Figure: ALD vs. Moss depth
plot1 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=ald_cm, y=MossDepth) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_vline(xintercept = 140, linetype="longdash") +
  xlab("ALD (cm)") + ylab("Moss depth (cm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Figure: Soil temperature vs. Moss depth
plot2 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=MossDepth, y=Temperature) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("Moss depth (cm)") + ylab(expression(paste("Soil temperature (",degree,"C)"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Figure: Soil temperature vs. ALD
plot3 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=ald_cm, y=Temperature) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_vline(xintercept = 140, linetype="longdash") +
  xlab("ALD (cm)") + ylab(expression(paste("Soil temperature (",degree,"C)"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(plot1, plot2, plot3, ncol=3, heights=c(3,3,3), widths=c(3,3,3))

# -----------------------------------------------------------------------------
## Figure: C:N vs. soil temperature
plot4 <- ggplot(subset(sla_all_soil.nosla2)) +
  aes(x=Temperature, y=CN) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab(expression(paste("Soil temperature (",degree,"C)"))) + ylab("Soil C:N") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Figure: C:N vs. Moss-organic layer depth
plot5 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=(MossDepth), y=(CN)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("Moss depth (cm)") + ylab("Soil C:N") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Figure: C:N vs. soil moisture
plot7 <- ggplot(subset(sla_all_soil.nosla2)) +
  aes(x=(Ow.mean), y=(CN)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("Soil moisture") + ylab("Soil C:N") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Figure: plot4, plot5
grid.arrange(plot4, plot5, ncol=2, heights=c(3,3), widths=c(3,3))

# Figure: plot4, plot5, plot7
fig <- grid.arrange(plot7, plot4, plot5, ncol=3, heights=c(3,3,3), widths=c(3,3,3))

pdf("soilCN_panel.pdf", width = 10, height = 10, useDingbats=FALSE) # Open a new pdf file
grid.arrange(plot7, plot4, plot5, ncol=3, heights=c(3,3,3), widths=c(3,3,3)) # Write the grid.arrange in the file
dev.off() # Close the file

# -----------------------------------------------------------------------------
## Figure: ALD vs. soil moisture
plot6 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=(ald_cm), y=(Ow.mean)) +
  geom_point() +
#  geom_smooth(method="lm") +
  geom_vline(xintercept = 140, linetype="longdash") +
  xlab("ALD (cm)") + ylab("Soil moisture") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all_soil, Depth_cm %in% 6)) +
  aes(x=ald_class2, y=Ow.mean) +
  geom_boxplot() +
  ggtitle("Soil Moisture by ALD Class") +
  xlab("ALD Class") + ylab("Soil Moisture") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Figure: plot6, plot7
grid.arrange(plot6, plot7, ncol=2, heights=c(3,3), widths=c(3,3))

# Figure: plot6 (soil moisture vs. ALD)
#pdf("./moisture_ald.pdf", width=3, height=3)
#print(plot6)
#graphics.off()

# -----------------------------------------------------------------------------
## Figure: SLA vs. C:N
ggplot(subset(sla_soil.ag.alder)) +
  aes(x=(CN), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. C:N (alder)") +
  xlab("Soil C:N") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_soil.ag.spruce)) +
  aes(x=(CN), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. C:N (spruce)") +
  xlab("Soil C:N") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Figure: SLA vs. soil moisture
ggplot(subset(sla_soil.ag.alder)) +
  aes(x=(Ow.mean), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. soil moisture (alder)") +
  xlab("Soil Moisture") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_soil.ag.spruce)) +
  aes(x=(Ow.mean), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. soil moisture (spruce)") +
  xlab("Soil Moisture") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


########
# rbind sla_soil.ag.alder and sla_soil.ag.spruce
sla_soil.ag.both <- rbind(sla_soil.ag.alder, sla_soil.ag.spruce)
########

plot8 <- ggplot(sla_soil.ag.both) +
  aes(x=Ow.mean, y=SLA.mean) +
  geom_point() +
  #geom_smooth(method="lm") +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1) + #http://statistics.ats.ucla.edu/stat/r/faq/smooths.htm
  facet_wrap(~species2, scales="free") +
  xlab("Soil moisture") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18))

plot9 <- ggplot(sla_soil.ag.both) +
  aes(x=CN, y=SLA.mean) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~species2, scales="free") +
  xlab("Soil C:N") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18))

grid.newpage()
grid.draw(rbind(ggplotGrob(plot8), ggplotGrob(plot9), size = "last"))

# -----------------------------------------------------------------------------
## Figure: SLA vs. ALD
ggplot(subset(sla_all, Species %in% "ALSP")) + #Species %in% c("PIMA","PIGL")
  aes(x=(ald_cm), y=(SLA)) + #, color=ald_class2
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. ALD (alder)") +
  xlab("ALD (cm)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all, Species %in% "ALSP")) + Species %in% c("PIMA","PIGL") +
  aes(x=ald_class2, y=SLA) +
  geom_boxplot() +
  ggtitle("SLA by ALD Class (alder)") +
  xlab("ALD Class") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.ag.alder)) +
  aes(x=(ald_cm), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. ALD (alder)") +
  xlab("ALD (cm)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.ag.spruce)) +
  aes(x=(ald_cm), y=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. ALD (spruce)") +
  xlab("ALD (cm)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## Figure: SLA vs. Slope
ggplot(subset(sla_all.ag.alder)) +
  aes(x=SlopePercent, y=SLA.mean) +
  geom_point() +
  geom_smooth() +
  ggtitle("SLA vs. Slope (alder)") +
  xlab("Slope (%)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.ag.spruce)) +
  aes(x=SlopePercent, y=SLA.mean) +
  geom_point() +
  geom_smooth() +
  ggtitle("SLA vs. Slope (spruce)") +
  xlab("Slope (%)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


########
# rbind sla_all.ag.alder and sla_all.ag.spruce
sla_all.ag.both <- rbind(sla_all.ag.alder, sla_all.ag.spruce)
########

plot10 <- ggplot(sla_all.ag.both) + ###4 x 7.29inches for pdf export
  aes(x=ald_cm, y=SLA.mean) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_vline(xintercept = 140, linetype="longdash") +
  facet_wrap(~species2, scales="free") +
  xlab("ALD (cm)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18))

plot11 <- ggplot(sla_all.ag.both) +
  aes(x=SlopePercent, y=SLA.mean) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~species2, scales="free") +
  xlab("Slope (%)") + ylab(expression(SLA~(cm^2~g^{-1}))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18))

grid.newpage()
grid.draw(rbind(ggplotGrob(plot10), ggplotGrob(plot11), size = "last"))
# -----------------------------------------------------------------------------
## Figure: SLA vs. NPP
ggplot(subset(sla_all, Species %in% "ALSP")) + #Species %in% c("PIMA","PIGL")
  aes(x=(SLA), y=(npp_gC)) +
  geom_point(aes(color=ald_class2)) +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. NPP (alder)") +
  xlab("Log(SLA)") + ylab("Log(NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_discrete(name="ALD Class")

ggplot(subset(sla_all.ag.alder)) +
  aes(x=log(SLA.mean), y=log(npp_gC)) +
  geom_point(aes()) +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. NPP (alder)") +
  xlab("Log(SLA)") + ylab("Log(NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.ag.spruce)) +
  aes(x=log(SLA.mean), y=log(npp_gC)) +
  geom_point(aes()) + #color=ald_class2
  geom_smooth(method="lm") +
  ggtitle("SLA vs. NPP (spruce)") +
  xlab("Log(SLA)") + ylab("Log(NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(sla_all.ag)) + #SLA vs. NPP, for Ben (log transformed)
  aes(x=log(SLA.mean), y=log(npp_gC)) +
  geom_point(aes()) +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. NPP") +
  xlab("Log(SLA)") + ylab("Log(NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18)) +
  facet_wrap(~species2, scales="free")

ggplot(subset(sla_all.ag)) + #SLA vs. NPP, for Ben (not transformed)
  aes(x=(SLA.mean), y=(npp_gC)) +
  geom_point(aes()) +
  geom_smooth(method="lm") +
  ggtitle("SLA vs. NPP") +
  xlab("SLA") + ylab("NPP") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18)) +
  facet_wrap(~species2, scales="free")

# Panel plot of NPP
ggplot(sla_all.ag.both) +
  aes(x=log(SLA.mean), y=log(npp_gC)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~species2, scales="free") +
  xlab("Log(SLA)") + ylab("Log(Total NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=18))

# -----------------------------------------------------------------------------
# NPP for black spruce (PIMA)
# -----------------------------------------------------------------------------
sla.spruce.interp <- read.csv("CPCRW_Data_2014/CPCRW_NPP/CPCRW_2014_NPP_Linear_Interp.pima.csv")
sla.spruce.55 <- subset(sla_all.ag.spruce, Position.E.to.W_m < 60)
sla.spruce.55.merge <- merge(sla.spruce.55, sla.spruce.interp)

sla.spruce.soil.55 <- subset(sla_soil.ag.spruce, Position.E.to.W_m < 60)
sla.spruce.55.soil.merge <- merge(sla.spruce.soil.55, sla.spruce.interp)

ggplot(subset(sla.spruce.55.merge, Species %in% "PIMA")) +
  aes(y=log(npp_gC), x=log(SLA.mean)) +
  geom_point() +
  geom_smooth() +
  ggtitle("PIMA SLA vs. total NPP") +
  xlab("ln(SLA)") + ylab("ln(NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_discrete(name="ALD Class")

plot12 <- ggplot(subset(sla.spruce.55.merge, Species %in% "PIMA")) +
  aes(y=log(npp_gC_m2), x=log(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
#  ggtitle("Spruce only") +
  xlab("ln(Spruce SLA)") + ylab("ln(Spruce NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot12b <- ggplot(subset(sla.spruce.55.merge, Species %in% "PIMA")) +
  aes(y=(npp_gC_m2), x=(SLA.mean)) +
  geom_point() +
  geom_smooth(method="lm") +
  #  ggtitle("Spruce only") +
#  xlab("Spruce SLA") + ylab("Spruce NPP") +
  xlab(expression(Spruce~SLA~(cm^2~g^{-1}))) + ylab(expression(Spruce~NPP~(g~C~m^2~yr^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#R-sq:
summary(lm(log(SLA.mean) ~ log(npp_gC_m2), data=sla.spruce.55.merge))
summary(lm((SLA.mean) ~ (npp_gC_m2), data=sla.spruce.55.merge))

# -----------------------------------------------------------------------------
# NPP for broadleaf (Alder, Birch)
# -----------------------------------------------------------------------------
sla.bepa.interp <- read.csv("CPCRW_Data_2014/CPCRW_NPP/CPCRW_2014_NPP_Linear_Interp.bepa.csv")
sla.bepa.interp$Position.E.to.W_m[sla.bepa.interp$Position.E.to.W_m==70] <- 72
sla.bepa.merge <- merge(sla_all.ag.alder, sla.bepa.interp)

ggplot(subset(sla.bepa.merge)) +
  aes(x=log(npp_gC), y=log(SLA.mean)) +
  geom_point() +
  geom_smooth()

ggplot(subset(sla.bepa.merge)) +
  aes(x=log(npp_gC_m2), y=log(SLA.mean)) +
  geom_point() +
  geom_smooth()

# -----------------------------------------------------------------------------
# NPP vs. Soil C:N
plot13 <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=log(CN), y=log(npp_gC)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("ln(Soil C:N)") + ylab("ln(Total NPP)") +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot13b <- ggplot(subset(sla_all_soil.nosla2, Depth_cm %in% 6)) +
  aes(x=(CN), y=(npp_gC)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("Soil C:N") + ylab(expression(Total~NPP~(g~C~m^2~yr^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


sla_all_soil.nosla2b <- subset(sla_all_soil.nosla2, Depth_cm %in% 6)

summary(lm(log(CN) ~ log(npp_gC), data=sla_all_soil.nosla2b))
summary(lm(CN ~ npp_gC, data=sla_all_soil.nosla2b))

grid.arrange(plot12b, plot13b, ncol=2, heights=c(3,3), widths=c(3,3))


ag.soil.depths <- do.call(data.frame, aggregate(cbind(CN, npp_gC) ~ Transect + Position.E.to.W_m, sla_all_soil.nosla2, function(x) c(
  min=min(x), 
  max=max(x), 
  mean=mean(x))))

ggplot(ag.soil.depths) +
  aes(x=CN.mean, y=npp_gC.mean) +
  geom_point() +
#  geom_smooth() +
  geom_vline(xintercept = 30.5, linetype="longdash")

# -----------------------------------------------------------------------------
## Summary tables
# -----------------------------------------------------------------------------
# Summary table: SLA (mean, sd) by species and ALD class.
ag.summary.all <- do.call(data.frame, aggregate(SLA ~ species2, sla_all, function(x) c(mean=mean(x), sd = sd(x))))

ag.summary.all2 <- do.call(data.frame, aggregate(SLA ~ species2, sla_all, function(x) c(
  min=min(x), 
  max=max(x), 
  mean=mean(x), 
  sd=sd(x), 
  CV=sd(x)/mean(x))))

# Summary table: SLA (mean, sd) by species and ALD class.
ag.summary.ald <- do.call(data.frame, aggregate(SLA ~ species2 + ald_class2, sla_all, function(x) c(mean=mean(x), sd = sd(x))))

# Summary tables for soil temp, C:N, moss depth, etc.
colnames(sla_all_soil.nosla2.6cm)
tmp <- lapply(sla_all_soil.nosla2.6cm[,c(8:18)], function(x) rbind(
  mean=mean(x),
  minimum=min(x),
  maximum=max(x)))
sla_soil.ag.summary <- data.frame(tmp)

ag.summary.soil.all <- do.call(data.frame, aggregate(cbind(BD.mean, Ow.mean, CN, Temperature) ~  
                                              Depth_cm,
                                            sla_all_soil.nosla2, function(x) c(
                                              mean=mean(x),
                                              minimum=min(x),
                                              maximum=max(x))))

# Summary tables for site-wide data (e.g. ALD, tree density, NPP, slope, etc.)
colnames(sla_all.nosla)
tmp2 <- lapply(sla_all.nosla[,c(4:10)], function(x) rbind(
  mean=mean(x),
  minimum=min(x),
  maximum=max(x),
  sd=sd(x)))
sla_all.ag.summary <- data.frame(tmp2)

# -----------------------------------------------------------------------------
## 1. ASSESS VARIATION WITHIN AND AMONG INDIVIDUAL TREES (NESTED LME WITH LEAVES)
# -----------------------------------------------------------------------------
# Assess variation within individual trees (nesting with lme, e.g. Messier et al. 2010)
varcomp.LME.alder <- varcomp(lme(log(SLA) ~ 1, random = ~1|Tree, data=sla_alder, na.action=na.omit),1)
varcomp.LME.alder

varcomp.LME.spruce <- varcomp(lme(log(SLA) ~ 1, random = ~1|Tree, data=sla_spruce, na.action=na.omit),1)
varcomp.LME.spruce

#varcomp.LME.all <- varcomp(lme(log(SLA) ~ 1, random = ~1|species2/Tree, data=sla_all_soil.red, na.action=na.omit),1)
#varcomp.LME.all

# -----------------------------------------------------------------------------
## 2. LINEAR MODELS WITH NO NESTING (USES AVERAGES OF REPLICATES)
# -----------------------------------------------------------------------------
# Linear models according to SLA schematic. For each species.
# QUESTION: Is there a relationship between slope and ALD? Use full dataset (sla_all.ag.spruce).

plot(ald_cm ~ SlopePercent, data=sla_all.nosla)
lm1 <- lm(ald_cm ~ SlopePercent, data = sla_all.nosla)
summary(lm1)
par(mfrow = c(2, 2))
plot(lm1)

##### Segmented regression
seg.mod.ald <- segmented(lm1, seg.Z = ~SlopePercent, psi=20)
summary(seg.mod.ald)
seg.mod.ald$psi
slope(seg.mod.ald)
my.fitted3 <- fitted(seg.mod.ald)
my.model3 <- data.frame(Slope = sla_all.nosla$SlopePercent, ALD = my.fitted3)
ggplot(my.model3, aes(x=Slope, y=ALD)) + geom_line()




cor.test(sla_all.nosla$ald_cm, sla_all.nosla$SlopePercent, method="spearman", exact=FALSE)
spearman_test(sla_all.nosla$ald_cm ~ sla_all.nosla$SlopePercent, distribution="asymptotic", ties.method="mid-ranks")

#Subset with slope < 23 to look at ald < 151 (threshold)
sla_all.ag.unique.lowslope <- subset(sla_all.nosla, SlopePercent < 23)
plot(ald_cm ~ SlopePercent, data=sla_all.ag.unique.lowslope)
summary(lm(ald_cm ~ SlopePercent, data = sla_all.ag.unique.lowslope))

cor.test(sla_all.ag.unique.lowslope$ald_cm, sla_all.ag.unique.lowslope$SlopePercent, method="spearman", exact=FALSE)
spearman_test(sla_all.ag.unique.lowslope$ald_cm ~ sla_all.ag.unique.lowslope$SlopePercent, distribution="asymptotic", ties.method="mid-ranks")
summary(lm(ald_cm ~ SlopePercent, data=sla_all.ag.unique.lowslope))

#Since ALD has an upper lmit of 150cm (>150 keyed to 151 in data), need to to a Tobit regression to reflect this upper limit.
#http://www.ats.ucla.edu/stat/r/dae/tobit.htm
#http://www1.karlin.mff.cuni.cz/~pesta/NMFM404/tobit.html#Tobit_model
#http://www.rdocumentation.org/packages/VGAM/functions/tobit
#http://stackoverflow.com/questions/19014122/pseudo-r2-using-vglm
#https://cran.r-project.org/web/packages/censReg/vignettes/censReg.pdf
ald.slope.tobit <- vglm(ald_cm ~ SlopePercent, tobit(Upper=151), data=sla_all.nosla)
summary(ald.slope.tobit)
ctable <- coef(summary(ald.slope.tobit))
pvals <- 2 * pt(abs(ctable[, "z value"]), df.residual(ald.slope.tobit), lower.tail=FALSE)
cbind(ctable, pvals)  
#To examine how well our model fits the data; plot of residuals to assess their absolute & relative (pearson) values and assumptions
sla_all.nosla$yhat <- fitted(ald.slope.tobit)[,1]
sla_all.nosla$rr <- resid(ald.slope.tobit, type = "response")
sla_all.nosla$rp <- resid(ald.slope.tobit, type = "pearson")[,1]

par(mfcol = c(2, 3))
with(sla_all.nosla, {
  plot(yhat, rr, main = "Fitted vs Residuals")
  qqnorm(rr)
  plot(yhat, rp, main = "Fitted vs Pearson Residuals")
  qqnorm(rp)
  plot(ald_cm, rp, main = "Actual vs Pearson Residuals")
  plot(ald_cm, yhat, main = "Actual vs Fitted")
})

summary(lm(yhat ~ ald_cm, data=sla_all.nosla))
(r <- with(sla_all.nosla, cor(yhat, ald_cm)))
r^2
plot(yhat ~ ald_cm, data=sla_all.nosla)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between ALD and the depth of the moss-organic layer?
sla_all_soil.nosla2.6cm <- subset(sla_all_soil.nosla2, Depth_cm %in% 6)

plot(MossDepth ~ ald_cm, data=sla_all_soil.nosla2.6cm)
lm2 <- (lm(MossDepth ~ ald_cm, sla_all_soil.nosla2.6cm))
summary(lm2)
par(mfrow = c(2, 2))
plot(lm2)

shapiro.test(sla_all_soil.nosla2.6cm$ald_cm) #ald not normal distribution
shapiro.test(sla_all_soil.nosla2.6cm$MossDepth)

t.test(sla_all_soil.nosla2.6cm$MossDepth ~ sla_all_soil.nosla2.6cm$ald_class2)
#var.equal = TRUE option to specify equal variances and a pooled variance estimate.
#alternative="less" or alternative="greater" option to specify a one tailed test.

var.test(sla_all_soil.nosla2.6cm$MossDepth[sla_all_soil.nosla2.6cm$ald_class2=="Shallow"],
         sla_all_soil.nosla2.6cm$MossDepth[sla_all_soil.nosla2.6cm$ald_class2=="Deep"])

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between soil temperature and moss depth? Use soil data set (sla_soil.ag.alder/spruce).
plot(Temperature ~ MossDepth, data = sla_soil.ag.spruce)
lm4 <- (lm(Temperature ~ MossDepth, data = sla_soil.ag.spruce))
summary(lm4)
par(mfrow = c(2, 2))
plot(lm4)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between soil temperature and ALD? Use soil data set (sla_soil.ag.alder/spruce).
plot(Temperature ~ ald_cm, data = sla_soil.ag.spruce)
lm3 <- (lm(Temperature ~ ald_cm, data = sla_soil.ag.spruce))
summary(lm3)
par(mfrow = c(2, 2))
plot(lm3)

shapiro.test(sla_all_soil.nosla2.6cm$Temperature)

t.test(sla_all_soil.nosla2.6cm$Temperature ~ sla_all_soil.nosla2.6cm$ald_class2)
#var.equal = TRUE option to specify equal variances and a pooled variance estimate.
#alternative="less" or alternative="greater" option to specify a one tailed test.

var.test(sla_all_soil.nosla2.6cm$Temperature[sla_all_soil.nosla2.6cm$ald_class2=="Shallow"],
         sla_all_soil.nosla2.6cm$Temperature[sla_all_soil.nosla2.6cm$ald_class2=="Deep"])

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between soil C:N and soil temperature?
plot(CN ~ Temperature, data = sla_all_soil.nosla2)
lm5 <- (lm(CN ~ Temperature, data = sla_all_soil.nosla2))
summary(lm5)
par(mfrow = c(2, 2))
plot(lm5)

# -----------------------------------------------------------------------------
# QUESTION: Indirect: Is there a relationship between soil C:N and the depth of the moss-organic layer? Use soil data set (sla_soil.ag.alder/spruce).
plot(CN ~ MossDepth, data = sla_all_soil.nosla2.6cm)
lm6 <- (lm(CN ~ MossDepth, data = sla_all_soil.nosla2.6cm))
summary(lm6)
par(mfrow = c(2, 2))
plot(lm6)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between soil moisture and ALD? Use soil data set (sla_soil.ag.alder/spruce).
plot(Ow.mean ~ ald_cm, data = sla_all_soil.nosla2.6cm)

lm7 <- (lm(Ow.mean ~ ald_cm, data = sla_all_soil.nosla2.6cm))
summary(lm7)
par(mfrow = c(2, 2))
plot(lm7)

shapiro.test(sla_all_soil.nosla2.6cm$Ow.mean)

t.test(sla_all_soil.nosla2.6cm$Ow.mean ~ sla_all_soil.nosla2.6cm$ald_class2)
#var.equal = TRUE option to specify equal variances and a pooled variance estimate.
#alternative="less" or alternative="greater" option to specify a one tailed test.

var.test(sla_all_soil.nosla2.6cm$Ow.mean[sla_all_soil.nosla2.6cm$ald_class2=="Shallow"],
         sla_all_soil.nosla2.6cm$Ow.mean[sla_all_soil.nosla2.6cm$ald_class2=="Deep"])

##### Testing segmented regression
seg.mod.ald2 <- segmented(lm7, seg.Z = ~ald_cm, psi=105)
summary(seg.mod.ald2)
seg.mod.ald2$psi
slope(seg.mod.ald2)
my.fitted3 <- fitted(seg.mod.ald2)
my.model3 <- data.frame(Slope = sla_all_soil.nosla2.6cm$ald_cm, ALD = my.fitted3)
ggplot(my.model3, aes(x=Slope, y=ALD)) + geom_line()

davies.test(lm7, seg.Z = ~ald_cm, k=3)




# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between soil C:N and soil moisture?
plot(CN ~ Ow.mean, data = sla_soil.ag2[sla_soil.ag2$species2 == "Spruce",])
lm8 <- (lm(CN ~ Ow.mean, data = sla_soil.ag2[sla_soil.ag2$species2 == "Spruce",]))
summary(lm8)
par(mfrow = c(2, 2))
plot(lm8)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between SLA and nutrient availability (via soil C:N)? Use soil data set (sla_soil.ag.alder/spruce).
plot(SLA.mean ~ CN, data = sla_soil.ag.alder)
lm9 <- (lm(SLA.mean ~ CN, data = sla_soil.ag.alder))
summary(lm9)
par(mfrow = c(2, 2))
plot(lm9)

plot(SLA.mean ~ CN, data = sla_soil.ag.spruce)
lm10 <- (lm(SLA.mean ~ CN, data = sla_soil.ag.spruce))
summary(lm10)
par(mfrow = c(2, 2))
plot(lm10)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between SLA and soil moisture? Use soil data set (sla_soil.ag.alder/spruce).
plot(SLA.mean ~ Ow.mean, data = sla_soil.ag.alder)
lm11 <- (lm(SLA.mean ~ Ow.mean, data = sla_soil.ag.alder))
summary(lm11)
par(mfrow = c(2, 2))
plot(lm11)
AIC(lm11)
AICc(lm11) #AICc (small sample size corrected)

p1 <- ggplot(lm11, aes(.fitted, .resid)) +
  geom_point() + stat_smooth(method="loess") + geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
#  ggtitle("Residual vs Fitted Plot")+theme_bw() +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(SLA.mean ~ Ow.mean, data = sla_soil.ag.spruce)
lm12 <- (lm(SLA.mean ~ Ow.mean, data = sla_soil.ag.spruce))
summary(lm12)
par(mfrow = c(2, 2))
plot(lm12)
AIC(lm12)
AICc(lm12) #AICc (small sample size corrected)

p2 <- ggplot(lm12, aes(.fitted, .resid)) +
  geom_point() + stat_smooth(method="loess") + geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  #  ggtitle("Residual vs Fitted Plot")+theme_bw() +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(p1, p2, ncol=2, heights=c(3,3), widths=c(3,3))


##### Testing segmented regression for breakpoints
# https://rpubs.com/MarkusLoew/12164
# https://climateecology.wordpress.com/2012/08/19/r-for-ecologists-putting-together-a-piecewise-regression/
# davies.test tests for a non-zero difference-in-slope parameter of a segmented relationship. 

lm11 <- (lm(SLA.mean ~ Ow.mean, data = sla_soil.ag.alder))
lm12 <- (lm(SLA.mean ~ Ow.mean, data = sla_soil.ag.spruce))

# alder
seg.mod.alder <- segmented(lm11, seg.Z = ~Ow.mean, psi=0.7)
summary(seg.mod.alder)
seg.mod.alder$psi
slope(seg.mod.alder)
my.fitted <- fitted(seg.mod.alder)
my.model <- data.frame(Ow = sla_soil.ag.alder$Ow.mean, SLA = my.fitted)
ggplot(my.model, aes(x=Ow, y=SLA)) + geom_line()
davies.test(lm11, seg.Z = ~Ow.mean, k=2)

# spruce
seg.mod.spruce <- segmented(lm12, seg.Z = ~Ow.mean, psi=0.7)
summary(seg.mod.spruce)
seg.mod.spruce$psi
slope(seg.mod.spruce)
my.fitted2 <- fitted(seg.mod.spruce)
my.model2 <- data.frame(Ow = sla_soil.ag.spruce$Ow.mean, SLA = my.fitted2)
ggplot(my.model2, aes(x=Ow, y=SLA)) + geom_line()
davies.test(lm12, seg.Z = ~Ow.mean, k=2)




##### Quadratic regression
# alder
Ow.mean2 <- sla_soil.ag.alder$Ow.mean^2
quad.mod1 <- (lm(SLA.mean ~ Ow.mean + Ow.mean2, data = sla_soil.ag.alder))
summary(quad.mod1)
AIC(quad.mod1)
AICc(quad.mod1)

moisture.values <- seq(0, 1, 0.05)
predicted.sla.alder <- predict(quad.mod1, list(Ow.mean = moisture.values, Ow.mean2 = moisture.values^2, data=sla_soil.ag.alder))
plot(SLA.mean ~ Ow.mean, pch=16, xlab = "Ow.mean", ylab = "SLA.mean", cex.lab = 1.3, col = "blue", data=sla_soil.ag.alder)
lines(moisture.values, predicted.sla.alder, col = "darkgreen", lwd = 3)


# spruce
Ow.mean2b <- sla_soil.ag.spruce$Ow.mean^2
quad.mod2 <- (lm(SLA.mean ~ Ow.mean + Ow.mean2b, data = sla_soil.ag.spruce))
summary(quad.mod2)
AIC(quad.mod2)
AICc(quad.mod2)

moisture.values <- seq(0, 1, 0.05)
predicted.sla.spruce <- predict(quad.mod2, list(Ow.mean = moisture.values, Ow.mean2b = moisture.values^2, data=sla_soil.ag.spruce))
plot(SLA.mean ~ Ow.mean, pch=16, xlab = "Ow.mean", ylab = "SLA.mean", cex.lab = 1.3, col = "blue", data=sla_soil.ag.spruce)
lines(moisture.values, predicted.sla.spruce, col = "darkgreen", lwd = 3)




# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between SLA and ALD?  Use full dataset (sla_all.ag.alder/spruce).
plot(SLA.mean ~ ald_cm, data = sla_all.ag.alder)
lm17 <- (lm(SLA.mean ~ ald_cm, data = sla_all.ag.alder))
summary(lm17)
par(mfrow = c(2, 2))
plot(lm17)

shapiro.test(sla_all.ag.alder$SLA.mean)

t.test(sla_all.ag.alder$SLA.mean ~ sla_all.ag.alder$ald_class2)
#var.equal = TRUE option to specify equal variances and a pooled variance estimate.
#alternative="less" or alternative="greater" option to specify a one tailed test.

var.test(sla_all.ag.alder$SLA.mean[sla_all.ag.alder$ald_class2=="Shallow"],
         sla_all.ag.alder$SLA.mean[sla_all.ag.alder$ald_class2=="Deep"])

#Spearman, SLA vs. ALD (alder)
cor.test(sla_all.ag.alder$SLA.mean, sla_all.ag.alder$ald_cm, method="spearman", exact=FALSE)



plot(SLA.mean ~ ald_cm, data = sla_all.ag.spruce)
lm18 <- (lm(SLA.mean ~ ald_cm, data = sla_all.ag.spruce))
summary(lm18)
par(mfrow = c(2, 2))
plot(lm18)

shapiro.test(sla_all.ag.spruce$SLA.mean)

t.test(sla_all.ag.spruce$SLA.mean ~ sla_all.ag.spruce$ald_class2)
#var.equal = TRUE option to specify equal variances and a pooled variance estimate.
#alternative="less" or alternative="greater" option to specify a one tailed test.

var.test(sla_all.ag.spruce$SLA.mean[sla_all.ag.spruce$ald_class2=="Shallow"],
         sla_all.ag.spruce$SLA.mean[sla_all.ag.spruce$ald_class2=="Deep"])


#Spearman, SLA vs. ALD (spruce)
cor.test(sla_all.ag.spruce$SLA.mean, sla_all.ag.spruce$ald_cm, method="spearman", exact=FALSE)

# -----------------------------------------------------------------------------
# QUESTION: Is there a relationship between SLA and NPP? Use full dataset (sla_all.ag.alder/spruce).
###LOG TRANSFORM?
plot(npp_gC ~ SLA.mean, data = sla_all.ag.alder)
lm15 <- (lm(log(npp_gC) ~ log(SLA.mean), data = sla_all.ag.alder))
summary(lm15)
par(mfrow = c(2, 2))
plot(lm15)

plot(npp_gC ~ SLA.mean, data = sla_all.ag.spruce)
lm16 <- (lm(log(npp_gC) ~ log(SLA.mean), data = sla_all.ag.spruce))
summary(lm16)
par(mfrow = c(2, 2))
plot(lm16)

# -----------------------------------------------------------------------------
# QUESTION: Indirect: Is there a relationship between slope and SLA? Use full dataset (sla_all.ag.alder/spruce).
plot(SLA.mean ~ SlopePercent, data = sla_all.ag.alder)
lm19 <- (lm(SLA.mean ~ SlopePercent, data = sla_all.ag.alder))
summary(lm19)
par(mfrow = c(2, 2))
plot(lm19)

plot(SLA.mean ~ SlopePercent, data = sla_all.ag.spruce)
lm20 <- (lm(SLA.mean ~ SlopePercent, data = sla_all.ag.spruce))
summary(lm20)
par(mfrow = c(2, 2))
plot(lm20)

# -----------------------------------------------------------------------------
## 3. PATH ANALYSIS
# -----------------------------------------------------------------------------
##### From Grace et al. 2015 (book chapter); calculating the number of samples per parameter (d); 2.5 is marginal; 
# reasoning: sample adequacy depends on model complexity
samples = 18 #or 36
parameter = 6 #or 7, 8, depending on full model
d = samples/parameter
d

##### piecewiseSEM
# https://jslefche.files.wordpress.com/2016/02/5_local-estimation3.pdf

library(devtools)
#install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)

# Data to use
path.alder <- sla_soil.ag.alder
path.spruce <- sla_soil.ag.spruce

# model with SLA, moisture, C:N, ALD
modelList = list(
  lm(SLA.mean ~ Ow.mean + CN + ald_cm, data=path.alder),
  lm(Ow.mean ~ ald_cm, data=path.alder),
  lm(CN ~ Ow.mean, data=path.alder)
)

sem.fit(modelList, data=path.alder)

# model assumptions
lapply(modelList, plot, which=1)
# individual fits
sem.model.fits(modelList)
# get coefficients
sem.coefs(modelList)
# scaled coefficients
sem.coefs(modelList, path.alder, standardize="scale")

# full model (Fig. 1)
modelList2 = list(
  lm(SLA.mean ~ Ow.mean + CN + ald_cm, data=path.alder),
  lm(Ow.mean ~ ald_cm, data=path.alder),
  lm(CN ~ Ow.mean + Temperature, data=path.alder),
  lm(Temperature ~ MossDepth + ald_cm, data=path.alder),
  lm(ald_cm ~ MossDepth + SlopePercent, data=path.alder)
)

sem.fit(modelList2, data=path.alder)

# model assumptions
lapply(modelList2, plot, which=1)
# individual fits
sem.model.fits(modelList2)
# get coefficients
sem.coefs(modelList2)
# scaled coefficients
sem.coefs(modelList2, path.alder, standardize="scale")


# full model with correlations (feedbacks?) (Fig. 1)
modelList3 = list(
  lm(SLA.mean ~ Ow.mean + CN + ald_cm, data=path.alder),
  lm(Ow.mean ~ ald_cm, data=path.alder),
  lm(CN ~ Ow.mean + Temperature, data=path.alder),
  lm(Temperature ~ MossDepth + ald_cm, data=path.alder),
  lm(ald_cm ~ MossDepth + SlopePercent, data=path.alder)
)

sem.fit(modelList3, path.alder, corr.errors="Temperature ~~ MossDepth")
sem.coefs(modelList3, path.alder, standardize="scale", corr.errors="Temperature ~~ MossDepth")

##### Using Lavaan (https://www.youtube.com/watch?v=-B37sK9NTfI)
library(lavaan)
library(corrplot)

# compare above piecewiseSEM with lavaan
modelList.lavaan = sem.lavaan(modelList, path.alder)
lavaan::summary(modelList.lavaan, rsq=TRUE, standardize=TRUE)
varTable(modelList.lavaan)

# Data to use
path.alder <- sla_soil.ag.alder
path.spruce <- sla_soil.ag.spruce

model.sla <- 'SLA.mean ~ CN + Ow.mean + ald_cm
              Ow.mean ~ ald_cm
              CN ~ Ow.mean'

model.sla2 <- 'SLA.mean ~ ald_cm + CN + Ow.mean
                CN ~ Ow.mean
                Ow.mean ~ ald_cm'

fit <- sem(model.sla, data=path.spruce)
fit <- sem(model.sla2, data=path.spruce)

summary(fit)
summary(fit, standardized=TRUE, fit.measures=TRUE, rsq=TRUE)

model.sla2 <- '0.489*SLA.mean ~ ald_cm + 0.345*CN + 0.435*Ow.mean
                0.345*CN ~ Ow.mean
                0.435*Ow.mean ~ ald_cm'
datsim <- simulateData(model.sla2, model.type='sem', sample.nobs=c(30,30))

# http://onlinelibrary.wiley.com/doi/10.1890/ES12-00048.1/epdf
# https://jonlefcheck.net/2014/07/06/piecewise-structural-equation-modeling-in-ecological-research/
# http://pareonline.net/getvn.asp?v=19&n=12
# http://web.stanford.edu/class/psych253/section/section_8/section8.html
plot_matrix <- function(matrix_toplot){
  corrplot(matrix_toplot, is.corr = FALSE, 
           type = 'lower', 
           order = "original", 
           tl.col='black', tl.cex=.75)
}

plot_matrix(path.alder)
alder.cov <- getCov(path.alder, names=labels1)
