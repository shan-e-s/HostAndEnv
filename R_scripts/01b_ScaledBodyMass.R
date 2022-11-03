# Scaled body mass, 16/4/2019, Shane Somers
## Calc scaled body mass for birds, need different methods for adult vs nestling birds

## follows procedure from Peig and Green (2009) https://onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0706.2009.17643.x
# 1) Plot and consider removing outliers
# 2) the best fit line to the remaining points is obtained by the standardised major axis (SMA) regression on ln‐transformed data. The slope is the bSMA value used in Eq. 2. 
# Alternatively, bSMA can be computed by dividing the regression coefficient bOLS from the OLS regression on ln‐variables by the Pearson's r correlation coefficient.
# 3) the scaled index is calculated for each individual (including outliers if desired) following Eq. 2 (also shown on the figure). 
# The arithmetic mean of L is a suitable value for L0, and Mi‐Li variables represent the raw data for each individual i. This index adjusts the mass of all
# individuals to that which they would have at length L0 (represented by the vertical dashed line).
## The scaled index can also be expected to adjust the whole body composition of each individual to that which it would have at the new length L0, according to allometry.

# Recommend to use the: L variable which has the strongest correlation with M on a log-log scale, since this is likely to be the L that best explains thatfraction of mass associated with structural size.

#install.packages("lmodel2")
library(lmodel2)
library(dplyr)
library(microbiome)

# rm(list=ls())
########### deprecating this section due to codependence w/ Data_setup
# need to read in raw biom as this is sourced in Data_setup and therefore would be codependent with Data_setup
#phylo.spring.sbm <- readRDS(file = "Data/phylo-spring.rds")
#metadata.sbm <- meta(phylo.spring.sbm) # doesnt maintain Date data-type

# subset to nestlings
#phylo.nestlings.sbm <- subset_samples(phylo.spring.sbm, ageBinned=="1week"|ageBinned=="2week")
#nestlings.meta.sbm <- meta(phylo.nestlings.sbm@sam_data)
# subset to adults
#phylo.adults.sbm <- subset_samples(phylo.spring.sbm, ageBinned=="adult")
#adults.meta.sbm <- meta(phylo.adults.sbm@sam_data)
##############

############## new data read-in
metadata.sbm <- read.csv2("Data/MasterBirdIDSheet.csv", sep = ",")

#filter out blanks ie. failed samples
metadata.sbm <- metadata.sbm[!is.na(metadata.sbm$BIOM.ID),]

# set rownames to match biom file ie. sample names
rownames(metadata.sbm) <- metadata.sbm$BIOM.ID

nestlings.meta.sbm <- subset(metadata.sbm, ageBinned=="1week"|ageBinned=="2week") 
adults.meta.sbm <- subset(metadata.sbm, ageBinned=="adult")
#############

######################## Adapted from GT_mb5
# Extract morphometric data for calculating Scaled BMI. Scaled mass index adjusts all individuals mass to that which they would have if they had the same body size. 
# Method from Peig and green 2009, uses linear regression of log-body mass against log-tarsus length estimated by type-2 (standardized major axis) regression.
# Unclear what type-2 means.. Hird calcs SMI as body mass × (mean(tarsus)/tarsus length)^1.87 (Peig and Green, 2009). 
# NOTE: need to use as.character before as.numeric to get correct values due to R storing factors internally as integers.
# Need separate calcs for adults (using wing) and chicks (using tarsus)

######################################################## 
## Overall scaled mass for both chicks and adults (probably not useful)
# type 2 regression, use standardized major axis
lmod2 <- lmodel2(log(as.numeric(as.character(metadata.sbm$Weight)))~log(as.numeric(as.character(metadata.sbm$Tarsus))))
sma <- lmod2$regression.results[[3]][3]
# careful,important to convert to character vector then numeric vector
mean.tarsus <- mean(as.numeric(as.character(metadata.sbm$Tarsus)),na.rm=T)
scaled.mass <- (as.numeric(as.character(metadata.sbm$Weight))) * (mean.tarsus/as.numeric(as.character(metadata.sbm$Tarsus)))^(as.numeric(as.character(sma)))
scaled.mass <- cbind(as.data.frame(scaled.mass),metadata.sbm$BIOM.ID) %>% dplyr::rename(BIOM.ID = "metadata.sbm$BIOM.ID")

######################################################## 
# check L variables
cor(log(as.numeric(as.character(adults.meta.sbm$Weight))),log(as.numeric(as.character(adults.meta.sbm$wing))))
cor(log(as.numeric(as.character(adults.meta.sbm$Weight))),log(as.numeric(as.character(adults.meta.sbm$Tarsus)))) # Tarsus has stronger correlation with mass

summary(lm(log(as.numeric(as.character(adults.meta.sbm$Weight)))~log(as.numeric(as.character(adults.meta.sbm$wing)))))
summary(lm(log(as.numeric(as.character(adults.meta.sbm$Weight)))~log(as.numeric(as.character(adults.meta.sbm$Tarsus)))))

##Scaled mass for adult birds, using tarsus as L variable
lmod2 <- lmodel2(log(as.numeric(as.character(adults.meta.sbm$Weight)))~log(as.numeric(as.character(adults.meta.sbm$Tarsus))))
sma <- lmod2$regression.results[[3]][3]
mean.tarsus.adult <- mean(as.numeric(as.character(adults.meta.sbm$Tarsus)),na.rm=T)
scaled.mass.tarsus.adult <- (as.numeric(as.character(adults.meta.sbm$Weight))) * (mean.tarsus.adult/as.numeric(as.character(adults.meta.sbm$Tarsus)))^(as.numeric(as.character(sma)))

scaled.mass.tarsus.adult <- cbind(as.data.frame(scaled.mass.tarsus.adult),adults.meta.sbm$BIOM.ID) %>% dplyr::rename(BIOM.ID = "adults.meta.sbm$BIOM.ID")

##Scaled mass for adult birds, using wing as L variable
lmod2 <- lmodel2(log(as.numeric(as.character(adults.meta.sbm$Weight)))~log(as.numeric(as.character(adults.meta.sbm$wing))))
sma <- lmod2$regression.results[[3]][3]
mean.wing.adult <- mean(as.numeric(as.character(adults.meta.sbm$wing)),na.rm=T)
scaled.mass.wing.adult <- (as.numeric(as.character(adults.meta.sbm$Weight))) * (mean.wing.adult/as.numeric(as.character(adults.meta.sbm$wing)))^(as.numeric(as.character(sma)))

scaled.mass.wing.adult <- cbind(as.data.frame(scaled.mass.wing.adult),adults.meta.sbm$BIOM.ID) %>% dplyr::rename(BIOM.ID = "adults.meta.sbm$BIOM.ID")

########################################################
# scaled body mass for adults calculating for each sex separately
#using wing as L variable
female.metrics <- adults.meta.sbm[adults.meta.sbm$IDLetter=="F",c("BIOM.ID","Weight", "wing", "Tarsus")]

lmod2 <- lmodel2(log(as.numeric(as.character(female.metrics$Weight)))~log(as.numeric(as.character(female.metrics$wing))))
sma <- lmod2$regression.results[[3]][3]
mean.wing.female <- mean(as.numeric(as.character(female.metrics$wing)),na.rm=T)
scaled.mass.wing.female <- (as.numeric(as.character(female.metrics$Weight))) * (mean.wing.female/as.numeric(as.character(female.metrics$wing)))^(as.numeric(as.character(sma)))

scaled.mass.wing.female.df <- cbind(female.metrics,scaled.mass.wing.female) %>% dplyr::rename("scaled.mass.wing" = "scaled.mass.wing.female")

male.metrics <- adults.meta.sbm[adults.meta.sbm$IDLetter=="M",c("BIOM.ID","Weight", "wing", "Tarsus")]

lmod2 <- lmodel2(log(as.numeric(as.character(male.metrics$Weight)))~log(as.numeric(as.character(male.metrics$wing))))
sma <- lmod2$regression.results[[3]][3]
mean.wing.male <- mean(as.numeric(as.character(male.metrics$wing)),na.rm=T)
scaled.mass.wing.male <- (as.numeric(as.character(male.metrics$Weight))) * (mean.wing.male/as.numeric(as.character(male.metrics$wing)))^(as.numeric(as.character(sma)))

scaled.mass.wing.male.df <- cbind(male.metrics,scaled.mass.wing.male) %>% dplyr::rename("scaled.mass.wing" = "scaled.mass.wing.male")
scaled.mass.wing.df <- rbind(scaled.mass.wing.female.df, scaled.mass.wing.male.df) %>% dplyr::select("BIOM.ID","scaled.mass.wing")

#########
# using Tarsus as L
##Scaled mass for adult birds, using tarsus as L variable
# calculate for males and females separately
lmod2 <- lmodel2(log(as.numeric(as.character(female.metrics$Weight)))~log(as.numeric(as.character(female.metrics$Tarsus))))
sma <- lmod2$regression.results[[3]][3]
mean.tarsus.female <- mean(as.numeric(as.character(female.metrics$Tarsus)),na.rm=T)
scaled.mass.tarsus.female <- (as.numeric(as.character(female.metrics$Weight))) * (mean.tarsus.female/as.numeric(as.character(female.metrics$Tarsus)))^(as.numeric(as.character(sma)))

scaled.mass.tarsus.female.df <- cbind(female.metrics,scaled.mass.tarsus.female) %>% dplyr::rename("scaled.mass.tarsus" = "scaled.mass.tarsus.female")

lmod2 <- lmodel2(log(as.numeric(as.character(male.metrics$Weight)))~log(as.numeric(as.character(male.metrics$Tarsus))))
sma <- lmod2$regression.results[[3]][3]
mean.tarsus.male <- mean(as.numeric(as.character(male.metrics$Tarsus)),na.rm=T)
scaled.mass.tarsus.male <- (as.numeric(as.character(male.metrics$Weight))) * (mean.tarsus.male/as.numeric(as.character(male.metrics$Tarsus)))^(as.numeric(as.character(sma)))

scaled.mass.tarsus.male.df <- cbind(male.metrics, scaled.mass.tarsus.male) %>% dplyr::rename("scaled.mass.tarsus" = "scaled.mass.tarsus.male")
scaled.mass.tarsus.df <- rbind(scaled.mass.tarsus.female.df, scaled.mass.tarsus.male.df) %>% dplyr::select("BIOM.ID","scaled.mass.tarsus")


######################################################## 
## Scaled mass for nestling birds using tarsus as L variable
lmod2 <- lmodel2(log(as.numeric(as.character(nestlings.meta.sbm$Weight)))~log(as.numeric(as.character(nestlings.meta.sbm$Tarsus))))
sma <- lmod2$regression.results[[3]][3]
mean.tarsus <- mean(as.numeric(as.character(nestlings.meta.sbm$Tarsus)),na.rm=T)
scaled.mass.chick <- (as.numeric(as.character(nestlings.meta.sbm$Weight))) * (mean.tarsus/as.numeric(as.character(nestlings.meta.sbm$Tarsus)))^(as.numeric(as.character(sma)))

#scaled.mass.chick <- as.data.frame(scaled.mass.chick, row.names = as.character(nestlings.meta.sbm$BIOM.ID))
scaled.mass.chick <- cbind(as.data.frame(scaled.mass.chick),nestlings.meta.sbm$BIOM.ID) %>% dplyr::rename(BIOM.ID = "nestlings.meta.sbm$BIOM.ID")

#################################  
# need to source and bind in Data_setup
scaled.mass.all <- merge(scaled.mass, scaled.mass.tarsus.adult, by="BIOM.ID", all.x=T)
scaled.mass.all <- merge(scaled.mass.all, scaled.mass.wing.adult, by="BIOM.ID", all.x=T)
scaled.mass.all <- merge(scaled.mass.all, scaled.mass.chick, by="BIOM.ID", all.x=T)
scaled.mass.all <- merge(scaled.mass.all, scaled.mass.tarsus.df, by="BIOM.ID", all.x=T)
scaled.mass.all <- merge(scaled.mass.all, scaled.mass.wing.df, by="BIOM.ID", all.x=T)


rm(list=c("adults.meta.sbm", "lmod2", "mean.tarsus", "mean.tarsus.adult", "mean.wing.adult", "metadata.sbm", 
          "nestlings.meta.sbm", "scaled.mass", "scaled.mass.wing.adult", "scaled.mass.tarsus.adult", "scaled.mass.chick", "scaled.mass.wing.df","scaled.mass.tarsus.df",
          "scaled.mass.tarsus.male","mean.tarsus.female","sma","female.metrics", "male.metrics", "mean.tarsus.male", "mean.wing.female", 
          "mean.wing.male", "scaled.mass.tarsus.female","scaled.mass.tarsus.female.df", "scaled.mass.tarsus.male.df", 
          "scaled.mass.wing.female", "scaled.mass.wing.female.df", "scaled.mass.wing.male","scaled.mass.wing.male.df"))
     