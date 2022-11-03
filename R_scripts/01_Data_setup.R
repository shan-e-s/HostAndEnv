# Great tit microbiome Data set up script, 4/4/2019, Shane Somers
## script that can be sourced to read in and create neccessary variables
## alternative to this is just RData object containing phyloseq object
### UPDATE: 30/8/2019 DADA2 data used

# Set-up

rm(list=ls()) ## run fresh session instead of rm-list
#set.seed(10)

#install.packages("car")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("metagenomeSeq")

# libraries
library(tidyverse)
library(phyloseq) 
library(metagenomeSeq)
library(microbiome)

# NOTE: in rmd wd is automatically set to directory of rmd file, but for script it will be in project directory
#getwd()
load("Data/dada_outputCOMPLETE.rda")
#biom.spring <- import_biom("Data/Bird_microbiome_Spring/table_tax.biom")

# read in tree (phylogeny)
## note this tree created by SS using FF's otus.fa file and possibly using different alignment and filtering options so might not be valid
#tree.spring <- read_tree("Data/Bird_microbiome_Spring/otu.tre")
#biom.spring <- merge_phyloseq(biom.spring, tree.spring) ## losing 143 taxa, this occurs at qiime align_seqs.py and are listed in pynast/otus_failures.fasta


#Why does tree drop samples? These are dropped at alignment step and listed in pynast/otus_failures.fasta
#taxa.tre <- taxa_names(tree.spring)
#taxa.biom <- taxa_names(biom.spring)
#setdiff(taxa.biom,taxa.tre)

#Read in metadata. Note 3 samples dropped from phylo.object when sample data merged with phylo.object because of lack of corresponding sample data.
metadata <- read.csv2("Data/MasterBirdIDSheet.csv", sep = ",")

#filter out blanks ie. failed samples (not in BIOM file)
metadata <- metadata[!is.na(metadata$BIOM.ID),] ## why does metadata have 2 less samples than BIOM file

# swap "." for "_" in sample names to match dada output
#metadata$BIOM.ID <- gsub("\\.","\\_",metadata$BIOM.ID)
sample_names(ps) <- gsub("\\_S","\\.S",sample_names(ps))

# set rownames to match biom file ie. sample names
rownames(metadata) <- metadata$BIOM.ID

#setdiff(metadata$BIOM.ID,sample_names(biom.spring))
#setdiff(sample_names(biom.spring),metadata$BIOM.ID) # "CB5-D-15-38.S38"  "CB60-B-15-51.S56" "KB905F-D15.S9" ## not in metadata for some reason

# filter out failures and controls 
controls <- c("DNAnegcontrol-R.S59","G-Neg.S19","PCR-6D-28910.S93","PCRnegcon2809186E.S60")
metadata <- metadata[metadata$pass.or.fail!="control",]

# FILTERING
#Note: removing chloroplasts and mitochondria BEFORE dropping low read samples otherwise ends up w/ low read samples in phylo.object
#Note: subset_taxa also removes NA's ie. "The OTUs retained in the dataset is equivalent to x[subset & !is.na(subset)]" unless is.na() given.
# remove chloroplasts
ps <- subset_taxa(ps, Order!="Chloroplast" | is.na(Order)) # 55478/56852 taxa remain
# remove mitochondria
ps <- subset_taxa(ps,Family!="Mitochondria"| is.na(Family)) # 55310/55478
# remove Eukaryota
ps <- subset_taxa(ps,Kingdom!="Eukaryota"| is.na(Kingdom)) # 55181/55310
# remove Archaea
ps <- subset_taxa(ps,Kingdom!="Archaea"| is.na(Kingdom)) # 55073/55181
# remove singletons?
#ps1 <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE) # prevalence filter, drop taxa which occur in one sample
#ps <- prune_taxa(taxa_sums(ps)>1, ps) #abundance filter, drop taxa that only occur once

# remove ASV's not assigned to genus level?
#ps.noNAgenus <- subset_taxa(ps, !is.na(Genus)) # 30115/55310
#write.csv(ps.noNAgenus@tax_table, "Data/tax-table-noNA.csv") # quick method of viewing tax table

# remove failures and controls only, keep sub 40k
#metadata <- metadata[metadata$pass.or.fail!="control",]
#metadata <- metadata[metadata$pass.or.fail!="failedAtSequencing",]

# add sequence read numbers to samples
seq_reads <- sample_sums(ps)
metadata <- merge(metadata, seq_reads, by = "row.names") %>% dplyr::rename("NumberOfReads"="y")
rownames(metadata) <- metadata$BIOM.ID
metadata <- metadata[,-1] # drop "Row.names" column introduced by merge

# remove low read samples <40,000
## maybe try 5000 or 10000, as there is natural separation there and it isnt so strict
metadata <- metadata[metadata$NumberOfReads>10000,] # drops 258-247=11 samples

# which samples being dropped? Due to not having metadata records?
#setdiff(sample_names(biom.spring),sample_names(phylo.spring))
#setdiff(sample_names(phylo.spring),sample_names(biom.spring))

# remove unneccessary columns
# dput(colnames(metadata))
metadata <- metadata[,c("NumberID", "date", "BIOM.ID", "SequenceID.dilutionSheet", 
  "nest", "habitat", "IDLetter", "IDRing", "ageDays", "ageBinned", 
  "duplicate", "dateDNAextracted", "extractedBy", "notes", "QubitDNA", 
  "SequencePlate", "pass.or.fail", "mislabelled", "Tarsus", "Weight", "wing", "broodSizeWhenSampled", 
  "broodSizeMax", "Survive", "surviveRinged", "survivedFledged", 
  "totalFledge", "clutchSize", "layDateFirst", "numberDeadPreRinged", 
  "preRingedDeadID", "numberDeadPostRinged", "postRingedDeadID", "adult.survival", "NumberOfReads")]
      
# fix IDLetters eg. C-neat
# if IDLEtter longer than 1 symbol replace with first letter?
metadata["DW25-C-9-65-neat.S94","IDLetter"] <- "C"
metadata["KB905-E-B-177-neat.S96","IDLetter"] <- "E"
metadata["KB922-D-10-196-neat.S95","IDLetter"] <- "D"
# give controls IDLetter==control?
# set KB924C-D8.S75 BroodSizeWhenSampled to mean=4.3
metadata["KB924C-D8.S75","broodSizeWhenSampled"] <- 4
# assign 95 as layDateFirst as it is the estimated layDate for this nest according to MasterBreeding
for (i in 1:nrow(metadata)){
  if(metadata[i,"nest"]=="DW42") {metadata[i,"layDateFirst"] <- 95}
}

# GD not certain of ID of below birds as marking came off: exclude or impute?
exclude <- FALSE
if (exclude==F) {
  metadata["CB15-A-010716.S52","IDLetter"] <- "A"
  metadata["CB15-B-D15.S51","IDLetter"] <- "B"
  metadata["CB15-D-D15.S25","IDLetter"] <- "D"
} else {
  metadata <- metadata[-c("CB15-A-010716.S52","CB15-B-D15.S51","CB15-D-D15.S25"),]
}

# fix dates (incorrect year)
metadata["KB909-D-9-402.S82","date"] <- "10/06/2016"
metadata["CB59-B-unk.S30","date"] <- "01/07/2016"

# Fix BP7A-D15 metadata
# no record of BP7A-D15 being sampled or extracted on datasheets but it was ringed and biometrics taken. Sample date corresponds to 1 week after Day-9 sample.
# Possible that sample taken but not entered on spreadsheet but only in GD notebook. Will proceed as if sampled.
#metadata <- metadata[which(rownames(metadata)!="BP7-A-D15.S59"),] # remove?
metadata["BP7-A-D15.S59","wing"] <- 47
metadata["BP7-A-D15.S59","Tarsus"] <- 19.9
metadata["BP7-A-D15.S59","Weight"] <- 15.5
metadata["BP7-A-D15.S59","ageDays"] <- 15
metadata["BP7-A-D9.S17","IDRing"] <- "VZ80084" # fix Day-9 ring

# convert adult Females IDLetter from "F" to "Fe to differentiate between nestlings and females
# ifelse(metadata[,"IDLetter"]=="F" & metadata[,"ageBinned"]=="adult",print(metadata[,"BIOM.ID"]),print("nope"))
levels(metadata$IDLetter) <- c(levels(metadata$IDLetter), "Fe") ## need to add "Fe" as a factor level first

for (i in 1:nrow(metadata)){
  if(metadata[i,"IDLetter"]=="F" & metadata[i,"ageBinned"]=="adult") {metadata[i,"IDLetter"] <- "Fe"}
}

# create Sex variable of Male, Female and Unsexed(nestling)
# empty "" included for controls
# ifelse works w/ vectors > 1
Sex <- as.vector(metadata$IDLetter)
Sex <- ifelse(!Sex %in% c("M","Fe",""), "Unsexed", Sex)


# create bird.ID column of nest + IDLetter
bird.ID <- paste(metadata$nest, metadata$IDLetter, sep="")
metadata <- cbind(bird.ID, metadata, Sex)

# fix broodSizeWhenSampled NA's, which are actually known 
# c(LS2-C-15-215.S28, KB905-E-B-177-neat.S96, KB27-B-8-150.S58). KB905 should be 6,LS2C should be 3, KB27 is 4
metadata["LS2-C-15-215.S28","broodSizeWhenSampled"] <- 3
metadata["KB905-E-B-177-neat.S96","broodSizeWhenSampled"] <- 6
metadata["KB27-B-8-150.S58","broodSizeWhenSampled"] <- 4

##############
# rarefy
## need to join metadata and ps in order to drop very low samples, where (min(samplesums)==10220)
rarefy.y.n <- FALSE
if(rarefy.y.n==TRUE){
  phylo.metadata <- sample_data(metadata) # labelling as metadata
  rownames(phylo.metadata) <- phylo.metadata$BIOM.ID
  ps <- merge_phyloseq(phylo.metadata, ps) 
  ps <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),
                    rngseed = 1189, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
}
##############

############################## RECALCULATE ALPHA DIVERSITY
## So, estimate_richness changes the row.names AND orders them differently so they no longer match the BIOM and metadata files, microbiome::diversity doesnt do this but 
# does order them in the same way BUT doesnt have all the diversity indices options i need so im calculating with both functions, cbinding them and dropping the extra shannon diversity
# recalculate alpha diversity having droppped failures, low read samples and controls
a.div.phyloseq <- estimate_richness(merge_phyloseq(ps,sample_data(metadata)), split = T, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
## ALSO function converts all rownames to having "." instead of "-" so just cbind but ENSURE ORDER OF SAMPLES IS CORRECT
a.div.vegan <- microbiome::diversity(merge_phyloseq(ps,sample_data(metadata)), index=c("shannon"))
### FIXED: diversity() command was failing (cannot coerce S4 to vector) when other scripts have been run in same session 
#### was using diversity from microbiome pkg not vegan
# correct row.names and order
a.div.correct.rownames <- cbind(a.div.vegan,a.div.phyloseq)
colnames(a.div.correct.rownames) <- tolower(colnames(a.div.correct.rownames)) # lower colnames so dont have to change upstream scripts
meta.alpha <- merge(metadata, a.div.correct.rownames[,-1], by = "row.names")
rownames(meta.alpha) <- meta.alpha$Row.names # fix rownames
meta.alpha <- meta.alpha[,-1] # drop Row.names column
############################

#Need to create and add site column to metadata. Fix control samples and typo "1N".
# substring of 
site <- str_sub(meta.alpha$BIOM.ID,1,2)
site[1] <- "IN"

which.control <- which(site %in% c("DN","G-","PC"))
site[which.control] <- "control"

meta.alpha <- cbind(meta.alpha,site)

#Add distance to edge variable for each nest. Calculating in separate script and sourcing. WHAT UNIT IS DISTANCE MEASURING??
source("R_scripts/01a_DistanceToEdge.R")
## FIND OUT WHAT UNIT THE DISTANCE IS MEASURING

# need to merge on nest.ID
meta.alpha.site <- merge.data.frame(meta.alpha, nestbox.edge[,c("nest","DistanceToEdge")], by = "nest", all.x=TRUE)

# Add scaled body mass variable for all birds and then nestlings and adults separately
source("R_scripts/01b_ScaledBodyMass.R")
meta.alpha.site <- merge.data.frame(meta.alpha.site, scaled.mass.all, by = "BIOM.ID", all.x=TRUE)


#Need to convert numeric variables from factors.
#dput(colnames(phylo.spring@sam_data))
#dput(colnames(metadata))

numeric.columns <- c( "Tarsus","Weight", "wing", "broodSizeWhenSampled", "broodSizeMax",
                      "totalFledge", "clutchSize", "layDateFirst", "numberDeadPreRinged", "numberDeadPostRinged","QubitDNA",
                      "chao1", "se.chao1", "shannon", "invsimpson","observed","DistanceToEdge","scaled.mass","scaled.mass.wing.adult",
                      "scaled.mass.tarsus.adult", "scaled.mass.chick","scaled.mass.tarsus","scaled.mass.wing", "NumberOfReads")

#"ageDays", "N7index", "S5index", "QubitDNA", 

# set columns as numeric vectors
metadata.numeric <- meta.alpha.site %>% mutate_at( numeric.columns, funs(as.numeric(as.character(.))))

# convert date variable to date type
meta.alpha.site$date <- as.Date(meta.alpha.site$date,"%d/%m/%Y")

# meta doesnt maintain Date data-type, is this because of NA's,blanks,'?-marks'?
#unique(meta.alpha.site$date)
#str(meta.alpha.site)
#str(meta(phylo.spring))

# merge w/ phylo object
phylo.metadata <- sample_data(metadata.numeric) # labelling as metadata
rownames(phylo.metadata) <- phylo.metadata$BIOM.ID
phylo.spring <- merge_phyloseq(phylo.metadata, ps) 

# Prune troublesome samples
# These possibly causing error in lmer. CB60F listed as duplicate but i think its confusion due to F indicating female and f-chick. 
# DW42 is listed as failing at sequencing- why still included in biom?== because didnt quite fail but had low reads (8060 reads).
phylo.spring <- prune_samples(!rownames(phylo.spring@sam_data) %in% c("CB60F-D15ii.S59","DW42-D534684.S68"), phylo.spring)

# remove ASV's which now have zero abundance because of removed samples
phylo.spring <- prune_taxa(taxa_sums(phylo.spring) > 0, phylo.spring) # 51602/55310 remain

############ Clean up ## should be mostly unneccessary because im not sourcing but saving an object
#rm(list = c("alpha.div", "biom.spring", "meta.alpha", "meta.alpha.site","metadata", "nestbox.edge", "numeric.columns", 
# "phylo.metadata", "site", "tree.spring", "which.control"))

############ Save outputs
# save data as R object
saveRDS(phylo.spring, file = "Data/phylo-spring.rds")
#saveRDS(phylo.spring, file = "Data/phylo-spring-rarefied.rds")

