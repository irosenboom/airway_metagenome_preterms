# Author: Ilona Rosenboom
# Last updated: 13 April 2022

# set working directory
setwd("/mnt/sfb900nfs/groups/tuemmler/ilona/preterm_manuscript")

# clean global R environment
rm(list = ls())

# required R packages
library(readr)
library(readxl)
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggnewscale)
library(Hmisc)
library(dplyr)
library(ggVennDiagram)
library(rcompanion)


########################################################################################################
# load final bacteria per human cell table

bphc <- read_delim("bacteria_per_human_cell_haybaler.csv", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)

bphc <- data.frame(bphc, check.names = FALSE)

#replace "_" between species names with space " "
bphc$species <- gsub(x = bphc$species, pattern = "\\_", replacement = " ")

rownames(bphc) <- bphc$species
bphc$species <- NULL
bphc$chr_length <- NULL
bphc$gc_ref <- NULL

bphc_t <- data.frame(t(bphc))

#exclude 19_037d
bphc_t_all <- bphc_t

bphc_t <- bphc_t[-c(38),]

#preterm and healthy controls

bphc_preterm_healthy <- read_delim("bphc_preterm_healthy.csv", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
bphc_preterm_healthy <- data.frame(bphc_preterm_healthy, check.names = FALSE)
rownames(bphc_preterm_healthy) <- bphc_preterm_healthy$species
bphc_preterm_healthy$species <- NULL
bphc_preterm_healthy$gc_ref <- NULL
bphc_preterm_healthy$chr_lngth <- NULL
bphc_preterm_healthy$`19_037d` <- NULL

bphc_preterm_healthy_t <- data.frame(t(bphc_preterm_healthy))

########################################################################################################
# load metadata

metadata <- read_excel("metadata_preterm.xlsx")
metadata <- data.frame(metadata, check.names = FALSE)

rownames(metadata) <- metadata$sample

#exclude 19_037d
metadata_all <- metadata

metadata <- metadata[-c(38),]


metadata_preterm_healthy <- read_excel("metadata_preterm_healthycontrols.xlsx")
metadata_preterm_healthy <- data.frame(metadata_preterm_healthy, check.names = FALSE)

rownames(metadata_preterm_healthy) <- metadata_preterm_healthy$sample_id

metadata_preterm_healthy <- metadata_preterm_healthy[-c(49),]

########################################################################################################
# order tables
bphc <- bphc[order(rownames(bphc)),]
bphc_t <- bphc_t[order(rownames(bphc_t)),]

metadata <- metadata[order(rownames(metadata)),]

bphc_preterm_healthy <- bphc_preterm_healthy[order(rownames(bphc_preterm_healthy)),]
bphc_preterm_healthy_t <- bphc_preterm_healthy_t[order(rownames(bphc_preterm_healthy_t)),]

metadata_preterm_healthy <- metadata_preterm_healthy[order(rownames(metadata_preterm_healthy)),]
metadata_healthy <- subset(metadata_preterm_healthy, metadata_preterm_healthy$state == "Healthy")

########################################################################################################
# create different datasets for each time point of sampling, non-BPD and BPD

#make different data sets for time points
metadata_a <- subset(metadata, timepoint_sampling == 1)
list_t_a <- rownames(metadata_a)
bphc_t_a <- bphc_t[rownames(bphc_t) %in% list_t_a,]

metadata_b <- subset(metadata, timepoint_sampling == 2)
list_t_b <- rownames(metadata_b)
bphc_t_b <- bphc_t[rownames(bphc_t) %in% list_t_b,]

metadata_c <- subset(metadata, timepoint_sampling == 3)
list_t_c <- rownames(metadata_c)
bphc_t_c <- bphc_t[rownames(bphc_t) %in% list_t_c,]

metadata_d <- subset(metadata, timepoint_sampling == 4)
list_t_d <- rownames(metadata_d)
bphc_t_d <- bphc_t[rownames(bphc_t) %in% list_t_d,]

metadata_c_d <- rbind(metadata_c, metadata_d)

#also divide different data sets into control and BPD
metadata_a_control <- subset(metadata_a, group == "control")
metadata_a_BPD <- subset(metadata_a, group == "BPD")
list_t_a_control <- rownames(metadata_a_control)
list_t_a_BPD <- rownames(metadata_a_BPD)
bphc_t_a_control <- bphc_t_a[rownames(bphc_t_a) %in% list_t_a_control,]
bphc_t_a_BPD <- bphc_t_a[rownames(bphc_t_a) %in% list_t_a_BPD,]

metadata_b_control <- subset(metadata_b, group == "control")
metadata_b_BPD <- subset(metadata_b, group == "BPD")
list_t_b_control <- rownames(metadata_b_control)
list_t_b_BPD <- rownames(metadata_b_BPD)
bphc_t_b_control <- bphc_t_b[rownames(bphc_t_b) %in% list_t_b_control,]
bphc_t_b_BPD <- bphc_t_b[rownames(bphc_t_b) %in% list_t_b_BPD,]

metadata_c_control <- subset(metadata_c, group == "control")
metadata_c_BPD <- subset(metadata_c, group == "BPD")
list_t_c_control <- rownames(metadata_c_control)
list_t_c_BPD <- rownames(metadata_c_BPD)
bphc_t_c_control <- bphc_t_c[rownames(bphc_t_c) %in% list_t_c_control,]
bphc_t_c_BPD <- bphc_t_c[rownames(bphc_t_c) %in% list_t_c_BPD,]

metadata_d_control <- subset(metadata_d, group == "control")
metadata_d_BPD <- subset(metadata_d, group == "BPD")
list_t_d_control <- rownames(metadata_d_control)
list_t_d_BPD <- rownames(metadata_d_BPD)
bphc_t_d_control <- bphc_t_d[rownames(bphc_t_d) %in% list_t_d_control,]
bphc_t_d_BPD <- bphc_t_d[rownames(bphc_t_d) %in% list_t_d_BPD,]

#######################################################################################################
#calculate core and rare species for time point c and d, threshold 95% vs 5%

#time point c
##control
###core species
bphc_c_control <- data.frame(t(bphc_t_c_control), check.names = FALSE)
bphc_c_control$rowsum <- rowSums(bphc_c_control)
sum_total_c_control= sum(bphc_c_control$rowsum)
bphc_c_control$abundance <- (bphc_c_control$rowsum / sum_total_c_control) * 100
# sort abundance column (decreasing)
bphc_c_control <- bphc_c_control[order(bphc_c_control$abundance, decreasing = TRUE),]  
# sum abundance values, make cutoff where cumsum <95%
bphc_c_control$cumsum <- cumsum(bphc_c_control$abundance)
bphc_c_control <- subset(bphc_c_control, abundance !=0)
core_c_control <- rownames(bphc_c_control[1:17,])
bphc_c_core_control <- bphc_c_control[rownames(bphc_c_control) %in% core_c_control ,]
species_c_core_control <- data.frame(t(bphc_c_core_control), check.names = FALSE)
species_c_core_control <- species_c_core_control[-c(13,14,15),]

###rare_species
rare_c_control <- bphc_c_control[!rownames(bphc_c_control) %in% core_c_control ,]
rare_c_control$rowsum <- NULL
rare_c_control$abundance <- NULL
rare_c_control$cumsum <- NULL
names_rare_c_control <- rownames(rare_c_control)
bphc_c_rare_control <- bphc_c_control[rownames(bphc_c_control) %in% names_rare_c_control ,]
species_c_rare_control <- data.frame(t(bphc_c_rare_control), check.names = FALSE)
species_c_rare_control <- species_c_rare_control[-c(13,14,15),]

##BPD
###core species
bphc_c_BPD <- data.frame(t(bphc_t_c_BPD), check.names = FALSE)
bphc_c_BPD$rowsum <- rowSums(bphc_c_BPD)
sum_total_c_BPD= sum(bphc_c_BPD$rowsum)
bphc_c_BPD$abundance <- (bphc_c_BPD$rowsum / sum_total_c_BPD) * 100
# sort abundance column (decreasing)
bphc_c_BPD <- bphc_c_BPD[order(bphc_c_BPD$abundance, decreasing = TRUE),]  
# sum abundance values, make cutoff where cumsum <95%
bphc_c_BPD$cumsum <- cumsum(bphc_c_BPD$abundance)
bphc_c_BPD <- subset(bphc_c_BPD, abundance !=0)
core_c_BPD <- rownames(bphc_c_BPD[1:18,])
bphc_c_core_BPD <- bphc_c_BPD[rownames(bphc_c_BPD) %in% core_c_BPD ,]
species_c_core_BPD <- data.frame(t(bphc_c_core_BPD), check.names = FALSE)
species_c_core_BPD <- species_c_core_BPD[-c(12,13,14),]

###rare species
rare_c_BPD <- bphc_c_BPD[!rownames(bphc_c_BPD) %in% core_c_BPD ,]
rare_c_BPD$rowsum <- NULL
rare_c_BPD$abundance <- NULL
rare_c_BPD$cumsum <- NULL
names_rare_c_BPD <- rownames(rare_c_BPD)
bphc_c_rare_BPD <- bphc_c_BPD[rownames(bphc_c_BPD) %in% names_rare_c_BPD ,]
species_c_rare_BPD <- data.frame(t(bphc_c_rare_BPD), check.names = FALSE)
species_c_rare_BPD <- species_c_rare_BPD[-c(12,13,14),]

#time point d
##control
###core species
bphc_d_control <- data.frame(t(bphc_t_d_control), check.names = FALSE)
bphc_d_control$rowsum <- rowSums(bphc_d_control)
sum_total_d_control= sum(bphc_d_control$rowsum)
bphc_d_control$abundance <- (bphc_d_control$rowsum / sum_total_d_control) * 100
# sort abundance column (decreasing)
bphc_d_control <- bphc_d_control[order(bphc_d_control$abundance, decreasing = TRUE),]  
# sum abundance values, make cutoff where cumsum <95%
bphc_d_control$cumsum <- cumsum(bphc_d_control$abundance)
bphc_d_control <- subset(bphc_d_control, abundance !=0)
core_d_control <- rownames(bphc_d_control[1:27,])
bphc_d_core_control <- bphc_d_control[rownames(bphc_d_control) %in% core_d_control ,]
species_d_core_control <- data.frame(t(bphc_d_core_control), check.names = FALSE)
species_d_core_control <- species_d_core_control[-c(9,10,11),]

####rare_species
rare_d_control <- bphc_d_control[!rownames(bphc_d_control) %in% core_d_control ,]
rare_d_control$rowsum <- NULL
rare_d_control$abundance <- NULL
rare_d_control$cumsum <- NULL
names_rare_d_control <- rownames(rare_d_control)
bphc_d_rare_control <- bphc_d_control[rownames(bphc_d_control) %in% names_rare_d_control ,]
species_d_rare_control <- data.frame(t(bphc_d_rare_control), check.names = FALSE)
species_d_rare_control <- species_d_rare_control[-c(9,10,11),]

##BPD
###core species
bphc_d_BPD <- data.frame(t(bphc_t_d_BPD), check.names = FALSE)
bphc_d_BPD$rowsum <- rowSums(bphc_d_BPD)
sum_total_d_BPD= sum(bphc_d_BPD$rowsum)
bphc_d_BPD$abundance <- (bphc_d_BPD$rowsum / sum_total_d_BPD) * 100
# sort abundance column (decreasing)
bphc_d_BPD <- bphc_d_BPD[order(bphc_d_BPD$abundance, decreasing = TRUE),]  
# sum abundance values, make cutoff where cumsum <95%
bphc_d_BPD$cumsum <- cumsum(bphc_d_BPD$abundance)
bphc_d_BPD <- subset(bphc_d_BPD, abundance !=0)
core_d_BPD <- rownames(bphc_d_BPD[1:24,])
bphc_d_core_BPD <- bphc_d_BPD[rownames(bphc_d_BPD) %in% core_d_BPD ,]
species_d_core_BPD <- data.frame(t(bphc_d_core_BPD), check.names = FALSE)
species_d_core_BPD <- species_d_core_BPD[-c(9,10,11),]

###rare_species
rare_d_BPD <- bphc_d_BPD[!rownames(bphc_d_BPD) %in% core_d_BPD ,]
rare_d_BPD$rowsum <- NULL
rare_d_BPD$abundance <- NULL
rare_d_BPD$cumsum <- NULL
names_rare_d_BPD <- rownames(rare_d_BPD)
bphc_d_rare_BPD <- bphc_d_BPD[rownames(bphc_d_BPD) %in% names_rare_d_BPD ,]
species_d_rare_BPD <- data.frame(t(bphc_d_rare_BPD), check.names = FALSE)
species_d_rare_BPD <- species_d_rare_BPD[-c(9,10,11),]

#######################################################################################################
#new datasets to work with

#time point sampling a
species_a <- subset(bphc_t_a[,colSums(bphc_t_a) !=0])
##filter, make rowSum column and decreasing
species_a_t <- data.frame(t(species_a), check.names = FALSE)
species_a_t$rowSum <- rowSums(species_a_t)
species_a_t <- species_a_t[order(species_a_t$rowSum, decreasing = TRUE),]
species_a_t$rowSum <- NULL
species_a <- data.frame(t(species_a_t), check.names = FALSE)

#time point sampling b
species_b <- subset(bphc_t_b[,colSums(bphc_t_b) !=0])
##filter, make rowSum column and decreasing
species_b_t <- data.frame(t(species_b), check.names = FALSE)
species_b_t$rowSum <- rowSums(species_b_t)
species_b_t <- species_b_t[order(species_b_t$rowSum, decreasing = TRUE),]
species_b_t$rowSum <- NULL
species_b <- data.frame(t(species_b_t), check.names = FALSE)

#time point sampling c: core_rare_species_c
#species_c <- core_rare_species_c
species_c_t <- data.frame(t(bphc_t_c), check.names = FALSE)
##filter, make rowSum column and decreasing
species_c_t$rowSum <- rowSums(species_c_t)
species_c_t <- species_c_t[order(species_c_t$rowSum, decreasing = TRUE),]
species_c_t$rowSum <- NULL
species_c <- data.frame(t(species_c_t), check.names = FALSE)

#time point sampling d: core_rare_species_d
species_d_t <- data.frame(t(bphc_t_d), check.names = FALSE)
##filter, make rowSum column and decreasing
species_d_t$rowSum <- rowSums(species_d_t)
species_d_t <- species_d_t[order(species_d_t$rowSum, decreasing = TRUE),]
species_d_t$rowSum <- NULL
species_d <- data.frame(t(species_d_t), check.names = FALSE)

#time point sampling c and d combined, no filtering for core and/or rare species
species_c_d <- bphc_t
species_c_d <- rbind(bphc_t[rownames(bphc_t) %in% list_t_c,], bphc_t[rownames(bphc_t) %in% list_t_d,])
##filter, make rowSum column and decreasing
species_c_d_t <- data.frame(t(species_c_d), check.names = FALSE)
species_c_d_t$rowSum <- rowSums(species_c_d_t)
species_c_d_t <- species_c_d_t[order(species_c_d_t$rowSum, decreasing = TRUE),]
species_c_d_t <- subset(species_c_d_t, rowSum!=0)
species_c_d_t$rowSum <- NULL
species_c_d <- data.frame(t(species_c_d_t), check.names = FALSE)
species_c_d <- species_c_d[order(rownames(species_c_d)),]

#healthy controls
species_healthycontrols <- bphc_preterm_healthy_t
species_healthycontrols$state <- metadata_preterm_healthy$state

species_healthycontrols <- subset(species_healthycontrols, state=="Healthy")
species_healthycontrols$state <- NULL


####################################################################################################
#manipulate metadata to obtain non-corrected age
#metadata$gestational_age_days = real age

metadata_c_d <- metadata_c_d[order(rownames(metadata_c_d)),]

metadata_preterm <- subset(metadata_preterm_healthy, state == "Preterm")
metadata_preterm$age_sampling <- metadata_c_d$age_sampling

metadata_healthy$age_sampling <- metadata_healthy$age_in_days

metadata_preterm_healthy <- rbind(metadata_preterm, metadata_healthy)


#######################################################################################################
#######################################################################################################
# calculate mean abundances for heat trees

species_a_delivery <- species_a
is.na(species_a_delivery) <- species_a_delivery==0
species_a_delivery$delivery <- metadata_a$delivery_details

species_a_electiveCS <- subset(species_a_delivery, delivery==1)
species_a_electiveCS$delivery <- NULL

species_a_VD <- subset(species_a_delivery, delivery==4)
species_a_VD$delivery <- NULL


##use mean for a, b, c, d and healthy
#species a
species_a_electiveCS_mean <- mapply(mean, species_a_electiveCS, na.rm = TRUE)
species_a_electiveCS_mean <- data.frame(species_a_electiveCS_mean)

species_a_VD_mean <- mapply(mean, species_a_VD, na.rm = TRUE)
species_a_VD_mean <- data.frame(species_a_VD_mean)

#species b
species_b_for_mean <- species_b

species_b_mean <- mapply(mean, species_b_for_mean, na.rm = TRUE)
species_b_mean <- data.frame(species_b_mean)

#species c
species_c_for_mean <- species_c

species_c_mean <- mapply(mean, species_c_for_mean, na.rm = TRUE)
species_c_mean <- data.frame(species_c_mean)

#species d
species_d_for_mean <- species_d

species_d_mean <- mapply(mean, species_d_for_mean, na.rm = TRUE)
species_d_mean <- data.frame(species_d_mean)

#species healthy controls
species_healthy_for_mean <- species_healthycontrols

species_healthy_mean <- mapply(mean, species_healthy_for_mean, na.rm = TRUE)
species_healthy_mean <- data.frame(species_healthy_mean)

##merge mean
species_mean_a <- merge(species_a_electiveCS_mean, species_a_VD_mean,
                        by = 'row.names', all = TRUE)
rownames(species_mean_a) <- species_mean_a$Row.names
species_mean_a$Row.names <- NULL

species_mean_a_b <- merge(species_mean_a, species_b_mean,
                          by = 'row.names', all = TRUE)
rownames(species_mean_a_b) <- species_mean_a_b$Row.names
species_mean_a_b$Row.names <- NULL

species_mean_all <- merge(species_c_mean, species_d_mean,
                          by = 'row.names', all = TRUE)
rownames(species_mean_all) <- species_mean_all$Row.names
species_mean_all$Row.names <- NULL

species_mean_all <- merge(species_mean_a_b, species_mean_all,
                          by = 'row.names', all = TRUE)
rownames(species_mean_all) <- species_mean_all$Row.names
species_mean_all$Row.names <- NULL

species_mean_all <- merge(species_mean_all, species_healthy_mean,
                          by = 'row.names', all = TRUE)
rownames(species_mean_all) <- species_mean_all$Row.names
species_mean_all$Row.names <- NULL

species_mean_all[is.na(species_mean_all)] <- 0
species_mean_all$rowsum <- rowSums(species_mean_all)
species_mean_all <- subset(species_mean_all, rowsum != 0)

species_mean_all <- species_mean_all[order(species_mean_all$rowsum, decreasing = TRUE),]
species_mean_all$rowsum <- NULL

write.table(species_mean_all,"species_mean_all.csv",
            row.names = TRUE, sep = "\t")

#######################################################################################################
#######################################################################################################
#######################################################################################################

#calculate alpha diversity
#order tables
metadata <- metadata[order(rownames(metadata)),]
bphc_t <- bphc_t[order(rownames(bphc_t)),]

#calculate alpha diversity
metadata$shannon_index <- diversity(bphc_t, index = "shannon")
metadata$specnumber <- specnumber(bphc_t)

metadata_healthy$shannon_index <- diversity(species_healthycontrols, index = "shannon")
metadata_healthy$specnumber <- specnumber(species_healthycontrols)

#time point sampling c

metadata_c$shannon_index <- diversity(bphc_t_c, index = "shannon")
metadata_c$simpson_index <- diversity(bphc_t_c, index = "simpson")

##control
metadata_c_control$shannon_core <- diversity(species_c_core_control, index = "shannon")
metadata_c_control$simpson_core <- diversity(species_c_core_control, index = "simpson")
metadata_c_control$specNumber_core <- specnumber(species_c_core_control)
metadata_c_control$richness_core <- metadata_c_control$specNumber_core / sqrt(length(metadata_c_control$specNumber_core))
metadata_c_control$bacterialBurden_core <- round(rowSums(species_c_core_control),0)

metadata_c_control$shannon_rare <- diversity(species_c_rare_control, index = "shannon")
metadata_c_control$simpson_rare <- diversity(species_c_rare_control, index = "simpson")
metadata_c_control$specNumber_rare <- specnumber(species_c_rare_control)
metadata_c_control$richness_rare <- metadata_c_control$specNumber_rare / sqrt(length(metadata_c$specNumber_rare))
metadata_c_control$bacterialBurden_rare <- round(rowSums(species_c_rare_control),0)

##BPD
metadata_c_BPD$shannon_core <- diversity(species_c_core_BPD, index = "shannon")
metadata_c_BPD$simpson_core <- diversity(species_c_core_BPD, index = "simpson")
metadata_c_BPD$specNumber_core <- specnumber(species_c_core_BPD)
metadata_c_BPD$richness_core <- metadata_c_BPD$specNumber_core / sqrt(length(metadata_c_BPD$specNumber_core))
metadata_c_BPD$bacterialBurden_core <- round(rowSums(species_c_core_BPD),0)

metadata_c_BPD$shannon_rare <- diversity(species_c_rare_BPD, index = "shannon")
metadata_c_BPD$simpson_rare <- diversity(species_c_rare_BPD, index = "simpson")
metadata_c_BPD$specNumber_rare <- specnumber(species_c_rare_BPD)
metadata_c_BPD$richness_rare <- metadata_c_BPD$specNumber_rare / sqrt(length(metadata_c$specNumber_rare))
metadata_c_BPD$bacterialBurden_rare <- round(rowSums(species_c_rare_BPD),0)

##combine metadata_c_control and _BPD
metadata_c_merge <- rbind(metadata_c_control, metadata_c_BPD)
metadata_c_merge <- metadata_c_merge[order(rownames(metadata_c_merge)),]
metadata_c_merge$shannon_index <- metadata_c$shannon_index
metadata_c_merge$simpson_index <- metadata_c$simpson_index

#time point sampling d

metadata_d$shannon_index <- diversity(bphc_t_d, index = "shannon")
metadata_d$simpson_index <- diversity(bphc_t_d, index = "simpson")

##control
metadata_d_control$shannon_core <- diversity(species_d_core_control, index = "shannon")
metadata_d_control$simpson_core <- diversity(species_d_core_control, index = "simpson")
metadata_d_control$specNumber_core <- specnumber(species_d_core_control)
metadata_d_control$richness_core <- metadata_d_control$specNumber_core / sqrt(length(metadata_d_control$specNumber_core))
metadata_d_control$bacterialBurden_core <- round(rowSums(species_d_core_control),0)

metadata_d_control$shannon_rare <- diversity(species_d_rare_control, index = "shannon")
metadata_d_control$simpson_rare <- diversity(species_d_rare_control, index = "simpson")
metadata_d_control$specNumber_rare <- specnumber(species_d_rare_control)
metadata_d_control$richness_rare <- metadata_d_control$specNumber_rare / sqrt(length(metadata_c$specNumber_rare))
metadata_d_control$bacterialBurden_rare <- round(rowSums(species_d_rare_control),0)

##BPD
metadata_d_BPD$shannon_core <- diversity(species_d_core_BPD, index = "shannon")
metadata_d_BPD$simpson_core <- diversity(species_d_core_BPD, index = "simpson")
metadata_d_BPD$specNumber_core <- specnumber(species_d_core_BPD)
metadata_d_BPD$richness_core <- metadata_d_BPD$specNumber_core / sqrt(length(metadata_d_BPD$specNumber_core))
metadata_d_BPD$bacterialBurden_core <- round(rowSums(species_d_core_BPD),0)

metadata_d_BPD$shannon_rare <- diversity(species_d_rare_BPD, index = "shannon")
metadata_d_BPD$simpson_rare <- diversity(species_d_rare_BPD, index = "simpson")
metadata_d_BPD$specNumber_rare <- specnumber(species_d_rare_BPD)
metadata_d_BPD$richness_rare <- metadata_d_BPD$specNumber_rare / sqrt(length(metadata_c$specNumber_rare))
metadata_d_BPD$bacterialBurden_rare <- round(rowSums(species_d_rare_BPD),0)

##combine metadata_d_control and _BPD
metadata_d_merge <- rbind(metadata_d_control, metadata_d_BPD)
metadata_d_merge <- metadata_d_merge[order(rownames(metadata_d_merge)),]
metadata_d_merge$shannon_index <- metadata_d$shannon_index
metadata_d_merge$simpson_index <- metadata_d$simpson_index


##combine metadata_c and _d
metadata_c_d <- rbind(metadata_c_merge, metadata_d_merge)
metadata_c_d <- metadata_c_d[order(rownames(metadata_c_d)),]


#calculation for preterm and healthy combined
metadata_preterm_healthy <- metadata_preterm_healthy[order(rownames(metadata_preterm_healthy)),]
bphc_preterm_healthy_t <- bphc_preterm_healthy_t[order(rownames(bphc_preterm_healthy_t)),]

metadata_preterm_healthy$shannon_index <- diversity(bphc_preterm_healthy_t, index = "shannon")
metadata_preterm_healthy$specnumber <- specnumber(bphc_preterm_healthy_t)

#######################################################################################################
#plot and test alpha diversity

#all time points
##ggline plot with error bars

#all time points
##ggline plot with error bars

metadata$delivery_details <- factor(metadata$delivery_details)


metadata$timepoint_sampling <- as.factor(metadata$timepoint_sampling)
metadata$timepoint_median <- with(metadata,
                                  ifelse(metadata$timepoint_sampling == 1, "4",
                                         ifelse(metadata$timepoint_sampling == 2, "33",
                                                ifelse(metadata$timepoint_sampling == 3, "272",
                                                       ifelse(metadata$timepoint_sampling == 4, "449", NA)))))

metadata$timepoint_median <- as.numeric(metadata$timepoint_median)

preterm <- 
  ggline(metadata, x = "timepoint_median", y = "shannon_index",
         color = "delivery_details", group = "delivery_details", size = 0.7,
         numeric.x.axis = TRUE,
         add = c("mean_se"), position = position_dodge(0.1), 
         add.params = list(alpha = 0.5)) +
  stat_compare_means(aes(group = delivery_details), label = "p.signif",
                     label.y = c(3), hide.ns = TRUE, size = 6) +
  scale_color_discrete(breaks = c(1,2,4),
                       labels = c("Elective CS", "Secondary CS", "Vaginal delivery")) +
  scale_y_continuous(limits = c(0,3.1), breaks = c(1,2,3)) +
  scale_x_continuous(limits = c(0, 510),
                     breaks = c(0,100,200,300,400,500)) +
  labs(y = "Shannon diversity", x = "Age (in days)", 
       color = "Mode of delivery") +
  geom_point(aes(x=age_sampling, y=shannon_index, color=delivery_details),
             alpha=0.5, size = 2.5) 

preterm

metadata %>%
  group_by(timepoint_sampling) %>%
  kruskal_test(shannon_index ~ delivery_details)


metadata_first_t <- metadata[metadata$timepoint_sampling==1,]
epsilonSquared(x = metadata_first_t$shannon_index, 
               g = metadata_first_t$delivery_details, ci=TRUE)


ggline <- preterm +
  new_scale_color() +
  geom_point(data=metadata_healthy, aes(x = age_sampling, y = shannon_index, color = state,
                                        shape = state), alpha = 0.8, size = 2.2) +
  scale_color_manual(values = c("Healthy" = "grey60"),
                     labels = c("Healthy full-term"),
                     name = "State") +
  scale_shape_manual(values = c("Healthy" = 15),
                     labels = c("Healthy full-term"),
                     name = "State")

ggline

compare_means(shannon_index ~ timepoint_sampling,  data = metadata, p.adjust.method = "BH")


##time point c and d - antimicrobial therapy in neonatal period
#time point c

core_272_antibiotics <-
  ggerrorplot(metadata_c_merge, 
              x = "neonatal_antibiotics", y = "shannon_core",
              main = "", xlab = "", 
              ylab = "\nShannon diversity",
              desc_stat = "median", 
              ggtheme = theme_pubr(base_size = 14, 
                                   base_family = "Helvetica",
                                   border = TRUE),
              add = c("boxplot", "point"), add.params = list(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Untreated", "Treated")) +
  labs(title = "", subtitle = "Core species (m9)", caption = "") + 
  stat_compare_means(method = "wilcox.test", label.y = 3.5, size = 4.5, label = "p.signif",
                     family = "Helvetica", label.x.npc = "left") 

#calculate effect size
wilcox_effsize(metadata_c_merge, shannon_core ~ neonatal_antibiotics,
               ci = TRUE)
wilcox_test(metadata_c_merge, shannon_core ~ neonatal_antibiotics)

rare_272_antibiotics <-
  ggerrorplot(metadata_c_merge, 
              x = "neonatal_antibiotics", y = "shannon_rare",
              main = "", xlab = "", 
              ylab = "\nShannon diversity",
              desc_stat = "median", 
              ggtheme = theme_pubr(base_size = 14, 
                                   base_family = "Helvetica",
                                   border = TRUE),
              add = c("boxplot", "point"), add.params = list(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Untreated", "Treated")) +
  labs(title = "", subtitle = "Rare species (m9)", caption = "") +
  stat_compare_means(method = "wilcox.test", label.y = 3.5, size = 4.5, label = "p.signif",
                     family = "Helvetica", label.x.npc = "left")

#calculate effect size
wilcox_effsize(metadata_c_merge, shannon_rare ~ neonatal_antibiotics,
               ci = TRUE)
wilcox_test(metadata_c_merge, shannon_rare ~ neonatal_antibiotics)

#shannon diversity core and rare combined
ggerrorplot(metadata_c_merge, 
            x = "neonatal_antibiotics", y = "shannon_index",
            main = "", xlab = "", 
            ylab = "\nShannon diversity",
            desc_stat = "median", 
            ggtheme = theme_pubr(base_size = 14, 
                                 base_family = "Helvetica",
                                 border = TRUE),
            add = c("boxplot", "point"), add.params = list(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Untreated", "Treated")) +
  stat_compare_means(method = "wilcox.test", label.y = 3.5, size = 4.5, label = "p.signif",
                     family = "Helvetica", label.x.npc = "left")

core_rare_272_antibiotics <- 
  ggarrange(core_272_antibiotics, rare_272_antibiotics)
core_rare_272_antibiotics
                              

#time point d

core_450_antibiotics <-
  ggerrorplot(metadata_d_merge, 
              x = "neonatal_antibiotics", y = "shannon_core",
              main = "", xlab = "", 
              ylab = "\nShannon diversity",
              desc_stat = "median", 
              ggtheme = theme_pubr(base_size = 14, 
                                   base_family = "Helvetica",
                                   border = TRUE),
              add = c("boxplot", "point"), add.params = list(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Untreated", "Treated")) +
  labs(title = "", subtitle = "Core species (m15)", caption = "") + 
  stat_compare_means(method = "wilcox.test", label.y = 3.5, size = 4.5, label = "p.signif",
                     family = "Helvetica", label.x.npc = "left", hide.ns = TRUE) 

rare_450_antibiotics <-
  ggerrorplot(metadata_d_merge, 
              x = "neonatal_antibiotics", y = "shannon_rare",
              main = "", xlab = "", 
              ylab = "\nShannon diversity",
              desc_stat = "median", 
              ggtheme = theme_pubr(base_size = 14, 
                                   base_family = "Helvetica",
                                   border = TRUE),
              add = c("boxplot", "point"), add.params = list(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3)) +
  scale_x_discrete(breaks = c("0", "1"), labels = c("Untreated", "Treated")) +
  labs(title = "", subtitle = "Rare species (m15)", caption = "") +
  stat_compare_means(method = "wilcox.test", label.y = 3.5, size = 4.5, label = "p.signif",
                     family = "Helvetica", label.x.npc = "left", hide.ns = TRUE)


#calculate effect size
wilcox_effsize(metadata_d_merge, shannon_core ~ neonatal_antibiotics,
               ci = TRUE)
wilcox_test(metadata_d_merge, shannon_core ~ neonatal_antibiotics)


wilcox_effsize(metadata_d_merge, shannon_rare ~ neonatal_antibiotics,
               ci = TRUE)
wilcox_test(metadata_d_merge, shannon_rare ~ neonatal_antibiotics)

core_rare_450_antibiotics <- 
  ggarrange(core_450_antibiotics, rare_450_antibiotics)


#######################################################################################################
#create alpha diversity figure (figure 2)

core_rare_antibiotics <-
  ggarrange(core_rare_272_antibiotics, core_rare_450_antibiotics, labels = c("B"),
            ncol = 2)
core_rare_antibiotics


ggline <- ggarrange(ggline, labels = c("A"))
ggline


ggarrange(ggline, core_rare_antibiotics, nrow = 2, heights = c(1,0.8))


#######################################################################################################
#######################################################################################################
#beta diversity

#beta diversity across all sampling time points in the preterm cohort
set.seed(2)

mds_preterm <- metaMDS((bphc_t+0.000001), distance = "bray", k = 2, trymax = 100, 
                       autotransform = TRUE)

stressplot(mds_preterm)

data.scores_preterm <- as.data.frame(scores(mds_preterm))

#order files
metadata <- metadata[order(rownames(metadata)),]

data.scores_preterm$gender <- metadata$gender
data.scores_preterm$timepoint_sampling <- metadata$timepoint_sampling
data.scores_preterm$gestational_age_group <- metadata$gestational_age_group
data.scores_preterm$delivery_details <- metadata$delivery_details
data.scores_preterm$group <- metadata$group
data.scores_preterm$surfactant <- metadata$surfactant
data.scores_preterm$ventilation <- metadata$ventilation
data.scores_preterm$previous_antibiotics <- metadata$previous_antibiotics
data.scores_preterm$current_antibiotics <- metadata$current_antibiotics
data.scores_preterm$neonatal_antibiotics <- metadata$neonatal_antibiotics
data.scores_preterm$shannon_index <- metadata$shannon_index
data.scores_preterm$specnumber <- metadata$specnumber

#permutate metadata or species
envit_preterm_metadata <- data.scores_preterm
envit_preterm_metadata$NMDS1 <- NULL
envit_preterm_metadata$NMDS2 <- NULL

envit_preterm_metadata <- data.frame(cbind(envit_preterm_metadata, bphc_t))

data.envit_preterm <- envfit(mds_preterm, envit_preterm_metadata, permutation = 999, 
                             na.rm = TRUE)
data.envit_preterm

nmds_preterm <-
  ggplot(data.scores_preterm) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint_sampling),
                 shape = as.factor(ventilation)),
             size = 3) +
  theme_pubr(border = TRUE) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint_sampling), 
                   fill=as.factor(timepoint_sampling)), 
               geom="polygon", alpha=0.1) +
  scale_y_continuous(limits = c(-1.7, 1.5),
                     breaks = c(-1.0, 0, 1.0),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(limits = c(-2, 1.5),
                     breaks = c(-2.0, -1.0, 0, 1.0),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_color_discrete(breaks = c(1,2,3,4), 
                       labels = c("w1", "m1", "m9", "m15"),
                       name ="Age") +
  scale_fill_discrete(labels= c("w1", "m1", "m9", "m15"),
                      name ="Agr") +
  scale_shape_discrete(breaks = c("0", "1", "2", "3", "4"),
                       labels = c("None","Intubation", "CPAP", "High-flow", "Low-flow"),
                       name = "Ventilation mode") +
  theme(legend.position = "right") +
  guides(fill = "none",
         color = guide_legend(order = 1),
         shape = guide_legend(order = 2))

nmds_preterm

#calculate distances from each point to group centroid
vegdist_object <- vegdist(bphc_t+0.000001, method="bray", binary=FALSE)
metadata$timepoint_sampling_fac <- as.factor(as.numeric(metadata$timepoint_sampling))
timepoint_vector <- factor(metadata$timepoint_sampling_fac, levels = c(1,2,3,4))

betadisper_result <-
  betadisper(vegdist_object, timepoint_vector, 
             type = c("centroid"), 
             bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

betadisper_df_distances <- 
  data.frame(value=betadisper_result$distances,
             group=timepoint_vector)

pairwise_comp <- list(#c(1,2),
  c(1,3),
  c(1,4),
  c(2,3),
  c(2,4))
#c(3,4))

stat.test <- betadisper_df_distances %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")
stat.test <- stat.test %>%
  add_xy_position(x = "group", fun = "max")
stat.test


epsilonSquared(x = betadisper_df_distances$value, g = betadisper_df_distances$group,
               ci = TRUE)


distance_preterm <-
  ggplot(betadisper_df_distances, aes(x=group, y=value, color=group)) +
  geom_boxplot() +
  geom_point() +
  theme_pubr(border = TRUE) +
  scale_y_continuous(limits = c(0.3,1)) +
  scale_x_discrete(breaks = c(1,2,3,4),
                   labels = c("w1", "m1", "m9", "m15")) +
  labs(x="Age", y="Distance to centroid") +
  stat_compare_means(comparisons = pairwise_comp, label = "p.signif",
                     hide.ns = TRUE) +
  guides(color = FALSE) 

distance_preterm

"anova"(betadisper_result)


######################################################################################################
#beta diversity preterm m15/sampling d and age-matched healthy controls
list_samples_c <- c("19_003c", "19_004c", "19_013c", "19_020c", "19_027c", "19_028c", "19_029c",
                    "19_030c", "19_033c", "19_037c", "19_041c", "19_046c", "19_050c", "20_004c",
                    "20_017c", "20_021c", "20_022c", "20_023c", "20_028c", "20_041c", "20_052c",
                    "20_058c", "20_061c", 
                    "KGCF56", "KGCF41", "KGCF49", "KGCF48", "KGCF44", "KGCF16", "KGCF53", "KGCF58", 
                    "KGCF57", "KGCF51", "KGCF45", "KGCF55", "KGCF04")

list_samples_d <- c("19_003d", "19_004d", "19_013d", "19_027d", "19_028d", "19_029d", "19_030d",
                    "19_033d", "19_037d", "19_041d", "19_046d", "19_050d", "20_004d", "20_017d",
                    "20_021d", "20_022d", "20_023d",
                    "KGCF16", "KGCF53", "KGCF58", "KGCF57", "KGCF51", "KGCF45", "KGCF55", "KGCF04",
                    "KGCF36", "KGCF11", "KGCF50")

bphc_selected_preterm_healthy_t <- 
  bphc_preterm_healthy_t[rownames(bphc_preterm_healthy_t) %in% list_samples_d,]

metadata_selected_preterm_healthy <-
  metadata_preterm_healthy[rownames(metadata_preterm_healthy) %in% list_samples_d,]


set.seed(2)

mds_selected_preterm_control <- metaMDS((bphc_selected_preterm_healthy_t), distance = "bray", k = 2, 
                                        trymax = 100, autotransform = TRUE)

stressplot(mds_selected_preterm_control)

data.scores_selected_preterm_control <- as.data.frame(scores(mds_selected_preterm_control))

#order files
metadata_selected_preterm_healthy <-
  metadata_selected_preterm_healthy[order(rownames(metadata_selected_preterm_healthy)),]


data.scores_selected_preterm_control$gender <- metadata_selected_preterm_healthy$gender
data.scores_selected_preterm_control$state <- metadata_selected_preterm_healthy$state
data.scores_selected_preterm_control$antimicrobial_therapy <-
  metadata_selected_preterm_healthy$antimicrobial_therapy
data.scores_selected_preterm_control$group <- metadata_selected_preterm_healthy$group
data.scores_selected_preterm_control$age_sampling <- metadata_selected_preterm_healthy$age_sampling
data.scores_selected_preterm_control$shannon_index <- metadata_selected_preterm_healthy$shannon_index
data.scores_selected_preterm_control$specnumber <- metadata_selected_preterm_healthy$specnumber

#permutate metadata or species
envit_selected_preterm_control_metadata <- data.scores_selected_preterm_control
envit_selected_preterm_control_metadata$NMDS1 <- NULL
envit_selected_preterm_control_metadata$NMDS2 <- NULL

envit_selected_preterm_control_metadata <- data.frame(cbind(envit_selected_preterm_control_metadata, 
                                                            bphc_selected_preterm_healthy_t))

data.envit_selected_preterm_control <- envfit(mds_selected_preterm_control,
                                              envit_selected_preterm_control_metadata,
                                              permutation = 999, na.rm = TRUE)
data.envit_selected_preterm_control

data.scores_selected_preterm_control$state <- factor(data.scores_selected_preterm_control$state,
                                                     levels = c("Preterm", "Healthy"))

levels(data.scores_selected_preterm_control$state)

nmds_selected_preterm_healthy <-
  ggplot(data.scores_selected_preterm_control,
         aes(x=NMDS1, y=NMDS2, color=as.factor(state))) +
  geom_point(size = 3) +
  theme_pubr(border = TRUE) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=as.factor(state)), 
               geom="polygon", alpha=0.1) +
  scale_y_continuous(limits = c(-1.7, 1.5),
                     breaks = c(-1, 0, 1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(limits = c(-2, 1.5),
                     breaks = c(-2, -1, 0, 1),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(breaks = c("Preterm", "Healthy"),
                     labels = c("Preterm", "Full-term"),
                     name = "State",
                     values = c("Preterm" = "slateblue4",
                                "Healthy" = "grey60")) +
  scale_fill_discrete(breaks = c("Preterm", "Healthy"),
                      labels= c("Preterm", "Full-term"),
                      name = "State") +
  scale_fill_manual(values = c("Preterm" = "slateblue4",
                               "Healthy" = "grey60")) +
  theme(legend.position = "right") +
  guides(fill = "none", color = guide_legend(order = 1))

nmds_selected_preterm_healthy

#calculate distances from each point to group centroid
vegdist_object_selected_preterm_healthy <- 
  vegdist(bphc_selected_preterm_healthy_t, method="bray", binary=FALSE)

betadisper_result_selected_preterm_healthy <-
  betadisper(vegdist_object_selected_preterm_healthy, metadata_selected_preterm_healthy$state, 
             type = c("centroid"), 
             bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

betadisper_df_distances_selected_preterm_healthy <- 
  data.frame(value=betadisper_result_selected_preterm_healthy$distances,
             group=metadata_selected_preterm_healthy$state)

betadisper_df_distances_selected_preterm_healthy$group <-
  factor(betadisper_df_distances_selected_preterm_healthy$group,
         levels = c("Preterm", "Healthy"))

distance_selected_preterm_healthy <-
  ggplot(betadisper_df_distances_selected_preterm_healthy, aes(x=group, y=value, color=group)) +
  geom_boxplot() +
  geom_point() +
  theme_pubr(border = TRUE) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_x_discrete(breaks = c("Healthy", "Preterm"),
                   labels = c("Healthy full-term", "Preterm")) +
  scale_color_manual(breaks = c("Healthy", "Preterm"),
                     labels = c("Healthy full-term", "Preterm"),
                     name ="State",
                     values = c("Healthy" = "grey60",
                                "Preterm" = "slateblue4")) +
  labs(x="State", y="Distance to centroid") +
  stat_compare_means(label = "p.signif", label.y = 0.85,
                     hide.ns = TRUE) +
  #stat_compare_means(label.y = 1) +
  guides(color = FALSE) 

distance_selected_preterm_healthy

"anova"(betadisper_result_preterm_healthy)

#calculate effect size
wilcox_effsize(betadisper_df_distances_selected_preterm_healthy, value ~ group,
               ci = TRUE)
wilcox_test(betadisper_df_distances_selected_preterm_healthy, value ~ group)

ggarrange(nmds_preterm, distance_preterm, 
          nmds_selected_preterm_healthy, distance_selected_preterm_healthy,
          ncol = 2, nrow = 2, widths = c(1, 0.4, 1, 0.4),
          labels = c("A", "B", "C", "D"))


######################################################################################################
######################################################################################################
######################################################################################################
#plot contamination/ negative controls

raspir_readcount <- read_delim("read_count_raspir_blanks_samples.csv", 
                               ";", escape_double = FALSE, trim_ws = TRUE)

raspir_readcount <- data.frame(raspir_readcount)
rownames(raspir_readcount) <- raspir_readcount$species
raspir_readcount$species <- NULL
raspir_readcount_t <- data.frame(t(raspir_readcount))
raspir_readcount_t$microbial <- as.numeric(raspir_readcount_t$microbial)
raspir_readcount_t$blank <- rownames(raspir_readcount_t)

readcount_blanks <- subset(raspir_readcount_t, type == "control")

#####################################################################################################
#########co-occurrence networks
# Preterms, time-point a
preterm_count_00_a <- species_a
names(preterm_count_00_a) <- gsub("\\.", "_", names(preterm_count_00_a))

# prevalence filtering
prev_filt_abund = 90
preterm_count_01_a <- data.frame(t(preterm_count_00_a))
n_preterm_a=round(length(colnames(preterm_count_01_a)) * prev_filt_abund / 100) 
preterm_count_01_a[preterm_count_01_a == 0] <- NA
preterm_count_01_a <- preterm_count_01_a[rowSums(is.na(preterm_count_01_a)) < n_preterm_a, ]
preterm_count_01_a[is.na(preterm_count_01_a)] <- 0
preterm_count_02_a <- compositions::clr(preterm_count_01_a)
preterm_count_03_a <- data.frame(t(preterm_count_02_a))

# perform Spearman's rank correlation analysis
spearman_a <- rcorr(as.matrix(preterm_count_03_a), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_a_00 <- spearman_a$P

# convert data to long format and stack columns into a single column of data for p-values
spearman_a_edges <- reshape2::melt(spearman_p_a_00)
spearman_a_edges <- na.omit(spearman_a_edges)
colnames(spearman_a_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_a <- spearman_a$r

# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_a_cor_edges <- reshape2::melt(spearman_r_a)
spearman_a_cor_edges <- na.omit(spearman_a_cor_edges)
spearman_a_cor_edges$value <- round(spearman_a_cor_edges$value, 5)
colnames(spearman_a_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_a <- spearman_a_edges %>% inner_join(spearman_a_cor_edges, by=c("node1","node2"))

spearman_table_a <- spearman_table_a[spearman_table_a$p<0.01,]
spearman_table_a$Id <- spearman_table_a$node1
# rename columns
colnames(spearman_table_a) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_a) <- NULL
# extract significant correlations
spearman_table_a <- subset(spearman_table_a, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_a$Correlation <- ifelse(spearman_table_a$Weight > 0, "pos", "neg")
spearman_table_a$Type <- "Undirected"

# generate edge list
spearman_edges_a <- select(spearman_table_a, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_a$Genus <- sapply(strsplit(as.character(spearman_edges_a$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_a, file="edges_a.csv", sep=";", col.names = TRUE, row.names = FALSE)

# generate node list
spearman_nodes_a <- select(spearman_edges_a, c("Target"))
spearman_nodes_a$Target2 <- spearman_nodes_a$Target
spearman_nodes_a$Correlation <- spearman_edges_a$Correlation
# remove duplicate entries
spearman_nodes_a = spearman_nodes_a[!duplicated(spearman_nodes_a$Target),]
# make data frame
spearman_nodes_a <- data.frame(spearman_nodes_a)
# re-index rows
rownames(spearman_nodes_a) <- NULL
# rename columns
colnames(spearman_nodes_a) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_a$Genus <- sapply(strsplit(as.character(spearman_nodes_a$Label),"_"), `[`, 1)

# spearman_nodes_a$mean_abundance <- mean_list_a
# export node list
write.table(spearman_nodes_a, file="nodes_a.csv", sep=";", col.names = TRUE, row.names = FALSE)



###############################################################################
# Preterms, time-point b
preterm_count_00_b <- species_b
names(preterm_count_00_b) <- gsub("\\.", "_", names(preterm_count_00_b))

# prevalence filtering
preterm_count_01_b <- data.frame(t(preterm_count_00_b))
n_preterm_b=round(length(colnames(preterm_count_01_b)) * prev_filt_abund / 100) 
preterm_count_01_b[preterm_count_01_b == 0] <- NA
preterm_count_01_b <- preterm_count_01_b[rowSums(is.na(preterm_count_01_b)) < n_preterm_b, ]
preterm_count_01_b[is.na(preterm_count_01_b)] <- 0
preterm_count_02_b <- compositions::clr(preterm_count_01_b)
preterm_count_03_b <- data.frame(t(preterm_count_02_b))

# perform Spearman's rank correlation analysis
spearman_b <- rcorr(as.matrix(preterm_count_03_b), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_b_00 <- spearman_b$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_b_edges <- reshape2::melt(spearman_p_b_00)
spearman_b_edges <- na.omit(spearman_b_edges)
colnames(spearman_b_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_b <- spearman_b$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_b_cor_edges <- reshape2::melt(spearman_r_b)
spearman_b_cor_edges <- na.omit(spearman_b_cor_edges)
spearman_b_cor_edges$value <- round(spearman_b_cor_edges$value, 5)
colnames(spearman_b_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_b <- spearman_b_edges %>% inner_join(spearman_b_cor_edges, by=c("node1","node2"))
spearman_table_b <- spearman_table_b[spearman_table_b$p<0.01,]
spearman_table_b$Id <- spearman_table_b$node1
# rename columns
colnames(spearman_table_b) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_b) <- NULL
# extract significant correlations
spearman_table_b <- subset(spearman_table_b, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_b$Correlation <- ifelse(spearman_table_b$Weight > 0, "pos", "neg")
spearman_table_b$Type <- "Undirected"

# generate edge list
spearman_edges_b <- select(spearman_table_b, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_b$Genus <- sapply(strsplit(as.character(spearman_edges_b$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_b, file="edges_b_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

# generate node list
spearman_nodes_b <- select(spearman_edges_b, c("Target"))
spearman_nodes_b$Target2 <- spearman_nodes_b$Target
spearman_nodes_b$Correlation <- spearman_edges_b$Correlation
# remove duplicate entries
spearman_nodes_b = spearman_nodes_b[!duplicated(spearman_nodes_b$Target),]
# make data frame
spearman_nodes_b <- data.frame(spearman_nodes_b)
# re-index rows
rownames(spearman_nodes_b) <- NULL
# rename columns
colnames(spearman_nodes_b) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_b$Genus <- sapply(strsplit(as.character(spearman_nodes_b$Label),"_"), `[`, 1)

# spearman_nodes_b$mean_abundance <- mean_list_b
# export node list
write.table(spearman_nodes_b, file="nodes_b_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)



###############################################################################
# Preterms, time-point c
preterm_count_00_c <- species_c
names(preterm_count_00_c) <- gsub("\\.", "_", names(preterm_count_00_c))

# prevalence filtering
preterm_count_01_c <- data.frame(t(preterm_count_00_c))
n_preterm_c=round(length(colnames(preterm_count_01_c)) * prev_filt_abund / 100) 
preterm_count_01_c[preterm_count_01_c == 0] <- NA
preterm_count_01_c <- preterm_count_01_c[rowSums(is.na(preterm_count_01_c)) < n_preterm_c, ]
preterm_count_01_c[is.na(preterm_count_01_c)] <- 0
preterm_count_02_c <- compositions::clr(preterm_count_01_c)
preterm_count_03_c <- data.frame(t(preterm_count_02_c))

# perform Spearman's rank correlation analysis
spearman_c <- rcorr(as.matrix(preterm_count_03_c), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_c_00 <- spearman_c$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_c_edges <- reshape2::melt(spearman_p_c_00)
spearman_c_edges <- na.omit(spearman_c_edges)
colnames(spearman_c_edges) <- c("node1", "node2", "p")


# extract and store correlation coefficients of analysis
spearman_r_c <- spearman_c$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_c_cor_edges <- reshape2::melt(spearman_r_c)
spearman_c_cor_edges <- na.omit(spearman_c_cor_edges)
spearman_c_cor_edges$value <- round(spearman_c_cor_edges$value, 5)
colnames(spearman_c_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_c <- spearman_c_edges %>% inner_join(spearman_c_cor_edges, by=c("node1","node2"))
spearman_table_c <- spearman_table_c[spearman_table_c$p<0.01,]
spearman_table_c$Id <- spearman_table_c$node1
# rename columns
colnames(spearman_table_c) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_c) <- NULL
# extract significant correlations
spearman_table_c <- subset(spearman_table_c, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_c$Correlation <- ifelse(spearman_table_c$Weight > 0, "pos", "neg")
spearman_table_c$Type <- "Undirected"

# generate edge list
spearman_edges_c <- select(spearman_table_c, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_c$Genus <- sapply(strsplit(as.character(spearman_edges_c$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_c, file="edges_c_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)

# generate node list
spearman_nodes_c <- select(spearman_edges_c, c("Target"))
spearman_nodes_c$Target2 <- spearman_nodes_c$Target
spearman_nodes_c$Correlation <- spearman_edges_c$Correlation
# remove duplicate entries
spearman_nodes_c = spearman_nodes_c[!duplicated(spearman_nodes_c$Target),]
# make data frame
spearman_nodes_c <- data.frame(spearman_nodes_c)
# re-index rows
rownames(spearman_nodes_c) <- NULL
# rename columns
colnames(spearman_nodes_c) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_c$Genus <- sapply(strsplit(as.character(spearman_nodes_c$Label),"_"), `[`, 1)
# spearman_nodes_c$mean_abundance <- mean_list_c
# export node list
write.table(spearman_nodes_c, file="nodes_c_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)



###############################################################################
# Preterms, time-point d
preterm_count_00_d <- species_d
names(preterm_count_00_d) <- gsub("\\.", "_", names(preterm_count_00_d))

# prevalence filtering
preterm_count_01_d <- data.frame(t(preterm_count_00_d))
n_preterm_d=round(length(colnames(preterm_count_01_d)) * prev_filt_abund / 100) 
preterm_count_01_d[preterm_count_01_d == 0] <- NA
preterm_count_01_d <- preterm_count_01_d[rowSums(is.na(preterm_count_01_d)) < n_preterm_d, ]
preterm_count_01_d[is.na(preterm_count_01_d)] <- 0
preterm_count_02_d <- compositions::clr(preterm_count_01_d)
preterm_count_03_d <- data.frame(t(preterm_count_02_d))

# perform Spearman's rank correlation analysis
spearman_d <- rcorr(as.matrix(preterm_count_03_d), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_d_00 <- spearman_d$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_d_edges <- reshape2::melt(spearman_p_d_00)
spearman_d_edges <- na.omit(spearman_d_edges)
colnames(spearman_d_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_d <- spearman_d$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_d_cor_edges <- reshape2::melt(spearman_r_d)
spearman_d_cor_edges <- na.omit(spearman_d_cor_edges)
spearman_d_cor_edges$value <- round(spearman_d_cor_edges$value, 5)
colnames(spearman_d_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_d <- spearman_d_edges %>% inner_join(spearman_d_cor_edges, by=c("node1","node2"))
spearman_table_d <- spearman_table_d[spearman_table_d$p<0.01,]
spearman_table_d$Id <- spearman_table_d$node1
# rename columns
colnames(spearman_table_d) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_d) <- NULL
# extract significant correlations
spearman_table_d <- subset(spearman_table_d, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_d$Correlation <- ifelse(spearman_table_d$Weight > 0, "pos", "neg")
spearman_table_d$Type <- "Undirected"

# generate edge list
spearman_edges_d <- select(spearman_table_d, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_d$Genus <- sapply(strsplit(as.character(spearman_edges_d$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_d, file="edges_d_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)

# generate node list
spearman_nodes_d <- select(spearman_edges_d, c("Target"))
spearman_nodes_d$Target2 <- spearman_nodes_d$Target
spearman_nodes_d$Correlation <- spearman_edges_d$Correlation
# remove duplicate entries
spearman_nodes_d = spearman_nodes_d[!duplicated(spearman_nodes_d$Target),]
# make data frame
spearman_nodes_d <- data.frame(spearman_nodes_d)
# re-index rows
rownames(spearman_nodes_d) <- NULL
# rename columns
colnames(spearman_nodes_d) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_d$Genus <- sapply(strsplit(as.character(spearman_nodes_d$Label),"_"), `[`, 1)
# spearman_nodes_d$mean_abundance <- mean_list_d
# export node list
write.table(spearman_nodes_d, file="nodes_d_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)



###############################################################################
# healthy
fullterm_count_00 <- bphc_preterm_healthy_t[grepl("KGCF", 
                                                  rownames(bphc_preterm_healthy_t)),]
names(fullterm_count_00) <- gsub("\\.", "_", names(fullterm_count_00))

# prevalence filtering
fullterm_count_01 <- data.frame(t(fullterm_count_00))
n_fullterm=round(length(colnames(fullterm_count_01)) * prev_filt_abund / 100) 
fullterm_count_01[fullterm_count_01 == 0] <- NA
fullterm_count_01 <- fullterm_count_01[rowSums(is.na(fullterm_count_01)) < n_fullterm, ]
fullterm_count_01[is.na(fullterm_count_01)] <- 0
fullterm_count_02 <- compositions::clr(fullterm_count_01)
fullterm_count_03 <- data.frame(t(fullterm_count_02))

# perform Spearman's rank correlation analysis
spearman_full <- rcorr(as.matrix(fullterm_count_03), type ='spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_full_00 <- spearman_full$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_full_edges <- reshape2::melt(spearman_p_full_00)
spearman_full_edges <- na.omit(spearman_full_edges)
colnames(spearman_full_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_full <- spearman_full$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_full_cor_edges <- reshape2::melt(spearman_r_full)
spearman_full_cor_edges <- na.omit(spearman_full_cor_edges)
spearman_full_cor_edges$value <- round(spearman_full_cor_edges$value, 5)
colnames(spearman_full_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_full <- spearman_full_edges %>% merge(spearman_full_cor_edges, by=c("node1","node2"))
spearman_table_full <- spearman_table_full[spearman_table_full$p<0.01,]
spearman_table_full$Id <- spearman_table_full$node1
# rename columns
colnames(spearman_table_full) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_full) <- NULL
# extract significant correlations
spearman_table_full <- subset(spearman_table_full, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_full$Correlation <- ifelse(spearman_table_full$Weight > 0, "pos", "neg")
spearman_table_full$Type <- "Undirected"

# generate edge list
spearman_edges_full <- select(spearman_table_full, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_full$Genus <- sapply(strsplit(as.character(spearman_edges_full$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_full, file="edges_full_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)

# generate node list
spearman_nodes_full <- select(spearman_edges_full, c("Target"))
spearman_nodes_full$Target2 <- spearman_nodes_full$Target
spearman_nodes_full$Correlation <- spearman_edges_full$Correlation
# remove duplicate entries
spearman_nodes_full = spearman_nodes_full[!duplicated(spearman_nodes_full$Target),]
# make data frame
spearman_nodes_full <- data.frame(spearman_nodes_full)
# re-index rows
rownames(spearman_nodes_full) <- NULL
# rename columns
colnames(spearman_nodes_full) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_full$Genus <- sapply(strsplit(as.character(spearman_nodes_full$Label),"_"), `[`, 1)

# spearman_nodes_full$mean_abundance <- mean_list_full
# export node list
write.table(spearman_nodes_full, file="nodes_full_pos.csv", sep=";", col.names = TRUE, row.names = FALSE)

#Healthy full-terms matched c (8m-12m)
fullterm_m9_count <- 
  bphc_preterm_healthy_t[rownames(bphc_preterm_healthy_t) %in% list_samples_c,]
fullterm_m9_count_00 <- fullterm_m9_count[grepl("KGCF", 
                                                rownames(fullterm_m9_count)),]
names(fullterm_m9_count_00) <- gsub("\\.", "_", names(fullterm_m9_count_00))

# prevalence filtering
fullterm_m9_count_01 <- data.frame(t(fullterm_m9_count_00))
n_fullterm_m9=round(length(colnames(fullterm_m9_count_01)) * prev_filt_abund / 100) 
fullterm_m9_count_01[fullterm_m9_count_01 == 0] <- NA
fullterm_m9_count_01 <- fullterm_m9_count_01[rowSums(is.na(fullterm_m9_count_01)) < n_fullterm_m9, ]
fullterm_m9_count_01[is.na(fullterm_m9_count_01)] <- 0
fullterm_m9_count_02 <- compositions::clr(fullterm_m9_count_01)
fullterm_m9_count_03 <- data.frame(t(fullterm_m9_count_02))

# perform Spearman's rank correlation analysis
spearman_full_m9 <- rcorr(as.matrix(fullterm_m9_count_03), type ='spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_full_m9_00 <- spearman_full_m9$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_full_m9_edges <- reshape2::melt(spearman_p_full_m9_00)
spearman_full_m9_edges <- na.omit(spearman_full_m9_edges)
colnames(spearman_full_m9_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_full_m9 <- spearman_full_m9$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_full_m9_cor_edges <- reshape2::melt(spearman_r_full_m9)
spearman_full_m9_cor_edges <- na.omit(spearman_full_m9_cor_edges)
spearman_full_m9_cor_edges$value <- round(spearman_full_m9_cor_edges$value, 5)
colnames(spearman_full_m9_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_full_m9 <- spearman_full_m9_edges %>% merge(spearman_full_m9_cor_edges, by=c("node1","node2"))
spearman_table_full_m9 <- spearman_table_full_m9[spearman_table_full_m9$p<0.01,]
spearman_table_full_m9$Id <- spearman_table_full_m9$node1
# rename columns
colnames(spearman_table_full_m9) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_full_m9) <- NULL
# extract significant correlations
spearman_table_full_m9 <- subset(spearman_table_full_m9, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_full_m9$Correlation <- ifelse(spearman_table_full_m9$Weight > 0, "pos", "neg")
spearman_table_full_m9$Type <- "Undirected"

# generate edge list
spearman_edges_full_m9 <- select(spearman_table_full_m9, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_full_m9$Genus <- sapply(strsplit(as.character(spearman_edges_full_m9$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_full_m9, file="edges_full_m9_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

# generate node list
spearman_nodes_full_m9 <- select(spearman_edges_full_m9, c("Target"))
spearman_nodes_full_m9$Target2 <- spearman_nodes_full_m9$Target
spearman_nodes_full_m9$Correlation <- spearman_edges_full_m9$Correlation
# remove duplicate entries
spearman_nodes_full_m9 = spearman_nodes_full_m9[!duplicated(spearman_nodes_full_m9$Target),]
# make data frame
spearman_nodes_full_m9 <- data.frame(spearman_nodes_full_m9)
# re-index rows
rownames(spearman_nodes_full_m9) <- NULL
# rename columns
colnames(spearman_nodes_full_m9) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_full_m9$Genus <- sapply(strsplit(as.character(spearman_nodes_full_m9$Label),"_"), `[`, 1)

# spearman_nodes_full_m9$mean_abundance <- mean_list_full_m9
# export node list
write.table(spearman_nodes_full_m9, file="nodes_full_m9_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

#Healthy full-terms matched d (10m-14m)
fullterm_m15_count <- 
  bphc_preterm_healthy_t[rownames(bphc_preterm_healthy_t) %in% list_samples_d,]
fullterm_m15_count_00 <- fullterm_m15_count[grepl("KGCF", 
                                                  rownames(fullterm_m15_count)),]
names(fullterm_m15_count_00) <- gsub("\\.", "_", names(fullterm_m15_count_00))

# prevalence filtering
fullterm_m15_count_01 <- data.frame(t(fullterm_m15_count_00))
n_fullterm_m15=round(length(colnames(fullterm_m15_count_01)) * prev_filt_abund / 100) 
fullterm_m15_count_01[fullterm_m15_count_01 == 0] <- NA
fullterm_m15_count_01 <- fullterm_m15_count_01[rowSums(is.na(fullterm_m15_count_01)) < n_fullterm_m15, ]
fullterm_m15_count_01[is.na(fullterm_m15_count_01)] <- 0
fullterm_m15_count_02 <- compositions::clr(fullterm_m15_count_01)
fullterm_m15_count_03 <- data.frame(t(fullterm_m15_count_02))

# perform Spearman's rank correlation analysis
spearman_full_m15 <- rcorr(as.matrix(fullterm_m15_count_03), type ='spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_full_m15_00 <- spearman_full_m15$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_full_m15_edges <- reshape2::melt(spearman_p_full_m15_00)
spearman_full_m15_edges <- na.omit(spearman_full_m15_edges)
colnames(spearman_full_m15_edges) <- c("node1", "node2", "p")

# extract and store correlation coefficients of analysis
spearman_r_full_m15 <- spearman_full_m15$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_full_m15_cor_edges <- reshape2::melt(spearman_r_full_m15)
spearman_full_m15_cor_edges <- na.omit(spearman_full_m15_cor_edges)
spearman_full_m15_cor_edges$value <- round(spearman_full_m15_cor_edges$value, 5)
colnames(spearman_full_m15_cor_edges) <- c("node1", "node2", "r")

# merge tables
spearman_table_full_m15 <- spearman_full_m15_edges %>% merge(spearman_full_m15_cor_edges, by=c("node1","node2"))
spearman_table_full_m15 <- spearman_table_full_m15[spearman_table_full_m15$p<0.01,]
spearman_table_full_m15$Id <- spearman_table_full_m15$node1
# rename columns
colnames(spearman_table_full_m15) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_full_m15) <- NULL
# extract significant correlations
spearman_table_full_m15 <- subset(spearman_table_full_m15, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_full_m15$Correlation <- ifelse(spearman_table_full_m15$Weight > 0, "pos", "neg")
spearman_table_full_m15$Type <- "Undirected"

# generate edge list
spearman_edges_full_m15 <- select(spearman_table_full_m15, c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_full_m15$Genus <- sapply(strsplit(as.character(spearman_edges_full_m15$Target),"_"), `[`, 1)
# export edge list
write.table(spearman_edges_full_m15, file="edges_full_m15_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

# generate node list
spearman_nodes_full_m15 <- select(spearman_edges_full_m15, c("Target"))
spearman_nodes_full_m15$Target2 <- spearman_nodes_full_m15$Target
spearman_nodes_full_m15$Correlation <- spearman_edges_full_m15$Correlation
# remove duplicate entries
spearman_nodes_full_m15 = spearman_nodes_full_m15[!duplicated(spearman_nodes_full_m15$Target),]
# make data frame
spearman_nodes_full_m15 <- data.frame(spearman_nodes_full_m15)
# re-index rows
rownames(spearman_nodes_full_m15) <- NULL
# rename columns
colnames(spearman_nodes_full_m15) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_full_m15$Genus <- sapply(strsplit(as.character(spearman_nodes_full_m15$Label),"_"), `[`, 1)

# spearman_nodes_full_m15$mean_abundance <- mean_list_full_m15
# export node list
write.table(spearman_nodes_full_m15, file="nodes_full_m15_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

##matched preterm sample sizes for m9 and m15

##m9

preterm_count_03_c_subset <- preterm_count_03_c[sample(nrow(preterm_count_03_c), 13), ]

# perform Spearman's rank correlation analysis
spearman_c_subset <- rcorr(as.matrix(preterm_count_03_c_subset), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_c_00_subset <- spearman_c_subset$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_c_edges_subset <- reshape2::melt(spearman_p_c_00_subset)
spearman_c_edges_subset <- na.omit(spearman_c_edges_subset)
colnames(spearman_c_edges_subset) <- c("node1", "node2", "p")


# extract and store correlation coefficients of analysis
spearman_r_c_subset <- spearman_c_subset$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_c_cor_edges_subset <- reshape2::melt(spearman_r_c_subset)
spearman_c_cor_edges_subset <- na.omit(spearman_c_cor_edges_subset)
spearman_c_cor_edges_subset$value <- round(spearman_c_cor_edges_subset$value, 5)
colnames(spearman_c_cor_edges_subset) <- c("node1", "node2", "r")

# merge tables
spearman_table_c_subset <- spearman_c_edges_subset %>% inner_join(spearman_c_cor_edges_subset, 
                                                                  by=c("node1","node2"))
spearman_table_c_subset <- spearman_table_c_subset[spearman_table_c_subset$p<0.01,]
spearman_table_c_subset$Id <- spearman_table_c_subset$node1
# rename columns
colnames(spearman_table_c_subset) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_c_subset) <- NULL
# extract significant correlations
spearman_table_c_subset <- subset(spearman_table_c_subset, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_c_subset$Correlation <- ifelse(spearman_table_c_subset$Weight > 0, "pos", "neg")
spearman_table_c_subset$Type <- "Undirected"

# generate edge list
spearman_edges_c_subset <- select(spearman_table_c_subset, 
                                  c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_c_subset$Genus <- sapply(strsplit(as.character(spearman_edges_c_subset$Target),"_"), 
                                        `[`, 1)
# export edge list
write.table(spearman_edges_c_subset, file="edges_c_subset_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

# generate node list
spearman_nodes_c_subset <- select(spearman_edges_c_subset, c("Target"))
spearman_nodes_c_subset$Target2 <- spearman_nodes_c_subset$Target
spearman_nodes_c_subset$Correlation <- spearman_edges_c_subset$Correlation
# remove duplicate entries
spearman_nodes_c_subset = spearman_nodes_c_subset[!duplicated(spearman_nodes_c_subset$Target),]
# make data frame
spearman_nodes_c_subset <- data.frame(spearman_nodes_c_subset)
# re-index rows
rownames(spearman_nodes_c_subset) <- NULL
# rename columns
colnames(spearman_nodes_c_subset) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_c_subset$Genus <- sapply(strsplit(as.character(spearman_nodes_c_subset$Label),"_"),
                                        `[`, 1)
# spearman_nodes_c$mean_abundance <- mean_list_c
# export node list
write.table(spearman_nodes_c_subset, file="nodes_c_subset_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

##m15
preterm_count_03_d_subset <- preterm_count_03_d[sample(nrow(preterm_count_03_d), 11), ]

# perform Spearman's rank correlation analysis
spearman_d_subset <- rcorr(as.matrix(preterm_count_03_d_subset), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_d_00_subset <- spearman_d_subset$P
# convert data to long format and stack columns into a single column of data for p-values
spearman_d_edges_subset <- reshape2::melt(spearman_p_d_00_subset)
spearman_d_edges_subset <- na.omit(spearman_d_edges_subset)
colnames(spearman_d_edges_subset) <- c("node1", "node2", "p")


# extract and store correlation coefficients of analysis
spearman_r_d_subset <- spearman_d_subset$r
# convert data to long format and stack columns into a single column of data for correlation coefficients
spearman_d_cor_edges_subset <- reshape2::melt(spearman_r_d_subset)
spearman_d_cor_edges_subset <- na.omit(spearman_d_cor_edges_subset)
spearman_d_cor_edges_subset$value <- round(spearman_d_cor_edges_subset$value, 5)
colnames(spearman_d_cor_edges_subset) <- c("node1", "node2", "r")

# merge tables
spearman_table_d_subset <- spearman_d_edges_subset %>% inner_join(spearman_d_cor_edges_subset, 
                                                                  by=c("node1","node2"))
spearman_table_d_subset <- spearman_table_d_subset[spearman_table_d_subset$p<0.01,]
spearman_table_d_subset$Id <- spearman_table_d_subset$node1
# rename columns
colnames(spearman_table_d_subset) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_table_d_subset) <- NULL
# extract significant correlations
spearman_table_d_subset <- subset(spearman_table_d_subset, Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_table_d_subset$Correlation <- ifelse(spearman_table_d_subset$Weight > 0, "pos", "neg")
spearman_table_d_subset$Type <- "Undirected"

# generate edge list
spearman_edges_d_subset <- select(spearman_table_d_subset, 
                                  c("Source", "Target", "Weight","Type", "Correlation"))
spearman_edges_d_subset$Genus <- sapply(strsplit(as.character(spearman_edges_d_subset$Target),"_"), 
                                        `[`, 1)
# export edge list
write.table(spearman_edges_d_subset, file="edges_d_subset_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

# generate node list
spearman_nodes_d_subset <- select(spearman_edges_d_subset, c("Target"))
spearman_nodes_d_subset$Target2 <- spearman_nodes_d_subset$Target
spearman_nodes_d_subset$Correlation <- spearman_edges_d_subset$Correlation
# remove duplicate entries
spearman_nodes_d_subset = spearman_nodes_d_subset[!duplicated(spearman_nodes_d_subset$Target),]
# make data frame
spearman_nodes_d_subset <- data.frame(spearman_nodes_d_subset)
# re-index rows
rownames(spearman_nodes_d_subset) <- NULL
# rename columns
colnames(spearman_nodes_d_subset) <- c("Id", "Label", "Correlation")
# add genus information
spearman_nodes_d_subset$Genus <- sapply(strsplit(as.character(spearman_nodes_d_subset$Label),"_"),
                                        `[`, 1)
# spearman_nodes_c$mean_abundance <- mean_list_c
# export node list
write.table(spearman_nodes_d_subset, file="nodes_d_subset_pos.csv", sep=";", col.names = TRUE, 
            row.names = FALSE)

################################################################
#import gephi statistics

gephi_out_b <- read_csv("gephi_out_b.csv")
gephi_out_b <- data.frame(gephi_out_b)
gephi_out_b$time_point <- "preterm m1"
gephi_out_b$Experiment <- "Subgroups"
gephi_out_b$group <- "Preterm"

gephi_out_c <- read_csv("gephi_out_c.csv")
gephi_out_c <- data.frame(gephi_out_c)
gephi_out_c$time_point <- "preterm m9"
gephi_out_c$Experiment <- "Subgroups"
gephi_out_c$group <- "Preterm"

gephi_out_d <- read_csv("gephi_out_d.csv")
gephi_out_d <- data.frame(gephi_out_d)
gephi_out_d$time_point <- "preterm m15"
gephi_out_d$Experiment <- "Subgroups"
gephi_out_d$group <- "Preterm"

gephi_out_full <- read_csv("gephi_out_fullterm.csv")
gephi_out_full <- data.frame(gephi_out_full)
gephi_out_full$time_point <- "full-term m1-m14"
gephi_out_full$Experiment <- "Reference"
gephi_out_full$group <- "Reference"

gephi_out_full_m9 <- read_csv("gephi_out_fullterm_m9.csv")
gephi_out_full_m9 <- data.frame(gephi_out_full_m9)
gephi_out_full_m9$time_point <- "full-term m8-m12"
gephi_out_full_m9$Experiment <- "Subgroups"
gephi_out_full_m9$group <- "Full-term"

gephi_out_full_m15 <- read_csv("gephi_out_fullterm_m15.csv")
gephi_out_full_m15 <- data.frame(gephi_out_full_m15)
gephi_out_full_m15$time_point <- "full-term m10-m14"
gephi_out_full_m15$Experiment <- "Subgroups"
gephi_out_full_m15$group <- "Full-term"

#gephi_out_c_subset


#gephi_out_d_subset

gephi_out <- data.frame(rbind(gephi_out_b, gephi_out_c, gephi_out_d, gephi_out_full,
                              gephi_out_full_m9, gephi_out_full_m15))
#gephi_out$time_point <- factor(gephi_out$Experiment, levels = c("m1", "m9", "m15", "Reference"))
compare_time <- list(c("preterm m1", "preterm m9"), 
                     c("preterm m9", "preterm m15"),
                     c("preterm m9", "full-term m8-m12"),
                     c("preterm m15", "full-term m10-m14"))
# c("full-term m8-m12", "full-term m10-m14"))

gephi_out$time_point <- factor(gephi_out$time_point, 
                               levels = c("preterm m1", "preterm m9", "preterm m15",
                                          "full-term m8-m12", "full-term m10-m14", 
                                          "full-term m1-m14"))
gephi_out$Experiment <- factor(gephi_out$Experiment, 
                               levels = c("Reference", "Subgroups"))

dg_plot <- 
  ggplot(gephi_out, aes(x=time_point, y=Degree, colour=group)) +
  geom_violin(width=0.5) + geom_point() + 
  stat_compare_means(comparisons = compare_time, label = "p.signif", 
                     label.y = c(40, 48, 56, 64)) +
  scale_x_discrete(labels = c("preterm m1" = "Preterm\nm1", "preterm m9" = "Preterm\nm9",
                              "preterm m15" = "Preterm\nm15",
                              "full-term m8-m12" = "Full-term\nm8-m12", 
                              "full-term m10-m14" = "Full-term\nm10-m14",
                              "full-term m1-m14" = "Full-term\nm1-m14")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.position = "none") + 
  xlab("") + ylab("Degree centrality") +
  facet_grid(~Experiment, scales = "free_x", space="free") + 
  scale_colour_manual(values = c("Preterm"="black", "Reference"="gold3", 
                                 "Full-term"="grey60")) +
  scale_y_continuous(limits = c(0,72), breaks = c(0,20,40,60)) +
  geom_pointrange(mapping = aes(x = time_point, y = Degree),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.3)


dg_plot

cc_plot <-
  ggplot(gephi_out, aes(x=time_point, y=closnesscentrality, colour=group)) +
  geom_violin(width=0.5) + geom_point() + 
  stat_compare_means(comparisons = compare_time, label = "p.signif", label.y = c(1.1, 1.25,
                                                                                 1.4, 1.55)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.position = "none") + 
  xlab("") + ylab("Closeness centrality") +
  facet_grid(~Experiment, scales = "free_x", space="free") + 
  scale_colour_manual(values=c("Preterm"="black", "Reference"="gold3", "Full-term"="grey60")) +
  scale_x_discrete(labels = c("preterm m1" = "Preterm\nm1", "preterm m9" = "Preterm\nm9",
                              "preterm m15" = "Preterm\nm15",
                              "full-term m8-m12" = "Full-term\nm8-m12", 
                              "full-term m10-m14" = "Full-term\nm10-m14",
                              "full-term m1-m14" = "Full-term\nm1-m14"))  +
  scale_y_continuous(limits = c(0,1.65)) +
  geom_pointrange(mapping = aes(x = time_point, y = closnesscentrality),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.3)

cc_plot

bc_plot <-
  ggplot(gephi_out, aes(x=time_point, y=betweenesscentrality, colour=group)) +
  geom_violin(width=0.5) + geom_point() + 
  stat_compare_means(comparisons = compare_time, label = "p.signif", label.y = c(0.17, 0.20,
                                                                                 0.23, 0.26)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.position = "none") + 
  xlab("") + ylab("Betweeness centrality") +
  facet_grid(~Experiment, scales = "free_x", space="free") + 
  scale_colour_manual(values=c("Preterm"="black", "Reference"="gold3", "Full-term"="grey60")) +
  scale_x_discrete(labels = c("preterm m1" = "Preterm\nm1", "preterm m9" = "Preterm\nm9",
                              "preterm m15" = "Preterm\nm15",
                              "full-term m8-m12" = "Full-term\nm8-m12", 
                              "full-term m10-m14" = "Full-term\nm10-m14",
                              "full-term m1-m14" = "Full-term\nm1-m14"))  +
  scale_y_continuous(limits = c(0,0.28)) +
  geom_pointrange(mapping = aes(x = time_point, y = betweenesscentrality),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.3) 

bc_plot

network_stats_plot <- ggarrange(dg_plot, cc_plot, bc_plot, nrow= 1, labels = c("A", "B", "C"))

network_stats_plot

# Statistics
#gephi_out_preterms <- gephi_out[gephi_out$Experiment=="Preterms",]
kruskal.test(gephi_out$Degree, g=gephi_out$time_point)
rcompanion::epsilonSquared(gephi_out$Degree, g=gephi_out$time_point, ci=TRUE) 

kruskal.test(gephi_out$betweenesscentrality, g=gephi_out$time_point) 
rcompanion::epsilonSquared(gephi_out$betweenesscentrality, g=gephi_out$time_point, ci=TRUE) 

kruskal.test(gephi_out$closnesscentrality, g=gephi_out$time_point) 
rcompanion::epsilonSquared(gephi_out$closnesscentrality, g=gephi_out$time_point, ci=TRUE)


# export images
pdf("network_stats.pdf", width=12, height=4)
network_stats_plot
dev.off()

