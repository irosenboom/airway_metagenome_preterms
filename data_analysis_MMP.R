# Author: Marie-Madlen Pust
# Last updated: 11 April 2022

# set working directory
setwd("C:/Users/marie/Desktop/Ilona_Rosenboom/R")

# clean global R environment
rm(list = ls())

# required R packages
library(ggplot2)
library(readr)
library(readxl)
library(randomForest)
library(Boruta)
library(plyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(tidylog)
library(RcmdrMisc)
library(Hmisc)
library(tidyverse)

# import and clean data tables
preterm_count_00 <- read_delim("bacteria_per_human_cell_haybaler.csv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
preterm_count_00 <- data.frame(preterm_count_00)
rownames(preterm_count_00) <- preterm_count_00$species
preterm_count_00$species <- NULL
preterm_count_00$rowsum <- NULL
preterm_count_01 <- data.frame(t(preterm_count_00))
preterm_count_01 <- preterm_count_01[order(rownames(preterm_count_01)),]

preterm_meta_00 <- read_excel("metadata_preterm.xlsx", 
                              col_types = c("text", "text", "numeric", "numeric", "numeric", "numeric", 
                                            "text", "text", "text", "text", "text", "numeric", "numeric", 
                                            "text", "text", "text", "text", "text", "text", "text", "text", 
                                            "numeric", "text", "text", "text", "text", "text", "numeric", 
                                            "text", "numeric", "text", "numeric", "text", "text", "text", 
                                            "text", "text"))
preterm_meta_00 <- data.frame(preterm_meta_00)
rownames(preterm_meta_00) <- preterm_meta_00$sample
preterm_meta_00$sample <- NULL
preterm_meta_00 <- preterm_meta_00[order(rownames(preterm_meta_00)),]
preterm_meta_00$Longi_names_0 <- str_sub(rownames(preterm_meta_00), end=-2)

preterm_meta_01 <- preterm_meta_00
preterm_meta_01$RowNAMES <- rownames(preterm_meta_01)
set.seed(1)
preterm_meta_02 <- preterm_meta_01 %>% group_by(Longi_names_0) %>% slice_sample(n = 1)
preterm_meta_03 <- data.frame(preterm_meta_02)
rownames(preterm_meta_03) <- preterm_meta_03$RowNAMES
preterm_meta_03$RowNAMES <- NULL
preterm_meta_03$Longi_names_0 <- NULL
preterm_meta_03 <- preterm_meta_03[order(rownames(preterm_meta_03)),]

final_df_rf <- preterm_count_01
rownames(final_df_rf) <- str_replace(rownames(final_df_rf), "X", "")
final_df_rf <- final_df_rf[rownames(final_df_rf)%in%rownames(preterm_meta_03),]
final_df_rf$group <- factor(preterm_meta_03$group)
final_df_rf$age <- preterm_meta_03$age_sampling
final_df_rf$gender <- factor(preterm_meta_03$gender, labels = c(0,1))
final_df_rf$gestational_age_group <- factor(preterm_meta_03$gestational_age_group, 
                                            levels = c("moderate", "very", "extremely"), labels = c(1,2,3))
final_df_rf$delivery_details <- factor(preterm_meta_03$delivery_details, levels = c(1,2,4))
final_df_rf$group <- factor(preterm_meta_03$group, levels = c("control", "BPD"), labels = c(1,2))
final_df_rf$surfactant <- factor(preterm_meta_03$surfactant, labels = c(0,1,2))
final_df_rf$ventilation <- factor(preterm_meta_03$ventilation, labels = c(0,1,2,4))
final_df_rf$KISS_binary <- factor(preterm_meta_03$KISS_binary, labels = c(0,1))
final_df_rf$previous_antibiotics <- factor(preterm_meta_03$previous_antibiotics, labels = c(0,1))
final_df_rf$neonatal_antibiotics <- factor(preterm_meta_03$neonatal_antibiotics, labels = c(0,1))
final_df_rf$current_antibiotics <- factor(preterm_meta_03$current_antibiotics, labels = c(0,1))
final_df_rf$lung_maturity <- factor(preterm_meta_03$lung_maturity, labels = c(1, 0, 0.5))
final_df_rf$nutrition_MM <- factor(preterm_meta_03$nutrition_MM, labels = c(0,1))
final_df_rf$age_mother <- preterm_meta_03$age_mother
final_df_rf$days_in_hospital <- preterm_meta_03$days_in_hospital



#######################################################################################
#######################################################################################
#######################################################################################
# run first random forest
set.seed(1)
# repeat random forest with 100 different seeds set
random_seeds <- sample(1:200, 100, replace = FALSE) 

imp_final = NULL
error_rate_final = NULL
final_df_rf$group

for (i in random_seeds){
  set.seed(i)
  final_rf <- randomForest(group ~ ., data=final_df_rf, na.action = na.omit, ntree=80, mtry=12, importance=TRUE)
  # store error rate locally
  error_rate_rf <- final_rf$err.rate
  # make data frame of error rate
  error_rate_rf <- as.data.frame(error_rate_rf)
  # add meta data
  error_rate_rf$Seed <- i
  # re-name columns
  colnames(error_rate_rf) <- c('OOB_all', 'non-BPD preterms', 'BPD preterms',  'Seed')
  error_rate_rf <- ddply(error_rate_rf, "Seed", numcolwise(mean))
  # transfer error rate data frame to global environment
  error_rate_final <- rbind(error_rate_final, error_rate_rf)
  
  # confirm random forest with boruta
  # run boruta
  final_df_rf_naFix <- na.roughfix(final_df_rf)
  boruta_rf <- Boruta(group~., data = final_df_rf_naFix, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_rf <- as.data.frame(boruta_rf$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a3 <- importance(final_rf)
  imp_osps_a3 <- data.frame(predictors = rownames(imp_osps_a3), imp_osps_a3)
  
  # Make a data frame with predictor names and their importance
  imp_rf <- importance(final_rf)
  imp_rf <- data.frame(predictors = rownames(imp_rf), imp_rf)
  # add meta data
  imp_rf$Boruta_name <- rownames(boruta_rf)
  imp_rf$Boruta_predict <- boruta_rf$`boruta_rf$finalDecision`
  imp_rf_sub <- subset(imp_rf, MeanDecreaseAccuracy > 0.0)
  # add meta data
  imp_rf_sub$Seed <- i
  # store table globally
  imp_final <- rbind(imp_final, imp_rf_sub)
}

imp_final_conf <- imp_final[imp_final$Boruta_predict=="Confirmed",]
bpd_nbpd_features <- c(table(imp_final_conf$Boruta_name))
bpd_nbpd_features <- bpd_nbpd_features[bpd_nbpd_features>50]
bpd_nbpd_features <- names(bpd_nbpd_features)

imp_final_conf <- imp_final_conf[imp_final_conf$Boruta_name %in%bpd_nbpd_features,]
error_rate_final_L <- gather(error_rate_final, key="OOB", value="value", -c("Seed"))
error_rate_final_L$OOB <- factor(error_rate_final_L$OOB, levels = c("OOB_all", "non-BPD preterms", "BPD preterms"))
error_rate_mad <- ddply(error_rate_final_L, "OOB", numcolwise(mad))
error_rate_med <- ddply(error_rate_final_L, "OOB", numcolwise(median))

imp_final_conf$Type <- "Other"
type_df_bpd_00 <- data.frame(table(imp_final_conf$Seed, imp_final_conf$Type))
type_df_bpd_01 <- ddply(type_df_bpd_00, "Var2", numcolwise(median))
type_df_bpd_01$Per <- round((type_df_bpd_01$Freq / sum(type_df_bpd_01$Freq)) * 100)
type_df_bpd_01$Label_pos <- type_df_bpd_01$Per-.3*type_df_bpd_01$Per+2

pie_plot_bpd <- ggplot(data = type_df_bpd_01, 
                       aes(x = 1, y = Per, fill = Var2))+
  geom_bar(stat = "identity")+
  coord_polar("y", start = 0) +
  geom_text(aes(y = Label_pos, label = paste(Per,"%", sep = "")), col = "white") +
  theme_void() + guides(fill=guide_legend(ncol=1)) +
  theme(legend.title = element_blank(), legend.position = "right") +
  scale_fill_manual(values=c("Other"="black")) +
  xlim(0.1,1.5)

mda_bpd_nbpd <-
  ggplot(imp_final_conf, aes(x=MeanDecreaseAccuracy, y=reorder(Boruta_name, MeanDecreaseAccuracy))) +
  geom_violin() +
  geom_point(size=.6) +
  theme_bw() + xlab("MeanDecreaseAccuracy") + ylab(" ") + theme(panel.grid = element_blank()) +
  scale_y_discrete(labels = c("days_in_hospital"="Days in hospital")) +
  geom_pointrange(mapping = aes(x = MeanDecreaseAccuracy, y = Boruta_name),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.3) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2,4))

mdg_bpd_nbpd <-
  ggplot(imp_final_conf, aes(x=MeanDecreaseGini, y=reorder(Boruta_name, MeanDecreaseGini))) +
  geom_violin() +
  geom_point(size=.6) +
  theme_bw() + xlab("MeanDecreaseGini") + ylab(" ") + theme(panel.grid = element_blank()) +
  scale_y_discrete(labels = c("days_in_hospital"="Days in hospital")) +
  geom_pointrange(mapping = aes(x = MeanDecreaseGini, y = Boruta_name),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median,
                  colour="red", size=0.3) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2,4))

err_bpd_nbpd <-
  ggplot(error_rate_final_L, aes(x=value, y=OOB)) +
  geom_jitter(size=.6) +
  geom_pointrange(mapping = aes(x = value, y = OOB),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median,
                  colour="red", size=0.3) +
  xlim(0,1) +
  theme_bw() + ylab(" ") + xlab("OOB estimate of error rate") + theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 0.5, colour="grey", linetype="dashed")

rf_bpd_nbpd_plots <- ggarrange(mda_bpd_nbpd, mdg_bpd_nbpd, err_bpd_nbpd, pie_plot_bpd, 
                               nrow=2, ncol=2, labels = c("A", "B", "C", "D"),
                               widths = c(1,1,1,1), heights = c(1,1))


#######################################################################################
#######################################################################################
#######################################################################################
# compare pre-term vs. full-term children
# prepare meta data
preterm_fullterm_meta_00 <- read_excel("metadata_preterm_healthycontrols.xlsx")
preterm_fullterm_meta_00 <- data.frame(preterm_fullterm_meta_00)
preterm_fullterm_meta_00$state <- factor(preterm_fullterm_meta_00$state, levels = c("Healthy", "Preterm"), labels = c("Full-term", "Preterm"))
preterm_fullterm_meta_00$group <- factor(preterm_fullterm_meta_00$group, levels = c("Healthy", "Control", "BPD"), 
                                         labels = c("Full-term", "Preterm-nBPD", "Preterm-BPD"))
preterm_fullterm_meta_00$antimicrobial_therapy <- factor(preterm_fullterm_meta_00$antimicrobial_therapy, levels = c("y", "n"), labels = c("Y", "N"))
preterm_fullterm_meta_00$gender <- factor(preterm_fullterm_meta_00$gender, levels = c("m", "f"), labels = c("M", "F"))
rownames(preterm_fullterm_meta_00) <- preterm_fullterm_meta_00$sample_id
preterm_fullterm_meta_00$sample_id <- NULL
preterm_fullterm_meta_00 <- preterm_fullterm_meta_00[order(rownames(preterm_fullterm_meta_00)),]
preterm_fullterm_meta_01 <- preterm_fullterm_meta_00[preterm_fullterm_meta_00$age_in_days>100,]
# remove longitudinal samples
preterm_fullterm_meta_01$Longi_names_0 <- str_sub(rownames(preterm_fullterm_meta_01), end=-2)
preterm_fullterm_meta_01$rowNAMES <- rownames(preterm_fullterm_meta_01)
preterm_fullterm_meta_01$Longi_names_1 <- ifelse(preterm_fullterm_meta_01$state == "Preterm", 
                                                 preterm_fullterm_meta_01$Longi_names, 
                                                 preterm_fullterm_meta_01$rowNAMES)
preterm_fullterm_meta_02 <- preterm_fullterm_meta_01 %>% group_by(Longi_names_1) %>% slice_sample(n = 1)
preterm_fullterm_meta_03 <- data.frame(preterm_fullterm_meta_02)
rownames(preterm_fullterm_meta_03) <- preterm_fullterm_meta_03$rowNAMES
preterm_fullterm_meta_03$Longi_names_1 <- NULL
preterm_fullterm_meta_03$rowNAMES <- NULL
preterm_fullterm_meta_03$Longi_names_0 <- NULL

# prepare count data
preterm_fullterm_count <- read_delim("bphc_preterm_healthy.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
preterm_fullterm_count <- data.frame(preterm_fullterm_count)
rownames(preterm_fullterm_count) <- preterm_fullterm_count$species
preterm_fullterm_count$species <- NULL
preterm_fullterm_count_t <- data.frame(t(preterm_fullterm_count))

rownames(preterm_fullterm_count_t) <- str_replace(rownames(preterm_fullterm_count_t), "X", "")
preterm_fullterm_count_t <- preterm_fullterm_count_t[order(rownames(preterm_fullterm_count_t)),]
preterm_fullterm_count_t <- preterm_fullterm_count_t[rownames(preterm_fullterm_count_t)%in%rownames(preterm_fullterm_meta_03),]

# merge count and meta data
preterm_fullterm_meta_03$Shannon_div <- vegan::diversity(preterm_fullterm_count_t, index = "shannon")
final_df_rf_pre_vs_full <- data.frame(cbind(preterm_fullterm_meta_03, preterm_fullterm_count_t))
final_df_rf_pre_vs_full$group <- NULL
preterm_fullterm_meta_03
# find information on core and rare species
preterm_fullterm_count_02 <- preterm_fullterm_count
preterm_fullterm_count_02$rowSums_00 <- rowSums(preterm_fullterm_count_02)
preterm_fullterm_count_02$rowSums_01 <- preterm_fullterm_count_02$rowSums_00 / sum(preterm_fullterm_count_02$rowSums_00)

preterm_fullterm_count_03 <- preterm_fullterm_count_02[order(preterm_fullterm_count_02$rowSums_01, decreasing = TRUE),]
preterm_fullterm_count_03$cumsum <- cumsum(preterm_fullterm_count_03$rowSums_01)
preterm_fullterm_count_03$species_type <- ifelse(preterm_fullterm_count_03$cumsum < 0.95, "Core species", "Rare species")
preterm_fullterm_core <- rownames(subset(preterm_fullterm_count_03, species_type=="Core species"))
preterm_fullterm_rare <- rownames(subset(preterm_fullterm_count_03, species_type=="Rare species"))


# run second random forest (preterm vs full term)
imp_final_pre_vs_full = NULL
error_rate_final_pre_vs_full = NULL

for (i in random_seeds){
  set.seed(i)
  final_rf <- randomForest(state ~ ., data=final_df_rf_pre_vs_full, na.action = na.omit, ntree=80, mtry=12, importance=TRUE)
  # store error rate locally
  error_rate_rf <- final_rf$err.rate
  # make data frame of error rate
  error_rate_rf <- as.data.frame(error_rate_rf)
  # add meta data
  error_rate_rf$Seed <- i
  # re-name columns
  colnames(error_rate_rf) <- c('OOB_all', 'Full-term infants', 'Pre-term infants',  'Seed')
  error_rate_rf <- ddply(error_rate_rf, "Seed", numcolwise(mean))
  # transfer error rate data frame to global environment
  error_rate_final_pre_vs_full <- rbind(error_rate_final_pre_vs_full, error_rate_rf)
  
  # confirm random forest with boruta
  # run boruta
  final_df_rf_naFix <- na.roughfix(final_df_rf_pre_vs_full)
  boruta_rf <- Boruta(state~., data = final_df_rf_naFix, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_rf <- as.data.frame(boruta_rf$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a3 <- importance(final_rf)
  imp_osps_a3 <- data.frame(predictors = rownames(imp_osps_a3), imp_osps_a3)
  
  # Make a data frame with predictor names and their importance
  imp_rf <- importance(final_rf)
  imp_rf <- data.frame(predictors = rownames(imp_rf), imp_rf)
  # add meta data
  imp_rf$Boruta_name <- rownames(boruta_rf)
  imp_rf$Boruta_predict <- boruta_rf$`boruta_rf$finalDecision`
  imp_rf_sub <- subset(imp_rf, MeanDecreaseAccuracy > 0.0)
  # add meta data
  imp_rf_sub$Seed <- i
  # store table globally
  imp_final_pre_vs_full <- rbind(imp_final_pre_vs_full, imp_rf_sub)
}

imp_final_conf_pre_vs_full <- imp_final_pre_vs_full[imp_final_pre_vs_full$Boruta_predict=="Confirmed",]
preterm_fullterm_features <- c(table(imp_final_conf_pre_vs_full$Boruta_name))
preterm_fullterm_features <- preterm_fullterm_features[preterm_fullterm_features>50]
preterm_fullterm_features <- names(preterm_fullterm_features)

imp_final_conf_pre_vs_full <- imp_final_conf_pre_vs_full[imp_final_conf_pre_vs_full$Boruta_name %in%preterm_fullterm_features,]
imp_final_conf_pre_vs_full$Type <- ifelse(imp_final_conf_pre_vs_full$Boruta_name %in% preterm_fullterm_core, "High-abundance taxa",
                                          ifelse(imp_final_conf_pre_vs_full$Boruta_name %in% preterm_fullterm_rare, "Low-abundance taxa", "Other"))

error_rate_final_L_pre_vs_full <- gather(error_rate_final_pre_vs_full, key="OOB", value="value", -c("Seed"))
error_rate_final_L_pre_vs_full$OOB <- factor(error_rate_final_L_pre_vs_full$OOB, levels = c("OOB_all", "Full-term infants", "Pre-term infants"))
error_rate_mad_preterm_fullterm <- ddply(error_rate_final_L_pre_vs_full, "OOB", numcolwise(mad))
error_rate_med_preterm_fullterm <- ddply(error_rate_final_L_pre_vs_full, "OOB", numcolwise(median))

mda_full_preterm <-
  ggplot(imp_final_conf_pre_vs_full, aes(x=MeanDecreaseAccuracy, y=Boruta_name)) +
  geom_violin() +
  geom_point(size=.6) +
  theme_bw() + xlab("MeanDecreaseAccuracy") + ylab(" ") + theme(panel.grid = element_blank()) +
  scale_y_discrete(labels = c("Campylobacter_concisus"=expression(italic("Campylobacter concisus")),
                              "Fusobacterium_periodonticum"=expression(italic("Fusobacterium periodonticum")),
                              "Fusobacterium_nucleatum"=expression(italic("Fusobacterium nucleatum")),
                              "Gemella_sanguinis"=expression(italic("Gemella sanguinis")),
                              "Haemophilus_parainfluenzae"=expression(italic("Haemophilus parainfluenzae")),
                              "Neisseria_subflava"=expression(italic("Neisseria subflava")),
                              "Streptococcus_agalactiae"=expression(italic("Streptococcus agalactiae")),
                              "Streptococcus_anginosus"=expression(italic("Streptococcus anginosus")),
                              "Streptococcus_constellatus"=expression(italic("Streptococcus constellatus")),
                              "Streptococcus_cristatus"=expression(italic("Streptococcus cristatus")),
                              "Streptococcus_himalayensis"=expression(italic("Streptococcus himalayensis")),
                              "Streptococcus_iniae"=expression(italic("Streptococcus iniae")),
                              "Streptococcus_intermedius"=expression(italic("Streptococcus intermedius")),
                              "Streptococcus_marmotae"=expression(italic("Streptococcus marmotae")),
                              "Streptococcus_infantarius"=expression(italic("Streptococcus infantarius")),
                              "Streptococcus_oralis"=expression(italic("Streptococcus oralis")),
                              "Streptococcus_pneumoniae"=expression(italic("Streptococcus pneumoniae")),
                              "Streptococcus_pseudopneumoniae"=expression(italic("Streptococcus pseudopneumoniae")),
                              "Streptococcus_suis"=expression(italic("Streptococcus suis")),
                              "Streptococcus_thermophilus"=expression(italic("Streptococcus thermophilus")),
                              "Rothia_dentocariosa"=expression(italic("Rothia dentocariosa")),
                              "antimicrobial_therapy"="Antimicrobial therapy")) +
  geom_pointrange(mapping = aes(x = MeanDecreaseAccuracy, y = Boruta_name),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.3) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2,4))

mdg_full_preterm <-
  ggplot(imp_final_conf_pre_vs_full, aes(x=MeanDecreaseGini, y=Boruta_name)) +
  geom_violin() +
  geom_point(size=.6) +
  theme_bw() + xlab("MeanDecreaseGini") + ylab(" ") + theme(panel.grid = element_blank()) +
  scale_y_discrete(labels = c("Campylobacter_concisus"=expression(italic("Campylobacter concisus")),
                              "Fusobacterium_periodonticum"=expression(italic("Fusobacterium periodonticum")),
                              "Fusobacterium_nucleatum"=expression(italic("Fusobacterium nucleatum")),
                              "Gemella_sanguinis"=expression(italic("Gemella sanguinis")),
                              "Haemophilus_parainfluenzae"=expression(italic("Haemophilus parainfluenzae")),
                              "Neisseria_subflava"=expression(italic("Neisseria subflava")),
                              "Streptococcus_agalactiae"=expression(italic("Streptococcus agalactiae")),
                              "Streptococcus_anginosus"=expression(italic("Streptococcus anginosus")),
                              "Streptococcus_constellatus"=expression(italic("Streptococcus constellatus")),
                              "Streptococcus_cristatus"=expression(italic("Streptococcus cristatus")),
                              "Streptococcus_himalayensis"=expression(italic("Streptococcus himalayensis")),
                              "Streptococcus_iniae"=expression(italic("Streptococcus iniae")),
                              "Streptococcus_intermedius"=expression(italic("Streptococcus intermedius")),
                              "Streptococcus_marmotae"=expression(italic("Streptococcus marmotae")),
                              "Streptococcus_infantarius"=expression(italic("Streptococcus infantarius")),
                              "Streptococcus_oralis"=expression(italic("Streptococcus oralis")),
                              "Streptococcus_pneumoniae"=expression(italic("Streptococcus pneumoniae")),
                              "Streptococcus_pseudopneumoniae"=expression(italic("Streptococcus pseudopneumoniae")),
                              "Streptococcus_suis"=expression(italic("Streptococcus suis")),
                              "Streptococcus_thermophilus"=expression(italic("Streptococcus thermophilus")),
                              "Rothia_dentocariosa"=expression(italic("Rothia dentocariosa")),
                              "antimicrobial_therapy"="Antimicrobial therapy")) +
  geom_pointrange(mapping = aes(x = MeanDecreaseGini, y = Boruta_name),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median,
                  colour="red", size=0.3) +
  scale_x_continuous(limits = c(0,5), breaks = c(0,2,4))

err_full_preterm <-
  ggplot(error_rate_final_L_pre_vs_full, aes(x=value, y=OOB)) +
  geom_jitter(size=.6) +
  geom_pointrange(mapping = aes(x = value, y = OOB),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, 
                  fun = median,
                  colour="red", size=0.3) +
  xlim(0,1) +
  theme_bw() + ylab(" ") + xlab("OOB estimate of error rate") + theme(panel.grid = element_blank()) +
  scale_y_discrete(labels=c("Pre-term infants"="    Preterm infants",
                            "Full-term infants"="    Full-term infants")) +
  geom_vline(xintercept = 0.5, colour="grey", linetype="dashed")

type_df_00 <- data.frame(table(imp_final_conf_pre_vs_full$Seed, imp_final_conf_pre_vs_full$Type))
type_df_01 <- ddply(type_df_00, "Var2", numcolwise(median))
type_df_01$Per <- round((type_df_01$Freq / sum(type_df_01$Freq)) * 100)
type_df_01$Label_pos <- type_df_01$Per-.3*type_df_01$Per+2

pie_plot <- ggplot(data = type_df_01, 
                   aes(x = 1, y = Per, fill = Var2))+
  geom_bar(stat = "identity")+
  coord_polar("y", start = 0) +
  geom_text(aes(y = Label_pos, label = paste(Per,"%", sep = "")), col = "white") +
  theme_void() + guides(fill=guide_legend(ncol=1)) +
  theme(legend.title = element_blank(), legend.position = "right") +
  scale_fill_manual(values=c("Other"="black", "High-abundance taxa"="goldenrod1", "Low-abundance taxa"="green4")) +
  xlim(0.1,1.5)

rf_full_preterm_plots <- ggarrange(mda_full_preterm, mdg_full_preterm, err_full_preterm, pie_plot, 
                                   nrow=2, ncol=2, labels = c("A", "B", "C", "D"),
                                   widths = c(1,1,1,1), heights = c(1,0.7))


# export pdf image (random forest)
pdf("random_forest_primal_01.pdf", width=10, height = 8)
rf_full_preterm_plots
dev.off()

# export pdf image (random forest for supplementary materials)
pdf("random_forest_primal_supp_01.pdf", width=8, height = 6)
rf_bpd_nbpd_plots
dev.off()

