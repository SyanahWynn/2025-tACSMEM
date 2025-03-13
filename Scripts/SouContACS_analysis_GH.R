#### LIBRARIES ####
library(tidyverse)
library(lme4)
library(glmmTMB)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)
library(reshape2)
library(effsize)
library(emmeans)
library(performance)
library(see)
library(png)
library(patchwork)
library(DHARMa)
library(olsrr)
library(conflicted)

#### VARIABLES ####
rm(list=ls()) # clear the global environment
dirroot <- dirname(rstudioapi::getSourceEditorContext()$path)
data_path_tbt <- paste(dirroot,'/Data_TbT.csv', sep="") # trial by trial
data_path_tbc <- paste(dirroot,'/Data_TbC.csv', sep="") # by condition
data_path_mani <- paste(dirroot,'/ManiCheck.csv', sep="") # stimulation check
outp_path <- paste(dirroot,'/Output/', sep="")
brpl_path <- paste(dirroot,'/BrainPlots/', sep="")
outl_thres <- 3 #SD
cp1 <- c('royalblue','turquoise','firebrick','plum')
cp2 <- c('forestgreen','purple3')
do_stats_print <- FALSE # to print or not print the stats/behavioral data
do_assum_check <- TRUE
do_mani_check <- FALSE
the_stim_freq <- 4
gam_stim_freq <- 50

#### FUNCTIONS ####
stndz <- function(x, na.rm = TRUE, outl_filt = 'None') {
  if (outl_filt[1]=='None') {
    (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
  } else {
    (x - mean(x[outl_filt==FALSE], na.rm = na.rm)) / sd(x[outl_filt==FALSE], na.rm)
  }
}
#### PREPARE THE DATA: ADJUST & CALCULATE NEW VARIABLES ####
# Read in the data
data_mat_tbt <- read_csv(data_path_tbt, show_col_types = FALSE)
data_mat_tbc <- read_csv(data_path_tbc, show_col_types = FALSE)
# Replace NaNs with NA since NaN means "Not a Number" and NA is "Not Available"
data_mat_tbt[sapply(data_mat_tbt, is.nan)] <- NA
data_mat_tbt$Source[data_mat_tbt$Source=="N/A"] <- NA
data_mat_tbc[sapply(data_mat_tbc, is.nan)] <- NA
data_mat_tbc$SourceAcc[data_mat_tbc$SourceAcc=="NaN"] <- NA
data_mat_tbc$SourceConf[data_mat_tbc$SourceConf=="NaN"] <- NA
# remove the rows that represent combinations that are not in the actual data
data_mat_tbc <- data_mat_tbc[complete.cases(data_mat_tbc$label),]
# add additional labels
data_mat_tbc <- mutate(data_mat_tbc, label_itemmem = NA)
data_mat_tbc <- mutate(data_mat_tbc, label_itemconf = NA)
data_mat_tbc <- mutate(data_mat_tbc, label_itemconfacc = NA)
data_mat_tbc <- mutate(data_mat_tbc, label_sourcemem = NA)
data_mat_tbc <- mutate(data_mat_tbc, label_sourceconf = NA)
data_mat_tbc$label_itemmem[select(data_mat_tbc,OldNew)=='Old' & select(data_mat_tbc,ItemAcc)=='iCor'] <- 'iHit'
data_mat_tbc$label_itemmem[select(data_mat_tbc,OldNew)=='Old' & select(data_mat_tbc,ItemAcc)=='iInc'] <- 'iMiss'
data_mat_tbc$label_itemmem[select(data_mat_tbc,OldNew)=='New' & select(data_mat_tbc,ItemAcc)=='iCor'] <- 'CR'
data_mat_tbc$label_itemmem[select(data_mat_tbc,OldNew)=='New' & select(data_mat_tbc,ItemAcc)=='iInc'] <- 'FA'
data_mat_tbc$label_itemconf[select(data_mat_tbc,ItemConf)=='iHC'] <- 'iHC'
data_mat_tbc$label_itemconf[select(data_mat_tbc,ItemConf)=='iLC'] <- 'iLC'
data_mat_tbc$label_itemconfacc[select(data_mat_tbc,ItemConf)=='iHC' & select(data_mat_tbc,ItemAcc)=='iCor'] <- 'iHCCor'
data_mat_tbc$label_sourcemem[select(data_mat_tbc,OldNew)=='Old' & select(data_mat_tbc,SourceAcc)=='sCor'] <- 'sHit'
data_mat_tbc$label_sourcemem[select(data_mat_tbc,OldNew)=='Old' & select(data_mat_tbc,SourceAcc)=='sInc'] <- 'sMiss'
data_mat_tbc$label_sourceconf[select(data_mat_tbc,SourceConf)=='sHC'] <- 'sHC'
data_mat_tbc$label_sourceconf[select(data_mat_tbc,SourceConf)=='sLC'] <- 'sLC'
# calculate the accuracy of the memory responses
data_mat_tbt <- mutate(data_mat_tbt, ONAcc = NA)
data_mat_tbt$ONAcc[select(data_mat_tbt,OldNew)=='Old' & select(data_mat_tbt,ONResp)<3] <- 1 # hits
data_mat_tbt$ONAcc[select(data_mat_tbt,OldNew)=='New' & select(data_mat_tbt,ONResp)>3] <- 1 # correct rejections
data_mat_tbt$ONAcc[select(data_mat_tbt,OldNew)=='Old' & select(data_mat_tbt,ONResp)>3] <- 0 # misses
data_mat_tbt$ONAcc[select(data_mat_tbt,OldNew)=='New' & select(data_mat_tbt,ONResp)<3] <- 0 # false alarms
data_mat_tbt$ONAcc[select(data_mat_tbt,ONResp)==3] = 0 # guesses, here classified as incorrect
data_mat_tbt <- mutate(data_mat_tbt, SourceAcc = NA)
data_mat_tbt$SourceAcc[select(data_mat_tbt,Source)=='Pleasant' & select(data_mat_tbt,SourceResp)<3] <- 1 # hits
data_mat_tbt$SourceAcc[select(data_mat_tbt,Source)=='Place' & select(data_mat_tbt,SourceResp)>3] <- 1 # hits
data_mat_tbt$SourceAcc[select(data_mat_tbt,Source)=='Pleasant' & select(data_mat_tbt,SourceResp)>3] <- 0 # misses
data_mat_tbt$SourceAcc[select(data_mat_tbt,Source)=='Place' & select(data_mat_tbt,SourceResp)<3] <- 0 # misses
data_mat_tbt$SourceAcc[(select(data_mat_tbt,Source)=='Pleasant' | select(data_mat_tbt,Source)=='Place') & (select(data_mat_tbt,SourceResp)==3)] <- 0 # guesses, here classified as incorrect
# calculate the confidence of the memory responses
data_mat_tbt <- mutate(data_mat_tbt, ONConf = NA)
data_mat_tbt$ONConf[select(data_mat_tbt,ONResp)==5 | select(data_mat_tbt,ONResp)==1] <- 1 # high-confident
data_mat_tbt$ONConf[select(data_mat_tbt,ONResp)<5 & select(data_mat_tbt,ONResp)>1] <- 0 # low-confident
data_mat_tbt <- mutate(data_mat_tbt, SourceConf = NA)
data_mat_tbt$SourceConf[select(data_mat_tbt,SourceResp)==5 | select(data_mat_tbt,SourceResp)==1] <- 1 # high-confident
data_mat_tbt$SourceConf[select(data_mat_tbt,SourceResp)<5 & select(data_mat_tbt,SourceResp)>1] <- 0 # low-confident
# Change the categorical variables into factors..
cur_cat_vars <- c('Subject','Sex','OldNew','Source','ONAcc','SourceAcc','ONConf','SourceConf')
data_mat_tbt[,cur_cat_vars] <- lapply(data_mat_tbt[,cur_cat_vars], factor)
cur_cat_vars <- c('Subject','OldNew','ItemAcc','SourceAcc','ItemConf','SourceConf', 'label', 'label_itemmem','label_itemconf','label_sourcemem','label_sourceconf')
data_mat_tbc[,cur_cat_vars] <- lapply(data_mat_tbc[,cur_cat_vars], factor)
# .. Or ordered factors
cur_cat_vars <- c('Session','EncResp','ONResp','SourceResp')
data_mat_tbt[,cur_cat_vars] <- lapply(data_mat_tbt[,cur_cat_vars], ordered)
# Adjust the order of the Stim factor levels
data_mat_tbt$Stimulation <- factor(data_mat_tbt$Stimulation, levels = c('None','Sham','Gamma','Theta'))
# Change the "pek" (peak amplitude) to the absolute deviation of the stimulation frequency: 4 & 50
data_mat_tbc$d_Enc_pek_theta_frontal  = abs(data_mat_tbc$Enc_pek_theta_frontal - the_stim_freq)
data_mat_tbc$d_Enc_pek_theta_parietal = abs(data_mat_tbc$Enc_pek_theta_parietal - the_stim_freq)
data_mat_tbc$d_Ret_pek_theta_frontal  = abs(data_mat_tbc$Ret_pek_theta_frontal - the_stim_freq)
data_mat_tbc$d_Ret_pek_theta_parietal = abs(data_mat_tbc$Ret_pek_theta_parietal - the_stim_freq)
data_mat_tbc$d_Enc_pek_gamma_frontal  = abs(data_mat_tbc$Enc_pek_gamma_frontal - gam_stim_freq)
data_mat_tbc$d_Enc_pek_gamma_parietal = abs(data_mat_tbc$Enc_pek_gamma_parietal - gam_stim_freq)
data_mat_tbc$d_Ret_pek_gamma_frontal  = abs(data_mat_tbc$Ret_pek_gamma_frontal - gam_stim_freq)
data_mat_tbc$d_Ret_pek_gamma_parietal = abs(data_mat_tbc$Ret_pek_gamma_parietal - gam_stim_freq)
# Making a filtering variable for the EEG outliers (ret pow) (threshold SD > outl_thres)
pow_thet_vars <- c('Ret_ThetaFro', 'Ret_ThetaPar')
pow_gamm_vars <- c('Ret_GammaFro', 'Ret_GammaPar')
pac_enc_vars  <- c('Enc_frontalPAC', 'Enc_parietalPAC')
pac_ret_vars  <- c('Ret_frontalPAC', 'Ret_parietalPAC')
data_mat_tbt <- mutate(data_mat_tbt, outlrs_powtf = FALSE)
data_mat_tbt <- mutate(data_mat_tbt, outlrs_powtp = FALSE)
data_mat_tbt <- mutate(data_mat_tbt, outlrs_powgf = FALSE)
data_mat_tbt <- mutate(data_mat_tbt, outlrs_powgp = FALSE)
z_scores <- mutate_at(data_mat_tbt, c('Ret_ThetaFro', 'Ret_ThetaPar', 'Ret_GammaFro', 'Ret_GammaPar'), stndz, outl_filt = 'None')
data_mat_tbt$outlrs_powtf[which(rowSums(select(z_scores,'Ret_ThetaFro')>outl_thres 
                                   | select(z_scores,'Ret_ThetaFro')<(-outl_thres),na.rm = TRUE)>=1)] = TRUE
data_mat_tbt$outlrs_powtp[which(rowSums(select(z_scores,'Ret_ThetaPar')>outl_thres 
                                   | select(z_scores,'Ret_ThetaPar')<(-outl_thres),na.rm = TRUE)>=1)] = TRUE
data_mat_tbt$outlrs_powgf[which(rowSums(select(z_scores,'Ret_GammaFro')>outl_thres 
                                   | select(z_scores,'Ret_GammaFro')<(-outl_thres),na.rm = TRUE)>=1)] = TRUE
data_mat_tbt$outlrs_powgp[which(rowSums(select(z_scores,'Ret_GammaPar')>outl_thres 
                                   | select(z_scores,'Ret_GammaPar')<(-outl_thres),na.rm = TRUE)>=1)] = TRUE
rm(z_scores)

#### PREPARE THE DATA: AGGREGATED DATA ####
# build the basic structure
data_agr = unique(data_mat_tbt[ , c('Subject','Session','Stimulation')])
vars <- c('Age','Sex','EncResp','EncRespPla','EncRespPle','EncRT','EncRTPla','EncRTPle',
    'ItemHitRate','ItemHitRatePla','ItemHitRatePle','ItemHitRateLC','ItemHitRateHC',
    'ItemFARate','ItemCRRate','ItemCRRateLC','ItemCRRateHC',
    'ItemDprime','ItemDprimePla','ItemDprimePle',
    'ItemHCRate','ItemHCHitRate','ItemHCHitRatePla','ItemHCHitRatePle','ItemHCCRRate',
    'ItemRT','ItemRTPla','ItemRTPle','ItemRTCor','ItemRTCorPla','ItemRTCorPle','ItemRTIncor','ItemRTIncorPla','ItemRTIncorPle',
    'SourceHitRate','SourceHitRateLC','SourceHitRateHC','SourceFARate','SourceDprime',
    'SourceHCRate','SourceHCHitRate','SourceHCHitRatePla','SourceHCHitRatePle',
    'SourceRT','SourceRTPla','SourceRTPle','SourceRTCor','SourceRTCorPla','SourceRTCorPle','SourceRTIncor','SourceRTIncorPla','SourceRTIncorPle',
    'Ret_ThetaFro', 'Ret_GammaFro', 'Ret_ThetaPar', 'Ret_GammaPar'
)
data_agr[, vars] <- NA
# get the data
cat('\n\nAGGREGATE DATA')
for (p in unique(data_mat_tbt$Subject)) {
  cat('\n\nPARTICIPANT: ',p)
  for (s in unique(data_mat_tbt$Session)) {
    # get the data of this participant and session
    cur_data <- dplyr::filter(data_mat_tbt,Subject==p,Session==s)
    data_agr$Age[data_agr$Subject==p & data_agr$Session==s] <- cur_data$Age[1]
    data_agr$Sex[data_agr$Subject==p & data_agr$Session==s] <- cur_data$Sex[1]
    # ENCODING
    data_agr$EncResp[data_agr$Subject==p & data_agr$Session==s] <- mean(as.numeric(levels(cur_data$EncResp)[cur_data$EncResp]), na.rm = TRUE)/3
    data_agr$EncRespPla[data_agr$Subject==p & data_agr$Session==s] <- with(dplyr::filter(cur_data, Source == 'Place'), mean(as.numeric(as.character(EncResp)), na.rm = TRUE) / 3)
    data_agr$EncRespPle[data_agr$Subject==p & data_agr$Session==s] <- with(dplyr::filter(cur_data, Source == 'Pleasant'), mean(as.numeric(as.character(EncResp)), na.rm = TRUE) / 3)
    data_agr$EncRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$EncRT, na.rm = TRUE)
    data_agr$EncRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Place')$EncRT, na.rm = TRUE)
    data_agr$EncRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Pleasant')$EncRT, na.rm = TRUE)
    # RETRIEVAL
    # hit rate & fa rate: item
    data_agr$ItemHitRate[data_agr$Subject==p & data_agr$Session==s] <-
        with(dplyr::filter(cur_data, OldNew == 'Old'), sum(as.numeric(as.character(ONResp)) < 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # (HitHC + HitLC) / old items
    data_agr$ItemHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Place'), sum(as.numeric(as.character(ONResp)) < 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # (HitHC + HitLC) / old place items
    data_agr$ItemHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Pleasant'), sum(as.numeric(as.character(ONResp)) < 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # (HitHC + HitLC) / old pleasant items
    data_agr$ItemHitRateLC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'Old', ONConf == 0), sum(as.numeric(as.character(ONResp)) == 2, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # HitLC / LC old items
    data_agr$ItemHitRateHC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'Old', ONConf == 1), sum(as.numeric(as.character(ONResp)) == 1, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # HitHC / HC old items
    data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'New'), sum(as.numeric(as.character(ONResp)) < 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # (FAHC + FALC) / new items
    data_agr$ItemCRRate[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'New'), sum(as.numeric(as.character(ONResp)) > 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # (CRHC + CRLC) / new items
    data_agr$ItemCRRateLC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'New', ONConf == 0), sum(as.numeric(as.character(ONResp)) == 4, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # CRLC / new items
    data_agr$ItemCRRateHC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, OldNew == 'New', ONConf == 1), sum(as.numeric(as.character(ONResp)) == 5, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(ONResp))))) # CRHC / new items
    # hit rate & fa rate: source
    data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Place'), sum(as.numeric(as.character(SourceResp)) > 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(SourceResp))))) # Place Hits / Place items
    data_agr$SourceHitRateLC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Place', SourceConf == 0), sum(as.numeric(as.character(SourceResp)) == 4, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(SourceResp))))) # Place LC Hits / Place items
    data_agr$SourceHitRateHC[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Place', SourceConf == 1), sum(as.numeric(as.character(SourceResp)) == 5, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(SourceResp))))) # Place HC Hits / place items
    data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] <- 
        with(dplyr::filter(cur_data, Source == 'Pleasant'), sum(as.numeric(as.character(SourceResp)) > 3, na.rm = TRUE) / 
        sum(!is.na(as.numeric(as.character(SourceResp))))) # Pleasant Misses / Pleasant items
    # d-prime: item
    # adjust the hit and fa rate if it is 1 or 0 to make d' not inf
    if (data_agr$ItemHitRate[data_agr$Subject==p & data_agr$Session==s] == 1) {data_agr$ItemHitRate[data_agr$Subject==p & data_agr$Session==s] <- .999}
    if (data_agr$ItemHitRatePla[data_agr$Subject==p & data_agr$Session==s] == 1) {data_agr$ItemHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- .999}
    if (data_agr$ItemHitRatePle[data_agr$Subject==p & data_agr$Session==s] == 1) {data_agr$ItemHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- .999}
    if (data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s] == 0) {data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s] <- .001}
    data_agr$ItemDprime[data_agr$Subject==p & data_agr$Session==s] <- qnorm(data_agr$ItemHitRate[data_agr$Subject==p & data_agr$Session==s])-
      qnorm(data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s])
    data_agr$ItemDprimePla[data_agr$Subject==p & data_agr$Session==s] <- qnorm(data_agr$ItemHitRatePla[data_agr$Subject==p & data_agr$Session==s])-
      qnorm(data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s])
    data_agr$ItemDprimePle[data_agr$Subject==p & data_agr$Session==s] <- qnorm(data_agr$ItemHitRatePle[data_agr$Subject==p & data_agr$Session==s])-
      qnorm(data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s])
    # d-prime: Source
    # adjust the hit and fa rate if it is 1 or 0 to make d' not inf
    if (data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] == 1) {data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] <- .999}
    if (data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] == 0) {data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] <- .001}
    data_agr$SourceDprime[data_agr$Subject==p & data_agr$Session==s] <- qnorm(data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s])-
      qnorm(data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s])
    # confidence: item
    data_agr$ItemHCRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(cur_data$ONResp)[cur_data$ONResp])==1 | as.numeric(levels(cur_data$ONResp)[cur_data$ONResp])==5, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(cur_data$ONResp)[cur_data$ONResp]))) # HC  / all responses
    data_agr$ItemHCHitRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old')$ONResp)[dplyr::filter(cur_data,OldNew=='Old')$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old')$ONResp)[dplyr::filter(cur_data,OldNew=='Old')$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old' & Source=='Place')$ONResp)[dplyr::filter(cur_data,OldNew=='Old' & Source=='Place')$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old' & Source=='Place')$ONResp)[dplyr::filter(cur_data,OldNew=='Old' & Source=='Place')$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old' & Source=='Pleasant')$ONResp)[dplyr::filter(cur_data,OldNew=='Old' & Source=='Pleasant')$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='Old' & Source=='Pleasant')$ONResp)[dplyr::filter(cur_data,OldNew=='Old' & Source=='Pleasant')$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCCRRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='New')$ONResp)[dplyr::filter(cur_data,OldNew=='New')$ONResp])==5, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,OldNew=='New')$ONResp)[dplyr::filter(cur_data,OldNew=='New')$ONResp])>3, na.rm = TRUE) # HCCRs  / CRs
    # confidence: source
    data_agr$SourceHCRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp])==1 | as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp])==5, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp]))) # HC  / all responses
    data_agr$SourceHCHitRate[data_agr$Subject==p & data_agr$Session==s] <- (sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Pleasant')$SourceResp)[dplyr::filter(cur_data,Source=='Pleasant')$SourceResp])==1, na.rm = TRUE) + sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Place')$SourceResp)[dplyr::filter(cur_data,Source=='Place')$SourceResp])==5, na.rm = TRUE))/
      (sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Pleasant')$SourceResp)[dplyr::filter(cur_data,Source=='Pleasant')$SourceResp])<3, na.rm = TRUE) + sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Place')$SourceResp)[dplyr::filter(cur_data,Source=='Place')$SourceResp])>3, na.rm = TRUE)) # HCHits  / Hits
    data_agr$SourceHCHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Place')$SourceResp)[dplyr::filter(cur_data,Source=='Place')$SourceResp])==5, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Place')$SourceResp)[dplyr::filter(cur_data,Source=='Place')$SourceResp])>3, na.rm = TRUE) # HCHits  / Hits
    data_agr$SourceHCHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Pleasant')$SourceResp)[dplyr::filter(cur_data,Source=='Pleasant')$SourceResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(dplyr::filter(cur_data,Source=='Pleasant')$SourceResp)[dplyr::filter(cur_data,Source=='Pleasant')$SourceResp])<3, na.rm = TRUE) # HCHits  / Hits
    # reaction time: item
    data_agr$ItemRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$ONRT, na.rm = TRUE)
    data_agr$ItemRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Place')$ONRT, na.rm = TRUE)
    data_agr$ItemRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Pleasant')$ONRT, na.rm = TRUE)
    data_agr$ItemRTCor[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==1)$ONRT, na.rm = TRUE)
    data_agr$ItemRTCorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==1 & Source=='Place')$ONRT, na.rm = TRUE)
    data_agr$ItemRTCorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==1 & Source=='Pleasant')$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncor[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==0)$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==0 & Source=='Place')$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,ONAcc==0 & Source=='Pleasant')$ONRT, na.rm = TRUE)
    # reaction time: Source
    data_agr$SourceRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$SourceRT, na.rm = TRUE)
    data_agr$SourceRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Place')$SourceRT, na.rm = TRUE)
    data_agr$SourceRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,Source=='Pleasant')$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCor[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==1)$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==1 & Source=='Place')$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==1 & Source=='Pleasant')$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncor[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==0)$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==0 & Source=='Place')$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,SourceAcc==0 & Source=='Pleasant')$SourceRT, na.rm = TRUE)
    # EEG
    data_agr$Ret_ThetaFro[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,outlrs_powtf!=TRUE)$Ret_ThetaFro, na.rm = TRUE)
    data_agr$Ret_GammaFro[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,outlrs_powgf!=TRUE)$Ret_GammaFro, na.rm = TRUE)
    data_agr$Ret_ThetaPar[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,outlrs_powtp!=TRUE)$Ret_ThetaPar, na.rm = TRUE)
    data_agr$Ret_GammaPar[data_agr$Subject==p & data_agr$Session==s] <- mean(dplyr::filter(cur_data,outlrs_powgp!=TRUE)$Ret_GammaPar, na.rm = TRUE)
  }
}
# adjust the string variable(s)
data_agr$Sex[data_agr$Sex==1] <- levels(cur_data$Sex)[1]
data_agr$Sex[data_agr$Sex==2] <- levels(cur_data$Sex)[2]
# PRINT THE DESCRIPTIVES PER SESSION
if (do_stats_print) {
  for (s in unique(data_mat_tbt$Session)) {
    cat('\n\nXXXXXXXXXXXXXXXXXX\nXX SESSION: ',s,' XX\nXXXXXXXXXXXXXXXXXX\n')
    cat('\nEncoding success (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncResp)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncResp)*100,2))
    cat('\nEncoding success (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncRespPla)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncRespPla)*100,2))
    cat('\nEncoding success (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncRespPle)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncRespPle)*100,2))
    x <- dplyr::filter(data_agr,Session==s)$EncRespPla
    y <- dplyr::filter(data_agr,Session==s)$EncRespPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nEncoding RT (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncRT),0))
    cat('\nEncoding RT (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncRTPla),0))
    cat('\nEncoding RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$EncRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$EncRTPle),0))
    x <- dplyr::filter(data_agr,Session==s)$EncRTPla
    y <- dplyr::filter(data_agr,Session==s)$EncRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Item Hitrate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHitRate),3))
    cat('\nRetrieval Item Hitrate (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHitRatePla),3))
    cat('\nRetrieval Item Hitrate (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHitRatePle),3))
    x <- dplyr::filter(data_agr,Session==s)$ItemHitRatePla
    y <- dplyr::filter(data_agr,Session==s)$ItemHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item FArate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemFARate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemFARate),3))
    cat('\nRetrieval Item d-prime (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemDprime),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemDprime),3))
    cat('\nRetrieval Item d-prime (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemDprimePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemDprimePla),3))
    cat('\nRetrieval Item d-prime (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemDprimePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemDprimePle),3))
    x <- dplyr::filter(data_agr,Session==s)$ItemDprimePla
    y <- dplyr::filter(data_agr,Session==s)$ItemDprimePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Item HC rate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHCRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHCRate),3))
    cat('\nRetrieval Item HC rate (Hits): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHCHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHCHitRate),3))
    cat('\nRetrieval Item HC rate (Hits, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHCHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHCHitRatePla),3))
    cat('\nRetrieval Item HC rate (Hits, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHCHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHCHitRatePle),3))
    x <- dplyr::filter(data_agr,Session==s)$ItemHCHitRatePla
    y <- dplyr::filter(data_agr,Session==s)$ItemHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item HC rate (CRs): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemHCCRRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemHCCRRate),3))
    
    cat('\n\nRetrieval Item RT (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRT),0))
    cat('\nRetrieval Item RT (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTPla),0))
    cat('\nRetrieval Item RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTPle),0))
    x <- dplyr::filter(data_agr,Session==s)$ItemRTPla
    y <- dplyr::filter(data_agr,Session==s)$ItemRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Correct): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTCor),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTCor),0))
    cat('\nRetrieval Item RT (Incorrect): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTIncor),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTIncor),0))
    x <- dplyr::filter(data_agr,Session==s)$ItemRTCor
    y <- dplyr::filter(data_agr,Session==s)$ItemRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nCorrect vs. Incorrect : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Correct, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTCorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTCorPla),0))
    cat('\nRetrieval Item RT (Correct, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTCorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTCorPle),0))
    x <- dplyr::filter(data_agr,Session==s)$ItemRTCorPla
    y <- dplyr::filter(data_agr,Session==s)$ItemRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Incorrect, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTIncorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTIncorPla),0))
    cat('\nRetrieval Item RT (Incorrect, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$ItemRTIncorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$ItemRTIncorPle),0))
    x <- dplyr::filter(data_agr,Session==s)$ItemRTIncorPla
    y <- dplyr::filter(data_agr,Session==s)$ItemRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Source Hitrate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceHitRate),3))
    cat('\nRetrieval Source FArate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceFARate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceFARate),3))
    cat('\nRetrieval Source d-prime (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceDprime),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceDprime),3))
    
    cat('\n\nRetrieval Source HC rate (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceHCRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceHCRate),3))
    cat('\nRetrieval Source HC rate (Hits): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceHCHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceHCHitRate),3))
    cat('\nRetrieval Source HC rate (Hits, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceHCHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceHCHitRatePla),3))
    cat('\nRetrieval Source HC rate (Hits, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceHCHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceHCHitRatePle),3))
    x <- dplyr::filter(data_agr,Session==s)$SourceHCHitRatePla
    y <- dplyr::filter(data_agr,Session==s)$SourceHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Source RT (all): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRT),0))
    cat('\nRetrieval Source RT (Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTPla),0))
    cat('\nRetrieval Source RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTPle),0))
    x <- dplyr::filter(data_agr,Session==s)$SourceRTPla
    y <- dplyr::filter(data_agr,Session==s)$SourceRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Correct): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTCor),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTCor),0))
    cat('\nRetrieval Source RT (Incorrect): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTIncor),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTIncor),0))
    x <- dplyr::filter(data_agr,Session==s)$SourceRTCor
    y <- dplyr::filter(data_agr,Session==s)$SourceRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nCorrect vs. Incorrect : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Correct, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTCorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTCorPla),0))
    cat('\nRetrieval Source RT (Correct, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTCorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTCorPle),0))
    x <- dplyr::filter(data_agr,Session==s)$SourceRTCorPla
    y <- dplyr::filter(data_agr,Session==s)$SourceRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Incorrect, Place): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTIncorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTIncorPla),0))
    cat('\nRetrieval Source RT (Incorrect, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Session==s)$SourceRTIncorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Session==s)$SourceRTIncorPle),0))
    x <- dplyr::filter(data_agr,Session==s)$SourceRTIncorPla
    y <- dplyr::filter(data_agr,Session==s)$SourceRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
  }
  # PRINT THE DESCRIPTIVES PER STIMULATION
  for (s in unique(data_mat_tbt$Stimulation)) {
    cat('\n\nXXXXXXXXXXXXXXXXXX\nXX STIMULATION: ',s,' XX\nXXXXXXXXXXXXXXXXXX\n')
    cat('\nEncoding success (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncResp)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncResp)*100,2))
    cat('\nEncoding success (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncRespPla)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncRespPla)*100,2))
    cat('\nEncoding success (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncRespPle)*100,2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncRespPle)*100,2))
    x <- dplyr::filter(data_agr,Stimulation==s)$EncRespPla
    y <- dplyr::filter(data_agr,Stimulation==s)$EncRespPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nEncoding RT (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncRT),0))
    cat('\nEncoding RT (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncRTPla),0))
    cat('\nEncoding RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$EncRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$EncRTPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$EncRTPla
    y <- dplyr::filter(data_agr,Stimulation==s)$EncRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Item Hitrate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHitRate),3))
    cat('\nRetrieval Item Hitrate (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePla),3))
    cat('\nRetrieval Item Hitrate (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePle),3))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item FArate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemFARate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemFARate),3))
    cat('\nRetrieval Item d-prime (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemDprime),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemDprime),3))
    cat('\nRetrieval Item d-prime (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemDprimePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemDprimePla),3))
    cat('\nRetrieval Item d-prime (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemDprimePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemDprimePle),3))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemDprimePla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemDprimePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Item HC rate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHCRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHCRate),3))
    cat('\nRetrieval Item HC rate (Hits): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRate),3))
    cat('\nRetrieval Item HC rate (Hits, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePla),3))
    cat('\nRetrieval Item HC rate (Hits, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePle),3))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item HC rate (CRs): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemHCCRRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemHCCRRate),3))
    
    cat('\n\nRetrieval Item RT (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRT),0))
    cat('\nRetrieval Item RT (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTPla),0))
    cat('\nRetrieval Item RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemRTPla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Correct): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTCor),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTCor),0))
    cat('\nRetrieval Item RT (Incorrect): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncor),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncor),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemRTCor
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nCorrect vs. Incorrect : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Correct, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPla),0))
    cat('\nRetrieval Item RT (Correct, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Item RT (Incorrect, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPla),0))
    cat('\nRetrieval Item RT (Incorrect, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPla
    y <- dplyr::filter(data_agr,Stimulation==s)$ItemRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Source Hitrate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceHitRate),3))
    cat('\nRetrieval Source FArate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceFARate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceFARate),3))
    cat('\nRetrieval Source d-prime (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceDprime),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceDprime),3))
    
    cat('\n\nRetrieval Source HC rate (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceHCRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceHCRate),3))
    cat('\nRetrieval Source HC rate (Hits): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRate),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRate),3))
    cat('\nRetrieval Source HC rate (Hits, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePla),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePla),3))
    cat('\nRetrieval Source HC rate (Hits, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePle),2),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePle),3))
    x <- dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePla
    y <- dplyr::filter(data_agr,Stimulation==s)$SourceHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    
    cat('\n\nRetrieval Source RT (all): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRT),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRT),0))
    cat('\nRetrieval Source RT (Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTPla),0))
    cat('\nRetrieval Source RT (Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$SourceRTPla
    y <- dplyr::filter(data_agr,Stimulation==s)$SourceRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Correct): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTCor),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTCor),0))
    cat('\nRetrieval Source RT (Incorrect): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncor),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncor),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$SourceRTCor
    y <- dplyr::filter(data_agr,Stimulation==s)$SourceRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nCorrect vs. Incorrect : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Correct, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPla),0))
    cat('\nRetrieval Source RT (Correct, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPla
    y <- dplyr::filter(data_agr,Stimulation==s)$SourceRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
    cat('\nRetrieval Source RT (Incorrect, Place): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPla),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPla),0))
    cat('\nRetrieval Source RT (Incorrect, Pleasantness): M = ',round(mean(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPle),0),', SD = ',round(sd(dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPle),0))
    x <- dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPla
    y <- dplyr::filter(data_agr,Stimulation==s)$SourceRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat('\nPlace vs. Pleasant : t(', res$parameter, ') = ', round(res$statistic,2),', p-value:', round(res$p.value,3), D$method, ':', round(D$estimate,3))
  }
}
# Add the variables we need for the models to data_mat_tbc
vars = c('ItemDprime', 'ItemHCRate', 'ItemHCHitRate', 'ItemHCCRRate')
for (v in vars){
  for (s in unique(data_mat_tbt$Stimulation)){
    var_name = paste('Beh_',v,'_',s ,sep='')
    data_mat_tbc[var_name] <- NA
    for (p in unique(data_mat_tbt$Subject)){
    data_mat_tbc[var_name][data_mat_tbc$Subject==p,] <- dplyr::filter(data_agr, Subject==p & Stimulation==s)[[v]]
    }
  }
}

#### EEG AND MEMORY CORRELATION ####
# get the data from the first session only
cur_data <- dplyr::filter(data_agr, Session==1)
xs = c("ItemDprime", "ItemHCRate", "SourceDprime", "SourceHCRate")
ys = c("Ret_ThetaFro", "Ret_GammaFro", "Ret_ThetaPar", "Ret_GammaPar")
for (x in xs){
    cat("\n")
    for (y in ys){
        cur_c = cor.test(cur_data[[x]], cur_data[[y]], use = "complete.obs", method = "spearman")
        cat("\nCorrelation:", x, "&", y, "=", "r:", round(cur_c$estimate,3), "S:", round(cur_c$statistic,3), "p:", round(cur_c$p.value,3))
        plot(cur_data[[x]], cur_data[[y]], xlab = x, ylab = y)
        title(main = c(x,y))
    }
}

#### BEHAVIORAL DESCRIPTIVES ####
# stimulation blinding check 
if (do_mani_check) {
  # load in the data
  data_mat_mani <- read_csv(data_path_mani, show_col_types = FALSE)
  # check the accuracy
  data_mat_mani$TargetAcc <- ifelse(data_mat_mani$GammaReal == data_mat_mani$TargetSes | data_mat_mani$ThetaReal == data_mat_mani$TargetSes, 1, 0)
  data_mat_mani$TargetAccGamma <- ifelse(data_mat_mani$GammaReal == data_mat_mani$TargetSes, 1, 0)
  data_mat_mani$TargetAccTheta <- ifelse(data_mat_mani$ThetaReal == data_mat_mani$TargetSes, 1, 0)
  data_mat_mani$ShamAcc <- ifelse(data_mat_mani$ShamReal == data_mat_mani$ShamSes, 1, 0)
  mean(data_mat_mani$TargetAcc)
  mean(data_mat_mani$TargetAccGamma)
  mean(data_mat_mani$TargetAccTheta)
  mean(data_mat_mani$ShamAcc)
  # standardize confidence ratings
  data_mat_mani <- mutate(data_mat_mani, across(ends_with('Conf'), .fns = stndz, .names = "s_{.col}"))
  # change the categorical variables into factors
  cur_cat_vars <- c('ppn','GammaReal','ThetaReal','ShamReal','TargetSes','ShamSes')
  lapply(data_mat_mani[,cur_cat_vars], factor)
  print('ACTIVE')
  # create a new variable showing whether the guess matches an Active stimulation session
  data_mat_mani$correct_match <- ifelse(data_mat_mani$TargetSes == data_mat_mani$GammaReal | data_mat_mani$TargetSes == data_mat_mani$ThetaReal, "Active", "Sham")
  # create a contingency table for the guesses vs correct matches
  contingency_table <- table(data_mat_mani$TargetSes, data_mat_mani$correct_match)
  # perform the Chi-square test (not the most reliable given the small count)
  chi_square_result <- chisq.test(contingency_table)
  print(chi_square_result)
  print('SHAM')
  # create a contingency_table
  contingency_table <- table(data_mat_mani$ShamReal, data_mat_mani$ShamSes)
  # perform the Chi-square test (not the most reliable given the small count)
  chi_square_result <- chisq.test(contingency_table)
  print(chi_square_result)
  # perform a logistic regression ACTIVE
  print('ACTIVE')
  chance_level = 2/3 
  logistic_model <- glm(TargetAcc ~ s_TargetConf, 
                      data = data_mat_mani, family = "binomial")
  logistic_model.sum <- summary(logistic_model)
  logistic_model.sum$coefficients[,1]
  logistic_model.sum$coefficients[,1] <- plogis(logistic_model.sum$coefficients[,1]) # in probablities
  print(logistic_model.sum)
  # test if the predicted accuracy is higher than chance
  intercept_log_odds <- coef(logistic_model)[1]
  # Calculate the Log-Odds of a 1/3 Probability
  log_odds_chance_level <- qlogis(chance_level)
  # Standard error of the intercept
  intercept_se <- summary(logistic_model)$coefficients[1, "Std. Error"]
  # Calculate the z-value for one-sided test
  z_value <- (intercept_log_odds - log_odds_chance_level) / intercept_se
  # p-value for a two-sided test
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  cat("Z-Statistic:", z_value, "\n")
  cat("p-Value:", p_value, "\n")

  # perform a logistic regression SHAM
  print('SHAM')
  chance_level = 1/3 
  logistic_model <- glm(ShamAcc ~ s_ShamConf, 
                      data = data_mat_mani, family = "binomial")
  logistic_model.sum <- summary(logistic_model)
  logistic_model.sum$coefficients[,1]
  logistic_model.sum$coefficients[,1] <- plogis(logistic_model.sum$coefficients[,1]) # in probablities
  print(logistic_model.sum)
  # test if the predicted accuracy is higher than chance
  intercept_log_odds <- coef(logistic_model)[1]
  # Calculate the Log-Odds of a 1/3 Probability
  log_odds_chance_level <- qlogis(chance_level)
  # Standard error of the intercept
  intercept_se <- summary(logistic_model)$coefficients[1, "Std. Error"]
  # Calculate the z-value for one-sided test
  z_value <- (intercept_log_odds - log_odds_chance_level) / intercept_se
  # p-value for a two-sided test
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  cat("Z-Statistic:", z_value, "\n")
  cat("p-Value:", p_value, "\n")

}

# get the means and standard deviations
for (cur_var in colnames(data_agr)[6:length(colnames(data_agr))]){
  cur_mean <- aggregate(data_agr[[cur_var]], list(data_agr$Stimulation), FUN=mean) 
  cur_sd <- aggregate(data_agr[[cur_var]], list(data_agr$Stimulation), FUN=sd)
  cur_sum <- data.frame(cur_mean[1],cur_mean[2],cur_sd[2])
  colnames(cur_sum) <- c('stim', 'mean', 'sd')
  cat('\n\n',cur_var,':')
  cat('\n', as.character(cur_sum[1,1]), ':', format(round(cur_sum[1,2],3), nsmall = 3), paste0('(', format(round(cur_sum[1,3],3), nsmall = 3),')'))
  cat('\n', as.character(cur_sum[2,1]), ':', format(round(cur_sum[2,2],3), nsmall = 3), paste0('(', format(round(cur_sum[2,3],3), nsmall = 3),')'))
  cat('\n', as.character(cur_sum[3,1]), ':', format(round(cur_sum[3,2],3), nsmall = 3), paste0('(', format(round(cur_sum[3,3],3), nsmall = 3),')'))
  cat('\n', as.character(cur_sum[4,1]), ':', format(round(cur_sum[4,2],3), nsmall = 3), paste0('(', format(round(cur_sum[4,3],3), nsmall = 3),')'))
  cat('\n', 'average', ':', format(round(sapply(cur_sum[2],mean),3), nsmall = 3), paste0('(', format(round(sapply(cur_sum[3],mean),3), nsmall = 3),')'))
}

#### PREPARE THE DATA: ADJUST THE DATA FOR MODELS ####
# https://marissabarlaz.github.io/portfolio/contrastcoding/
# Sum coding ('centering') of the categorical variables used in the models
data_mat_tbt <- mutate(data_mat_tbt, s_OldNew = OldNew)
contrasts(data_mat_tbt$s_OldNew) <- contr.sum(2)            # New = 1, Old = -1
data_mat_tbt <- mutate(data_mat_tbt, s_Source = Source)
contrasts(data_mat_tbt$s_Source) <- contr.sum(2)            # New = 1, Old = -1
data_mat_tbt <- mutate(data_mat_tbt, s_ONAcc = ONAcc) 
contrasts(data_mat_tbt$s_ONAcc) <- contr.sum(2)           # Incorrect (0) = 1, Correct (1) = -1
data_mat_tbt <- mutate(data_mat_tbt, s_SourceAcc = SourceAcc) 
contrasts(data_mat_tbt$s_SourceAcc) <- contr.sum(2)         # Incorrect (0) = 1, Correct (1) = -1
data_mat_tbt <- mutate(data_mat_tbt, s_ONConf = ONConf) 
contrasts(data_mat_tbt$s_ONConf) <- contr.sum(2)          # LC = 1, HC = -1
data_mat_tbt <- mutate(data_mat_tbt, s_SourceConf = SourceConf) 
contrasts(data_mat_tbt$s_SourceConf) <- contr.sum(2)        # LC = 1, HC = -1
# Make the stimulation difference scores
vars = c('ItemDprime', 'ItemHCRate', 'ItemHCHitRate', 'ItemHCCRRate')
data_mat_tbc$Beh_ItemDprime_GamSha <- data_mat_tbc$Beh_ItemDprime_Gamma - data_mat_tbc$Beh_ItemDprime_Sham
data_mat_tbc$Beh_ItemDprime_TheSha <- data_mat_tbc$Beh_ItemDprime_Theta - data_mat_tbc$Beh_ItemDprime_Sham
data_mat_tbc$Beh_ItemHCRate_GamSha <- data_mat_tbc$Beh_ItemHCRate_Gamma - data_mat_tbc$Beh_ItemHCRate_Sham
data_mat_tbc$Beh_ItemHCRate_TheSha <- data_mat_tbc$Beh_ItemHCRate_Theta - data_mat_tbc$Beh_ItemHCRate_Sham
data_mat_tbc$Beh_ItemHCHitRate_GamSha <- data_mat_tbc$Beh_ItemHCHitRate_Gamma - data_mat_tbc$Beh_ItemHCHitRate_Sham
data_mat_tbc$Beh_ItemHCHitRate_TheSha <- data_mat_tbc$Beh_ItemHCHitRate_Theta - data_mat_tbc$Beh_ItemHCHitRate_Sham
data_mat_tbc$Beh_ItemHCCRRate_GamSha <- data_mat_tbc$Beh_ItemHCCRRate_Gamma - data_mat_tbc$Beh_ItemHCCRRate_Sham
data_mat_tbc$Beh_ItemHCCRRate_TheSha <- data_mat_tbc$Beh_ItemHCCRRate_Theta - data_mat_tbc$Beh_ItemHCCRRate_Sham
# Standardizing of the EEG and behav data. The mean of the data is now 0 and the sd is 1.
data_mat_tbc <- mutate(data_mat_tbc, across(starts_with('Enc') | starts_with('Ret') | starts_with('Beh'), .fns = stndz, .names = "s_{.col}"))
# Transform the data to a wide format
EEG_vars = colnames(select(data_mat_tbc, starts_with('s_Enc') | starts_with('s_Ret') | starts_with('d_Enc') | starts_with('d_Ret')))
beh_vars = colnames(select(data_mat_tbc, starts_with('s_Beh') | starts_with('Beh')))
label_vars = colnames(select(data_mat_tbc, starts_with('label_')))
data_mat_tbc_wide = data.frame()
frst <- TRUE
for (l in label_vars){
  tmp <- select(data_mat_tbc,c('Subject',all_of(l),all_of(beh_vars),all_of(EEG_vars)))
  tmp <- aggregate(tmp[, EEG_vars], tmp[, 1:(2+length(beh_vars))], FUN = mean, na.rm = TRUE)
  tmp2 <- pivot_wider(data = tmp, 
                     names_from = l, 
                     values_from = EEG_vars
                     )
  if (frst){
    data_mat_tbc_wide <- tmp2
    frst <- FALSE
  } else {
    data_mat_tbc_wide <- merge(data_mat_tbc_wide,tmp2,by=c('Subject',all_of(beh_vars)))
  }
}

#### BEHAVIORAL ANALYSIS: Stimulation Effects ####
# https://www.ssc.wisc.edu/sscc/pubs/MM/MM_TestEffects.html
# https://stats.stackexchange.com/questions/120768/different-p-values-for-fixed-effects-in-summary-of-glmer-and-likelihood-rati
# https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#contrasts
# https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html
# https://rdrr.io/cran/emmeans/man/eff_size.html
# https://medium.com/analytics-vidhya/odds-ratio-and-effect-size-a59c968ddda6
# https://stackoverflow.com/questions/73536308/how-to-get-emmeans-to-print-degrees-of-freedom-for-glmer-class
# https://stats.stackexchange.com/questions/524376/testing-glmer-model-assumptions-optionally-in-r
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

# Step 1: select the model, by removing non-significant levels of interactions until none are left
# Step 2: do (relevant) pairwise post-hoc t-tests on the estimated marginal means (predicted values) to interpret the model

# Display the warnings when they occur (to "hide" fixed issues)
options(warn=1)

### MODEL SPECIFICATION
cur_model.name <- ''          # name of model
cur_model.formula.high <- ''  # complex model
cur_model.formula.low <- ''   # simple model
cur_model.formula.pcomp <- '' # pairwise comparison formula for estimated marginal means

# Item Accuracy (iacc) ##
cur_model.name[1]         <- 'model.iacc_test'
cur_model.formula.high[1] <- 'ONAcc ~ 
                              Stimulation * s_OldNew * s_ONConf 
                            + (1 | Subject)'
cur_model.formula.low[1] <- 'ONAcc ~ 
                              Stimulation * s_OldNew 
                            + Stimulation * s_ONConf 
                            + s_OldNew    * s_ONConf 
                            + (1 | Subject)'
cur_model.formula.pcomp[1] <- '~ Stimulation |
                              s_OldNew
                            + s_ONConf'

# Item Confidence (icnf) ##
cur_model.name[2]         <- 'model.icnf_test'
cur_model.formula.high[2] <- 'ONConf ~ 
                              Stimulation * s_OldNew * s_ONAcc 
                            + (1 | Subject)'
cur_model.formula.low[2] <- 'ONConf ~ 
                              Stimulation * s_OldNew 
                            + Stimulation * s_ONAcc 
                            + s_OldNew    * s_ONAcc 
                            + (1 | Subject)'
cur_model.formula.pcomp[2] <- '~ Stimulation |
                              s_OldNew
                            + s_ONAcc'

# Source Accuracy (sacc) ##
cur_model.name[3]         <- 'model.sacc_test'
cur_model.formula.high[3] <- 'SourceAcc ~ 
                              Stimulation * s_SourceConf
                            + s_Source
                            + (1 | Subject)'
cur_model.formula.low[3] <- 'SourceAcc ~ 
                              Stimulation 
                            + s_SourceConf 
                            + s_Source 
                            + (1 | Subject)'
cur_model.formula.pcomp[3] <- '~ Stimulation |
                              s_SourceConf'

# Source Confidence (scnf) ##
cur_model.name[4]         <- 'model.scnf_test'
cur_model.formula.high[4] <- 'SourceConf ~ 
                              Stimulation * s_SourceAcc
                            + s_Source
                            + (1 | Subject)'
cur_model.formula.low[4] <- 'SourceConf ~ 
                              Stimulation 
                            + s_SourceAcc 
                            + s_Source 
                            + (1 | Subject)'
cur_model.formula.pcomp[4] <- '~ Stimulation |
                              s_SourceAcc'


### RUNNING THE MODELS
for (m in 1:length(cur_model.name)) {
  cat(cur_model.name[m])
  # Check if a 3-way interaction is significantly better than a 2-way interaction
  cur_data <- dplyr::filter(data_mat_tbt,Session!=1)
  cur_model.high <- glmer(str_replace_all(cur_model.formula.high[m], "\n", ""),
                          data = cur_data,
                          family = "binomial",
                          control = glmerControl(optimizer = "bobyqa")
                          )
  cur_model.low <- glmer(str_replace_all(cur_model.formula.low[m], "\n", ""),
                          data = cur_data,
                          family = "binomial",
                          control = glmerControl(optimizer = "bobyqa")
                          )
  anova_test <- anova(cur_model.high,cur_model.low)
  if (anova_test$`Pr(>Chisq)`[2] < .05) {
    cat('The anova Likelihood ratio test shows that there a significant effect of the higher level interaction, \nso we will use the more complex model.') 
    cur_model <- cur_model.high
  } else {
    cat('The anova Likelihood ratio test shows that there a no significant effect of the higher level interaction, \nso we will use the less complex model.')
    cur_model <- cur_model.low
  }
  # summarize the model
  cur_model.sum <- as.data.frame(round(coef(summary(cur_model)),3))
  cur_model.sum$Estimate <- plogis( cur_model.sum$Estimate) # Since the estimates are in log odds, we can transform that to probabilities
  cur_model.sum <- subset( cur_model.sum, select = -c(`Std. Error`)) # remove the Std. Err as they are on the log-odds scale (I assume)
  cat('\n***', cur_model.name[m],'***')
  print(summary(cur_model))
  # also save a copy with the current model name
  assign(cur_model.name[m], cur_model)
  assign(paste(cur_model.name[m],'.sum',sep=''), cur_model.sum)
  # check the assumptions
  if (do_assum_check){
    cat('\n*********************\n* Model Performance *\n*********************\n')
    print(model_performance(cur_model))
    cat('\n******************************************************\n* Multicollinearity: Variance Inflation Factor (VIF) *\n******************************************************\n')
    cat('    VIF equal to 1 = variables are not correlated\n')
    cat('    VIF between 1 and 5 = variables are moderately correlated\n')
    cat('    VIF greater than 5 = variables are highly correlated\n')
    print(car::vif(cur_model))
    if (all(car::vif(cur_model) < 5)){
      cat('All variables are max moderately correlated\n')
    } else if (any(car::vif(cur_model) >= 5)){
      warning('\n!!!\nAt least one variable is highly correlated, this needs to be handeled!\n!!!')
    }
    cat('\n*********************************************************************\n* Normality of Residuals: *\n*********************************************************************\n')
    # https://stats.stackexchange.com/questions/232011/ties-should-not-be-present-in-one-sample-kolmgorov-smirnov-test-in-r
    # https://stats.stackexchange.com/questions/92394/checking-residuals-for-normality-in-generalised-linear-models
    cat('Residuals are not expected to have a normal distribution for GLM, but it could not hurt to look at them anyways\n')
    layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
    layout(layout.matrix, heights = c(1, 1), widths = c(1, 1))
    h <- hist(residuals(cur_model),breaks = 10, density = 10,
              col = "gray")
    g <- residuals(cur_model)
    xfit <- seq(min(g), max(g), length = 40) 
    yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
    yfit <- yfit * diff(h$mids[1:2]) * length(g) 
    lines(xfit, yfit, col = "black", lwd = 2)
    qqnorm(residuals(cur_model));qqline(residuals(cur_model))
    plot(fitted(cur_model), residuals(cur_model));abline(h=0)
  }
  ## POST-HOC TESTS ##
  # Now let's look at all (relevant) pairwise comparisons
  em   <- emmeans(cur_model, formula(str_replace_all(paste('pairwise', cur_model.formula.pcomp[m]), "\n", "")))
  emef <- eff_size(em, sigma = sigma(cur_model), edf = df.residual(cur_model) ) # calculate the Cohen's D
  em_means          <- summary(em)$emmean
  em_means$emmean   <- plogis(em_means$emmean) # Since the predicted values are in log odds, we can transform that to probabilities
  em_diffs          <- summary(em)$contrasts
  em_diffs$p.value  <- round(em_diffs$p.value,4)
  em_effsz          <- summary(emef)
  em_means$cond <- NA
  em_means$cond[em_means$s_OldNew=='Old' & em_means$s_ONConf==1] <- 'ItemOldHC'
  em_means$cond[em_means$s_OldNew=='Old' & em_means$s_ONConf==0] <- 'ItemOldLC'
  em_means$cond[em_means$s_OldNew=='New' & em_means$s_ONConf==1] <- 'ItemNewHC'
  em_means$cond[em_means$s_OldNew=='New' & em_means$s_ONConf==0] <- 'ItemNewLC'
  em_means$cond[em_means$s_OldNew=='Old' & em_means$s_ONAcc==1]  <- 'ItemHit'
  em_means$cond[em_means$s_OldNew=='Old' & em_means$s_ONAcc==0]  <- 'ItemMiss'
  em_means$cond[em_means$s_OldNew=='New' & em_means$s_ONAcc==1]  <- 'ItemCR'
  em_means$cond[em_means$s_OldNew=='New' & em_means$s_ONAcc==0]  <- 'ItemFA'
  em_means$cond[em_means$s_SourceConf==1]                        <- 'SourceHC'
  em_means$cond[em_means$s_SourceConf==0]                        <- 'SourceLC'
  em_means$cond[em_means$s_SourceAcc==1]                         <- 'SourceHit'
  em_means$cond[em_means$s_SourceAcc==0]                         <- 'SourceMiss'
  em_diffs$cond <- NA
  em_diffs$cond[em_diffs$s_OldNew=='Old' & em_diffs$s_ONConf==1] <- '1ItemOldHC'
  em_diffs$cond[em_diffs$s_OldNew=='Old' & em_diffs$s_ONConf==0] <- '2ItemOldLC'
  em_diffs$cond[em_diffs$s_OldNew=='New' & em_diffs$s_ONConf==1] <- '3ItemNewHC'
  em_diffs$cond[em_diffs$s_OldNew=='New' & em_diffs$s_ONConf==0] <- '4ItemNewLC'
  em_diffs$cond[em_diffs$s_OldNew=='Old' & em_diffs$s_ONAcc==1]  <- '1ItemHit'
  em_diffs$cond[em_diffs$s_OldNew=='Old' & em_diffs$s_ONAcc==0]  <- '2ItemMiss'
  em_diffs$cond[em_diffs$s_OldNew=='New' & em_diffs$s_ONAcc==1]  <- '3ItemCR'
  em_diffs$cond[em_diffs$s_OldNew=='New' & em_diffs$s_ONAcc==0]  <- '4ItemFA'
  em_diffs$cond[em_diffs$s_SourceConf==1]                        <- '1SourceHC'
  em_diffs$cond[em_diffs$s_SourceConf==0]                        <- '2SourceLC'
  em_diffs$cond[em_diffs$s_SourceAcc==1]                         <- '1SourceHit'
  em_diffs$cond[em_diffs$s_SourceAcc==0]                         <- '2SourceMiss'
  em_diffs$diff[em_diffs$contrast=='Sham - Theta']  <- '1Theta - Sham'
  em_diffs$diff[em_diffs$contrast=='Sham - Gamma']  <- '2Gamma - Sham'
  em_diffs$diff[em_diffs$contrast=='Gamma - Theta'] <- '3Theta - Gamma'
  em_diffs$effsz <- -1*em_effsz$effect.size # flip the sign of the z-ratio as we report the opposite differences in the paper
  # as the difference is the difference in log odds, it cannot directly be transferred to probabilities, therefore, we add the separate probabilities here
  em_diffs$est1[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1ItemOldHC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemOldHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2ItemOldLC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemOldLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='3ItemNewHC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemNewHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='4ItemNewLC' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='ItemNewLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1ItemHit' & em_diffs$diff=='1Theta - Sham']        <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='1ItemHit' & em_diffs$diff=='2Gamma - Sham']        <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1ItemHit' & em_diffs$diff=='3Theta - Gamma']       <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='1ItemHit' & em_diffs$diff=='1Theta - Sham']        <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1ItemHit' & em_diffs$diff=='2Gamma - Sham']        <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1ItemHit' & em_diffs$diff=='3Theta - Gamma']       <- em_means$emmean[em_means$cond=='ItemHit' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2ItemMiss' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='ItemMiss' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='3ItemCR' & em_diffs$diff=='1Theta - Sham']         <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='3ItemCR' & em_diffs$diff=='2Gamma - Sham']         <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='3ItemCR' & em_diffs$diff=='3Theta - Gamma']        <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='3ItemCR' & em_diffs$diff=='1Theta - Sham']         <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='3ItemCR' & em_diffs$diff=='2Gamma - Sham']         <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='3ItemCR' & em_diffs$diff=='3Theta - Gamma']        <- em_means$emmean[em_means$cond=='ItemCR' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='4ItemFA' & em_diffs$diff=='1Theta - Sham']         <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='4ItemFA' & em_diffs$diff=='2Gamma - Sham']         <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='4ItemFA' & em_diffs$diff=='3Theta - Gamma']        <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='4ItemFA' & em_diffs$diff=='1Theta - Sham']         <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='4ItemFA' & em_diffs$diff=='2Gamma - Sham']         <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='4ItemFA' & em_diffs$diff=='3Theta - Gamma']        <- em_means$emmean[em_means$cond=='ItemFA' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1SourceHC' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='1SourceHC' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1SourceHC' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='1SourceHC' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1SourceHC' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1SourceHC' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='SourceHC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2SourceLC' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='2SourceLC' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2SourceLC' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='2SourceLC' & em_diffs$diff=='1Theta - Sham']       <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2SourceLC' & em_diffs$diff=='2Gamma - Sham']       <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2SourceLC' & em_diffs$diff=='3Theta - Gamma']      <- em_means$emmean[em_means$cond=='SourceLC' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1SourceHit' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='1SourceHit' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='1SourceHit' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='1SourceHit' & em_diffs$diff=='1Theta - Sham']      <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1SourceHit' & em_diffs$diff=='2Gamma - Sham']      <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='1SourceHit' & em_diffs$diff=='3Theta - Gamma']     <- em_means$emmean[em_means$cond=='SourceHit' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='1Theta - Sham']     <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Theta']
  em_diffs$est1[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='2Gamma - Sham']     <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Gamma']
  em_diffs$est1[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='3Theta - Gamma']    <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Theta']
  em_diffs$est2[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='1Theta - Sham']     <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='2Gamma - Sham']     <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Sham']
  em_diffs$est2[em_diffs$cond=='2SourceMiss' & em_diffs$diff=='3Theta - Gamma']    <- em_means$emmean[em_means$cond=='SourceMiss' & em_means$Stimulation=='Gamma']
  em_diffs$z.ratio <- -1*em_diffs$z.ratio # flip the sign of the z-ratio as we report the opposite differences in the paper
  em_diffs <- select(em_diffs,'cond','diff','est1','est2','z.ratio','p.value','effsz')
  em_diffs <- em_diffs[order(em_diffs$cond,em_diffs$diff),]
  # store the estimated probabilities, the pairwise differences and the effect size (Cohen's d)
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep=''), em_means)
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.diffs',sep=''), em_diffs)
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.effsz',sep=''), em_effsz)
  # get the estimated probabilities per stimulation, averaged over the conditions
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemNewLC',sep=''),   round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemNewLC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemNewHC',sep=''),   round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemNewHC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemOldLC',sep=''),   round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemOldLC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemOldHC',sep=''),   round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemOldHC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemHit',sep=''),     round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemHit')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemMiss',sep=''),    round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemMiss')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemCR',sep=''),      round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemCR')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.ItemFA',sep=''),      round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='ItemFA')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.SourceHC',sep=''),    round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='SourceHC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.SourceLC',sep=''),    round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='SourceLC')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.SourceHit',sep=''),   round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='SourceHit')$emmean),2))
  assign(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.SourceMiss',sep=''),  round(mean(dplyr::filter(get(paste(unlist(strsplit(cur_model.name[m], split='.', fixed=TRUE))[2],'.means',sep='')),cond=='SourceMiss')$emmean),2))
}

# reset the warnings
options(warn=0)

### Participant stimulation effect count
# difference in high-confidence ihit rate between active and sham
cnf.HCHit <- {}
cnf.HCHit$Subject <- dplyr::filter(data_agr, Session==1)$Subject
cnf.HCHit$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$ItemHCHitRate
cnf.HCHit$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$ItemHCHitRate
cnf.HCHit$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$ItemHCHitRate
cnf.HCHit$TheSha <- cnf.HCHit$Theta - cnf.HCHit$Sham 
cnf.HCHit$GamSha <- cnf.HCHit$Gamma - cnf.HCHit$Sham
mean(cnf.HCHit$TheSha, na.rm = TRUE)
sum(cnf.HCHit$TheSha>0, na.rm = TRUE)
mean(cnf.HCHit$GamSha, na.rm = TRUE)
sum(cnf.HCHit$GamSha>0, na.rm = TRUE)
# difference in high-confidence correct rejection rate between active and sham
cnf.HCCR <- {}
cnf.HCCR$Subject <- dplyr::filter(data_agr, Session==1)$Subject
cnf.HCCR$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$ItemHCCRRate
cnf.HCCR$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$ItemHCCRRate
cnf.HCCR$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$ItemHCCRRate
cnf.HCCR$TheSha <- cnf.HCCR$Theta - cnf.HCCR$Sham 
cnf.HCCR$GamSha <- cnf.HCCR$Gamma - cnf.HCCR$Sham 
mean(cnf.HCCR$TheSha, na.rm = TRUE)
sum(cnf.HCCR$TheSha<0, na.rm = TRUE)
mean(cnf.HCCR$GamSha, na.rm = TRUE)
sum(cnf.HCCR$GamSha>0, na.rm = TRUE)
# difference in hit rate for HC between active and sham
acc.HitHC <- {}
acc.HitHC$Subject <- dplyr::filter(data_agr, Session==1)$Subject
acc.HitHC$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$ItemHitRateHC
acc.HitHC$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$ItemHitRateHC
acc.HitHC$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$ItemHitRateHC
acc.HitHC$TheSha <- acc.HitHC$Theta - acc.HitHC$Sham 
acc.HitHC$GamSha <- acc.HitHC$Gamma - acc.HitHC$Sham
mean(acc.HitHC$TheSha, na.rm = TRUE)
sum(acc.HitHC$TheSha>0, na.rm = TRUE)
mean(acc.HitHC$GamSha, na.rm = TRUE)
sum(acc.HitHC$GamSha>0, na.rm = TRUE)
# difference in correct rejection rate for HC between active and sham
acc.CRHC <- {}
acc.CRHC$Subject <- dplyr::filter(data_agr, Session==1)$Subject
acc.CRHC$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$ItemCRRateHC
acc.CRHC$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$ItemCRRateHC
acc.CRHC$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$ItemCRRateHC
acc.CRHC$TheSha <- acc.CRHC$Theta - acc.CRHC$Sham 
acc.CRHC$GamSha <- acc.CRHC$Gamma - acc.CRHC$Sham
mean(acc.CRHC$TheSha, na.rm = TRUE)
sum(acc.CRHC$TheSha<0, na.rm = TRUE)
mean(acc.CRHC$GamSha, na.rm = TRUE)
sum(acc.CRHC$GamSha>0, na.rm = TRUE)
# difference in high-confidence shit rate between active and sham
cnf.HCsHit <- {}
cnf.HCsHit$Subject <- dplyr::filter(data_agr, Session==1)$Subject
cnf.HCsHit$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$SourceHCHitRate
cnf.HCsHit$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$SourceHCHitRate
cnf.HCsHit$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$SourceHCHitRate
cnf.HCsHit$TheSha <- cnf.HCsHit$Theta - cnf.HCsHit$Sham 
cnf.HCsHit$GamSha <- cnf.HCsHit$Gamma - cnf.HCsHit$Sham
mean(cnf.HCsHit$TheSha, na.rm = TRUE)
sum(cnf.HCsHit$TheSha<0, na.rm = TRUE)
mean(cnf.HCsHit$GamSha, na.rm = TRUE)
sum(cnf.HCsHit$GamSha>0, na.rm = TRUE)
# difference in source hit rate for HC between active and sham
acc.sHitHC <- {}
acc.sHitHC$Subject <- dplyr::filter(data_agr, Session==1)$Subject
acc.sHitHC$Sham <- dplyr::filter(data_agr, Stimulation=='Sham')$SourceHitRateHC
acc.sHitHC$Theta <- dplyr::filter(data_agr, Stimulation=='Theta')$SourceHitRateHC
acc.sHitHC$Gamma <- dplyr::filter(data_agr, Stimulation=='Gamma')$SourceHitRateHC
acc.sHitHC$TheSha <- acc.sHitHC$Theta - acc.sHitHC$Sham 
acc.sHitHC$GamSha <- acc.sHitHC$Gamma - acc.sHitHC$Sham
mean(acc.sHitHC$TheSha, na.rm = TRUE)
sum(acc.sHitHC$TheSha<0, na.rm = TRUE)
mean(acc.sHitHC$GamSha, na.rm = TRUE)
sum(acc.sHitHC$GamSha>0, na.rm = TRUE)

# plot the results
plotnames = c(
    'StimEffs_iAcc',
    'StimEffs_iConf',
    'StimEffs_sAcc',
    'StimEffs_sConf'
    )
y1s = list(
    c(acc.HitHC$TheSha, acc.HitHC$GamSha),
    c(cnf.HCHit$TheSha, cnf.HCHit$GamSha),
    c(acc.sHitHC$TheSha, acc.sHitHC$GamSha),
    c(cnf.HCsHit$TheSha, cnf.HCsHit$GamSha)
    )
y2s = list(
    c(acc.CRHC$TheSha, acc.CRHC$GamSha),
    c(cnf.HCCR$TheSha, cnf.HCCR$GamSha),
    '',
    ''
    )
ylabs = c(
    'Difference in Proportion of Correct Responses',
    'Difference in Proportion of High-Confidenct Responses',
    'Difference in Proportion of Correct Responses',
    'Difference in Proportion of High-Confidenct Responses'
    )
y1tit = c(
    'High-Confident Old',
    'Hits',
    'High-Confident Source',
    'Source Hits'
    )
y2tit = c(
    'High-Confident New',
    'Correct Rejections',
    "",
    ""
    )

for (i in seq_along(plotnames)) {
    stmeffects <- data.frame(
        Subject = rep(dplyr::filter(data_agr, Session == 1)$Subject, times = 2),
        Stim = rep(c('Theta', 'Gamma'), each = length(dplyr::filter(data_agr, Session == 1)$Subject)),
        y1 = y1s[[i]],
        y2 = y2s[[i]]
    )
    plotname            <- plotnames[i]
    y1_mean_theta = mean(dplyr::filter(stmeffects, Stim == "Theta")$y1, na.rm = TRUE)
    y1_mean_gamma = mean(dplyr::filter(stmeffects, Stim == "Gamma")$y1, na.rm = TRUE)
    if (y1_mean_theta > 0) {
        y1_ypos_theta   = .6
        y1_label_theta  = paste("N>0 = ",sum(dplyr::filter(stmeffects, Stim == "Theta")$y1>0, na.rm = TRUE))
    } else if (y1_mean_theta < 0) {
        y1_ypos_theta   = -.6
        y1_label_theta  = paste("N<0 = ",sum(dplyr::filter(stmeffects, Stim == "Theta")$y1<0, na.rm = TRUE))
    }
    if (y1_mean_gamma > 0) {
        y1_ypos_gamma   = .6
        y1_label_gamma  = paste("N>0 = ",sum(dplyr::filter(stmeffects, Stim == "Gamma")$y1>0, na.rm = TRUE))
    } else if (y1_mean_gamma < 0) {
        y1_ypos_gamma   = -.6
        y1_label_gamma  = paste("N<0 = ",sum(dplyr::filter(stmeffects, Stim == "Gamma")$y1<0, na.rm = TRUE))
    }
    p1 <- ggplot(stmeffects, aes(x=Stim, y=y1, fill=Stim)) + 
        geom_violin() + 
        ylab(ylabs[i]) +
        xlab(element_blank()) +
        ggtitle(y1tit[i]) +
        scale_x_discrete(limits=c("Theta", "Gamma")) +
        ylim(-1,1) +
        geom_hline(yintercept=0, linetype="dashed", color = "black") +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=2, binwidth=.01, position=position_jitter(width = .04, height = 0, seed = 11)) +
        stat_summary(fun = "mean", na.rm = TRUE, geom = "crossbar", width = 0.5, colour = "black") +
        geom_label(x=1, y=y1_ypos_theta, label=y1_label_theta, fill="white") +
        geom_label(x=2, y=y1_ypos_gamma, label=y1_label_gamma, fill="white") +
        theme_minimal() +
        theme(panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.position="none",
                plot.title = element_text(hjust=0.5,size=14),
                axis.text.x = element_text(size=14, color = "black"),
                axis.title = element_text(size=14)) + 
        scale_fill_brewer(palette="Accent")
    if (length(y2s[[i]]) != 1) {
        y2_mean_theta = mean(dplyr::filter(stmeffects, Stim == "Theta")$y2, na.rm = TRUE)
        y2_mean_gamma = mean(dplyr::filter(stmeffects, Stim == "Gamma")$y2, na.rm = TRUE)
        if (y2_mean_theta > 0) {
            y2_ypos_theta   = .6
            y2_label_theta  = paste("N>0 = ",sum(dplyr::filter(stmeffects, Stim == "Theta")$y2>0, na.rm = TRUE))
        } else if (y2_mean_theta < 0) {
            y2_ypos_theta   = -.6
            y2_label_theta  = paste("N<0 = ",sum(dplyr::filter(stmeffects, Stim == "Theta")$y2<0, na.rm = TRUE))
        }
        if (y2_mean_gamma > 0) {
            y2_ypos_gamma   = .6
            y2_label_gamma  = paste("N>0 = ",sum(dplyr::filter(stmeffects, Stim == "Gamma")$y2>0, na.rm = TRUE))
        } else if (y2_mean_gamma < 0) {
            y2_ypos_gamma   = -.6
            y2_label_gamma  = paste("N<0 = ",sum(dplyr::filter(stmeffects, Stim == "Gamma")$y2<0, na.rm = TRUE))
        }
        p2 <- ggplot(stmeffects, aes(x=Stim, y=y2, fill=Stim)) + 
            geom_violin() + 
            ylab(ylabs[i]) +
            xlab(element_blank()) +
            ggtitle(y2tit[i]) +
            scale_x_discrete(limits=c("Theta", "Gamma")) +
            ylim(-1,1) +
            geom_hline(yintercept=0, linetype="dashed", color = "black") +
            geom_dotplot(binaxis='y', stackdir='center', dotsize=2, binwidth=.01, position=position_jitter(width = .04, height = 0, seed = 11)) +
            stat_summary(fun = "mean", na.rm = TRUE, geom = "crossbar", width = 0.5, colour = "black") +
            geom_label(x=1, y=y2_ypos_theta, label=y2_label_theta, fill="white") +
            geom_label(x=2, y=y2_ypos_gamma, label=y2_label_gamma, fill="white") +
            theme_minimal() +
            theme(panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    legend.position="none",
                    plot.title = element_text(hjust=0.5,size=14),
                    axis.text.x = element_text(size=14, color = "black"),
                    axis.title = element_text(size=14)) + 
            scale_fill_brewer(palette="Accent")
        p <- arrangeGrob(p1, p2, ncol = 2)
    } else {
        p = p1
    }
    ggsave(paste(outp_path, plotname, '.pdf',sep=''), width = 2000, height = 2000, units = 'px', p)
    ggsave(paste(outp_path, plotname, '.png',sep=''), width = 2000, height = 2000, units = 'px', p)
}

#### EEG ANALYSIS: Brain predicting behavior ####
# https://cran.r-project.org/web/packages/emmeans/vignettes/models.html
# https://stackoverflow.com/questions/53034261/warning-lme4-model-failed-to-converge-with-maxgrad
# https://stats.stackexchange.com/questions/352285/paired-data-comparison-regression-or-paired-t-test
# https://www.investopedia.com/terms/v/variance-inflation-factor.asp
# https://rpubs.com/DragonflyStats/Cooks-Distance
# https://statisticsbyjim.com/regression/multicollinearity-in-regression-analysis
# http://www.sthda.com/english/wiki/normality-test-in-r

# Display the warnings when they occur (to "hide" fixed issues)
options(warn=1)

### MODEL SPECIFICATION
cur_model.name <- ''
cur_model.formula <- ''
# Theta HC Hit (TheHCHit) ##
cur_model.name[1] <- 'model.thehchit'
cur_model.formula[1] <- 'Beh_ItemHCHitRate_TheSha ~ 
                        s_Enc_pow_theta_frontal_iHCCor
                      + s_Enc_pow_theta_parietal_iHCCor
                      + d_Enc_pek_theta_frontal_iHCCor
                      + d_Enc_pek_theta_parietal_iHCCor
                      + s_Enc_phc_theta_iHCCor
                      + s_Enc_pac_frontal_iHCCor
                      + s_Enc_pac_parietal_iHCCor
                      + s_Ret_pow_theta_frontal_iHCCor
                      + s_Ret_pow_theta_parietal_iHCCor
                      + d_Ret_pek_theta_frontal_iHCCor
                      + d_Ret_pek_theta_parietal_iHCCor
                      + s_Ret_phc_theta_iHCCor
                      + s_Ret_pac_frontal_iHCCor
                      + s_Ret_pac_parietal_iHCCor'
# Theta HC CR  (TheHCCR) ##
cur_model.name[2] <- 'model.thehccr'
cur_model.formula[2] <- 'Beh_ItemHCCRRate_TheSha ~ 
                        s_Ret_pow_theta_frontal_iHCCor
                      + s_Ret_pow_theta_parietal_iHCCor
                      + d_Ret_pek_theta_frontal_iHCCor
                      + d_Ret_pek_theta_parietal_iHCCor
                      + s_Ret_phc_theta_iHCCor
                      + s_Ret_pac_frontal_iHCCor
                      + s_Ret_pac_parietal_iHCCor'
# Gamma HC Hit (GamHCHit) ##
cur_model.name[3] <- 'model.gamhchit'
cur_model.formula[3] <- 'Beh_ItemHCHitRate_GamSha ~ 
                        s_Enc_pow_theta_frontal_iHCCor
                      + s_Enc_pow_theta_parietal_iHCCor
                      + d_Enc_pek_theta_frontal_iHCCor
                      + d_Enc_pek_theta_parietal_iHCCor
                      + s_Enc_phc_theta_iHCCor
                      + s_Enc_pac_frontal_iHCCor
                      + s_Enc_pac_parietal_iHCCor
                      + s_Ret_pow_theta_frontal_iHCCor
                      + s_Ret_pow_theta_parietal_iHCCor
                      + d_Ret_pek_theta_frontal_iHCCor
                      + d_Ret_pek_theta_parietal_iHCCor
                      + s_Ret_phc_theta_iHCCor
                      + s_Ret_pac_frontal_iHCCor
                      + s_Ret_pac_parietal_iHCCor'
# Gamma HC CR  (GamHCCR) ##
cur_model.name[4] <- 'model.gamhccr'
cur_model.formula[4] <- 'Beh_ItemHCCRRate_GamSha ~ 
                        s_Ret_pow_gamma_frontal_iHCCor
                      + s_Ret_pow_gamma_parietal_iHCCor
                      + d_Ret_pek_gamma_frontal_iHCCor
                      + d_Ret_pek_gamma_parietal_iHCCor
                      + s_Ret_phc_gamma_iHCCor
                      + s_Ret_pac_frontal_iHCCor
                      + s_Ret_pac_parietal_iHCCor'

### MODEL CHECKS
# For each model, we have looked at the checks below and decided if it was (significantly) better to exclude outliers or not.
# When running new/changed models, set them all to TRUE and evaluate again which ones should be set to FALSE.
# TRUE: remove the outliers based on the Cook's distance
# FALSE: do not remove outliers / there are no outliers
cur_model.check = c(TRUE,TRUE,TRUE,TRUE, # TheHCHit, TheHCCR, GamHCHit, GamHCCR
                    )

### RUNNING THE MODELS
for (m in 5:length(cur_model.name)) {
  outlrs = c()
  outlrs_log = TRUE
  cnt = 0
  while (outlrs_log == TRUE && cnt < 2) {
    # Redo the model once if there are any outliers
    cur_data <- data_mat_tbc_wide[!data_mat_tbc_wide$Subject %in% outlrs,]
    cur_model <- lm(str_replace_all(cur_model.formula[m], "\n", ""),
                    data = cur_data)
    # also save a copy with the current model name
    assign(cur_model.name[m], cur_model)
    assign(paste(cur_model.name[m],'.sum',sep=''), as.data.frame(round(coef(summary(cur_model)),3)))
    cat('\n***', cur_model.name[m],'***')
    print(summary(cur_model))
    if (do_assum_check){
      cat('\n*********************\n* Model Performance *\n*********************\n')
      print(model_performance(cur_model))
      cat('\n******************************************************\n* Multicollinearity: Variance Inflation Factor (VIF) *\n******************************************************\n')
      cat('    VIF equal to 1 = variables are not correlated\n')
      cat('    VIF between 1 and 5 = variables are moderately correlated\n')
      cat('    VIF greater than 5 = variables are highly correlated\n')
      print(car::vif(cur_model))
      if (all(car::vif(cur_model) < 5)){
        cat('All variables are max moderately correlated\n')
      } else if (any(car::vif(cur_model) >= 5)){
        warning('\n!!!\nAt least one variable is highly correlated, this needs to be handeled!\n!!!')
      }
      cat('\n*********************************************************************\n* Normality of Residuals: Shapiro-Wilk’s & Kolmogorov–Smirnov tests *\n*********************************************************************\n')
      s <- shapiro.test(residuals(cur_model))
      cat('SW: statistic:',s$statistic,', p-value:',s$p.value,'\n')
      k <- ks.test(residuals(cur_model), "pnorm")
      cat('KS: statistic:',k$statistic,', p-value:',k$p.value,'\n')
      if (s$p.value >= .05 && k$p.value >= .05){
        cat('Residuals are normally distributed\n')
      } else {
        warning('\n!!!\nThe residuals may not be normally distributed, this needs to be checked!\n!!!')
      }
      layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
      layout(layout.matrix, heights = c(1, 1), widths = c(1, 1))
      h <- hist(residuals(cur_model),breaks = 10, density = 10,
                col = "gray")
      g <- residuals(cur_model)
      xfit <- seq(min(g), max(g), length = 40) 
      yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
      yfit <- yfit * diff(h$mids[1:2]) * length(g) 
      lines(xfit, yfit, col = "black", lwd = 2)
      qqnorm(residuals(cur_model));qqline(residuals(cur_model))
      plot(fitted(cur_model), residuals(cur_model));abline(h=0)
      cat("\n**************************************\n* Outlier Detection: Cook's Distance & *\n**************************************\n")
      cooks.distance(cur_model)
      plot(cur_model,which=4)
      # for the cut-off: "4/(number of observations - number of explanatory variables - 1)"
      cook_cutoff = 4/(length(cur_model$effects)-(length(cur_model$coefficients)-1)-1)
      abline(h=cook_cutoff,lty=2)
      cat('Cook cutoff:\n',cook_cutoff,'\n')
      outlpps = as.numeric(names(cooks.distance(cur_model))[cooks.distance(cur_model)>cook_cutoff])+100
      cat('participants above cutoff:',outlpps)
      cat('\n')
      # check if there are any outliers
      if (length(outlpps)>0 && cnt==0 && cur_model.check[m]==TRUE) {
        outlrs_log = TRUE
        outlrs = c(outlrs, outlpps)
      } else {
        outlrs_log = FALSE
      }
      par(mfrow=c(1,1))
      mtext(cur_model.name[m], side = 3, line = -1, outer = TRUE)
    }
    cnt = cnt+1
  }
}

# reset the warnings
# options(warn=0)

#### PLOTTING: EEG RESULTS ####
freqs = c('Theta', 'Gamma')
locs = c('Frontal', 'Parietal')
name_bases = c('Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_', 
                'Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_')
name_conds = c('iHit-CR', 'iHC-iLC', 
                'iHitHC-iHitLC', 'iCRHC-iCRLC')
name_ext = '.png'
title_conds = c('Hit - Correct rejections', 'High-confidence - Low-confidence', 
                'Hit HC - Hit LC', 'Correct rejections HC - LC')

for (name_cond in seq_along(name_conds)){
    p <- list()
    plotname <- name_conds[name_cond]
    for (freq in seq_along(freqs)){
        brn_img = readPNG(paste(brpl_path, name_bases[2], tolower(freqs[freq]), '_', name_conds[name_cond], name_ext, sep=''), native = TRUE)
        p[[freqs[freq]]][['topo']] <- ggplot() + 
        geom_blank() +
        ggtitle(paste(freqs[freq], '\n', title_conds[name_cond], sep = "")) +
        theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
                panel.background = element_blank()) +
        inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = 'full')
        for (loc in 1:length(locs)){
            brn_img = readPNG(paste(brpl_path, name_bases[1], tolower(freqs[freq]), '_', tolower(locs[loc]), '_', name_conds[name_cond], name_ext, sep=''), native = TRUE)
            p[[freqs[freq]]][[locs[loc]]] <- ggplot() + 
                    geom_blank() +
                    ggtitle(paste(locs[loc], '\n', title_conds[name_cond], sep = "")) +
                    theme(plot.title = element_text(size = 9,hjust = 0.55 ,vjust = 8),
                            panel.background = element_blank()) +
                    inset_element(brn_img,
                                    left = 0,
                                    bottom = 0,
                                    right = 1,
                                    top = .93,
                                    on_top = TRUE,
                                    align_to = 'full')
        }
    }
    ## Combine ##
    ggdraw() +
    draw_plot(p[[freqs[1]]][['topo']],  x = 0.00,  y = .50, width = .50, height = .48) +
    draw_plot(p[[freqs[2]]][['topo']],  x = 0.50,  y = .50, width = .50, height = .48) +
    draw_plot(p[[freqs[1]]][[locs[1]]], x = 0.02,  y = .00, width = .25, height = .46) +
    draw_plot(p[[freqs[1]]][[locs[2]]], x = 0.23,  y = .00, width = .25, height = .46) +
    draw_plot(p[[freqs[2]]][[locs[1]]], x = 0.52,  y = .00, width = .25, height = .46) +
    draw_plot(p[[freqs[2]]][[locs[2]]], x = 0.73, y = .00, width = .25, height = .46) +
    draw_plot_label(label = c("A" , "B" , "C" ,  "D" , "E" , "F"), size = 15,
                    x = c(0.10,  0.60, 0.04, 0.25, 0.54, 0.75), 
                    y = c(1.00, 1.00, 0.49, 0.49, 0.49, 0.49))
    ## Save the plots ##
    ggsave(paste(outp_path, plotname, '.pdf',sep=''), width = 4000, height = 3000, units = 'px')
    ggsave(paste(outp_path, plotname, '.png',sep=''), width = 4000, height = 3000, units = 'px')
    ggsave(paste(outp_path, plotname, '.jpg',sep=''), width = 4000, height = 3000, units = 'px')
}

