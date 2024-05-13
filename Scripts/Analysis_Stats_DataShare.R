#### LIBRARIES ####
library(plyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(reshape2)
library(effsize)
library(emmeans)
library(performance)
library(see)
library(png)
library(patchwork)
library(data.table)
library(Rmisc)

#### VARIABLES ####
rm(list = ls()) # clear the global environment
dirroot <- dirname(rstudioapi::getSourceEditorContext()$path)
data_path <- paste(dirroot, "/Data_TbT.csv", sep = "")
outp_path <- paste(dirroot, "/Output/", sep = "")
brpl_path <- paste(dirroot, "/BrainPlots/", sep = "")
outl_thres <- 3 #SD
cp1 <- c("royalblue", "turquoise", "firebrick", "plum")
cp2 <- c("forestgreen", "purple3")
do_stats_print <- FALSE # to print or not print the stats/behavioral data

#### FUNCTIONS ####
stndz <- function(x, na_rm = TRUE, outl_filt = "None") {
  if (outl_filt[1] == "None") {
    (x - mean(x, na.rm = na_rm)) / sd(x, na.rm = na_rm)
  } else {
    (x - mean(x[outl_filt == FALSE], na.rm = na_rm)) / sd(x[outl_filt == FALSE], na.rm = na_rm)
  }
}
#### PREPARE THE DATA: CALCULATE NEW VARIABLES ####
# Read in the data
data_mat <- read_csv(data_path, show_col_types = FALSE)
# Replace NaNs with NA since NaN means "Not a Number" and NA is "Not Available"
data_mat[sapply(data_mat, is.nan)] <- NA
data_mat$Source[data_mat$Source == "N/A"] <- NA
# calculate the accuracy of the memory responses
data_mat <- mutate(data_mat, ONAcc = NA)
data_mat$ONAcc[select(data_mat, OldNew) == "Old" & select(data_mat,ONResp)<3] = 1 # hits
data_mat$ONAcc[select(data_mat, OldNew) == "New" & select(data_mat,ONResp)>3] = 1 # correct rejections
data_mat$ONAcc[select(data_mat, OldNew) == "Old" & select(data_mat,ONResp)>3] = 0 # misses
data_mat$ONAcc[select(data_mat, OldNew) == "New" & select(data_mat,ONResp)<3] = 0 # false alarms
data_mat$ONAcc[select(data_mat, ONResp) == 3] = 0 # guesses, here classified as incorrect
data_mat <- mutate(data_mat, SourceAcc = NA)
data_mat$SourceAcc[select(data_mat, Source) == "Pleasant" & select(data_mat,SourceResp)<3] = 1 # hits
data_mat$SourceAcc[select(data_mat, Source) == "Place" & select(data_mat,SourceResp)>3] = 1 # hits
data_mat$SourceAcc[select(data_mat, Source) == "Pleasant" & select(data_mat,SourceResp)>3] = 0 # misses
data_mat$SourceAcc[select(data_mat, Source) == "Place" & select(data_mat,SourceResp)<3] = 0 # misses
data_mat$SourceAcc[(select(data_mat, Source) == "Pleasant" | select(data_mat,Source)=="Place") & (select(data_mat,SourceResp)==3)] = 0 # guesses, here classified as incorrect
# calculate the confidence of the memory responses
data_mat <- mutate(data_mat, ONConf = NA)
data_mat$ONConf[select(data_mat, ONResp) == 5 | select(data_mat,ONResp) == 1] = 1 # high-confident
data_mat$ONConf[select(data_mat, ONResp) < 5 & select(data_mat,ONResp) > 1] = 0 # low-confident
data_mat <- mutate(data_mat, SourceConf = NA)
data_mat$SourceConf[select(data_mat, SourceResp) == 5 | select(data_mat,SourceResp) == 1] = 1 # high-confident
data_mat$SourceConf[select(data_mat, SourceResp) < 5 & select(data_mat,SourceResp) > 1] = 0 # low-confident
# Change the categorical variables into factors..
cur_cat_vars <- c("Subject", "Sex", "OldNew", "Source", "ONAcc", "SourceAcc", "ONConf", "SourceConf")
data_mat[,cur_cat_vars] <- lapply(data_mat[, cur_cat_vars], factor)
# .. Or ordered factors
cur_cat_vars <- c("Session", "EncResp", "ONResp", "SourceResp")
data_mat[,cur_cat_vars] <- lapply(data_mat[, cur_cat_vars], ordered)
# Adjust the order of the Stim factor levels
data_mat$Stimulation <- factor(data_mat$Stimulation, levels = c("None","Sham","Gamma","Theta"))

#### PREPARE THE DATA: ADJUST THE DATA FOR MODELS ####
# https://marissabarlaz.github.io/portfolio/contrastcoding/
# Sum coding ("centering") of the categorical variables used in the models
data_mat <- mutate(data_mat, s_OldNew = OldNew)
contrasts(data_mat$s_OldNew) <- contr.sum(2)            # New = 1, Old = -1
data_mat <- mutate(data_mat, s_Source = Source)
contrasts(data_mat$s_Source) <- contr.sum(2)            # New = 1, Old = -1
data_mat <- mutate(data_mat, s_ONAcc = ONAcc)
contrasts(data_mat$s_ONAcc) <- contr.sum(2)           # Incorrect (0) = 1, Correct (1) = -1
data_mat <- mutate(data_mat, s_SourceAcc = SourceAcc)
contrasts(data_mat$s_SourceAcc) <- contr.sum(2)         # Incorrect (0) = 1, Correct (1) = -1
data_mat <- mutate(data_mat, s_ONConf = ONConf)
contrasts(data_mat$s_ONConf) <- contr.sum(2)          # LC = 1, HC = -1
data_mat <- mutate(data_mat, s_SourceConf = SourceConf)
contrasts(data_mat$s_SourceConf) <- contr.sum(2)        # LC = 1, HC = -1
# Making a filtering variable for the EEG outliers (pow & pac) (threshold SD > outl_thres)
pow_thet_vars <- c("Ret_ThetaFro", "Ret_ThetaPar")
pow_gamm_vars <- c("Ret_GammaFro", "Ret_GammaPar")
pac_enc_vars  <- c("Enc_frontalPAC", "Enc_parietalPAC")
pac_ret_vars  <- c("Ret_frontalPAC", "Ret_parietalPAC")
data_mat <- mutate(data_mat, outlrs_powt = FALSE)
data_mat <- mutate(data_mat, outlrs_powg = FALSE)
data_mat <- mutate(data_mat, outlrs_pace = FALSE)
data_mat <- mutate(data_mat, outlrs_pacr = FALSE)
z_scores <- mutate_at(data_mat, c(pow_thet_vars, ow_gamm_vars, pac_enc_vars, pac_ret_vars), stndz, outl_filt = "None")
data_mat$outlrs_powt[which(rowSums(select(z_scores, all_of(pow_thet_vars))>outl_thres 
                                   | select(z_scores, all_of(pow_thet_vars))<(-outl_thres), na_rm = TRUE) >= 1)] = TRUE
data_mat$outlrs_powg[which(rowSums(select(z_scores, all_of(pow_gamm_vars))>outl_thres 
                                   | select(z_scores, all_of(pow_gamm_vars))<(-outl_thres), na_rm = TRUE) >= 1)] = TRUE
data_mat$outlrs_pace[which(rowSums(select(z_scores, all_of(pac_enc_vars))>outl_thres 
                                   | select(z_scores, all_of(pac_enc_vars))<(-outl_thres), na_rm = TRUE) >= 1)] = TRUE
data_mat$outlrs_pacr[which(rowSums(select(z_scores, all_of(pac_ret_vars))>outl_thres 
                                   | select(z_scores, all_of(pac_ret_vars))<(-outl_thres), na_rm = TRUE) >= 1)] = TRUE
rm(z_scores)
# Standardizing of the EEG data (based on the data without the outliers). The mean of the non-outlier data is now 0 and the sd is 1. You do need to use the outlr filtering variable to achieve this.
data_mat <- mutate(data_mat, s_Ret_ThetaFro    = Ret_ThetaFro,
                  s_Ret_ThetaPar    = Ret_ThetaPar,
                  s_Ret_GammaFro    = Ret_GammaFro,
                  s_Ret_GammaPar    = Ret_GammaPar,
                  s_Enc_frontalPAC  = Enc_frontalPAC,
                  s_Enc_parietalPAC = Enc_parietalPAC,
                  s_Ret_frontalPAC  = Ret_frontalPAC,
                  s_Ret_parietalPAC = Ret_parietalPAC)
pow_thet_vars_s <- c("s_Ret_ThetaFro", "s_Ret_ThetaPar")
pow_gamm_vars_s <- c("s_Ret_GammaFro", "s_Ret_GammaPar")
pac_enc_vars_s  <- c("s_Enc_frontalPAC", "s_Enc_parietalPAC")
pac_ret_vars_s  <- c("s_Ret_frontalPAC", "s_Ret_parietalPAC")
data_mat <- mutate_at(data_mat, pow_thet_vars_s, stndz, outl_filt = data_mat$outlrs_powt)
data_mat <- mutate_at(data_mat, pow_gamm_vars_s, stndz, outl_filt = data_mat$outlrs_powg)
data_mat <- mutate_at(data_mat, pac_enc_vars_s, stndz, outl_filt = data_mat$outlrs_pace)
data_mat <- mutate_at(data_mat, pac_ret_vars_s, stndz, outl_filt = data_mat$outlrs_pacr)
# Making longer data
data_mat_long <- reshape(data_mat, direction = "long",
                varying = list(c("Ret_ThetaFro", "Ret_ThetaPar"), c("Ret_GammaFro", "Ret_GammaPar"), c("Enc_frontalPAC", "Enc_parietalPAC"), c("Ret_frontalPAC", "Ret_parietalPAC"),
                                 c("s_Ret_ThetaFro", "s_Ret_ThetaPar"), c("s_Ret_GammaFro", "s_Ret_GammaPar"), c("s_Enc_frontalPAC", "s_Enc_parietalPAC"), c("s_Ret_frontalPAC", "s_Ret_parietalPAC")
                            ),
                timevar = c("BrainRegion"),
                times   = c("Frontal", "Parietal"),
                v.names = c("Ret_Theta", "Ret_Gamma", "Enc_PAC", "Ret_PAC", "s_Ret_Theta", "s_Ret_Gamma", "s_Enc_PAC", "s_Ret_PAC"),
                idvar = c("Subject", "Session", "RetTrial")
                )
data_mat_long <- arrange(data_mat_long, Subject, Session, RetTrial)
# Change the categorical variable into a factor and make a sum-coded variant
data_mat_long$BrainRegion <- factor(data_mat_long$BrainRegion)
data_mat_long <- mutate(data_mat_long, s_BrainRegion = BrainRegion)
contrasts(data_mat_long$s_BrainRegion) <- contr.sum(2)            # Frontal = 1, Parietal = -1

#### AGGREGATED DATA ####
# build the basic structure
data_agr = unique(data_mat[ , c("Subject", "Session", "Stimulation")])
vars <- c("Age","Sex","EncResp","EncRespPla","EncRespPle","EncRT","EncRTPla","EncRTPle",
          "ItemHitRate","ItemHitRatePla","ItemHitRatePle","ItemFARate","ItemDprime","ItemDprimePla","ItemDprimePle",
          "ItemHCRate","ItemHCHitRate","ItemHCHitRatePla","ItemHCHitRatePle","ItemHCCRRate",
          "ItemRT","ItemRTPla","ItemRTPle","ItemRTCor","ItemRTCorPla","ItemRTCorPle","ItemRTIncor","ItemRTIncorPla","ItemRTIncorPle",
          "SourceHitRate","SourceFARate","SourceDprime",
          "SourceHCRate","SourceHCHitRate","SourceHCHitRatePla","SourceHCHitRatePle",
          "SourceRT","SourceRTPla","SourceRTPle","SourceRTCor","SourceRTCorPla","SourceRTCorPle","SourceRTIncor","SourceRTIncorPla","SourceRTIncorPle"
          )
data_agr[, vars] <- NA
# get the data
for (p in unique(data_mat$Subject)) {
  cat("\n\nPARTICIPANT: ",p)
  for (s in unique(data_mat$Session)) {
    # get the data of this participant and session
    cur_data <- filter(data_mat,Subject==p,Session==s)
    data_agr$Age[data_agr$Subject==p & data_agr$Session==s] <- cur_data$Age[1]
    data_agr$Sex[data_agr$Subject==p & data_agr$Session==s] <- cur_data$Sex[1]
    # ENCODING
    data_agr$EncResp[data_agr$Subject==p & data_agr$Session==s] <- mean(as.numeric(levels(cur_data$EncResp)[cur_data$EncResp]), na.rm = TRUE)/3
    data_agr$EncRespPla[data_agr$Subject==p & data_agr$Session==s] <- mean(as.numeric(levels(filter(cur_data,Source=="Place")$EncResp)[filter(cur_data,Source=="Place")$EncResp]), na.rm = TRUE)/3
    data_agr$EncRespPle[data_agr$Subject==p & data_agr$Session==s] <- mean(as.numeric(levels(filter(cur_data,Source=="Pleasant")$EncResp)[filter(cur_data,Source=="Pleasant")$EncResp]), na.rm = TRUE)/3
    data_agr$EncRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$EncRT, na.rm = TRUE)
    data_agr$EncRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Place")$EncRT, na.rm = TRUE)
    data_agr$EncRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Pleasant")$EncRT, na.rm = TRUE)
    # RETRIEVAL
    # hit rate & fa rate: item
    data_agr$ItemHitRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="Old")$ONResp)[filter(cur_data,OldNew=="Old")$ONResp])<3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,OldNew=="Old")$ONResp)[filter(cur_data,OldNew=="Old")$ONResp]))) # (HitHC + HitLC) / old responses
    data_agr$ItemHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Place")$ONResp)[filter(cur_data,Source=="Place")$ONResp])<3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Place")$ONResp)[filter(cur_data,Source=="Place")$ONResp]))) # (HitHC + HitLC) / old responses
    data_agr$ItemHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$ONResp)[filter(cur_data,Source=="Pleasant")$ONResp])<3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Pleasant")$ONResp)[filter(cur_data,Source=="Pleasant")$ONResp]))) # (HitHC + HitLC) / old responses
    data_agr$ItemFARate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="New")$ONResp)[filter(cur_data,OldNew=="New")$ONResp])<3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,OldNew=="New")$ONResp)[filter(cur_data,OldNew=="New")$ONResp]))) # (FAHC + FALC) / new responses
    # hit rate & fa rate: source
    data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp])>3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp]))) # Place Hits / (Place Hits + Place Misses)
    data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp])>3, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp]))) # Pleasant Hits / (Pleasant Hits + Pleasant Misses)
    # d-prime: item
      # adjust the hit and fa rate if it is 1 or 0 to make d" not inf
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
      # adjust the hit and fa rate if it is 1 or 0 to make d" not inf
    if (data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] == 1) {data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s] <- .999}
    if (data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] == 0) {data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s] <- .001}
    data_agr$SourceDprime[data_agr$Subject==p & data_agr$Session==s] <- qnorm(data_agr$SourceHitRate[data_agr$Subject==p & data_agr$Session==s])-
      qnorm(data_agr$SourceFARate[data_agr$Subject==p & data_agr$Session==s])
    # confidence: item
    data_agr$ItemHCRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(cur_data$ONResp)[cur_data$ONResp])==1 | as.numeric(levels(cur_data$ONResp)[cur_data$ONResp])==5, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(cur_data$ONResp)[cur_data$ONResp]))) # HC  / all responses
    data_agr$ItemHCHitRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="Old")$ONResp)[filter(cur_data,OldNew=="Old")$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,OldNew=="Old")$ONResp)[filter(cur_data,OldNew=="Old")$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Place")$ONResp)[filter(cur_data,OldNew=="Old" & Source=="Place")$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Place")$ONResp)[filter(cur_data,OldNew=="Old" & Source=="Place")$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Pleasant")$ONResp)[filter(cur_data,OldNew=="Old" & Source=="Pleasant")$ONResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,OldNew=="Old" & Source=="Pleasant")$ONResp)[filter(cur_data,OldNew=="Old" & Source=="Pleasant")$ONResp])<3, na.rm = TRUE) # HCHits  / Hits
    data_agr$ItemHCCRRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,OldNew=="New")$ONResp)[filter(cur_data,OldNew=="New")$ONResp])==5, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,OldNew=="New")$ONResp)[filter(cur_data,OldNew=="New")$ONResp])>3, na.rm = TRUE) # HCCRs  / CRs
    # confidence: source
    data_agr$SourceHCRate[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp])==1 | as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp])==5, na.rm = TRUE)/
      sum(!is.na(as.numeric(levels(cur_data$SourceResp)[cur_data$SourceResp]))) # HC  / all responses
    data_agr$SourceHCHitRate[data_agr$Subject==p & data_agr$Session==s] <- (sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp])==1, na.rm = TRUE) + sum(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp])==5, na.rm = TRUE))/
      (sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp])<3, na.rm = TRUE) + sum(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp])>3, na.rm = TRUE)) # HCHits  / Hits
    data_agr$SourceHCHitRatePla[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp])==5, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,Source=="Place")$SourceResp)[filter(cur_data,Source=="Place")$SourceResp])>3, na.rm = TRUE) # HCHits  / Hits
    data_agr$SourceHCHitRatePle[data_agr$Subject==p & data_agr$Session==s] <- sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp])==1, na.rm = TRUE)/
      sum(as.numeric(levels(filter(cur_data,Source=="Pleasant")$SourceResp)[filter(cur_data,Source=="Pleasant")$SourceResp])<3, na.rm = TRUE) # HCHits  / Hits
    # reaction time: item
    data_agr$ItemRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$ONRT, na.rm = TRUE)
    data_agr$ItemRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Place")$ONRT, na.rm = TRUE)
    data_agr$ItemRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Pleasant")$ONRT, na.rm = TRUE)
    data_agr$ItemRTCor[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==1)$ONRT, na.rm = TRUE)
    data_agr$ItemRTCorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==1 & Source=="Place")$ONRT, na.rm = TRUE)
    data_agr$ItemRTCorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==1 & Source=="Pleasant")$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncor[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==0)$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==0 & Source=="Place")$ONRT, na.rm = TRUE)
    data_agr$ItemRTIncorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,ONAcc==0 & Source=="Pleasant")$ONRT, na.rm = TRUE)
    # reaction time: Source
    data_agr$SourceRT[data_agr$Subject==p & data_agr$Session==s] <- mean(cur_data$SourceRT, na.rm = TRUE)
    data_agr$SourceRTPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Place")$SourceRT, na.rm = TRUE)
    data_agr$SourceRTPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,Source=="Pleasant")$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCor[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==1)$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==1 & Source=="Place")$SourceRT, na.rm = TRUE)
    data_agr$SourceRTCorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==1 & Source=="Pleasant")$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncor[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==0)$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncorPla[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==0 & Source=="Place")$SourceRT, na.rm = TRUE)
    data_agr$SourceRTIncorPle[data_agr$Subject==p & data_agr$Session==s] <- mean(filter(cur_data,SourceAcc==0 & Source=="Pleasant")$SourceRT, na.rm = TRUE)
  }
}
# adjust the string variable(s)
data_agr$Sex[data_agr$Sex==1] <- levels(cur_data$Sex)[1]
data_agr$Sex[data_agr$Sex==2] <- levels(cur_data$Sex)[2]
# PRINT THE DESCRIPTIVES PER SESSION
if (do_stats_print) {
  for (s in unique(data_mat$Session)) {
    cat("\n\nXXXXXXXXXXXXXXXXXX\nXX SESSION: ",s," XX\nXXXXXXXXXXXXXXXXXX\n")
    cat("\nEncoding success (all): M = ",round(mean(filter(data_agr,Session==s)$EncResp)*100,2),", SD = ",round(sd(filter(data_agr,Session==s)$EncResp)*100,2))
    cat("\nEncoding success (Place): M = ",round(mean(filter(data_agr,Session==s)$EncRespPla)*100,2),", SD = ",round(sd(filter(data_agr,Session==s)$EncRespPla)*100,2))
    cat("\nEncoding success (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$EncRespPle)*100,2),", SD = ",round(sd(filter(data_agr,Session==s)$EncRespPle)*100,2))
    x <- filter(data_agr,Session==s)$EncRespPla
    y <- filter(data_agr,Session==s)$EncRespPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    
    cat("\n\nEncoding RT (all): M = ",round(mean(filter(data_agr,Session==s)$EncRT),0),", SD = ",round(sd(filter(data_agr,Session==s)$EncRT),0))
    cat("\nEncoding RT (Place): M = ",round(mean(filter(data_agr,Session==s)$EncRTPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$EncRTPla),0))
    cat("\nEncoding RT (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$EncRTPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$EncRTPle),0))
    x <- filter(data_agr,Session==s)$EncRTPla
    y <- filter(data_agr,Session==s)$EncRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    
    cat("\n\nRetrieval Item Hitrate (all): M = ",round(mean(filter(data_agr,Session==s)$ItemHitRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHitRate),3))
    cat("\nRetrieval Item Hitrate (Place): M = ",round(mean(filter(data_agr,Session==s)$ItemHitRatePla),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHitRatePla),3))
    cat("\nRetrieval Item Hitrate (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemHitRatePle),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHitRatePle),3))
    x <- filter(data_agr,Session==s)$ItemHitRatePla
    y <- filter(data_agr,Session==s)$ItemHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Item FArate (all): M = ",round(mean(filter(data_agr,Session==s)$ItemFARate),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemFARate),3))
    cat("\nRetrieval Item d-prime (all): M = ",round(mean(filter(data_agr,Session==s)$ItemDprime),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemDprime),3))
    cat("\nRetrieval Item d-prime (Place): M = ",round(mean(filter(data_agr,Session==s)$ItemDprimePla),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemDprimePla),3))
    cat("\nRetrieval Item d-prime (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemDprimePle),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemDprimePle),3))
    x <- filter(data_agr,Session==s)$ItemDprimePla
    y <- filter(data_agr,Session==s)$ItemDprimePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    
    cat("\n\nRetrieval Item HC rate (all): M = ",round(mean(filter(data_agr,Session==s)$ItemHCRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHCRate),3))
    cat("\nRetrieval Item HC rate (Hits): M = ",round(mean(filter(data_agr,Session==s)$ItemHCHitRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHCHitRate),3))
    cat("\nRetrieval Item HC rate (Hits, Place): M = ",round(mean(filter(data_agr,Session==s)$ItemHCHitRatePla),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHCHitRatePla),3))
    cat("\nRetrieval Item HC rate (Hits, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemHCHitRatePle),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHCHitRatePle),3))
    x <- filter(data_agr,Session==s)$ItemHCHitRatePla
    y <- filter(data_agr,Session==s)$ItemHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Item HC rate (CRs): M = ",round(mean(filter(data_agr,Session==s)$ItemHCCRRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$ItemHCCRRate),3))
    
    cat("\n\nRetrieval Item RT (all): M = ",round(mean(filter(data_agr,Session==s)$ItemRT),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRT),0))
    cat("\nRetrieval Item RT (Place): M = ",round(mean(filter(data_agr,Session==s)$ItemRTPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTPla),0))
    cat("\nRetrieval Item RT (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemRTPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTPle),0))
    x <- filter(data_agr,Session==s)$ItemRTPla
    y <- filter(data_agr,Session==s)$ItemRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Item RT (Correct): M = ",round(mean(filter(data_agr,Session==s)$ItemRTCor),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTCor),0))
    cat("\nRetrieval Item RT (Incorrect): M = ",round(mean(filter(data_agr,Session==s)$ItemRTIncor),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTIncor),0))
    x <- filter(data_agr,Session==s)$ItemRTCor
    y <- filter(data_agr,Session==s)$ItemRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nCorrect vs. Incorrect : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Item RT (Correct, Place): M = ",round(mean(filter(data_agr,Session==s)$ItemRTCorPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTCorPla),0))
    cat("\nRetrieval Item RT (Correct, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemRTCorPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTCorPle),0))
    x <- filter(data_agr,Session==s)$ItemRTCorPla
    y <- filter(data_agr,Session==s)$ItemRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Item RT (Incorrect, Place): M = ",round(mean(filter(data_agr,Session==s)$ItemRTIncorPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTIncorPla),0))
    cat("\nRetrieval Item RT (Incorrect, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$ItemRTIncorPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$ItemRTIncorPle),0))
    x <- filter(data_agr,Session==s)$ItemRTIncorPla
    y <- filter(data_agr,Session==s)$ItemRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    
    cat("\n\nRetrieval Source Hitrate (all): M = ",round(mean(filter(data_agr,Session==s)$SourceHitRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceHitRate),3))
    cat("\nRetrieval Source FArate (all): M = ",round(mean(filter(data_agr,Session==s)$SourceFARate),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceFARate),3))
    cat("\nRetrieval Source d-prime (all): M = ",round(mean(filter(data_agr,Session==s)$SourceDprime),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceDprime),3))
    
    cat("\n\nRetrieval Source HC rate (all): M = ",round(mean(filter(data_agr,Session==s)$SourceHCRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceHCRate),3))
    cat("\nRetrieval Source HC rate (Hits): M = ",round(mean(filter(data_agr,Session==s)$SourceHCHitRate),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceHCHitRate),3))
    cat("\nRetrieval Source HC rate (Hits, Place): M = ",round(mean(filter(data_agr,Session==s)$SourceHCHitRatePla),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceHCHitRatePla),3))
    cat("\nRetrieval Source HC rate (Hits, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$SourceHCHitRatePle),2),", SD = ",round(sd(filter(data_agr,Session==s)$SourceHCHitRatePle),3))
    x <- filter(data_agr,Session==s)$SourceHCHitRatePla
    y <- filter(data_agr,Session==s)$SourceHCHitRatePle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    
    cat("\n\nRetrieval Source RT (all): M = ",round(mean(filter(data_agr,Session==s)$SourceRT),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRT),0))
    cat("\nRetrieval Source RT (Place): M = ",round(mean(filter(data_agr,Session==s)$SourceRTPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTPla),0))
    cat("\nRetrieval Source RT (Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$SourceRTPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTPle),0))
    x <- filter(data_agr,Session==s)$SourceRTPla
    y <- filter(data_agr,Session==s)$SourceRTPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Source RT (Correct): M = ",round(mean(filter(data_agr,Session==s)$SourceRTCor),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTCor),0))
    cat("\nRetrieval Source RT (Incorrect): M = ",round(mean(filter(data_agr,Session==s)$SourceRTIncor),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTIncor),0))
    x <- filter(data_agr,Session==s)$SourceRTCor
    y <- filter(data_agr,Session==s)$SourceRTIncor
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nCorrect vs. Incorrect : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Source RT (Correct, Place): M = ",round(mean(filter(data_agr,Session==s)$SourceRTCorPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTCorPla),0))
    cat("\nRetrieval Source RT (Correct, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$SourceRTCorPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTCorPle),0))
    x <- filter(data_agr,Session==s)$SourceRTCorPla
    y <- filter(data_agr,Session==s)$SourceRTCorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
    cat("\nRetrieval Source RT (Incorrect, Place): M = ",round(mean(filter(data_agr,Session==s)$SourceRTIncorPla),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTIncorPla),0))
    cat("\nRetrieval Source RT (Incorrect, Pleasantness): M = ",round(mean(filter(data_agr,Session==s)$SourceRTIncorPle),0),", SD = ",round(sd(filter(data_agr,Session==s)$SourceRTIncorPle),0))
    x <- filter(data_agr,Session==s)$SourceRTIncorPla
    y <- filter(data_agr,Session==s)$SourceRTIncorPle
    res = t.test(x, y, paired = TRUE, alternative = "two.sided")
    D = cohen.d(x, y)
    cat("\nPlace vs. Pleasant : t(", res$parameter, ") = ", round(res$statistic,2),", p-value:", round(res$p.value,3), D$method, ":", round(D$estimate,3))
  }
}
#### EEG ANALYSIS: Behavior predicting brain ####
# https://cran.r-project.org/web/packages/emmeans/vignettes/models.html

### THETA ~ ITEM MEMORY ###
## MODEL ##
model.thei   <- lmer(s_Ret_Theta ~ s_BrainRegion*(s_OldNew*s_ONAcc+s_ONConf) + (1 | Subject), 
                           data = filter(data_mat_long,outlrs_powt==FALSE), control = lmerControl(optimizer = "bobyqa"))
model.thei.sum <- as.data.frame(round(coef(summary(model.thei)),3))
model.thei.sum$df <- round(model.thei.sum$df,0)
print(car::vif(model.thei))
# Let"s explore this model with post-hoc tests.
## POST-HOC TESTS ##
# Now let"s look at all (relevant) pairwise comparisons
em1   <- emmeans(model.thei, pairwise ~ (s_OldNew*s_ONAcc)|s_BrainRegion, lmer.df = "asymp", weights = "cells")
em2   <- emmeans(model.thei, pairwise ~ s_ONConf|s_BrainRegion, lmer.df = "asymp", weights = "cells")
emef1 <- eff_size(em1, sigma = sigma(model.thei), edf = df.residual(model.thei) ) # calculate the Cohen's D
emef2 <- eff_size(em2, sigma = sigma(model.thei), edf = df.residual(model.thei) ) # calculate the Cohen's D
em_means          <- bind_rows(summary(em1)$emmean,summary(em2)$emmean)
em_means          <- em_means[,c(1,9,2:8)] # Just reordering the columns
em_diffs          <- bind_rows(summary(em1)$contrasts,summary(em2)$contrasts)
em_diffs$p.value  <- round(em_diffs$p.value,4)
em_effsz          <- bind_rows(summary(emef1),summary(emef2))
em_means$cond <- NA
em_means$cond[em_means$s_OldNew=="Old" & em_means$s_ONAcc==1] <- "Hits"
em_means$cond[em_means$s_OldNew=="Old" & em_means$s_ONAcc==0] <- "Misses"
em_means$cond[em_means$s_OldNew=="New" & em_means$s_ONAcc==1] <- "CRs"
em_means$cond[em_means$s_OldNew=="New" & em_means$s_ONAcc==0] <- "FAs"
em_means$cond[em_means$s_ONConf==0]                           <- "LC"
em_means$cond[em_means$s_ONConf==1]                           <- "HC"
em_diffs$cond <- NA
em_diffs$cond[em_diffs$contrast=="New s_ONAcc1 - Old s_ONAcc1"] <- "CRs-Hits"
em_diffs$cond[em_diffs$contrast=="Old s_ONAcc0 - Old s_ONAcc1"] <- "Misses-Hits"
em_diffs$cond[em_diffs$contrast=="New s_ONAcc0 - New s_ONAcc1"] <- "FAs-CRs"
em_diffs$cond[em_diffs$contrast=="s_ONConf0 - s_ONConf1"]       <- "LC-HC"
em_diffs <- em_diffs[!is.na(em_diffs$cond),]
em_effsz$cond <- NA
em_effsz$cond[em_effsz$contrast=="(New s_ONAcc1 - Old s_ONAcc1)"] <- "CRs-Hits"
em_effsz$cond[em_effsz$contrast=="(Old s_ONAcc0 - Old s_ONAcc1)"] <- "Misses-Hits"
em_effsz$cond[em_effsz$contrast=="(New s_ONAcc0 - New s_ONAcc1)"] <- "FAs-CRs"
em_effsz$cond[em_effsz$contrast=="(s_ONConf0 - s_ONConf1)"]       <- "LC-HC"
em_effsz <- em_effsz[!is.na(em_effsz$cond),]
# store the estimated probabilities, the pairwise differences and the effect size (Cohen's d)
thei.means <- em_means
thei.diffs <- em_diffs
thei.effsz <- em_effsz
thei.pcomp <- data.frame(em_diffs$s_BrainRegion, em_diffs$cond, round(em_diffs$estimate,3), round(em_diffs$SE,3), em_diffs$df, round(em_diffs$z.ratio,3), round(em_diffs$p.value,3), round(em_effsz$effect.size,3))
thei.pcomp <- plyr::rename(thei.pcomp,
             c("em_diffs.s_BrainRegion" = "BrainRegion",
               "em_diffs.cond" = "Comparison",
               "round.em_diffs.estimate..3." = "Estimate",
               "round.em_diffs.SE..3." = "Std. Error",
               "em_diffs.df" = "df",
               "round.em_diffs.z.ratio..3." = "Z-ratio",
               "round.em_diffs.p.value..3." = "P-value",
               "round.em_effsz.effect.size..3." = "Cohen's D")
)
thei.pcomp <- thei.pcomp[order(c(3,2,1,7,6,5,4,8)),]
# get the estimated probabilities per condition, averaged over the brain regions
thei.Hits     <- round(mean(filter(thei.means,cond=="Hits")$emmean),2)
thei.Misses   <- round(mean(filter(thei.means,cond=="Misses")$emmean),2)
thei.CRs      <- round(mean(filter(thei.means,cond=="CRs")$emmean),2)
thei.FAs      <- round(mean(filter(thei.means,cond=="FAs")$emmean),2)
thei.LC       <- round(mean(filter(thei.means,cond=="LC")$emmean),2)
thei.HC       <- round(mean(filter(thei.means,cond=="HC")$emmean),2)
# get the predicted values out of the models
model.thei.preds        <- na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_OldNew,s_ONAcc,s_ONConf))
model.thei.preds        <- distinct(na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_OldNew,s_ONAcc,s_ONConf)))
model.thei.preds$MemSta <- NA
model.thei.preds$Conf   <- NA
model.thei.preds$MemSta[model.thei.preds$s_OldNew=="Old" & model.thei.preds$s_ONAcc==1] <- "Hits"
model.thei.preds$MemSta[model.thei.preds$s_OldNew=="Old" & model.thei.preds$s_ONAcc==0] <- "Misses"
model.thei.preds$MemSta[model.thei.preds$s_OldNew=="New" & model.thei.preds$s_ONAcc==1] <- "CRs"
model.thei.preds$MemSta[model.thei.preds$s_OldNew=="New" & model.thei.preds$s_ONAcc==0] <- "FAs"
model.thei.preds$Conf[model.thei.preds$s_ONConf==0]                                     <- "LC"
model.thei.preds$Conf[model.thei.preds$s_ONConf==1]                                     <- "HC"
model.thei.preds$fit    <- predict(model.thei,model.thei.preds)

### THETA ~ SOURCE MEMORY ###
## MODEL ##
model.thes   <- lmer(s_Ret_Theta ~ s_BrainRegion*(s_SourceAcc+s_SourceConf) + (1 | Subject), 
                    data = filter(data_mat_long,outlrs_powt==FALSE), control = lmerControl(optimizer = "bobyqa"))
model.thes.sum <- as.data.frame(round(coef(summary(model.thes)),3))
model.thes.sum$df <- round(model.thes.sum$df,0)
print(car::vif(model.thes))
# Let"s explore this model with post-hoc tests.
## POST-HOC TESTS ##
# Now let"s look at all (relevant) pairwise comparisons
em0   <- emmeans(model.thes, pairwise ~ s_SourceConf, lmer.df = "asymp", weights = "cells")
em1   <- emmeans(model.thes, pairwise ~ s_SourceAcc|s_BrainRegion, lmer.df = "asymp", weights = "cells")
em2   <- emmeans(model.thes, pairwise ~ s_SourceConf|s_BrainRegion, lmer.df = "asymp", weights = "cells")
emef1 <- eff_size(em1, sigma = sigma(model.thes), edf = df.residual(model.thes) ) # calculate the Cohen"s D
emef2 <- eff_size(em2, sigma = sigma(model.thes), edf = df.residual(model.thes) ) # calculate the Cohen"s D
em_means          <- bind_rows(summary(em1)$emmean,summary(em2)$emmean)
em_means          <- em_means[,c(1,8,2:7)] # Just reordering the columns
em_diffs          <- bind_rows(summary(em1)$contrasts,summary(em2)$contrasts)
em_diffs$p.value  <- round(em_diffs$p.value,4)
em_effsz          <- bind_rows(summary(emef1),summary(emef2))
em_means$cond <- NA
em_means$cond[em_means$s_SourceAcc==1]  <- "SourceHit"
em_means$cond[em_means$s_SourceAcc==0]  <- "SourceMiss"
em_means$cond[em_means$s_SourceConf==1] <- "HC"
em_means$cond[em_means$s_SourceConf==0] <- "LC"
em_diffs$cond <- NA
em_diffs$cond[em_diffs$contrast=="s_SourceAcc0 - s_SourceAcc1"]     <- "SourceMiss-SourceHit"
em_diffs$cond[em_diffs$contrast=="s_SourceConf0 - s_SourceConf1"]   <- "LC-HC"
em_effsz$cond <- NA
em_effsz$cond[em_effsz$contrast=="(s_SourceAcc0 - s_SourceAcc1)"]   <- "SourceMiss-SourceHit"
em_effsz$cond[em_effsz$contrast=="(s_SourceConf0 - s_SourceConf1)"] <- "LC-HC"
# store the estimated probabilities, the pairwise differences and the effect size (Cohen"s d)
thes.means <- em_means
thes.diffs <- em_diffs
thes.effsz <- em_effsz
thes.pcomp <- data.frame(em_diffs$s_BrainRegion, em_diffs$cond, round(em_diffs$estimate,3), round(em_diffs$SE,3), em_diffs$df, round(em_diffs$z.ratio,3), round(em_diffs$p.value,3), round(em_effsz$effect.size,3))
thes.pcomp <- plyr::rename(thes.pcomp,
                           c("em_diffs.s_BrainRegion" = "BrainRegion",
                             "em_diffs.cond" = "Comparison",
                             "round.em_diffs.estimate..3." = "Estimate",
                             "round.em_diffs.SE..3." = "Std. Error",
                             "em_diffs.df" = "df",
                             "round.em_diffs.z.ratio..3." = "Z-ratio",
                             "round.em_diffs.p.value..3." = "P-value",
                             "round.em_effsz.effect.size..3." = "Cohen's D")
)
thes.pcomp <- thes.pcomp[order(thes.pcomp$BrainRegion),]
# get the estimated probabilities per condition, averaged over the brain regions
thes.SourceHits     <- round(mean(filter(thes.means,cond=="SourceHit")$emmean),2)
thes.SourceMisses   <- round(mean(filter(thes.means,cond=="SourceMiss")$emmean),2)
thes.LC             <- round(mean(filter(thes.means,cond=="LC")$emmean),2)
thes.HC             <- round(mean(filter(thes.means,cond=="HC")$emmean),2)
# get the predicted values out of the models
model.thes.preds        <- na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_SourceAcc,s_SourceConf))
model.thes.preds        <- distinct(na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_SourceAcc,s_SourceConf)))
model.thes.preds$MemSta <- NA
model.thes.preds$Conf   <- NA
model.thes.preds$MemSta[model.thes.preds$s_SourceAcc==1]  <- "SourceHit"
model.thes.preds$MemSta[model.thes.preds$s_SourceAcc==0]  <- "SourceMiss"
model.thes.preds$Conf[model.thes.preds$s_SourceConf==0]   <- "LC"
model.thes.preds$Conf[model.thes.preds$s_SourceConf==1]   <- "HC"
model.thes.preds$fit    <- predict(model.thes,model.thes.preds)

### GAMMA ~ ITEM MEMORY ###
## MODEL ##
model.gami   <- lmer(s_Ret_Gamma ~ s_BrainRegion*(s_OldNew*s_ONAcc+s_ONConf) + (1 | Subject), 
                     data = filter(data_mat_long,outlrs_powg==FALSE), control = lmerControl(optimizer = "bobyqa"))
model.gami.sum <- as.data.frame(round(coef(summary(model.gami)),3))
model.gami.sum$df <- round(model.gami.sum$df,0)
print(car::vif(model.gami))
# Let"s explore this model with post-hoc tests.
## POST-HOC TESTS ##
# Now let"s look at all (relevant) pairwise comparisons
em1   <- emmeans(model.gami, pairwise ~ (s_OldNew*s_ONAcc)|s_BrainRegion, lmer.df = "asymp", weights = "cells")
em2   <- emmeans(model.gami, pairwise ~ s_ONConf|s_BrainRegion, lmer.df = "asymp", weights = "cells")
emef1 <- eff_size(em1, sigma = sigma(model.gami), edf = df.residual(model.gami) ) # calculate the Cohen"s D
emef2 <- eff_size(em2, sigma = sigma(model.gami), edf = df.residual(model.gami) ) # calculate the Cohen"s D
em_means          <- bind_rows(summary(em1)$emmean,summary(em2)$emmean)
em_means          <- em_means[,c(1,9,2:8)] # Just reordering the columns
em_diffs          <- bind_rows(summary(em1)$contrasts,summary(em2)$contrasts)
em_diffs$p.value  <- round(em_diffs$p.value,4)
em_effsz          <- bind_rows(summary(emef1),summary(emef2))
em_means$cond <- NA
em_means$cond[em_means$s_OldNew=="Old" & em_means$s_ONAcc==1] <- "Hits"
em_means$cond[em_means$s_OldNew=="Old" & em_means$s_ONAcc==0] <- "Misses"
em_means$cond[em_means$s_OldNew=="New" & em_means$s_ONAcc==1] <- "CRs"
em_means$cond[em_means$s_OldNew=="New" & em_means$s_ONAcc==0] <- "FAs"
em_means$cond[em_means$s_ONConf==0]                           <- "LC"
em_means$cond[em_means$s_ONConf==1]                           <- "HC"
em_diffs$cond <- NA
em_diffs$cond[em_diffs$contrast=="New s_ONAcc1 - Old s_ONAcc1"] <- "CRs-Hits"
em_diffs$cond[em_diffs$contrast=="Old s_ONAcc0 - Old s_ONAcc1"] <- "Misses-Hits"
em_diffs$cond[em_diffs$contrast=="New s_ONAcc0 - New s_ONAcc1"] <- "FAs-CRs"
em_diffs$cond[em_diffs$contrast=="s_ONConf0 - s_ONConf1"]       <- "LC-HC"
em_diffs <- em_diffs[!is.na(em_diffs$cond),]
em_effsz$cond <- NA
em_effsz$cond[em_effsz$contrast=="(New s_ONAcc1 - Old s_ONAcc1)"] <- "CRs-Hits"
em_effsz$cond[em_effsz$contrast=="(Old s_ONAcc0 - Old s_ONAcc1)"] <- "Misses-Hits"
em_effsz$cond[em_effsz$contrast=="(New s_ONAcc0 - New s_ONAcc1)"] <- "FAs-CRs"
em_effsz$cond[em_effsz$contrast=="(s_ONConf0 - s_ONConf1)"]       <- "LC-HC"
em_effsz <- em_effsz[!is.na(em_effsz$cond),]
# store the estimated probabilities, the pairwise differences and the effect size (Cohen"s d)
gami.means <- em_means
gami.diffs <- em_diffs
gami.effsz <- em_effsz
gami.pcomp <- data.frame(em_diffs$s_BrainRegion, em_diffs$cond, round(em_diffs$estimate,3), round(em_diffs$SE,3), em_diffs$df, round(em_diffs$z.ratio,3), round(em_diffs$p.value,3), round(em_effsz$effect.size,3))
gami.pcomp <- plyr::rename(gami.pcomp,
                           c("em_diffs.s_BrainRegion" = "BrainRegion",
                             "em_diffs.cond" = "Comparison",
                             "round.em_diffs.estimate..3." = "Estimate",
                             "round.em_diffs.SE..3." = "Std. Error",
                             "em_diffs.df" = "df",
                             "round.em_diffs.z.ratio..3." = "Z-ratio",
                             "round.em_diffs.p.value..3." = "P-value",
                             "round.em_effsz.effect.size..3." = "Cohen's D")
)
gami.pcomp <- gami.pcomp[order(c(3,2,1,7,6,5,4,8)),]
# get the estimated probabilities per condition, averaged over the brain regions
gami.Hits     <- round(mean(filter(gami.means,cond=="Hits")$emmean),2)
gami.Misses   <- round(mean(filter(gami.means,cond=="Misses")$emmean),2)
gami.CRs      <- round(mean(filter(gami.means,cond=="CRs")$emmean),2)
gami.FAs      <- round(mean(filter(gami.means,cond=="FAs")$emmean),2)
gami.LC       <- round(mean(filter(gami.means,cond=="LC")$emmean),2)
gami.HC       <- round(mean(filter(gami.means,cond=="HC")$emmean),2)
# get the predicted values out of the models
model.gami.preds        <- na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_OldNew,s_ONAcc,s_ONConf))
model.gami.preds        <- distinct(na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_OldNew,s_ONAcc,s_ONConf)))
model.gami.preds$MemSta <- NA
model.gami.preds$Conf   <- NA
model.gami.preds$MemSta[model.gami.preds$s_OldNew=="Old" & model.gami.preds$s_ONAcc==1] <- "Hits"
model.gami.preds$MemSta[model.gami.preds$s_OldNew=="Old" & model.gami.preds$s_ONAcc==0] <- "Misses"
model.gami.preds$MemSta[model.gami.preds$s_OldNew=="New" & model.gami.preds$s_ONAcc==1] <- "CRs"
model.gami.preds$MemSta[model.gami.preds$s_OldNew=="New" & model.gami.preds$s_ONAcc==0] <- "FAs"
model.gami.preds$Conf[model.gami.preds$s_ONConf==0]                                     <- "LC"
model.gami.preds$Conf[model.gami.preds$s_ONConf==1]                                     <- "HC"
model.gami.preds$fit    <- predict(model.gami,model.gami.preds)

### GAMMA ~ SOURCE MEMORY ###
## MODEL ##
model.gams   <- lmer(s_Ret_Gamma ~ s_BrainRegion*(s_SourceAcc+s_SourceConf) + (1 | Subject), 
                     data = filter(data_mat_long,outlrs_powg==FALSE), control = lmerControl(optimizer = "bobyqa"))
model.gams.sum <- as.data.frame(round(coef(summary(model.gams)),3))
model.gams.sum$df <- round(model.gams.sum$df,0)
print(car::vif(model.gams))
# Let"s explore this model with post-hoc tests.
## POST-HOC TESTS ##
# Now let"s look at all (relevant) pairwise comparisons
em0   <- emmeans(model.gams, pairwise ~ s_SourceAcc, lmer.df = "asymp", weights = "cells")
em1   <- emmeans(model.gams, pairwise ~ s_SourceAcc|s_BrainRegion, lmer.df = "asymp", weights = "cells")
em2   <- emmeans(model.gams, pairwise ~ s_SourceConf|s_BrainRegion, lmer.df = "asymp", weights = "cells")
emef1 <- eff_size(em1, sigma = sigma(model.gams), edf = df.residual(model.gams) ) # calculate the Cohen"s D
emef2 <- eff_size(em2, sigma = sigma(model.gams), edf = df.residual(model.gams) ) # calculate the Cohen"s D
em_means          <- bind_rows(summary(em1)$emmean,summary(em2)$emmean)
em_means          <- em_means[,c(1,8,2:7)] # Just reordering the columns
em_diffs          <- bind_rows(summary(em1)$contrasts,summary(em2)$contrasts)
em_diffs$p.value  <- round(em_diffs$p.value,4)
em_effsz          <- bind_rows(summary(emef1),summary(emef2))
em_means$cond <- NA
em_means$cond[em_means$s_SourceAcc==1]  <- "SourceHit"
em_means$cond[em_means$s_SourceAcc==0]  <- "SourceMiss"
em_means$cond[em_means$s_SourceConf==1] <- "HC"
em_means$cond[em_means$s_SourceConf==0] <- "LC"
em_diffs$cond <- NA
em_diffs$cond[em_diffs$contrast=="s_SourceAcc0 - s_SourceAcc1"]     <- "SourceMiss-SourceHit"
em_diffs$cond[em_diffs$contrast=="s_SourceConf0 - s_SourceConf1"]   <- "LC-HC"
em_effsz$cond <- NA
em_effsz$cond[em_effsz$contrast=="(s_SourceAcc0 - s_SourceAcc1)"]   <- "SourceMiss-SourceHit"
em_effsz$cond[em_effsz$contrast=="(s_SourceConf0 - s_SourceConf1)"] <- "LC-HC"
# store the estimated probabilities, the pairwise differences and the effect size (Cohen"s d)
gams.means <- em_means
gams.diffs <- em_diffs
gams.effsz <- em_effsz
gams.pcomp <- data.frame(em_diffs$s_BrainRegion, em_diffs$cond, round(em_diffs$estimate,3), round(em_diffs$SE,3), em_diffs$df, round(em_diffs$z.ratio,3), round(em_diffs$p.value,3), round(em_effsz$effect.size,3))
gams.pcomp <- plyr::rename(gams.pcomp,
                           c("em_diffs.s_BrainRegion" = "BrainRegion",
                             "em_diffs.cond" = "Comparison",
                             "round.em_diffs.estimate..3." = "Estimate",
                             "round.em_diffs.SE..3." = "Std. Error",
                             "em_diffs.df" = "df",
                             "round.em_diffs.z.ratio..3." = "Z-ratio",
                             "round.em_diffs.p.value..3." = "P-value",
                             "round.em_effsz.effect.size..3." = "Cohen's D")
)
gams.pcomp <- gams.pcomp[order(gams.pcomp$BrainRegion),]
# get the estimated probabilities per condition, averaged over the brain regions
gams.SourceHits     <- round(mean(filter(gams.means,cond=="SourceHit")$emmean),2)
gams.SourceMisses   <- round(mean(filter(gams.means,cond=="SourceMiss")$emmean),2)
gams.LC             <- round(mean(filter(gams.means,cond=="LC")$emmean),2)
gams.HC             <- round(mean(filter(gams.means,cond=="HC")$emmean),2)
# get the predicted values out of the models
model.gams.preds        <- na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_SourceAcc,s_SourceConf))
model.gams.preds        <- distinct(na.omit(select(filter(data_mat_long,Session==1),Subject,s_BrainRegion,s_SourceAcc,s_SourceConf)))
model.gams.preds$MemSta <- NA
model.gams.preds$Conf   <- NA
model.gams.preds$MemSta[model.gams.preds$s_SourceAcc==1]  <- "SourceHit"
model.gams.preds$MemSta[model.gams.preds$s_SourceAcc==0]  <- "SourceMiss"
model.gams.preds$Conf[model.gams.preds$s_SourceConf==0]   <- "LC"
model.gams.preds$Conf[model.gams.preds$s_SourceConf==1]   <- "HC"
model.gams.preds$fit    <- predict(model.gams,model.gams.preds)

#### PLOTTING: EEG RESULTS ####
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# https://statisticsglobe.com/add-image-to-plot-in-r
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
# https://stackoverflow.com/questions/56234806/how-to-reduce-the-legend-symbol-thickness-for-a-bar-chart-in-ggplot2
ylim <- c(-.05,.2)
### THETA ~ ITEM MEMORY ###
plotname <- "pred_thet_item"
## Memory accuracy ##
data <- thei.means
levels <- c("Hits","Misses","CRs","FAs")
pltcolors <- cp1
xlab <- "Memory accuracy"
ylab <- "Predicted standardized theta power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_frontal_iHit-CR.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_parietal_iHit-CR.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_theta_iHit-CR.png"
brn_tit1 = "Frontal\n Hits - Correct rejections"
brn_tit2 = "Parietal\n Hits - Correct rejections"
brn_tit3 = "Hits - Correct rejections"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p1.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .9) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p1.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55,vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p1.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55,vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p1.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 8),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Memory confidence ##
data <- thei.means
levels <- c("HC","LC")
pltcolors <- cp2
xlab <- "Memory confidence"
ylab <- "Predicted standardized theta power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_frontal_iHC-iLC.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_parietal_iHC-iLC.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_theta_iHC-iLC.png"
brn_tit1 = "Frontal\n High-confidence - Low-confidence"
brn_tit2 = "Parietal\n High-confidence - Low-confidence"
brn_tit3 = "High-confidence - Low-confidence"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p2.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p2.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55,vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p2.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p2.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Combine ##
ggdraw() +
  draw_plot(p1.1, x = 0,     y = .66, width = .48, height = .33) +
  draw_plot(p2.1, x = .52,   y = .66, width = .48, height = .33) +
  draw_plot(p1.2, x = 0.02,  y = .30, width = .25, height = .30) +
  draw_plot(p1.3, x = 0.245, y = .30, width = .25, height = .30) +
  draw_plot(p2.2, x = 0.54,  y = .30, width = .25, height = .30) +
  draw_plot(p2.3, x = 0.765, y = .30, width = .25, height = .30) +
  draw_plot(p1.4, x = 0.14,  y = 0,   width = .25, height = .27) +
  draw_plot(p2.4, x = 0.66,  y = 0,   width = .25, height = .27) +
  draw_plot_label(label = c("A" , "B" , "C" ,  "D" , "E" , "F",  "G",  "H"), size = 15,
                  x = c(0.00, 0.52, 0.04, 0.26, 0.56, 0.78, 0.17,  0.69), 
                  y = c(1.00, 1.00, 0.62, 0.62, 0.62, 0.62, 0.290, 0.290))
## Save the plots ##
ggsave(paste(outp_path, plotname, ".pdf",sep=""), width = 4000, height = 3000, units = "px")
ggsave(paste(outp_path, plotname, ".png",sep=""), width = 4000, height = 3000, units = "px")

### THETA ~ SOURCE MEMORY ###
plotname <- "pred_thet_source"
## Memory accuracy ##
data <- thes.means
levels <- c("SourceHit","SourceMiss")
pltcolors <- cp1
xlab <- "Memory accuracy"
ylab <- "Predicted standardized theta power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_frontal_sHit-sMiss.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_parietal_sHit-sMiss.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_theta_sHit-sMiss.png"
brn_tit1 = "Frontal\n Source hits - Source misses"
brn_tit2 = "Parietal\n Source hits - Source misses"
brn_tit3 = "Source hits - Source misses"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p1.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p1.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9,hjust = 0.55 ,vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p1.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9,hjust = 0.55 ,vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p1.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Memory confidence ##
data <- thes.means
levels <- c("HC","LC")
pltcolors <- cp2
xlab <- "Memory confidence"
ylab <- "Predicted standardized theta power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_frontal_sHC-sLC.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_theta_parietal_sHC-sLC.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_theta_sHC-sLC.png"
brn_tit1 = "Frontal\n High-confidence - Low-confidence"
brn_tit2 = "Parietal\n High-confidence - Low-confidence"
brn_tit3 = "High-confidence - Low-confidence"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p2.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p2.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p2.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p2.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Combine ##
ggdraw() +
  draw_plot(p1.1, x = 0,     y = .66, width = .48, height = .33) +
  draw_plot(p2.1, x = .52,   y = .66, width = .48, height = .33) +
  draw_plot(p1.2, x = 0.02,  y = .30, width = .25, height = .30) +
  draw_plot(p1.3, x = 0.245, y = .30, width = .25, height = .30) +
  draw_plot(p2.2, x = 0.54,  y = .30, width = .25, height = .30) +
  draw_plot(p2.3, x = 0.765, y = .30, width = .25, height = .30) +
  draw_plot(p1.4, x = 0.14,  y = 0,   width = .25, height = .27) +
  draw_plot(p2.4, x = 0.66,  y = 0,   width = .25, height = .27) +
  draw_plot_label(label = c("A" , "B" , "C" ,  "D" , "E" , "F",  "G",  "H"), size = 15,
                  x = c(0.00, 0.52, 0.04, 0.26, 0.56, 0.78, 0.17,  0.69), 
                  y = c(1.00, 1.00, 0.62, 0.62, 0.62, 0.62, 0.290, 0.290))
## Save the plots ##
ggsave(paste(outp_path, plotname, ".pdf",sep=""), width = 4000, height = 3000, units = "px")
ggsave(paste(outp_path, plotname, ".png",sep=""), width = 4000, height = 3000, units = "px")

### GAMMA ~ ITEM MEMORY ###
plotname <- "pred_gamm_item"
## Memory accuracy ##
data <- gami.means
levels <- c("Hits","Misses","CRs","FAs")
pltcolors <- cp1
xlab <- "Memory accuracy"
ylab <- "Predicted standardized gamma power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_frontal_CR-FA.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_parietal_CR-FA.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_gamma_CR-FA.png"
brn_tit1 = "Frontal\n Correct rejections - False alarms"
brn_tit2 = "Parietal\n Correct rejections - False alarms"
brn_tit3 = "Correct rejections - False alarms"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p1.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .9) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p1.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
        brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p1.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p1.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Memory confidence ##
data <- gami.means
levels <- c("HC","LC")
pltcolors <- cp2
xlab <- "Memory confidence"
ylab <- "Predicted standardized gamma power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_frontal_iHC-iLC.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_parietal_iHC-iLC.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_gamma_iHC-iLC.png"
brn_tit1 = "Frontal\n High-confidence - Low-confidence"
brn_tit2 = "Parietal\n High-confidence - Low-confidence"
brn_tit3 = "High-confidence - Low-confidence"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p2.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p2.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p2.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p2.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Combine ##
ggdraw() +
  draw_plot(p1.1, x = 0,     y = .66, width = .48, height = .33) +
  draw_plot(p2.1, x = .52,   y = .66, width = .48, height = .33) +
  draw_plot(p1.2, x = 0.02,  y = .30, width = .25, height = .30) +
  draw_plot(p1.3, x = 0.245, y = .30, width = .25, height = .30) +
  draw_plot(p2.2, x = 0.54,  y = .30, width = .25, height = .30) +
  draw_plot(p2.3, x = 0.765, y = .30, width = .25, height = .30) +
  draw_plot(p1.4, x = 0.14,  y = 0,   width = .25, height = .27) +
  draw_plot(p2.4, x = 0.66,  y = 0,   width = .25, height = .27) +
  draw_plot_label(label = c("A" , "B" , "C" ,  "D" , "E" , "F",  "G",  "H"), size = 15,
                  x = c(0.00, 0.52, 0.04, 0.26, 0.56, 0.78, 0.17,  0.69), 
                  y = c(1.00, 1.00, 0.62, 0.62, 0.62, 0.62, 0.290, 0.290))
## Save the plots ##
ggsave(paste(outp_path, plotname, ".pdf",sep=""), width = 4000, height = 3000, units = "px")
ggsave(paste(outp_path, plotname, ".png",sep=""), width = 4000, height = 3000, units = "px")


### GAMMA ~ SOURCE MEMORY ###
plotname <- "pred_gamm_source"
## Memory accuracy ##
data <- gams.means
levels <- c("SourceHit","SourceMiss")
pltcolors <- cp1
xlab <- "Memory accuracy"
ylab <- "Predicted standardized gamma power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_frontal_sHit-sMiss.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_parietal_sHit-sMiss.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_gamma_sHit-sMiss.png"
brn_tit1 = "Frontal\n Source hits - Source misses"
brn_tit2 = "Parietal\n Source hits - Source misses"
brn_tit3 = "Source hits - Source misses"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p1.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p1.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p1.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p1.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Memory confidence ##
data <- gams.means
levels <- c("HC","LC")
pltcolors <- cp2
xlab <- "Memory confidence"
ylab <- "Predicted standardized gamma power"
brainplotname1 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_frontal_sHC-sLC.png"
brainplotname2 <- "Session1_Ret_TF_DiffSing_GrandAverage_tvals_RetMem_gamma_parietal_sHC-sLC.png"
brainplotname3 <- "Session1_Ret_TF_DiffTopo_GrandAverage_tvals_RetMem_gamma_sHC-sLC.png"
brn_tit1 = "Frontal\n High-confidence - Low-confidence"
brn_tit2 = "Parietal\n High-confidence - Low-confidence"
brn_tit3 = "High-confidence - Low-confidence"
df <- filter(data, cond %in% levels)
df$cond <- factor(df$cond, levels = levels)
p2.1 <- ggplot(df, aes(x=cond, y=emmean, fill=cond, colour=cond)) + 
  geom_col(show.legend = FALSE, width = .45) +
  geom_point(size = 1, shape=15) + # just for the legend
  geom_hline(yintercept=0, alpha=.5) +
  xlab(xlab) +
  ylim(ylim) +
  ylab(ylab) +
  scale_x_discrete(position = "top") +
  facet_grid(~s_BrainRegion,scales="free") + 
  labs(labels=levels) +
  scale_fill_manual(values = pltcolors) +
  scale_color_manual(values = pltcolors) +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(size = 3) )) +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = c(.5,0),
        legend.title=element_blank(),
        legend.text = element_text(size=7))
brn_img = readPNG(paste(brpl_path,brainplotname1,sep=""), native = TRUE)
p2.2 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit1) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname2,sep=""), native = TRUE)
p2.3 <- ggplot() + 
          geom_blank() +
          ggtitle(brn_tit2) +
          theme(plot.title = element_text(size = 9, hjust = 0.55, vjust = 8),
                panel.background = element_blank()) +
          inset_element(brn_img,
                        left = 0,
                        bottom = 0,
                        right = 1,
                        top = .93,
                        on_top = TRUE,
                        align_to = "full")
brn_img = readPNG(paste(brpl_path,brainplotname3,sep=""), native = TRUE)
p2.4 <- ggplot() + 
  geom_blank() +
  ggtitle(brn_tit3) +
  theme(plot.title = element_text(size = 9, hjust = 0.5,vjust = 4),
        panel.background = element_blank()) +
  inset_element(brn_img,
                left = 0,
                bottom = 0,
                right = 1,
                top = .93,
                on_top = TRUE,
                align_to = "full")
## Combine ##
ggdraw() +
  draw_plot(p1.1, x = 0,     y = .66, width = .48, height = .33) +
  draw_plot(p2.1, x = .52,   y = .66, width = .48, height = .33) +
  draw_plot(p1.2, x = 0.02,  y = .30, width = .25, height = .30) +
  draw_plot(p1.3, x = 0.245, y = .30, width = .25, height = .30) +
  draw_plot(p2.2, x = 0.54,  y = .30, width = .25, height = .30) +
  draw_plot(p2.3, x = 0.765, y = .30, width = .25, height = .30) +
  draw_plot(p1.4, x = 0.14,  y = 0,   width = .25, height = .27) +
  draw_plot(p2.4, x = 0.66,  y = 0,   width = .25, height = .27) +
  draw_plot_label(label = c("A" , "B" , "C" ,  "D" , "E" , "F",  "G",  "H"), size = 15,
                      x = c(0.00, 0.52, 0.04, 0.26, 0.56, 0.78, 0.17,  0.69), 
                      y = c(1.00, 1.00, 0.62, 0.62, 0.62, 0.62, 0.290, 0.290))
## Save the plots ##
ggsave(paste(outp_path, plotname, ".pdf",sep=""), width = 4000, height = 3000, units = "px")
ggsave(paste(outp_path, plotname, ".png",sep=""), width = 4000, height = 3000, units = "px")

