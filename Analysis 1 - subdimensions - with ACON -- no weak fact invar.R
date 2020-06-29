
# Install packages and load libraries #
if (!require("tidyverse")) install.packages("tidyverse") # install this package first (once)
if (!require("lavaan")) install.packages("lavaan") # install this package first (once)
library(tidyverse)
library(lavaan)
#
if (!require("restriktor")) install.packages("restriktor") # install this package first (once)
library(restriktor) # for goric function


# Read in data file #
data <- read.table("04-TwoTimepoints_CLPA_12032021.dat", header = T)
colnames(data)
colnames(data)[1] <- "ID"
data <- replace(data , data == -999.00, NA)

data_SS <- select(data, 
                  THT1_SS,
                  TBT1_SS,
                  ACOMT1_SS,
                  SATT1_SS,
                  ACONT1_SS,
                  AB_T1_SS,
                  DE_T1_SS,
                  VI_T1_SS,
                  SLT1_SS,
                  TH_T2_SS,
                  TB_T2_SS,
                  ACOMT2SS,
                  SAT_T2SS,
                  ACONT2SS,
                  ABT2_SS,
                  DET2_SS,
                  VIT2_SS,
                  SLT2SS) # use tidyverse



## Fitting a CLPM ##

# Split dimension RQ into two subdimensions #

# Model 1: Configural Invariance 
CLPM_M1 <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ THT1_SS  + TBT1_SS  
  RQ12 =~ TH_T2_SS + TB_T2_SS
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ ACOMT1_SS + SATT1_SS + ACONT1_SS
  RQ22 =~ ACOMT2SS  + SAT_T2SS + ACONT2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ AB_T1_SS + DE_T1_SS + VI_T1_SS
  SE2 =~ ABT2_SS  + DET2_SS  + VIT2_SS

  
  ############
  # DYNAMICS #
  ############
  
  # Specify the lagged effects between the latent variables.
  RQ12 + RQ22 + SE2 + SLT2SS ~ RQ11 + RQ21 + SE1 + SLT1_SS
  
  # Estimate the correlations within the same wave.
  # T1
  RQ11 ~~ RQ21 + SE1 + SLT1_SS
  RQ21 ~~ SE1 + SLT1_SS
  SE1 ~~ SLT1_SS
  # T2
  RQ12 ~~ RQ22 + SE2 + SLT2SS
  RQ22 ~~ SE2 + SLT2SS
  SE2 ~~ SLT2SS

'
CLPM_M1.fit <- sem(CLPM_M1, data = data_SS, missing = 'ML')
#Warning message:
#  In lav_object_post_check(object) :
#  lavaan WARNING: covariance matrix of latent variables
#is not positive definite;
#use lavInspect(fit, "cov.lv") to investigate.
lavInspect(CLPM_M1.fit, "cov.lv")
lavInspect(CLPM_M1.fit, "cor.lv")
# Correlations between RQ11 & RQ21 and between RQ22 and RQ12 is very high, which makes sense since these two subdimensions belong to one dimension.


summary(CLPM_M1.fit, standardized = T, fit.measures=TRUE)[1]$FIT[c("chisq","df")]
#chisq      df 
#865.4532 109.0000 



# Model 2: weak factorial invariance 
CLPM_M2 <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ L1 * THT1_SS  + L2 * TBT1_SS  
  RQ12 =~ L1 * TH_T2_SS + L2 * TB_T2_SS
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ L3 * ACOMT1_SS + L4 * SATT1_SS  + L5 * ACONT1_SS
  RQ22 =~ L3 * ACOMT2SS  + L4 * SAT_T2SS  + L5 * ACONT2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ L6 * AB_T1_SS + L7 * DE_T1_SS + L8 * VI_T1_SS
  SE2 =~ L6 * ABT2_SS  + L7 * DET2_SS  + L8 * VIT2_SS

  
  ############
  # DYNAMICS #
  ############
  
  # Specify the lagged effects between the latent variables.
  RQ12 + RQ22 + SE2 + SLT2SS ~ RQ11 + RQ21 + SE1 + SLT1_SS
  
  # Estimate the correlations within the same wave.
  # T1
  RQ11 ~~ RQ21 + SE1 + SLT1_SS
  RQ21 ~~ SE1 + SLT1_SS
  SE1 ~~ SLT1_SS
  # T2
  RQ12 ~~ RQ22 + SE2 + SLT2SS
  RQ22 ~~ SE2 + SLT2SS
  SE2 ~~ SLT2SS

'
CLPM_M2.fit <- sem(CLPM_M2, data = data_SS, missing = 'ML')
#Warning message:
#  In lav_object_post_check(object) :
#  lavaan WARNING: covariance matrix of latent variables
#is not positive definite;
#use lavInspect(fit, "cov.lv") to investigate.
lavInspect(CLPM_M2.fit, "cov.lv")
lavInspect(CLPM_M2.fit, "cor.lv")
# Correlations between RQ11 & RQ21 and betweem RQ22 and RQ12 is very high, which makes sense since these two subdimensions belong to one dimension.


summary(CLPM_M2.fit, standardized = T, fit.measures=TRUE)[1]$FIT[c("chisq","df")]
#chisq      df 
#885.5345 114.0000 

# Chi-square difference test
#Df = 114 - 109 = 5
#Check on number of constrained factor loadings (note, first is set to 1):
#1 + 2 + 2 = 5
#Two models nested, so chi-square difference test with this df.
#Diff = 885.5345 - 865.4532 = 20.0813
#https://www.socscistatistics.com/pvalues/chidistribution.aspx
#The P-Value is .001207. The result is significant at p < .05.
#
#When the chi-square test is nonsignificant, this implies the factor loadings are not significantly different from each other over time. In other words, we can assume weak factorial invariance holds.
#If however the test is significant, this implies that the factor loadings cannot be constrained over time. This makes further comparisons between the latent variables very problematic or even impossible.
#
#Hence, we cannot assume weak factorial invariance… 
# When looking at the estimates of the factor loadings from Model 1, 
summary(CLPM_M1.fit, standardized = T, fit.measures=TRUE)
#it seems that the factor loading of ACON changes over time.
# Therefore, we do the analysis once more (in another file) without ACON.
