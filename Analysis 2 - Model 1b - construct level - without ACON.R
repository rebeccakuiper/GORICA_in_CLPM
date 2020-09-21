
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
                  #ACONT1_SS,
                  AB_T1_SS,
                  DE_T1_SS,
                  VI_T1_SS,
                  SLT1_SS,
                  TH_T2_SS,
                  TB_T2_SS,
                  ACOMT2SS,
                  SAT_T2SS,
                  #ACONT2SS,
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
  #RQ1  =~ 1 * RQ11 + 1 * RQ12
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ ACOMT1_SS + SATT1_SS
  RQ22 =~ ACOMT2SS  + SAT_T2SS
  #RQ2  =~ 1 * RQ21 + 1 * RQ22
  
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
#715.867  78.000 



# Model 2: weak factorial invariance 
CLPM_M2 <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ L1 * THT1_SS  + L2 * TBT1_SS  
  RQ12 =~ L1 * TH_T2_SS + L2 * TB_T2_SS 
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ L3 * ACOMT1_SS + L4 * SATT1_SS
  RQ22 =~ L3 * ACOMT2SS  + L4 * SAT_T2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ L5 * AB_T1_SS + L6 * DE_T1_SS + L7 * VI_T1_SS
  SE2 =~ L5 * ABT2_SS  + L6 * DET2_SS  + L7 * VIT2_SS

  
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
# Correlations between RQ11 & RQ21 and between RQ22 and RQ12 is very high, which makes sense since these two subdimensions belong to one dimension.


summary(CLPM_M2.fit, standardized = T, fit.measures=TRUE)[1]$FIT[c("chisq","df")]
#chisq      df 
#721.021  82.000  

# Chi-square difference test
#Df = 82 - 78 = 4
#Check on number of constrained factor loadings (note, first is set to 1):
#1 + 1 + 2 = 4
#Two models nested, so chi-square difference test with this df.
#Diff = 721.021 - 715.867 = 5.154
#https://www.socscistatistics.com/pvalues/chidistribution.aspx
#The P-Value is .271858. The result is not significant at p < .05.
# 1-pchisq(5.154, df = 4)
##[1] 0.271858
#
#When the chi-square test is nonsignificant, this implies the factor loadings are not significantly different from each other over time. In other words, we can assume weak factorial invariance holds.


# Model 3: strong factorial invariance 
CLPM_M3 <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ L1 * THT1_SS  + L2 * TBT1_SS  
  RQ12 =~ L1 * TH_T2_SS + L2 * TB_T2_SS 
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ L3 * ACOMT1_SS + L4 * SATT1_SS
  RQ22 =~ L3 * ACOMT2SS  + L4 * SAT_T2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ L5 * AB_T1_SS + L6 * DE_T1_SS + L7 * VI_T1_SS
  SE2 =~ L5 * ABT2_SS  + L6 * DET2_SS  + L7 * VIT2_SS


  # Constrained intercepts over time
  THT1_SS ~ int_th*1 
  TH_T2_SS ~ int_th*1
  TBT1_SS ~ int_tb*1 
  TB_T2_SS ~ int_tb*1 
  ACOMT1_SS ~ int_acom*1 
  ACOMT2SS ~ int_acom*1 
  SATT1_SS ~ int_sat*1 
  SAT_T2SS ~ int_sat*1 
  #
  AB_T1_SS ~ int_ab*1 
  ABT2_SS  ~ int_ab*1 
  DE_T1_SS ~ int_de*1 
  DET2_SS ~ int_de*1 
  VI_T1_SS ~ int_vi*1 
  VIT2_SS ~ int_vi*1 
  #
  SLT1_SS ~ int_sl*1 
  SLT2SS ~ int_sl*1 
  
  
  # Free latent means on t=2
  RQ12 + RQ22 + SE2 + RQ11 + RQ21 + SE1 ~ 1

  
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
CLPM_M3.fit <- sem(CLPM_M3, data = data_SS, missing = 'ML')
#Warning message:
#  In lav_object_post_check(object) :
#  lavaan WARNING: covariance matrix of latent variables
#is not positive definite;
#use lavInspect(fit, "cov.lv") to investigate.
lavInspect(CLPM_M3.fit, "cov.lv")
lavInspect(CLPM_M3.fit, "cor.lv")
# Correlations between RQ11 & RQ21 and between RQ22 and RQ12 is very high, which makes sense since these two subdimensions belong to one dimension.


summary(CLPM_M3.fit, standardized = T, fit.measures=TRUE)[1]$FIT[c("chisq","df")]
#chisq      df 
#725.4913  84.0000  

# Chi-square difference test
#Df = 84 - 82 = 2
#Check on number of constrained parameters (intercepts and substracting the means):
#8 - 6 = 2
#Two models nested, so chi-square difference test with this df.
#Diff = 725.4913 - 721.021 = 4.4703
#https://www.socscistatistics.com/pvalues/chidistribution.aspx
#The P-Value is .106976. The result is not significant at p < .05.
# 1-pchisq(4.4703, df = 2)
#
#If this chi-square difference test is nonsignicant, this means we can assume that strong factorial invariance holds over time. 
#In that case we could consider investigating whether the means change over time. This is just optional. 


# Model 4: strong factorial invariance, without free latent means (so, constraining them over time)
CLPM_M4 <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ L1 * THT1_SS  + L2 * TBT1_SS  
  RQ12 =~ L1 * TH_T2_SS + L2 * TB_T2_SS 
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ L3 * ACOMT1_SS + L4 * SATT1_SS
  RQ22 =~ L3 * ACOMT2SS  + L4 * SAT_T2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ L5 * AB_T1_SS + L6 * DE_T1_SS + L7 * VI_T1_SS
  SE2 =~ L5 * ABT2_SS  + L6 * DET2_SS  + L7 * VIT2_SS


  # Constrained intercepts over time
  THT1_SS ~ int_th*1 
  TH_T2_SS ~ int_th*1
  TBT1_SS ~ int_tb*1 
  TB_T2_SS ~ int_tb*1 
  ACOMT1_SS ~ int_acom*1 
  ACOMT2SS ~ int_acom*1 
  SATT1_SS ~ int_sat*1 
  SAT_T2SS ~ int_sat*1 
  #
  AB_T1_SS ~ int_ab*1 
  ABT2_SS  ~ int_ab*1 
  DE_T1_SS ~ int_de*1 
  DET2_SS ~ int_de*1 
  VI_T1_SS ~ int_vi*1 
  VIT2_SS ~ int_vi*1 
  #
  SLT1_SS ~ int_sl*1 
  SLT2SS ~ int_sl*1 
  
  
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
CLPM_M4.fit <- sem(CLPM_M4, data = data_SS, missing = 'ML')
#Warning message:
#  In lav_object_post_check(object) :
#  lavaan WARNING: covariance matrix of latent variables
#is not positive definite;
#use lavInspect(fit, "cov.lv") to investigate.
lavInspect(CLPM_M4.fit, "cov.lv")
lavInspect(CLPM_M4.fit, "cor.lv")
# Correlations between RQ11 & RQ21 and between RQ22 and RQ12 is very high, which makes sense since these two subdimensions belong to one dimension.


summary(CLPM_M4.fit, standardized = T, fit.measures=TRUE)[1]$FIT[c("chisq","df")]
#chisq      df 
#757.1568  90.0000 

# Chi-square difference test
#Df = 90 - 84 = 6
#Check on number of constrained / freed means:
#6
#Two models nested, so chi-square difference test with this df.
#Diff = 757.1568 - 725.4913 = 31.6655
#https://www.socscistatistics.com/pvalues/chidistribution.aspx
#The P-Value is .000019. The result is significant at p < .05.
# 1-pchisq(31.6655, df = 6)
#
#Hence, we proceed with Model 3 (strong factorial invariance - with freeing the means; CLPM_M3.fit).


##############################################################################################


# Model 3: strong factorial invariance #

# Now, give names to lagged relations
clpmModel <- '
  
  #####################
  # MEASUREMENT MODEL #
  #####################
  
  # Factor models for RQ1 at 2 waves.
  RQ11 =~ L1 * THT1_SS  + L2 * TBT1_SS  
  RQ12 =~ L1 * TH_T2_SS + L2 * TB_T2_SS 
  
  # Factor models for RQ2 at 2 waves.
  RQ21 =~ L3 * ACOMT1_SS + L4 * SATT1_SS
  RQ22 =~ L3 * ACOMT2SS  + L4 * SAT_T2SS
  
  # Factor models for SE at 2 waves.
  SE1 =~ L5 * AB_T1_SS + L6 * DE_T1_SS + L7 * VI_T1_SS
  SE2 =~ L5 * ABT2_SS  + L6 * DET2_SS  + L7 * VIT2_SS


  # Constrained intercepts over time
  THT1_SS ~ int_th*1 
  TH_T2_SS ~ int_th*1
  TBT1_SS ~ int_tb*1 
  TB_T2_SS ~ int_tb*1 
  ACOMT1_SS ~ int_acom*1 
  ACOMT2SS ~ int_acom*1 
  SATT1_SS ~ int_sat*1 
  SAT_T2SS ~ int_sat*1 
  #
  AB_T1_SS ~ int_ab*1 
  ABT2_SS  ~ int_ab*1 
  DE_T1_SS ~ int_de*1 
  DET2_SS ~ int_de*1 
  VI_T1_SS ~ int_vi*1 
  VIT2_SS ~ int_vi*1 
  #
  SLT1_SS ~ int_sl*1 
  SLT2SS ~ int_sl*1 
  
  
  # Free latent means on t=2
  RQ12 + RQ22 + SE2 + RQ11 + RQ21 + SE1 ~ 1

  
  ############
  # DYNAMICS #
  ############
  
  # Specify the lagged effects between the latent variables.
  RQ12 ~ Phi11 * RQ11 + Phi12 * RQ21 + Phi13 * SE1 + Phi14 * SLT1_SS
  RQ22 ~ Phi21 * RQ11 + Phi22 * RQ21 + Phi23 * SE1 + Phi24 * SLT1_SS
  #
  SE2 ~ Phi31 * RQ11 + Phi32 * RQ21 + Phi33 * SE1 + Phi34 * SLT1_SS
  SLT2SS ~ Phi41 * RQ11 + Phi42 * RQ21 + Phi43 * SE1 + Phi44 * SLT1_SS
  
  
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
clpmUnc <- sem(clpmModel, data = data_SS, missing = 'ML')
stdClpmUnc <- standardizedsolution(clpmUnc, type = "std.all", se = TRUE, zstat = TRUE, 
                                   pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE, 
                                   remove.eq = TRUE, remove.ineq = TRUE, remove.def = FALSE, 
                                   partable = NULL, GLIST = NULL, est = NULL)

# Model fit and estimates etc
summary(clpmUnc, standardized = T, fit.measures=TRUE)
stdClpmUnc # p-values of standardized effects

# substract values of interest
#summary(clpmUnc, standardized = T)$PE
indices <- 37:52
#summary(clpmUnc, standardized = T)$PE[indices,'std.all']
stdClpmUnc[indices, 4] # Substracts estimates from the column 'Std.all' in summary above.
lavInspect(clpmUnc, "vcov.std.all")[indices-6, indices-6] 


# GORICA values and weights
est <- stdClpmUnc[indices, 4]
names(est) <- c("RQ12_RQ11", "RQ12_RQ21", "RQ12_SE1", "RQ12_SL1", 
                "RQ22_RQ11", "RQ22_RQ21", "RQ22_SE1", "RQ22_SL1", 
                "SE2_RQ11", "SE2_RQ21", "SE2_SE1", "SE2_SL1",
                "SL2_RQ11", "SL2_RQ21", "SL2_SE1", "SL2_SL1"
)
vcov <- lavInspect(clpmUnc, "vcov.std.all")[indices-6, indices-6]

# Since we cannot evaluate absolute values,
# first check whether there are negative estimates.
# If so, then adjust sign in hypothesis accordingly.
#
#
# Q1: Phi_21 > Phi_12
H1_Q1 <- "
RQ22_RQ11 > RQ12_RQ21 
"
est[c("RQ22_RQ11", "RQ12_RQ21")]
#
H1_Q1 <- "
RQ22_RQ11 > RQ12_RQ21 
"
#
# Q2
H1_Q2 <- "
SE2_RQ11 > RQ12_SE1;
SL2_RQ11 > RQ12_SL1;
SE2_RQ21 > RQ22_SE1;
SL2_RQ21 > RQ22_SL1
"
est[c("SE2_RQ11", "RQ12_SE1", "SL2_RQ11", "RQ12_SL1", "SE2_RQ21", "RQ22_SE1", "SL2_RQ21", "RQ22_SL1")]
#
H1_Q2 <- "
-SE2_RQ11 > RQ12_SE1;
SL2_RQ11 > RQ12_SL1;
SE2_RQ21 > RQ22_SE1;
-SL2_RQ21 > RQ22_SL1
"
#
# Gorica
#
# Q1
set.seed(123)
goricaResults_Q1 <- goric(est, VCOV = vcov, H1_Q1, comparison = "complement", type = "gorica")
summary(goricaResults_Q1)
#The order-restricted hypothesis ‘H1_Q1’ has 1.7 times more support than its complement.
#
# Q2 
set.seed(123)
goricaResults_Q2 <- goric(est, VCOV = vcov, H1_Q2, comparison = "complement", type = "gorica")
summary(goricaResults_Q2)
#The order-restricted hypothesis ‘H1_Q2’ has 1.4 times more support than its complement.


#####

# time-interval dependency #

# Install and load packages
library(devtools)
install_github("rebeccakuiper/CTmeta")
library(CTmeta)
#?PhiPlot
if (!require("expm")) install.packages("expm") # install this package first (once)
library(expm)

# Create Phi matrix from this and make Phi plot
##est <- stdClpmUnc[indices, 4]
#Phi <- matrix(est, byrow=T, ncol = sqrt(length(est)))
#Drift <- logm(Phi)/1 # the drift matrix = the continuous-time equivalent of the discrete-time CLPM lagged effects matrix

# Explanation why it does not work here:
#> eigen(Phi)
#eigen() decomposition
#$values
#[1] -2.2083625  0.8957727  0.6961380  0.3052565
# So, there is no continuous-time equivalent... 
# and it is even an exploding process: |ev_1| > 1
