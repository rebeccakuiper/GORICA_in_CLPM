
# Install packages and load libraries #
#
if (!require("tidyverse")) install.packages("tidyverse") # install this package first (once)
if (!require("lavaan")) install.packages("lavaan") # install this package first (once)
library(tidyverse)
library(lavaan)
#
if (!require("restriktor")) install.packages("restriktor") # install this package first (once)
library(restriktor) # for goric function
#
# Temporary, until this is on CRAN
#install.packages("devtools")
library(devtools)
devtools::install_github("yrosseel/lavaan")
devtools::install_github("leonardv/restriktor")
library(lavaan)
library(restriktor)


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
# On sum scores

clpmModel <- '
  
  ############
  # DYNAMICS #
  ############
  
  # Specify the lagged effects between the latent variables.
  TH_T2_SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  TB_T2_SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  ACOMT2SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  SAT_T2SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  #
  ABT2_SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  DET2_SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  VIT2_SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  #
  SLT2SS ~ THT1_SS + TBT1_SS + ACOMT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
  
  
  # Estimate the correlations within the same wave.
  # T1
   ACOMT1_SS ~~ THT1_SS + TBT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
   THT1_SS ~~ TBT1_SS + SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
   TBT1_SS ~~ SATT1_SS + AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
   SATT1_SS ~~ AB_T1_SS + DE_T1_SS + VI_T1_SS + SLT1_SS
   AB_T1_SS ~~ DE_T1_SS + VI_T1_SS + SLT1_SS
   DE_T1_SS ~~ VI_T1_SS + SLT1_SS
   VI_T1_SS ~~ SLT1_SS
   # T2
   TH_T2_SS ~~ TB_T2_SS + SAT_T2SS + ACOMT2SS + ABT2_SS + DET2_SS + VIT2_SS + SLT2SS
   TB_T2_SS ~~ SAT_T2SS + ACOMT2SS + ABT2_SS + DET2_SS + VIT2_SS + SLT2SS
   SAT_T2SS ~~ ACOMT2SS + ABT2_SS + DET2_SS + VIT2_SS + SLT2SS
   ACOMT2SS ~~ ABT2_SS + DET2_SS + VIT2_SS + SLT2SS
   ABT2_SS ~~ DET2_SS + VIT2_SS + SLT2SS
   DET2_SS ~~ VIT2_SS + SLT2SS
   VIT2_SS ~~ SLT2SS
'
clpmUnc <- sem(clpmModel, data = data_SS, missing = 'ML')
stdClpmUnc <- standardizedsolution(clpmUnc, type = "std.all", se = TRUE, zstat = TRUE, 
                                   pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE, 
                                   remove.eq = TRUE, remove.ineq = TRUE, remove.def = FALSE, 
                                   partable = NULL, GLIST = NULL, est = NULL)

# Model fit and estimates etc
summary(clpmUnc, standardized = T, fit.measures=TRUE)
# In this case, 'perfect' fit since degrees of freedom is 0.
stdClpmUnc # p-values of standardized effects
stdClpmUnc[112:dim(stdClpmUnc)[1],] # p-values of standardized effects


# substract values of interest
#summary(clpmUnc, standardized = T)$PE
indices <- 1:64
#summary(clpmUnc, standardized = T)$PE[indices,'std.all']
stdClpmUnc[indices, 'est.std'] # Substracts estimates from the column 'Std.all' in summary above.
lavInspect(clpmUnc, "vcov.std.all")[indices, indices] 


# GORICA values and weights
est <- stdClpmUnc[indices, 'est.std']
names(est) <- c("TH2_TH1", "TH2_TB1", "TH2_ACOM1", "TH2_SAT1", "TH2_AB1", "TH2_DE1", "TH2_VI1", "TH2_SL1",
                "TB2_TH1", "TB2_TB1", "TB2_ACOM1", "TB2_SAT1", "TB2_AB1", "TB2_DE1", "TB2_VI1", "TB2_SL1",
                "ACOM2_TH1", "ACOM2_TB1", "ACOM2_ACOM1", "ACOM2_SAT1", "ACOM2_AB1", "ACOM2_DE1", "ACOM2_VI1", "ACOM2_SL1",
                "SAT2_TH1", "SAT2_TB1", "SAT2_ACOM1", "SATM2_SAT1", "SAT2_AB1", "SAT2_DE1", "SAT2_VI1", "SAT2_SL1",
                #
                "AB2_TH1", "AB2_TB1", "AB2_ACOM1", "AB2_SAT1", "AB2_AB1", "AB2_DE1", "AB2_VI1", "AB2_SL1",
                "DE2_TH1", "DE2_TB1", "DE2_ACOM1", "DE2_SAT1", "DE2_AB1", "DE2_DE1", "DE2_VI1", "DE2_SL1",
                "VI2_TH1", "VI2_TB1", "VI2_ACOM1", "VI2_SAT1", "VI2_AB1", "VI2_DE1", "VI2_VI1", "VI2_SL1",
                #
                "SL2_TH1", "SL2_TB1", "SL2_ACOM1", "SL2_SAT1", "SL2_AB1", "SL2_DE1", "SL2_VI1", "SL2_SL1"
)
vcov <- lavInspect(clpmUnc, "vcov.std.all")[indices, indices]

# Specify hypotheses (using abs() notation)
#
# Q1
H1_Q1 <- "
abs(ACOM2_TH1) > abs(TH2_ACOM1); abs(SAT2_TH1) > abs(TH2_SAT1);
abs(ACOM2_TB1) > abs(TB2_ACOM1); abs(SAT2_TB1) > abs(TB2_SAT1)
"
# Q2
H1_Q2 <- "
abs(AB2_TH1) > abs(TH2_AB1); abs(DE2_TH1) > abs(TH2_DE1); abs(VI2_TH1) > abs(TH2_VI1); abs(SL2_TH1) > abs(TH2_SL1);
abs(AB2_TB1) > abs(TB2_AB1); abs(DE2_TB1) > abs(TB2_DE1); abs(VI2_TB1) > abs(TB2_VI1); abs(SL2_TB1) > abs(TB2_SL1);
abs(AB2_ACOM1) > abs(ACOM2_AB1); abs(DE2_ACOM1) > abs(ACOM2_DE1); abs(VI2_ACOM1) > abs(ACOM2_VI1); abs(SL2_ACOM1) > abs(ACOM2_SL1);
abs(AB2_SAT1) > abs(SAT2_AB1); abs(DE2_SAT1) > abs(SAT2_DE1); abs(VI2_SAT1) > abs(SAT2_VI1); abs(SL2_SAT1) > abs(SAT2_SL1)
"
#
#
# GORICA
#
# Q1
set.seed(123)
goricaResults_Q1 <- goric(est, VCOV = vcov, H1_Q1, comparison = "complement", type = "gorica")
summary(goricaResults_Q1)
#The order-restricted hypothesis ‘H1_Q1’ has  2.6 times more support than its complement.
#
# Q2
set.seed(123)
#goricaResults_Q2 <- goric(est, VCOV = vcov, H1_Q2, comparison = "complement", type = "gorica")
#summary(goricaResults_Q2)
#Since the default method takes too long, we will use the bootstrap method (to calculate the penalty of the GORICA):
if (!require("parallel")) install.packages("parallel") # install this package first (once)
library(parallel)
nrCPUcores <- detectCores(all.tests = FALSE, logical = TRUE) # 4
# Note: many restrictions, so takes long(er) to run
goricaResults_Q2_b <- goric(est, VCOV = vcov, H1_Q2, comparison = "complement", type = "gorica", 
                            mix.weights = "boot", parallel = "snow", ncpus = nrCPUcores, mix.bootstrap = 99999)
summary(goricaResults_Q2_b)
#The order-restricted hypothesis ‘H1_Q2’ has 14 times more support than its complement.


#####

# time-interval dependency #

# Install and load packages
library(devtools)
if (!require("CTmeta")) install_github("rebeccakuiper/CTmeta") #install_github("rebeccakuiper/CTmeta", force = TRUE)
library(CTmeta)
#?PhiPlot

# Create Phi matrix from this and make Phi plot
#est <- stdClpmUnc[indices, 'est.std']
Phi <- matrix(est, byrow=T, ncol = sqrt(length(est)))
Phi
PhiPlot(DeltaT = 1, Phi, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL)
#ggPhiPlot(DeltaT = 1, Phi, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL)

# Alternative
#if (!require("expm")) install.packages("expm") # install this package first (once)
#library(expm)
#Drift <- logm(Phi)/1 # the drift matrix = the continuous-time equivalent of the discrete-time CLPM lagged effects matrix
#PhiPlot(DeltaT = 1, Drift = Drift, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL)

# Denote which of the Phi elements we want to plot based on Q1
WhichEl_Q1 <- matrix(c(
  0, 0, 1, 1, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0
), ncol = sqrt(length(est))) 
PhiPlot(DeltaT = 1, Phi, Min = 0, Max = 10, Step = 0.05, WhichElements = WhichEl_Q1, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL)

# Denote which of the Phi elements we want to plot based on Q2
WhichEl_Q2 <- matrix(c(
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1, 1,
  1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0, 0
), ncol = sqrt(length(est))) 
PhiPlot(DeltaT = 1, Phi, Min = 0, Max = 10, Step = 0.05, WhichElements = WhichEl_Q2, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL)
