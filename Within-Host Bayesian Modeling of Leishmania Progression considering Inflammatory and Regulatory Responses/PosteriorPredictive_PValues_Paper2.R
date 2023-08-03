########################################################################
### Felix Pabon-Rodriguez
### R Code (NIMBLE) for Paper 2 
### Within-Host Bayesian Joint Modeling of Longitudinal and 
### Time-to-Event Data of Leishmania Infection
########################################################################

########################################################################
# Loading libraries
########################################################################

# Loading libraries
options(repos="https://cran.rstudio.com")
install_load <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

my_packages <- c("readxl", "dplyr", "coda", "ggplot2",
                 "lattice", "bayesplot", "BayesPostEst", "ggmcmc", 
                 "RCurl", "truncnorm", "kableExtra", "mvtnorm", "rlist", 
                 "extraDistr", "msm", "tmvtnorm", "runjags", "plotrix",
                 "lubridate", "ggpubr", "stringr", "nimble",
                 "igraph", "parallel", "doParallel", "MCMCvis", "Matrix")
invisible(install_load(my_packages))


########################################################################
# Reading Data and Helper Functions
########################################################################

# Reading Tables of Pathogen, Disease Status, Antibodies and Immune Responses
leish <- read.csv(file = "./Data/Updated Data/leish.csv", header = TRUE, 
                  na.strings = "", stringsAsFactors = TRUE)
dt_P <- read.csv(file = "./Data/Updated Data/dt_P.csv", 
                 header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_A <- read.csv(file = "./Data/Updated Data/dt_A.csv", 
                 header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_I1 <- read.csv(file = "./Data/Updated Data/dt_I1.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_I2 <- read.csv(file = "./Data/Updated Data/dt_I2.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_I3 <- read.csv(file = "./Data/Updated Data/dt_I3.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_R1 <- read.csv(file = "./Data/Updated Data/dt_R1.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_R2 <- read.csv(file = "./Data/Updated Data/dt_R2.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_R3 <- read.csv(file = "./Data/Updated Data/dt_R3.csv", 
                  header = TRUE, na.strings = "", stringsAsFactors = TRUE)

# Covariates 
levels(leish$Age_group)[levels(leish$Age_group)=="9-11 years old"] <-"6-8 years old"
levels(leish$Age_group)[levels(leish$Age_group)=="6-8 years old"] <-"6+ years old"
dt_age <- as.numeric(factor(leish$Age_group))

levels(leish$Ever.SNAP.Positive) <- c("0","1")
dt_snap <- as.numeric(as.vector(leish$Ever.SNAP.Positive))

levels(leish$treatment_group) <- c("0","1")
dt_trt <- as.numeric(as.vector(leish$treatment_group))

dt_dpp <- leish$DPP_Enroll

# New variable: disease_status per time point, based on leishvet score
# Healthy = 1 if leishvet score is 0 or 1
# Asymptomatic (inflammatory) = 2 if leishvet score is 1 or 2  
# Symptomatic (regulatory) = 3 if leishvet score is 3 or 4  
leish$DStatus_TP1 <- ifelse(leish$LeishVet.Score.TP1 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP1 == 2, 2, 3))
leish$DStatus_TP2 <- rep(NA,nrow(leish)) 
leish$DStatus_TP3 <- ifelse(leish$LeishVet.Score.TP3 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP3 == 2, 2, 3))
leish$DStatus_TP4 <- ifelse(leish$LeishVet.Score.TP4 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP4 == 2, 2, 3))
leish$DStatus_TP5 <- ifelse(leish$LeishVet.Score.TP5 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP5 == 2, 2, 3))
leish$DStatus_TP6 <- ifelse(leish$LeishVet.Score.TP6 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP6 == 2, 2, 3))
leish$DStatus_TP7 <- ifelse(leish$LeishVet.Score.TP7 %in% c(0,1), 1, 
                            ifelse(leish$LeishVet.Score.TP7 == 2, 2, 3))

dt_D <- as.matrix(leish[,c("DStatus_TP1","DStatus_TP2", 
                           "DStatus_TP3", "DStatus_TP4",
                           "DStatus_TP5", "DStatus_TP6",
                           "DStatus_TP7")])

# Removing column names
dt_P <- as.matrix(dt_P); colnames(dt_P) <- NULL
dt_A <- as.matrix(dt_A); colnames(dt_A) <- NULL 
dt_D <- as.matrix(dt_D); colnames(dt_D) <- NULL
dt_I1 <- matrix(as.numeric(as.matrix(dt_I1)), nrow = nrow(dt_I1))
dt_I2 <- matrix(as.numeric(as.matrix(dt_I2)), nrow = nrow(dt_I2))
dt_I3 <- matrix(as.numeric(as.matrix(dt_I3)), nrow = nrow(dt_I3))
dt_R1 <- matrix(as.numeric(as.matrix(dt_R1)), nrow = nrow(dt_R1))
dt_R2 <- matrix(as.numeric(as.matrix(dt_R2)), nrow = nrow(dt_R2))
dt_R3 <- matrix(as.numeric(as.matrix(dt_R3)), nrow = nrow(dt_R3))

# Remove time point 4 from PD1 responses for both CD4 and CD8
dt_R2[,4] <- rep(NA, nrow(dt_R2)) 
dt_R3[,4] <- rep(NA, nrow(dt_R3)) 

# Dimensions
N <- nrow(leish)
Time <- ncol(dt_P)
total_entries <- N*Time
ncolsM <- 24
ncolsX <- 6
dimIR <- 6

## Using Censoring for PCR load below detection limit (DL=10)
DL <- 10
Censored <- ifelse(dt_P <= DL, 0, 1)  # 0 represents the censored
dt_Pc <- ifelse(dt_P <= DL, NA, dt_P)    
censored_idxP <- which(Censored == 0)

# TRANSFORMATIONS
# Scale Pathogen load (divide by 20,000) 
# Scale Antibodies (divide by 10)
# Scale DPP, in Covariates (divide by 200)
generalmaxP <- 20000
generalmaxA <- 10
generalmaxDPP <- 200
dt_Pc <- dt_Pc / generalmaxP
dt_A <- dt_A / generalmaxA
dt_dpp <- dt_dpp / generalmaxDPP

DLt <- DL/generalmaxP
logit_DLt <- log(DLt/(1-DLt))

## Using Logit Transform
dt_Pc <- log(dt_Pc / (1 - dt_Pc))
dt_A <- log(dt_A / (1 - dt_A))

dt_I1 <- ifelse(dt_I1 == 0, 0.0001, dt_I1)  
dt_I2 <- ifelse(dt_I2 == 0, 0.0001, dt_I2) 
dt_I3 <- ifelse(dt_I3 == 0, 0.0001, dt_I3) 
dt_R1 <- ifelse(dt_R1 == 0, 0.0001, dt_R1) 
dt_R2 <- ifelse(dt_R2 == 0, 0.0001, dt_R2) 
dt_R3 <- ifelse(dt_R3 == 0, 0.0001, dt_R3) 

dt_I1 <- log(dt_I1 / (1 - dt_I1))
dt_I2 <- log(dt_I2 / (1 - dt_I2))
dt_I3 <- log(dt_I3 / (1 - dt_I3))
dt_R1 <- log(dt_R1 / (1 - dt_R1))
dt_R2 <- log(dt_R2 / (1 - dt_R2))
dt_R3 <- log(dt_R3 / (1 - dt_R3))

# Creating dummy variables for age group
# Since there is only one subject in 8+ years category, the last
# two categories are going to be combined and used as baseline 
age_0to2 <- age_3to5 <- numeric(length(dt_age))
for(i in 1:length(dt_age)){
  age_0to2[i] <- ifelse(dt_age[i]==1,1,0) 
  age_3to5[i] <- ifelse(dt_age[i]==2,1,0)
}

# X matrix
dt_X <- as.matrix(cbind(rep(1,50),age_0to2,age_3to5,
                        dt_snap,dt_dpp,dt_trt)) 
colnames(dt_X) <- NULL

# Looking at time indexes for removed subjects (LTFU)
first_timepoint <- rep(1,N)
last_timepoint <- rep(7,N)
# up to TP7
last_timepoint[which(leish$ID=="91")] <- 5   #10 removed, non-Leish related (ALIVE)
last_timepoint[which(leish$ID=="737")] <- 4  #13 removed, non-Leish related (DEATH)
last_timepoint[which(leish$ID=="743")] <- 6  #18 removed, severe Leish (DEATH)
last_timepoint[which(leish$ID=="750")] <- 2  #25 removed, non-Leish related (DEATH)
last_timepoint[which(leish$ID=="758")] <- 6  #33 removed, severe Leish (DEATH)
last_timepoint[which(leish$ID=="761")] <- 4  #36 removed, severe Leish (DEATH)
last_timepoint[which(leish$ID=="770")] <- 5  #44 removed, non-Leish related (ALIVE)
last_timepoint[which(leish$ID=="773")] <- 2  #47 removed, non-Leish related (DEATH)

# death status (including TP8 and TP9)
# status=1(dead) or 0 otherwise
status <- matrix(0,nrow = N, ncol = Time+2)
status[which(leish$ID=="743"),7:9] <- c(1,1,1)
status[which(leish$ID=="758"),7:9] <- c(1,1,1)
status[which(leish$ID=="761"),5:9] <- c(1,1,1,1,1)

status[which(leish$ID=="45"),9] <- 1
status[which(leish$ID=="55"),9] <- 1
status[which(leish$ID=="74"),9] <- 1
status[which(leish$ID=="77"),8:9] <- c(1,1)
status[which(leish$ID=="89"),9] <- 1
status[which(leish$ID=="119"),9] <- 1
status[which(leish$ID=="738"),9] <- 1
status[which(leish$ID=="747"),8:9] <- c(1,1)
status[which(leish$ID=="748"),9] <- 1
status[which(leish$ID=="753"),8:9] <- c(1,1)
status[which(leish$ID=="756"),9] <- 1
status[which(leish$ID=="759"),8:9] <- c(1,1)
status[which(leish$ID=="760"),9] <- 1
status[which(leish$ID=="763"),9] <- 1
status[which(leish$ID=="766"),8:9] <- c(1,1)

# death indicator as vector 
death <- status[,9]

# time until death or censoring
time <- rep(9,N)
time[which(leish$ID=="91")] <- 6   #dead, non-Leish related (ALIVE)
time[which(leish$ID=="737")] <- 5  #censored, non-Leish related (DEATH)
time[which(leish$ID=="743")] <- 7  #dead, severe Leish (DEATH)
time[which(leish$ID=="750")] <- 3  #censored, non-Leish related (DEATH)
time[which(leish$ID=="758")] <- 7  #dead, severe Leish (DEATH)
time[which(leish$ID=="761")] <- 5  #dead, severe Leish (DEATH)
time[which(leish$ID=="770")] <- 6  #censored, non-Leish related (ALIVE)
time[which(leish$ID=="773")] <- 2  #censored, non-Leish related (DEATH)
time[which(leish$ID=="45")] <- 9   #dead, severe Leish (DEATH)
time[which(leish$ID=="55")] <- 9   #dead, severe Leish (DEATH)
time[which(leish$ID=="74")] <- 9   #dead, severe Leish (DEATH)
time[which(leish$ID=="77")] <- 8   #dead, severe Leish (DEATH)
time[which(leish$ID=="89")] <- 9   #dead, severe Leish (DEATH)
time[which(leish$ID=="119")] <- 9  #dead, severe Leish (DEATH)
time[which(leish$ID=="738")] <- 9  #dead, severe Leish (DEATH) 
time[which(leish$ID=="747")] <- 8  #dead, severe Leish (DEATH)
time[which(leish$ID=="748")] <- 9  #dead, severe Leish (DEATH)
time[which(leish$ID=="753")] <- 8  #dead, severe Leish (DEATH) 
time[which(leish$ID=="756")] <- 9  #dead, severe Leish (DEATH)
time[which(leish$ID=="759")] <- 8  #dead, severe Leish (DEATH)
time[which(leish$ID=="760")] <- 9  #dead, severe Leish (DEATH)
time[which(leish$ID=="763")] <- 9  #dead, severe Leish (DEATH)
time[which(leish$ID=="766")] <- 8  #dead, severe Leish (DEATH) 
time[which(leish$ID=="766")] <- 8  #dead, severe Leish (DEATH) 
time[which(leish$ID=="64")] <- 7   #censored, (ALIVE) 
time[which(leish$ID=="754")] <- 7  #censored, (ALIVE) 
time[which(leish$ID=="755")] <- 7  #censored, (ALIVE) 


# Creating Multivariate Array
# Each row is a subject and column is each outcome
MV_IR_Object_at_time <- function(j,I1,I2,I3,R1,R2,R3){
  cbind(I1[,j], I2[,j], I3[,j], R1[,j], R2[,j], R3[,j]) 
}

# Sample Covariance of IR Data
MIR <- rbind(MV_IR_Object_at_time(2,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
             MV_IR_Object_at_time(3,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
             MV_IR_Object_at_time(4,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
             MV_IR_Object_at_time(5,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
             MV_IR_Object_at_time(6,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
             MV_IR_Object_at_time(7,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3))
covSampleIR <- cov(na.omit(MIR))

# Multivariate array of IR data
dt_IR <- array(c(MV_IR_Object_at_time(1,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(2,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(3,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(4,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(5,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(6,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3),
                 MV_IR_Object_at_time(7,dt_I1,dt_I2,dt_I3,dt_R1,dt_R2,dt_R3)),
               c(N,dimIR,Time))

# Fixing IR values 
dt_IR[25,,2] <- NA
dt_IR[28,,7] <- NA
dt_IR[47,,2] <- NA

# Missing indexes
missing_ind_D <- which(is.na(dt_D))
missing_ind_A <- which(is.na(dt_A))
missing_ind_P <- which(is.na(dt_P))
missing_ind_IR <- which(is.na(dt_IR))

# Entries to ignore due to lost-to-follow-up 
# For P, A, and D
ignore_mat1 <- matrix(data = 1,N,Time)
for(i in 1:N){
  ignore_mat1[i,(first_timepoint[i]):last_timepoint[i]] <- 0
}
ignore_entries <- which(ignore_mat1==1)

# For IR
ignore_mat2 <- nimArray(value = 1, dim = c(N,dimIR,Time))
for(i in 1:N){
  for(k in 1:6){
    ignore_mat2[,k,][i,(first_timepoint[i]):last_timepoint[i]] <- 0
  }
}
ignore_entries_IR <- which(ignore_mat2==1)

# Indexes of those missing latent quantities that needs to be monitored
idx_impute_D <- which(is.na(dt_D))[!(which(is.na(dt_D)) %in% ignore_entries)]
idx_impute_A <- which(is.na(dt_A))[!(which(is.na(dt_A)) %in% ignore_entries)]
idx_impute_IR <- which(is.na(dt_IR))[!(which(is.na(dt_D)) %in% ignore_entries_IR)]
idx_impute_P <- sort(c(which(is.na(dt_P))[!(which(is.na(dt_P)) %in% ignore_entries)],
                       censored_idxP))

## Missing indexes i,j,k ====================
idx_miss_D <- matrix(0,N,Time)
idx_miss_D[idx_impute_D] <- 1
vec_idx_miss_D <- c()
for (i in 1:nrow(dt_D)){
  for(j in first_timepoint[i]:last_timepoint[i]){
    if(idx_miss_D[i,j]==1){
      vec_idx_miss_D <- c(vec_idx_miss_D,paste0("D[",i,", ",j,"]"))             
    }
  }
}

idx_miss_P <- matrix(0,N,Time)
idx_miss_P[idx_impute_P] <- 1
vec_idx_miss_P <- c()
for (i in 1:nrow(dt_P)){
  for(j in first_timepoint[i]:last_timepoint[i]){
    if(idx_miss_P[i,j]==1){
      vec_idx_miss_P <- c(vec_idx_miss_P,paste0("P[",i,", ",j,"]"))             
    }
  }
}

idx_miss_A <- matrix(0,N,Time)
idx_miss_A[idx_impute_A] <- 1
vec_idx_miss_A <- c()
for (i in 1:nrow(dt_A)){
  for(j in first_timepoint[i]:last_timepoint[i]){
    if(idx_miss_A[i,j]==1){
      vec_idx_miss_A <- c(vec_idx_miss_A,paste0("A[",i,", ",j,"]"))             
    }
  }
}

idx_miss_IR <- nimArray(value = 0, dim = c(N,dimIR,Time))
idx_miss_IR[idx_impute_IR] <- 1
vec_idx_miss_IR <- c()
for(i in 1:nrow(dt_P)){
  for(j in first_timepoint[i]:last_timepoint[i]){
    for(k in 1:6){
      if(idx_miss_IR[i,k,j]==1){
        vec_idx_miss_IR <- c(vec_idx_miss_IR,paste0("IR[",i,", ",k,", ",j,"]"))             
      }
    }
  }
}


########################################################################
# Save Results as RDS File
########################################################################

pred_P_subj1 <- readRDS("pred_P_subj1_V3.rds")
pred_P_subj2 <- readRDS("pred_P_subj2_V3.rds")
pred_P_subj3 <- readRDS("pred_P_subj3_V3.rds")
pred_P_subj4 <- readRDS("pred_P_subj4_V3.rds")

pred_A_subj1 <- readRDS("pred_A_subj1_V3.rds")
pred_A_subj2 <- readRDS("pred_A_subj2_V3.rds")
pred_A_subj3 <- readRDS("pred_A_subj3_V3.rds")
pred_A_subj4 <- readRDS("pred_A_subj4_V3.rds")

pred_IR_subj1 <- readRDS("pred_IR_subj1_V3.rds")
pred_IR_subj2 <- readRDS("pred_IR_subj2_V3.rds")
pred_IR_subj3 <- readRDS("pred_IR_subj3_V3.rds")
pred_IR_subj4 <- readRDS("pred_IR_subj4_V3.rds")


########################################################################
# Posterior P-Values
########################################################################

keep <- c(which(leish$ID =="761"), which(leish$ID =="772"),
          which(leish$ID =="769"), which(leish$ID =="775"))

ppd_samplesP <- list(s1 = pred_P_subj1,
                     s2 = pred_P_subj2,
                     s3 = pred_P_subj3,
                     s4 = pred_P_subj4)

ppd_samplesA <- list(s1 = pred_A_subj1,
                     s2 = pred_A_subj2,
                     s3 = pred_A_subj3,
                     s4 = pred_A_subj4)

ppd_samplesIR1 <- list(s1 = pred_IR_subj1[,1,],
                     s2 = pred_IR_subj2[,1,],
                     s3 = pred_IR_subj3[,1,],
                     s4 = pred_IR_subj4[,1,])

ppd_samplesIR2 <- list(s1 = pred_IR_subj1[,2,],
                       s2 = pred_IR_subj2[,2,],
                       s3 = pred_IR_subj3[,2,],
                       s4 = pred_IR_subj4[,2,])

ppd_samplesIR3 <- list(s1 = pred_IR_subj1[,3,],
                       s2 = pred_IR_subj2[,3,],
                       s3 = pred_IR_subj3[,3,],
                       s4 = pred_IR_subj4[,3,])

ppd_samplesIR4 <- list(s1 = pred_IR_subj1[,4,],
                       s2 = pred_IR_subj2[,4,],
                       s3 = pred_IR_subj3[,4,],
                       s4 = pred_IR_subj4[,4,])

ppd_samplesIR5 <- list(s1 = pred_IR_subj1[,5,],
                       s2 = pred_IR_subj2[,5,],
                       s3 = pred_IR_subj3[,5,],
                       s4 = pred_IR_subj4[,5,])

ppd_samplesIR6 <- list(s1 = pred_IR_subj1[,6,],
                       s2 = pred_IR_subj2[,6,],
                       s3 = pred_IR_subj3[,6,],
                       s4 = pred_IR_subj4[,6,])



dt_A = (1/(1+exp(-1*dt_A)))*generalmaxA
dt_P = (1/(1+exp(-1*dt_Pc)))*generalmaxP
dt_IR1 = 1/(1+exp(-1*dt_I1))
dt_IR2 = 1/(1+exp(-1*dt_I2))
dt_IR3 = 1/(1+exp(-1*dt_I3))
dt_IR4 = 1/(1+exp(-1*dt_R1))
dt_IR5 = 1/(1+exp(-1*dt_R2))
dt_IR6 = 1/(1+exp(-1*dt_R3))


posterior.pvalues <- function(samples_P = ppd_samplesP,
                              samples_A = ppd_samplesA,
                              samples_IR1 = ppd_samplesIR1,
                              samples_IR2 = ppd_samplesIR2,
                              samples_IR3 = ppd_samplesIR3,
                              samples_IR4 = ppd_samplesIR4,
                              samples_IR5 = ppd_samplesIR5,
                              samples_IR6 = ppd_samplesIR6){
  
  # Empty objects
  pvalues_P <- rep(NA,4)
  pvalues_A <- rep(NA,4)
  pvalues_IR1 <- rep(NA,4)
  pvalues_IR2 <- rep(NA,4)
  pvalues_IR3 <- rep(NA,4)
  pvalues_IR4 <- rep(NA,4)
  pvalues_IR5 <- rep(NA,4)
  pvalues_IR6 <- rep(NA,4)
  
  for(subj in 1:4){
  
  # Getting data objects
  if(subj == 1){
    pathogen = samples_P$s1
    antibodies = samples_A$s1
    IRresp1 = samples_IR1$s1
    IRresp2 = samples_IR2$s1
    IRresp3 = samples_IR3$s1
    IRresp4 = samples_IR4$s1
    IRresp5 = samples_IR5$s1
    IRresp6 = samples_IR6$s1
    
    pathogen <- 1/(1+exp(-1*pathogen))
    antibodies <- 1/(1+exp(-1*antibodies))
    pathogen = pathogen*generalmaxP
    antibodies = antibodies*generalmaxA
    IRresp1 <- 1/(1+exp(-1*IRresp1))
    IRresp2 <- 1/(1+exp(-1*IRresp2))
    IRresp3 <- 1/(1+exp(-1*IRresp3))
    IRresp4 <- 1/(1+exp(-1*IRresp4))
    IRresp5 <- 1/(1+exp(-1*IRresp5))
    IRresp6 <- 1/(1+exp(-1*IRresp6))
    
    id = 761
    row = which(leish$ID==id)
    
    # Compute test statistic on observed data and on each simulated dataset
    mean_obs <- c(mean(dt_P[row,2:4]),
                  mean(dt_A[row,2:4]),
                  mean(dt_IR1[row,2:4]),
                  mean(dt_IR2[row,2:4]),
                  mean(dt_IR3[row,2:4]),
                  mean(dt_IR4[row,2:4]),
                  mean(dt_IR5[row,2:3]),
                  mean(dt_IR6[row,2:3]))
    mse_obs <- c(mean((dt_P[row,2:4] - mean_obs[1])^2), 
                 mean((dt_A[row,2:4] - mean_obs[2])^2),
                 mean((dt_IR1[row,2:4] - mean_obs[3])^2),
                 mean((dt_IR2[row,2:4] - mean_obs[4])^2),
                 mean((dt_IR3[row,2:4] - mean_obs[5])^2),
                 mean((dt_IR4[row,2:4] - mean_obs[6])^2),
                 mean((dt_IR5[row,2:3] - mean_obs[7])^2),
                 mean((dt_IR6[row,2:3] - mean_obs[8])^2))
    
    mse_post_P <- apply(pathogen[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_A <- apply(antibodies[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR1 <- apply(IRresp1[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR2 <- apply(IRresp2[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR3 <- apply(IRresp3[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR4 <- apply(IRresp4[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR5 <- apply(IRresp5[,2:3], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR6 <- apply(IRresp6[,2:3], 1, 
                        function(x) mean((x - mean(x))^2))
    
    # Calculate p-values
    p_mse_P <- mean(mse_post_P >= mse_obs[1])
    p_mse_A <- mean(mse_post_A >= mse_obs[2])
    p_mse_IR1 <- mean(mse_post_IR1 >= mse_obs[3])
    p_mse_IR2 <- mean(mse_post_IR2 >= mse_obs[4])
    p_mse_IR3 <- mean(mse_post_IR3 >= mse_obs[5])
    p_mse_IR4 <- mean(mse_post_IR4 >= mse_obs[6])
    p_mse_IR5 <- mean(mse_post_IR5 >= mse_obs[7])
    p_mse_IR6 <- mean(mse_post_IR6 >= mse_obs[8])
    
  } else if(subj == 2){
    pathogen = samples_P$s2
    antibodies = samples_A$s2
    IRresp1 = samples_IR1$s2
    IRresp2 = samples_IR2$s2
    IRresp3 = samples_IR3$s2
    IRresp4 = samples_IR4$s2
    IRresp5 = samples_IR5$s2
    IRresp6 = samples_IR6$s2
    
    pathogen <- 1/(1+exp(-1*pathogen))
    antibodies <- 1/(1+exp(-1*antibodies))
    pathogen = pathogen*generalmaxP
    antibodies = antibodies*generalmaxA
    IRresp1 <- 1/(1+exp(-1*IRresp1))
    IRresp2 <- 1/(1+exp(-1*IRresp2))
    IRresp3 <- 1/(1+exp(-1*IRresp3))
    IRresp4 <- 1/(1+exp(-1*IRresp4))
    IRresp5 <- 1/(1+exp(-1*IRresp5))
    IRresp6 <- 1/(1+exp(-1*IRresp6))
    
    id = 772
    row = which(leish$ID==id)
    
    # Compute test statistic on observed data and on each simulated dataset
    mean_obs <- c(mean(dt_P[row,2:4]),
                  mean(dt_A[row,2:4]),
                  mean(dt_IR1[row,2:4]),
                  mean(dt_IR2[row,2:4]),
                  mean(dt_IR3[row,2:4]),
                  mean(dt_IR4[row,2:4]),
                  mean(dt_IR5[row,2:3]),
                  mean(dt_IR6[row,2:3]))
    mse_obs <- c(mean((dt_P[row,2:4] - mean_obs[1])^2), 
                 mean((dt_A[row,2:4] - mean_obs[2])^2),
                 mean((dt_IR1[row,2:4] - mean_obs[3])^2),
                 mean((dt_IR2[row,2:4] - mean_obs[4])^2),
                 mean((dt_IR3[row,2:4] - mean_obs[5])^2),
                 mean((dt_IR4[row,2:4] - mean_obs[6])^2),
                 mean((dt_IR5[row,2:3] - mean_obs[7])^2),
                 mean((dt_IR6[row,2:3] - mean_obs[8])^2))
    
    mse_post_P <- apply(pathogen[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_A <- apply(antibodies[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR1 <- apply(IRresp1[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR2 <- apply(IRresp2[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR3 <- apply(IRresp3[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR4 <- apply(IRresp4[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR5 <- apply(IRresp5[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR6 <- apply(IRresp6[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    
    # Calculate p-values
    p_mse_P <- mean(mse_post_P >= mse_obs[1])
    p_mse_A <- mean(mse_post_A >= mse_obs[2])
    p_mse_IR1 <- mean(mse_post_IR1 >= mse_obs[3])
    p_mse_IR2 <- mean(mse_post_IR2 >= mse_obs[4])
    p_mse_IR3 <- mean(mse_post_IR3 >= mse_obs[5])
    p_mse_IR4 <- mean(mse_post_IR4 >= mse_obs[6])
    p_mse_IR5 <- mean(mse_post_IR5 >= mse_obs[7])
    p_mse_IR6 <- mean(mse_post_IR6 >= mse_obs[8])
    
  } else if(subj == 3){
    pathogen = samples_P$s3
    antibodies = samples_A$s3
    IRresp1 = samples_IR1$s3
    IRresp2 = samples_IR2$s3
    IRresp3 = samples_IR3$s3
    IRresp4 = samples_IR4$s3
    IRresp5 = samples_IR5$s3
    IRresp6 = samples_IR6$s3
    
    pathogen <- 1/(1+exp(-1*pathogen))
    antibodies <- 1/(1+exp(-1*antibodies))
    pathogen = pathogen*generalmaxP
    antibodies = antibodies*generalmaxA
    IRresp1 <- 1/(1+exp(-1*IRresp1))
    IRresp2 <- 1/(1+exp(-1*IRresp2))
    IRresp3 <- 1/(1+exp(-1*IRresp3))
    IRresp4 <- 1/(1+exp(-1*IRresp4))
    IRresp5 <- 1/(1+exp(-1*IRresp5))
    IRresp6 <- 1/(1+exp(-1*IRresp6))
    
    id = 769
    row = which(leish$ID==id)
    
    # Compute test statistic on observed data and on each simulated dataset
    mean_obs <- c(mean(dt_P[row,2:4]),
                  mean(dt_A[row,2:4]),
                  mean(dt_IR1[row,2:4]),
                  mean(dt_IR2[row,2:4]),
                  mean(dt_IR3[row,2:4]),
                  mean(dt_IR4[row,2:4]),
                  mean(dt_IR5[row,2:3]),
                  mean(dt_IR6[row,2:3]))
    mse_obs <- c(mean((dt_P[row,2:4] - mean_obs[1])^2), 
                 mean((dt_A[row,2:4] - mean_obs[2])^2),
                 mean((dt_IR1[row,2:4] - mean_obs[3])^2),
                 mean((dt_IR2[row,2:4] - mean_obs[4])^2),
                 mean((dt_IR3[row,2:4] - mean_obs[5])^2),
                 mean((dt_IR4[row,2:4] - mean_obs[6])^2),
                 mean((dt_IR5[row,2:3] - mean_obs[7])^2),
                 mean((dt_IR6[row,2:3] - mean_obs[8])^2))
    
    mse_post_P <- apply(pathogen[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_A <- apply(antibodies[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR1 <- apply(IRresp1[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR2 <- apply(IRresp2[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR3 <- apply(IRresp3[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR4 <- apply(IRresp4[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR5 <- apply(IRresp5[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR6 <- apply(IRresp6[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    
    # Calculate p-values
    p_mse_P <- mean(mse_post_P >= mse_obs[1])
    p_mse_A <- mean(mse_post_A >= mse_obs[2])
    p_mse_IR1 <- mean(mse_post_IR1 >= mse_obs[3])
    p_mse_IR2 <- mean(mse_post_IR2 >= mse_obs[4])
    p_mse_IR3 <- mean(mse_post_IR3 >= mse_obs[5])
    p_mse_IR4 <- mean(mse_post_IR4 >= mse_obs[6])
    p_mse_IR5 <- mean(mse_post_IR5 >= mse_obs[7])
    p_mse_IR6 <- mean(mse_post_IR6 >= mse_obs[8])
    
  } else if(subj == 4){
    pathogen = samples_P$s4
    antibodies = samples_A$s4
    IRresp1 = samples_IR1$s4
    IRresp2 = samples_IR2$s4
    IRresp3 = samples_IR3$s4
    IRresp4 = samples_IR4$s4
    IRresp5 = samples_IR5$s4
    IRresp6 = samples_IR6$s4
    
    pathogen <- 1/(1+exp(-1*pathogen))
    antibodies <- 1/(1+exp(-1*antibodies))
    pathogen = pathogen*generalmaxP
    antibodies = antibodies*generalmaxA
    IRresp1 <- 1/(1+exp(-1*IRresp1))
    IRresp2 <- 1/(1+exp(-1*IRresp2))
    IRresp3 <- 1/(1+exp(-1*IRresp3))
    IRresp4 <- 1/(1+exp(-1*IRresp4))
    IRresp5 <- 1/(1+exp(-1*IRresp5))
    IRresp6 <- 1/(1+exp(-1*IRresp6))
    
    id = 775
    row = which(leish$ID==id)
    
    # Compute test statistic on observed data and on each simulated dataset
    mean_obs <- c(mean(dt_P[row,2:4]),
                  mean(dt_A[row,2:4]),
                  mean(dt_IR1[row,2:4]),
                  mean(dt_IR2[row,2:4]),
                  mean(dt_IR3[row,2:4]),
                  mean(dt_IR4[row,2:4]),
                  mean(dt_IR5[row,2:3]),
                  mean(dt_IR6[row,2:3]))
    mse_obs <- c(mean((dt_P[row,2:4] - mean_obs[1])^2), 
                 mean((dt_A[row,2:4] - mean_obs[2])^2),
                 mean((dt_IR1[row,2:4] - mean_obs[3])^2),
                 mean((dt_IR2[row,2:4] - mean_obs[4])^2),
                 mean((dt_IR3[row,2:4] - mean_obs[5])^2),
                 mean((dt_IR4[row,2:4] - mean_obs[6])^2),
                 mean((dt_IR5[row,2:3] - mean_obs[7])^2),
                 mean((dt_IR6[row,2:3] - mean_obs[8])^2))
    
    mse_post_P <- apply(pathogen[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_A <- apply(antibodies[,2:4], 1, 
                        function(x) mean((x - mean(x))^2))
    mse_post_IR1 <- apply(IRresp1[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR2 <- apply(IRresp2[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR3 <- apply(IRresp3[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR4 <- apply(IRresp4[,2:4], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR5 <- apply(IRresp5[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    mse_post_IR6 <- apply(IRresp6[,2:3], 1, 
                          function(x) mean((x - mean(x))^2))
    
    # Calculate p-values
    p_mse_P <- mean(mse_post_P >= mse_obs[1])
    p_mse_A <- mean(mse_post_A >= mse_obs[2])
    p_mse_IR1 <- mean(mse_post_IR1 >= mse_obs[3])
    p_mse_IR2 <- mean(mse_post_IR2 >= mse_obs[4])
    p_mse_IR3 <- mean(mse_post_IR3 >= mse_obs[5])
    p_mse_IR4 <- mean(mse_post_IR4 >= mse_obs[6])
    p_mse_IR5 <- mean(mse_post_IR5 >= mse_obs[7])
    p_mse_IR6 <- mean(mse_post_IR6 >= mse_obs[8])
    
    
  } else {cat("Error: No valid subject")}
  
    pvalues_P[subj] <- p_mse_P
    pvalues_A[subj] <- p_mse_A
    pvalues_IR1[subj] <- p_mse_IR1 
    pvalues_IR2[subj] <- p_mse_IR2 
    pvalues_IR3[subj] <- p_mse_IR3 
    pvalues_IR4[subj] <- p_mse_IR4 
    pvalues_IR5[subj] <- p_mse_IR5 
    pvalues_IR6[subj] <- p_mse_IR6 
    
  }
  
  pvalues <- data.frame(Post_P_Pathogen = mean(pvalues_P, na.rm=TRUE),
                        Post_P_Antibodies = mean(pvalues_A, na.rm=TRUE),
                        Post_P_IR1 = mean(pvalues_IR1, na.rm=TRUE),
                        Post_P_IR2 = mean(pvalues_IR2, na.rm=TRUE),
                        Post_P_IR3 = mean(pvalues_IR3, na.rm=TRUE),
                        Post_P_IR4 = mean(pvalues_IR4, na.rm=TRUE),
                        Post_P_IR5 = mean(pvalues_IR5, na.rm=TRUE),
                        Post_P_IR6 = mean(pvalues_IR6, na.rm=TRUE))
  
  return(t(pvalues))
  
}


posterior.pvalues()


# Results
# Post_P_Pathogen   0.2426667
# Post_P_Antibodies 0.4858333
# Post_P_IR1        0.3478333
# Post_P_IR2        0.8683333
# Post_P_IR3        0.5628333
# Post_P_IR4        0.8010000
# Post_P_IR5        0.4936667
# Post_P_IR6        0.8471667
