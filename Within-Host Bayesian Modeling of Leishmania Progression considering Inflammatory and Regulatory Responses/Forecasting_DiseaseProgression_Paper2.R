########################################################################
### Felix Pabon-Rodriguez
### R Code (NIMBLE) for Paper 2 
### Within-Host Bayesian Modeling of Leishmania Progression 
###   considering Inflammatory and Regulatory Responses
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
# Forecasting
########################################################################

# Combine MCMC chains
results_file <- readRDS(file = "resultsBJLS_model.rds")
results_mcmc <- as.mcmc.list(lapply(1:3, function(x){as.mcmc(results_file[[x]])}))
combined_list <- combine.MCMC(results_mcmc)

# Use examples from data set, to make predictions
# Most important IDs: 761 (will be predicting missing data), 769, 772, and 775
# Corresponding row number 36, 43, 46, 49

# Add mode or mean to missing values on time-points 1 or 2
dtD <- dt_D
dtP <- dt_P
dtPc <- dt_Pc
dtA <- dt_A
dtI1 <- dt_I1
dtI2 <- dt_I2
dtI3 <- dt_I3
dtR1 <- dt_R1
dtR2 <- dt_R2
dtR3 <- dt_R3

dtPc[43,1] <- median(combined_list[,"P[43, 1]"])
dtPc[46,1] <- median(combined_list[,"P[46, 1]"])
dtPc[49,1] <- median(combined_list[,"P[49, 1]"])

dtPc[43,2] <- median(combined_list[,"P[43, 2]"])
dtPc[46,2] <- median(combined_list[,"P[46, 2]"])
dtPc[49,2] <- median(combined_list[,"P[49, 2]"])

dtD[36,2] <- which.max(table(combined_list[,"D[36, 2]"]))
dtD[43,2] <- which.max(table(combined_list[,"D[43, 2]"]))
dtD[46,2] <- which.max(table(combined_list[,"D[46, 2]"]))
dtD[49,2] <- which.max(table(combined_list[,"D[49, 2]"]))

dtI1[36,1] <- median(combined_list[,"IR[36, 1, 1]"])
dtI1[43,1] <- median(combined_list[,"IR[43, 1, 1]"])
dtI1[46,1] <- median(combined_list[,"IR[46, 1, 1]"])
dtI1[49,1] <- median(combined_list[,"IR[49, 1, 1]"])

dtI2[36,1] <- median(combined_list[,"IR[36, 2, 1]"])
dtI2[43,1] <- median(combined_list[,"IR[43, 2, 1]"])
dtI2[46,1] <- median(combined_list[,"IR[46, 2, 1]"])
dtI2[49,1] <- median(combined_list[,"IR[49, 2, 1]"])

dtI3[36,1] <- median(combined_list[,"IR[36, 3, 1]"])
dtI3[43,1] <- median(combined_list[,"IR[43, 3, 1]"])
dtI3[46,1] <- median(combined_list[,"IR[46, 3, 1]"])
dtI3[49,1] <- median(combined_list[,"IR[49, 3, 1]"])

dtR1[36,1] <- median(combined_list[,"IR[36, 4, 1]"])
dtR1[43,1] <- median(combined_list[,"IR[43, 4, 1]"])
dtR1[46,1] <- median(combined_list[,"IR[46, 4, 1]"])
dtR1[49,1] <- median(combined_list[,"IR[49, 4, 1]"])

dtR2[36,1] <- median(combined_list[,"IR[36, 5, 1]"])
dtR2[43,1] <- median(combined_list[,"IR[43, 5, 1]"])
dtR2[46,1] <- median(combined_list[,"IR[46, 5, 1]"])
dtR2[49,1] <- median(combined_list[,"IR[49, 5, 1]"])

dtR3[36,1] <- median(combined_list[,"IR[36, 6, 1]"])
dtR3[43,1] <- median(combined_list[,"IR[43, 6, 1]"])
dtR3[46,1] <- median(combined_list[,"IR[46, 6, 1]"])
dtR3[49,1] <- median(combined_list[,"IR[49, 6, 1]"])

# Create dataframes for predictions 
keep <- c(which(leish$ID =="761"), which(leish$ID =="772"),
          which(leish$ID =="769"), which(leish$ID =="775"))

mX <- dt_X[keep,]
mD <- dtD[keep,]
mP <- dtP[keep,]
mPc <- dtPc[keep,]
mA <- dtA[keep,]
mI1 <- dtI1[keep,]
mI2 <- dtI2[keep,]
mI3 <- dtI3[keep,]
mR1 <- dtR1[keep,]
mR2 <- dtR2[keep,]
mR3 <- dtR3[keep,]
mIR <- array(c(MV_IR_Object_at_time(1,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(2,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(3,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(4,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(5,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(6,mI1,mI2,mI3,mR1,mR2,mR3),
               MV_IR_Object_at_time(7,mI1,mI2,mI3,mR1,mR2,mR3)),
             c(4,dimIR,Time))

# Empty objects to be filled in by predictions
nsim <- 1500
Nstar <- 4
time_pred <- 4

pred_P_subj1 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_P_subj2 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_P_subj3 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_P_subj4 <- matrix(data = NA, nrow = nsim, ncol = time_pred)

pred_A_subj1 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_A_subj2 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_A_subj3 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_A_subj4 <- matrix(data = NA, nrow = nsim, ncol = time_pred)

pred_D_subj1 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_D_subj2 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_D_subj3 <- matrix(data = NA, nrow = nsim, ncol = time_pred)
pred_D_subj4 <- matrix(data = NA, nrow = nsim, ncol = time_pred)

pred_IR_subj1 <- nimArray(value = NA, dim = c(nsim,dimIR,time_pred))
pred_IR_subj2 <- nimArray(value = NA, dim = c(nsim,dimIR,time_pred))
pred_IR_subj3 <- nimArray(value = NA, dim = c(nsim,dimIR,time_pred))
pred_IR_subj4 <- nimArray(value = NA, dim = c(nsim,dimIR,time_pred))

pred_P_subj1[,1:2] <- cbind(rep(mPc[1,1],nsim),rep(mPc[1,2],nsim))
pred_P_subj2[,1:2] <- cbind(rep(mPc[2,1],nsim),rep(mPc[2,2],nsim)) 
pred_P_subj3[,1:2] <- cbind(rep(mPc[3,1],nsim),rep(mPc[3,2],nsim)) 
pred_P_subj4[,1:2] <- cbind(rep(mPc[4,1],nsim),rep(mPc[4,2],nsim)) 
pred_A_subj1[,1:2] <- cbind(rep(mA[1,1],nsim),rep(mA[1,2],nsim)) 
pred_A_subj2[,1:2] <- cbind(rep(mA[2,1],nsim),rep(mA[2,2],nsim))  
pred_A_subj3[,1:2] <- cbind(rep(mA[3,1],nsim),rep(mA[3,2],nsim))  
pred_A_subj4[,1:2] <- cbind(rep(mA[4,1],nsim),rep(mA[4,2],nsim))  
pred_D_subj1[,1:2] <- cbind(rep(mD[1,1],nsim),rep(mD[1,2],nsim)) 
pred_D_subj2[,1:2] <- cbind(rep(mD[2,1],nsim),rep(mD[2,2],nsim))  
pred_D_subj3[,1:2] <- cbind(rep(mD[3,1],nsim),rep(mD[3,2],nsim))  
pred_D_subj4[,1:2] <- cbind(rep(mD[4,1],nsim),rep(mD[4,2],nsim))  
pred_IR_subj1[,,1] <- matrix(rep(mIR[1,1:6,1],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj1[,,2] <- matrix(rep(mIR[1,1:6,2],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj2[,,1] <- matrix(rep(mIR[2,1:6,1],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj2[,,2] <- matrix(rep(mIR[2,1:6,2],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj3[,,1] <- matrix(rep(mIR[3,1:6,1],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj3[,,2] <- matrix(rep(mIR[3,1:6,2],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj4[,,1] <- matrix(rep(mIR[4,1:6,1],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)
pred_IR_subj4[,,2] <- matrix(rep(mIR[4,1:6,2],nsim),byrow = TRUE,
                             nrow = nsim, ncol = dimIR)

get_subject_M <- function(iter, time, P, A, D, IR){
  D2 <- 1*(D[iter,time]==2)
  D3 <- 1*(D[iter,time]==3)
  refD <- (1-D2-D3)
  M <- c(refD*P[iter,time],
         refD*A[iter,time],
         refD*IR[iter,1,time],
         refD*IR[iter,2,time],
         refD*IR[iter,3,time],
         refD*IR[iter,4,time],
         refD*IR[iter,5,time],
         refD*IR[iter,6,time],
         D2*P[iter,time],
         D2*A[iter,time],
         D2*IR[iter,1,time],
         D2*IR[iter,2,time],
         D2*IR[iter,3,time],
         D2*IR[iter,4,time],
         D2*IR[iter,5,time],
         D2*IR[iter,6,time],
         D3*P[iter,time],
         D3*A[iter,time],
         D3*IR[iter,1,time],
         D3*IR[iter,2,time],
         D3*IR[iter,3,time],
         D3*IR[iter,4,time],
         D3*IR[iter,5,time],
         D3*IR[iter,6,time])#24
  return(M)
}


for(i in 1:nsim){
  id_sim <- sample(1:dim(combined_list)[1],size = 1)
  betaP <- combined_list[id_sim,grepl("betaP",colnames(combined_list))]
  betaD2 <- combined_list[id_sim,grepl("betaD2",colnames(combined_list))]
  betaD3 <- combined_list[id_sim,grepl("betaD3",colnames(combined_list))]
  betaA <- combined_list[id_sim,grepl("betaA",colnames(combined_list))]
  betaI1 <- combined_list[id_sim,grepl("betaI1",colnames(combined_list))]
  betaI2 <- combined_list[id_sim,grepl("betaI2",colnames(combined_list))]
  betaI3 <- combined_list[id_sim,grepl("betaI3",colnames(combined_list))]
  betaR1 <- combined_list[id_sim,grepl("betaR1",colnames(combined_list))]
  betaR2 <- combined_list[id_sim,grepl("betaR2",colnames(combined_list))]
  betaR3 <- combined_list[id_sim,grepl("betaR3",colnames(combined_list))]
  alphaP <- combined_list[id_sim,grepl("alphaP",colnames(combined_list))]
  alphaA <- combined_list[id_sim,grepl("alphaA",colnames(combined_list))]
  alphaD2 <- combined_list[id_sim,grepl("alphaD2",colnames(combined_list))]
  alphaD3 <- combined_list[id_sim,grepl("alphaD3",colnames(combined_list))]
  alphaI1 <- combined_list[id_sim,grepl("alphaI1",colnames(combined_list))]
  alphaI2 <- combined_list[id_sim,grepl("alphaI2",colnames(combined_list))]
  alphaI3 <- combined_list[id_sim,grepl("alphaI3",colnames(combined_list))]
  alphaR1 <- combined_list[id_sim,grepl("alphaR1",colnames(combined_list))]
  alphaR2 <- combined_list[id_sim,grepl("alphaR2",colnames(combined_list))]
  alphaR3 <- combined_list[id_sim,grepl("alphaR3",colnames(combined_list))]
  sigmaP <- combined_list[id_sim,"sigmaP"]
  sigmaA <- combined_list[id_sim,"sigmaA"]
  SigmaIR <- matrix(combined_list[id_sim,grepl("SigmaIR",colnames(combined_list))],6,6)  
  thetaP <- combined_list[id_sim,"thetaP"]
  thetaA <- combined_list[id_sim,"thetaA"]
  thetaIR <- matrix(combined_list[id_sim,grepl("thetaIR",colnames(combined_list))],6,6)
  
  for(t in 3:time_pred){
    # Parameters
    wPtm1 <- combined_list[id_sim,paste0("wP","[",c(36,43,46,49),", ",t-1,"]")]
    wAtm1 <- combined_list[id_sim,paste0("wA","[",c(36,43,46,49),", ",t-1,"]")]
    wIRtm1_1 <- combined_list[id_sim,paste0("wIR","[",36,", ",1:6,", ", t-1,"]")]
    wIRtm1_2 <- combined_list[id_sim,paste0("wIR","[",43,", ",1:6,", ", t-1,"]")]
    wIRtm1_3 <- combined_list[id_sim,paste0("wIR","[",46,", ",1:6,", ", t-1,"]")]
    wIRtm1_4 <- combined_list[id_sim,paste0("wIR","[",49,", ",1:6,", ", t-1,"]")]

    # Matrices
    M1 <- get_subject_M(iter = i, time = t-1, P = pred_P_subj1, A = pred_A_subj1, 
                        D = pred_D_subj1, IR = pred_IR_subj1)
    M2 <- get_subject_M(iter = i, time = t-1, P = pred_P_subj2, A = pred_A_subj2, 
                        D = pred_D_subj2, IR = pred_IR_subj2)
    M3 <- get_subject_M(iter = i, time = t-1, P = pred_P_subj3, A = pred_A_subj3, 
                        D = pred_D_subj3, IR = pred_IR_subj3)
    M4 <- get_subject_M(iter = i, time = t-1, P = pred_P_subj4, A = pred_A_subj4, 
                        D = pred_D_subj4, IR = pred_IR_subj4)
      
    # Pathogen load (P)
    pred_P_subj1[i,t] <- rnorm(1, mean = M1%*%betaP + mX[1,]%*%alphaP + thetaP*wPtm1[1],
                               sd = sigmaP)
    pred_P_subj2[i,t] <- rnorm(1, mean = M2%*%betaP + mX[2,]%*%alphaP + thetaP*wPtm1[2],
                               sd = sigmaP)
    pred_P_subj3[i,t] <- rnorm(1, mean = M3%*%betaP + mX[3,]%*%alphaP + thetaP*wPtm1[3],
                               sd = sigmaP)
    pred_P_subj4[i,t] <- rnorm(1, mean = M4%*%betaP + mX[4,]%*%alphaP + thetaP*wPtm1[4],
                               sd = sigmaP)
    
    # Antibody levels (A)
    pred_A_subj1[i,t] <- rnorm(1,mean = M1%*%betaA + mX[1,]%*%alphaA + thetaA*wAtm1[1], 
                               sd = sigmaA)
    pred_A_subj2[i,t] <- rnorm(1,mean = M2%*%betaA + mX[2,]%*%alphaA + thetaA*wAtm1[2], 
                               sd = sigmaA)
    pred_A_subj3[i,t] <- rnorm(1,mean = M3%*%betaA + mX[3,]%*%alphaA + thetaA*wAtm1[3], 
                               sd = sigmaA)
    pred_A_subj4[i,t] <- rnorm(1,mean = M4%*%betaA + mX[4,]%*%alphaA + thetaA*wAtm1[4], 
                               sd = sigmaA)
    
    # Disease status (D)
    exp_term1 <- c(exp(M1%*%betaD2 + mX[1,]%*%alphaD2),
                   exp(M1%*%betaD3 + mX[1,]%*%alphaD3))
    total_term1 <- 1 + sum(exp_term1)
    prob1 <- c(1-sum(exp_term1/total_term1),exp_term1/total_term1) 
    pred_D_subj1[i,t] <- rcat(1,prob = prob1)
    
    exp_term2 <- c(exp(M2%*%betaD2 + mX[2,]%*%alphaD2),
                   exp(M2%*%betaD3 + mX[2,]%*%alphaD3))
    total_term2 <- 1 + sum(exp_term2)
    prob2 <- c(1-sum(exp_term2/total_term2),exp_term2/total_term2) 
    pred_D_subj2[i,t] <- rcat(1,prob = prob2)
    
    exp_term3 <- c(exp(M3%*%betaD2 + mX[3,]%*%alphaD2),
                   exp(M3%*%betaD3 + mX[3,]%*%alphaD3))
    total_term3 <- 1 + sum(exp_term3)
    prob3 <- c(1-sum(exp_term3/total_term3),exp_term3/total_term3) 
    pred_D_subj3[i,t] <- rcat(1,prob = prob3)
    
    exp_term4 <- c(exp(M4%*%betaD2 + mX[4,]%*%alphaD2),
                   exp(M4%*%betaD3 + mX[4,]%*%alphaD3))
    total_term4 <- 1 + sum(exp_term4)
    prob4 <- c(1-sum(exp_term4/total_term4),exp_term4/total_term4) 
    pred_D_subj4[i,t] <- rcat(1,prob = prob4)
    
    
    # Inflammatory-Regulatory Responses (IR)
    muIR1 <- c(M1%*%betaI1 + mX[1,]%*%alphaI1, 
               M1%*%betaI2 + mX[1,]%*%alphaI2,
               M1%*%betaI3 + mX[1,]%*%alphaI3, 
               M1%*%betaR1 + mX[1,]%*%alphaR1,
               M1%*%betaR2 + mX[1,]%*%alphaR2, 
               M1%*%betaR3 + mX[1,]%*%alphaR3)
    pred_IR_subj1[i,,t] <- rmvnorm(n = 1, mean = muIR1 + (thetaIR %*% wIRtm1_1), 
                                   sigma = as.matrix(forceSymmetric(SigmaIR)))
    
    muIR2 <- c(M2%*%betaI1 + mX[2,]%*%alphaI1, 
               M2%*%betaI2 + mX[2,]%*%alphaI2,
               M2%*%betaI3 + mX[2,]%*%alphaI3, 
               M2%*%betaR1 + mX[2,]%*%alphaR1,
               M2%*%betaR2 + mX[2,]%*%alphaR2, 
               M2%*%betaR3 + mX[2,]%*%alphaR3)
    pred_IR_subj2[i,,t] <- rmvnorm(n = 1, mean = muIR2 + (thetaIR %*% wIRtm1_2), 
                                   sigma = as.matrix(forceSymmetric(SigmaIR)))
    
    muIR3 <- c(M3%*%betaI1 + mX[3,]%*%alphaI1, 
               M3%*%betaI2 + mX[3,]%*%alphaI2,
               M3%*%betaI3 + mX[3,]%*%alphaI3, 
               M3%*%betaR1 + mX[3,]%*%alphaR1,
               M3%*%betaR2 + mX[3,]%*%alphaR2, 
               M3%*%betaR3 + mX[3,]%*%alphaR3)
    pred_IR_subj3[i,,t] <- rmvnorm(n = 1, mean = muIR3 + (thetaIR %*% wIRtm1_3), 
                                   sigma = as.matrix(forceSymmetric(SigmaIR)))
    
    muIR4 <- c(M4%*%betaI1 + mX[4,]%*%alphaI1, 
               M4%*%betaI2 + mX[4,]%*%alphaI2,
               M4%*%betaI3 + mX[4,]%*%alphaI3, 
               M4%*%betaR1 + mX[4,]%*%alphaR1,
               M4%*%betaR2 + mX[4,]%*%alphaR2, 
               M4%*%betaR3 + mX[4,]%*%alphaR3)
    pred_IR_subj4[i,,t] <- rmvnorm(n = 1, mean = muIR4 + (thetaIR %*% wIRtm1_4), 
                                   sigma = as.matrix(forceSymmetric(SigmaIR)))
  }
}
  
########################################################################
# Save Results as RDS File
########################################################################

saveRDS(pred_P_subj1, paste0("pred_P_subj1_V3.rds")) 
saveRDS(pred_P_subj2, paste0("pred_P_subj2_V3.rds")) 
saveRDS(pred_P_subj3, paste0("pred_P_subj3_V3.rds")) 
saveRDS(pred_P_subj4, paste0("pred_P_subj4_V3.rds")) 

saveRDS(pred_A_subj1, paste0("pred_A_subj1_V3.rds")) 
saveRDS(pred_A_subj2, paste0("pred_A_subj2_V3.rds")) 
saveRDS(pred_A_subj3, paste0("pred_A_subj3_V3.rds")) 
saveRDS(pred_A_subj4, paste0("pred_A_subj4_V3.rds")) 

saveRDS(pred_D_subj1, paste0("pred_D_subj1_V3.rds")) 
saveRDS(pred_D_subj2, paste0("pred_D_subj2_V3.rds")) 
saveRDS(pred_D_subj3, paste0("pred_D_subj3_V3.rds")) 
saveRDS(pred_D_subj4, paste0("pred_D_subj4_V3.rds")) 

saveRDS(pred_IR_subj1, paste0("pred_IR_subj1_V3.rds")) 
saveRDS(pred_IR_subj2, paste0("pred_IR_subj2_V3.rds")) 
saveRDS(pred_IR_subj3, paste0("pred_IR_subj3_V3.rds")) 
saveRDS(pred_IR_subj4, paste0("pred_IR_subj4_V3.rds")) 


########################################################################
# Plot Predictions
########################################################################

pred.plot <- function(subj = 1, 
                      rescaled = TRUE,
                      pred_time = 4){
  
  # Getting data objects
  if(subj == 1){
    pathogen <- 1/(1+exp(-1*pred_P_subj1))
    antibodies <- 1/(1+exp(-1*pred_A_subj1))
    Dstatus <- pred_D_subj1
    IRresp <- 1/(1+exp(-1*pred_IR_subj1))
  } else if (subj == 2){
    pathogen <- 1/(1+exp(-1*pred_P_subj2))
    antibodies <- 1/(1+exp(-1*pred_A_subj2))
    Dstatus <- pred_D_subj2
    IRresp <- 1/(1+exp(-1*pred_IR_subj2))
  } else if(subj == 3){
    pathogen <- 1/(1+exp(-1*pred_P_subj3))
    antibodies <- 1/(1+exp(-1*pred_A_subj3))
    Dstatus <- pred_D_subj3
    IRresp <- 1/(1+exp(-1*pred_IR_subj3))
  } else if(subj == 4){
    pathogen <- 1/(1+exp(-1*pred_P_subj4))
    antibodies <- 1/(1+exp(-1*pred_A_subj4))
    Dstatus <- pred_D_subj4
    IRresp <- 1/(1+exp(-1*pred_IR_subj4))
  } else {
    cat("Error: No valid subject")
  }
  
  if(rescaled == TRUE){
    pathogen = pathogen*generalmaxP
    antibodies = antibodies*generalmaxA
  } 
  nsim <- nrow(pathogen)
  
  # Estimates and confidence limits
  Pestimate <- apply(pathogen,2,function(x)mean(x))
  Aestimate <- apply(antibodies,2,function(x)mean(x))
  IRestimate <- matrix(NA,6,pred_time)
  for(i in 1:6){
    IRestimate[i,] <- apply(IRresp[,i,],2,function(x)mean(x))
  }
  
  lowerp95 <- apply(pathogen, 2, function(x) quantile(x,0.025))
  upperp95 <- apply(pathogen, 2, function(x) quantile(x,0.975))
  
  lowera95 <- apply(antibodies, 2, function(x) quantile(x,0.025))
  uppera95 <- apply(antibodies, 2, function(x) quantile(x,0.975))

  lowerIR95 <- upperIR95 <-matrix(NA,6,pred_time)
  for(i in 1:6){
    lowerIR95[i,] <- apply(IRresp[,i,],2,function(x)quantile(x,0.025))
    upperIR95[i,] <- apply(IRresp[,i,],2,function(x)quantile(x,0.975))
  }
  
  trP <- 1/(1+exp(-1*mPc))
  trA <- 1/(1+exp(-1*mA))
  trIR <- 1/(1+exp(-1*mIR))
  
  # Plot: Pathogen Load
  pdf(file=paste0("P_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, pathogen[1,2:pred_time], 
       xlab = "Time Point",
       ylab = "Pathogen Load",
       main = paste0("Predicted Pathogen Load for Subject ", subj),
       col = "darkgray",
       type = "l",
       ylim = c(0,max(upperp95[2:pred_time])),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 2:nrow(pathogen)){
    lines(2:pred_time, pathogen[i,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerp95[2:pred_time],rev(upperp95[2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerp95[2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperp95[2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trP[subj,2:pred_time]*generalmaxP, col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, Pestimate[2:pred_time],col = "red",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted", "Observed", "95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("red", "black", "royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Antibody Levels
  pdf(file=paste0("A_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, antibodies[1,2:pred_time], 
       xlab = "Time Point",
       ylab = "Antibody Levels",
       main = paste0("Predicted Antibody Levels for Subject ", subj),
       col = "darkgray",
       type = "l",
       ylim = c(0,max(uppera95[2:pred_time])))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 2:nrow(antibodies)){
    lines(2:pred_time, antibodies[i,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowera95[2:pred_time],rev(uppera95[2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowera95[2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, uppera95[2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trA[subj,2:pred_time]*generalmaxA, col = "forestgreen",lwd=4, lty=1, type="b")
  lines(2:pred_time, Aestimate[2:pred_time], col = "black", lwd = 4)
  legend("topleft", legend=c("Predicted", "Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("black", "forestgreen","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Disease State
  df <- matrix(data = 0, 3, pred_time)
  if(Dstatus[1,1]==1) {
    df[,1] = c(1,0,0)
  } else if(Dstatus[1,1]==2) {
    df[,1] = c(0,1,0)
  } else {df[,1] = c(0,0,1)}
  
  if(Dstatus[2,1]==1) {
    df[,2] = c(1,0,0)
  } else if(Dstatus[2,1]==2) {
    df[,2] = c(0,1,0)
  } else {df[,2] = c(0,0,1)}
  
  for(i in 3:pred_time){
    tb <- table(Dstatus[,i])/nsim
    df[,i] <- c(1-sum(tb[2:3]),tb[2],tb[3])
  }
  
  pdf(file=paste0("D_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, df[1,2:pred_time], 
       xlab = "Time Point",
       ylab = "Probability",
       main = paste0("Marginal Probabilities of Disease Status for Subject ", subj),
       col = "darkorange",
       type = "l",
       ylim = c(0,1),
       lwd = 3, lty = 1)
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  lines(2:pred_time, df[1,2:pred_time], col = "darkorange1", lwd = 4, lty = 1)
  lines(2:pred_time, df[2,2:pred_time], col = "darkorange1", lwd = 4, lty = 2)
  lines(2:pred_time, df[3,2:pred_time], col = "darkorange2", lwd = 4, lty = 3)
  legend("topleft", legend=c("D=1","D=2","D=3"), 
         lty = c(1,2,3), lwd = 3, col = c("darkorange", "darkorange1", "darkorange2"),
         inset = 0.02, cex = 0.9, seg.len=2, bg='white')
  dev.off() 
  
  # Plot: Proportion of INF-gamma+ of CD4+
  pdf(file=paste0("I1_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,1,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop INF-gamma+ of CD4+",
       main = paste0("Predicted Proportion of IFN-gamma+ of CD4+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,1,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[1,2:pred_time],rev(upperIR95[1,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[1,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[1,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,1,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[1,2:pred_time], col = "blue",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted", "Observed", "95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("blue", "black", "royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Proportion of Proliferation of CD4+
  pdf(file=paste0("I2_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,2,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop Proliferation of CD4+",
       main = paste0("Predicted Proportion of Proliferation of CD4+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,2,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[2,2:pred_time],rev(upperIR95[2,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[2,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[2,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,2,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[2,2:pred_time] ,col = "blue",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted", "Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("blue", "black","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Proportion of Proliferation of CD8+
  pdf(file=paste0("I3_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,3,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop Proliferation of CD8+",
       main = paste0("Predicted Proportion of Proliferation of CD8+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,3,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[3,2:pred_time],rev(upperIR95[3,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[3,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[3,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,3,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[3,2:pred_time] ,col = "blue",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted", "Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("blue", "black","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Proportion of IL-10+ of CD4+
  pdf(file=paste0("R1_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,4,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop IL-10+ of CD4+",
       main = paste0("Predicted Proportion of IL-10+ of CD4+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,4,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[4,2:pred_time],rev(upperIR95[4,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[4,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[4,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,4,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[4,2:pred_time] ,col = "purple",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted", "Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("purple", "black","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Proportion of PD1+ of CD4+
  pdf(file=paste0("R2_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,5,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop PD1+ of CD4+",
       main = paste0("Predicted Proportion of PD1+ of CD4+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,5,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[5,2:pred_time],rev(upperIR95[5,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[5,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[5,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,5,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[5,2:pred_time] ,col = "purple",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted","Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("purple", "black","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
  # Plot: Proportion of PD1+ of CD8+
  pdf(file=paste0("R3_Subject",subj,".pdf"))
  par(mfrow=c(1,1))
  plot(2:pred_time, IRresp[1,6,2:pred_time], 
       xlab = "Time Point",
       ylab = "Prop PD1+ of CD8+",
       main = paste0("Predicted Proportion of PD1+ of CD8+ for Subject ", subj),
       col = "darkgray",
       type = "l",
       xlim = c(2,pred_time),
       ylim = c(0,1),
       panel.first = rgb(0.92, 0.92, 0.92))
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = "gray92") 
  grid(NULL, NULL, col = "white")
  for(i in 1:nrow(pathogen)){
    lines(2:pred_time, IRresp[i,6,2:pred_time], col = alpha(rgb(0.37, 0.37, 0.37), 0.5))
  }
  polygon(c(2:pred_time,rev(2:pred_time)),c(lowerIR95[6,2:pred_time],rev(upperIR95[6,2:pred_time])),
          col=adjustcolor("blue", alpha.f = 0.25),
          border = NA)
  lines(2:pred_time, lowerIR95[6,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, upperIR95[6,2:pred_time], col = "royalblue3", lwd = 4, lty = 2)
  lines(2:pred_time, trIR[subj,6,2:pred_time], col = "black",lwd = 4, lty=1, type="b")
  lines(2:pred_time, IRestimate[6,2:pred_time] ,col = "purple",lwd = 4)
  grid(NULL, NULL, col = "white")
  legend("topleft", legend=c("Predicted","Observed","95% CI Bounds"), 
         lty = c(1,2,2), lwd = c(4,4,4), col = c("purple", "black","royalblue3"),
         inset = 0.02, cex = 0.85, seg.len=3.5, pch = c(NA,19,NA), bg='white')
  dev.off() 
  
}

pred.plot(subj = 1, rescaled = TRUE, pred_time = 4)
pred.plot(subj = 2, rescaled = TRUE, pred_time = 4)
pred.plot(subj = 3, rescaled = TRUE, pred_time = 4)
pred.plot(subj = 4, rescaled = TRUE, pred_time = 4)

 