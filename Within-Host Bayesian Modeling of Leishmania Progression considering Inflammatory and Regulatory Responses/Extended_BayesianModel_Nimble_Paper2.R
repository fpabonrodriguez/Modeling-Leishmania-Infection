########################################################################
### Felix Pabon-Rodriguez
### Dissertation R Code (Paper 2 == Chapter 3)
### Bayesian Hierarchical Model of Immune Responses
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
                 "igraph", "parallel", "doParallel", "MCMCvis")
invisible(install_load(my_packages))


########################################################################
# Creating Parallel Environment
########################################################################

# Creating clusters
ncore <- 3    
cl <- makeCluster(ncore, outfile = "log_file_Updated.txt")
clusterSetRNGStream(cl)
registerDoParallel(cl)

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
dt_D <- read.csv(file = "./Data/Updated Data/dt_D.csv", 
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
dt_age <- read.csv(file = "./Data/Updated Data/dt_age.csv", 
                   header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_snap <- read.csv(file = "./Data/Updated Data/dt_snap.csv", 
                    header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_dpp <- read.csv(file = "./Data/Updated Data/dt_dpp.csv", 
                   header = TRUE, na.strings = "", stringsAsFactors = TRUE)
dt_trt <- read.csv(file = "./Data/Updated Data/dt_trt.csv", 
                   header = TRUE, na.strings = "", stringsAsFactors = TRUE)

# Removing column names
dt_P <- as.matrix(dt_P); colnames(dt_P) <- NULL
dt_A <- as.matrix(dt_A); colnames(dt_A) <- NULL 
dt_D <- as.matrix(dt_D); colnames(dt_D) <- NULL
dt_I1 <- as.matrix(dt_I1); colnames(dt_I1) <- NULL
dt_I2 <- as.matrix(dt_I2); colnames(dt_I2) <- NULL
dt_I3 <- as.matrix(dt_I3); colnames(dt_I3) <- NULL
dt_R1 <- as.matrix(dt_R1); colnames(dt_R1) <- NULL
dt_R2 <- as.matrix(dt_R2); colnames(dt_R2) <- NULL
dt_R3 <- as.matrix(dt_R3); colnames(dt_R3) <- NULL

# Dimensions
N <- nrow(leish)
Time <- ncol(dt_P)
total_entries <- N*Time
ncolsM <- 32
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
age_0to2 <- age_3to5 <- numeric(nrow(dt_age))
for(i in 1:nrow(dt_age)){
  age_0to2[i] <- ifelse(dt_age[i,1]==1,1,0) 
  age_3to5[i] <- ifelse(dt_age[i,1]==2,1,0)
}

# X matrix
dt_X <- as.matrix(cbind(rep(1,50),age_0to2,age_3to5,
                        dt_snap[,1],dt_dpp[,1],dt_trt[,1])) 
colnames(dt_X) <- NULL

# Looking at time indexes for removed subjects (LTFU)
first_timepoint <- rep(1,N)
last_timepoint <- rep(7,N)
#last_timepoint[10] <- 5  # removed, non-Leish related (ALIVE)
last_timepoint[13] <- 4  # removed, non-Leish related (DEATH)
last_timepoint[18] <- 6  # removed, severe Leish (DEATH)
last_timepoint[25] <- 2  # removed, non-Leish related (DEATH)
last_timepoint[33] <- 6  # removed, severe Leish (DEATH)
last_timepoint[36] <- 4  # removed, severe Leish (DEATH)
#last_timepoint[44] <- 5  # removed, non-Leish related (ALIVE)
last_timepoint[47] <- 2  # removed, non-Leish related (DEATH)


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

# Entries to ignore due to lost-to-follow-up (reasons unknown)
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
# Nimble Code
########################################################################

BayesianModel <- nimbleCode({
  
  # PRIORS =============================================
  
  # Coefficients for Autoregressive Outcomes
  for(i in 1:ncolsM){
    betaP[i] ~ dnorm(0, sd = sigma.betaP)
    betaA[i] ~ dnorm(0, sd = sigma.betaA)
    betaD1[i] ~ dnorm(0, sd = sigma.betaD1)
    betaD2[i] ~ dnorm(0, sd = sigma.betaD2)
    betaD3[i] ~ dnorm(0, sd = sigma.betaD3)
    betaI1[i] ~ dnorm(0, sd = sigma.betaI1)
    betaI2[i] ~ dnorm(0, sd = sigma.betaI2)
    betaI3[i] ~ dnorm(0, sd = sigma.betaI3)
    betaR1[i] ~ dnorm(0, sd = sigma.betaR1)
    betaR2[i] ~ dnorm(0, sd = sigma.betaR2)
    betaR3[i] ~ dnorm(0, sd = sigma.betaR3)
  }
  
  # Coefficients for Non-time-varying Predictors
  for(i in 1:ncolsX){
    alphaP[i] ~ dnorm(0, sd = 1)
    alphaA[i] ~ dnorm(0, sd = 1)
    alphaD1[i] ~ dnorm(0, sd = 1)
    alphaD2[i] ~ dnorm(0, sd = 1)
    alphaD3[i] ~ dnorm(0, sd = 1)
    alphaI1[i] ~ dnorm(0, sd = 1)
    alphaI2[i] ~ dnorm(0, sd = 1)
    alphaI3[i] ~ dnorm(0, sd = 1)
    alphaR1[i] ~ dnorm(0, sd = 1)
    alphaR2[i] ~ dnorm(0, sd = 1)
    alphaR3[i] ~ dnorm(0, sd = 1)
  }
  
  # Variances of Coefficients of longitudinal outcomes
  sigma2.betaP ~ dgamma(1,1)
  sigma2.betaA ~ dgamma(1,1)
  sigma2.betaD1 ~ dgamma(1,1)
  sigma2.betaD2 ~ dgamma(1,1)
  sigma2.betaD3 ~ dgamma(1,1)
  sigma2.betaI1 ~ dgamma(1,1)
  sigma2.betaI2 ~ dgamma(1,1)
  sigma2.betaI3 ~ dgamma(1,1)
  sigma2.betaR1 ~ dgamma(1,1)
  sigma2.betaR2 ~ dgamma(1,1)
  sigma2.betaR3 ~ dgamma(1,1)
  sigma.betaP <- sqrt(sigma2.betaP)
  sigma.betaA <- sqrt(sigma2.betaA)
  sigma.betaD1 <- sqrt(sigma2.betaD1)
  sigma.betaD2 <- sqrt(sigma2.betaD2)
  sigma.betaD3 <- sqrt(sigma2.betaD3)
  sigma.betaI1 <- sqrt(sigma2.betaI1)
  sigma.betaI2 <- sqrt(sigma2.betaI2)
  sigma.betaI3 <- sqrt(sigma2.betaI3)
  sigma.betaR1 <- sqrt(sigma2.betaR1)
  sigma.betaR2 <- sqrt(sigma2.betaR2)
  sigma.betaR3 <- sqrt(sigma2.betaR3)
  
  # Variances for Coefficients of non-time-varying predictors
  sigma2.alphaP ~ dgamma(1,1)
  sigma2.alphaA ~ dgamma(1,1)
  sigma2.alphaD1 ~ dgamma(1,1)
  sigma2.alphaD2 ~ dgamma(1,1)
  sigma2.alphaD3 ~ dgamma(1,1)
  sigma2.alphaI1 ~ dgamma(1,1)
  sigma2.alphaI2 ~ dgamma(1,1)
  sigma2.alphaI3 ~ dgamma(1,1)
  sigma2.alphaR1 ~ dgamma(1,1)
  sigma2.alphaR2 ~ dgamma(1,1)
  sigma2.alphaR3 ~ dgamma(1,1)
  sigma.alphaP <- sqrt(sigma2.alphaP)
  sigma.alphaA <- sqrt(sigma2.alphaA)
  sigma.alphaD1 <- sqrt(sigma2.alphaD1)
  sigma.alphaD2 <- sqrt(sigma2.alphaD2)
  sigma.alphaD3 <- sqrt(sigma2.alphaD3)
  sigma.alphaI1 <- sqrt(sigma2.alphaI1)
  sigma.alphaI2 <- sqrt(sigma2.alphaI2)
  sigma.alphaI3 <- sqrt(sigma2.alphaI3)
  sigma.alphaR1 <- sqrt(sigma2.alphaR1)
  sigma.alphaR2 <- sqrt(sigma2.alphaR2)
  sigma.alphaR3 <- sqrt(sigma2.alphaR3)
  
  # Variances and Covariance Matrix
  sigma2P ~ dgamma(1,1)
  sigma2A ~ dgamma(1,1)
  sigmaP <- sqrt(sigma2P)
  sigmaA <- sqrt(sigma2A)
  PrecIR[1:dimIR,1:dimIR] ~ dwish(R[1:dimIR,1:dimIR], df = N)
  R[1:dimIR,1:dimIR] <- (N-dimIR-1)*covSampleIR[1:dimIR,1:dimIR]
  SigmaIR[1:dimIR,1:dimIR] <- inverse(PrecIR[1:dimIR,1:dimIR])
  
  # Latent Variables 
  for(i in 1:N){
    P[i,1] ~ dnorm(0, sd = 1)
    A[i,1] ~ dnorm(0, sd = 1)
    D[i,1] ~ dcat(probD[1:3])
  }
  probD[1:3] <- c(1/3, 1/3, 1/3)
  for(i in 1:N){
    for(k in 1:dimIR){
      IR[i,k,1] ~ dnorm(0, sd = 1)
    }
  }
  
  # Moving Average Components
  Prec.wIR[1:dimIR,1:dimIR] ~ dwish(R.wIR[1:dimIR,1:dimIR], df = N)
  Sigma.wIR[1:dimIR,1:dimIR] <- inverse(Prec.wIR[1:dimIR,1:dimIR])
  sigma2.wP ~ dgamma(1,1)
  sigma2.wA ~ dgamma(1,1)
  sigma2.wD ~ dgamma(1,1)
  sigma.wP <- sqrt(sigma2.wP) 
  sigma.wA <- sqrt(sigma2.wA) 
  sigma.wD <- sqrt(sigma2.wD) 
  thetaP ~ dnorm(0, sd = 1)
  thetaA ~ dnorm(0, sd = 1)
  thetaD ~ dnorm(0, sd = 1) 
  for(k in 1:(dimIR*dimIR)){
    thIR[k] ~ dnorm(0, sd = 1) 
  }
  thetaIR[1:dimIR,1:dimIR] <- nimMatrix(value = thIR[1:(dimIR*dimIR)],
                                        nrow = dimIR, ncol = dimIR)
  for(i in 1:N){
    for(t in 1:Time){
      wP[i,t] ~ dnorm(0, sd = sigma.wP)
      wA[i,t] ~ dnorm(0, sd = sigma.wA)
      wD[i,t] ~ dnorm(0, sd = sigma.wD)
      wIR[i,1:dimIR,t] ~ dmnorm(zero.wIR[1:dimIR], cov = Sigma.wIR[1:dimIR,1:dimIR])
    }
  }
  zero.wIR[1:dimIR] <- c(nimRep(0,dimIR))[1:dimIR]

  # MODEL COMPONENTS ===================================
  
  # Likelihood (Model)
  for(i in 1:N){
    for(t in (first_timepoint[i]+1):last_timepoint[i]){
      # Pathogen load (P)
      Censored[i,t] ~ dinterval(P[i,t], LOD)
      P[i,t] ~ dnorm(mean = muP[i,t], sd = sigmaP)
      
      # Antibody levels (A)
      A[i,t] ~ dnorm(mean = muA[i,t], sd = sigmaA)
      
      # Disease status (D)
      exp_term[i,t,1] <- exp(muD1[i,t])
      exp_term[i,t,2] <- exp(muD2[i,t])
      exp_term[i,t,3] <- exp(muD3[i,t])
      prob[i,t,1:3] <- exp_term[i,t,1:3] / (1 + sum(exp_term[i,t,1:3]))
      prob[i,t,4] <- 1 - sum(prob[i,t,1:3])
      D[i,t] ~ dcat(prob[i,t,1:4])
      
      # Inflammatory-Regulatory Responses (IR)
      IR[i,1:dimIR,t] ~ dmnorm(IRmean[i,1:dimIR,t-1], cov = SigmaIR[1:dimIR,1:dimIR])
      
    }
  }
  
  # Mean Components
  for(i in 1:N){
    for(t in (first_timepoint[i]+1):last_timepoint[i]){
      # Pathogen load (P)
      muP[i,t] <- inprod(M[i,1:ncolsM,t-1],betaP[1:ncolsM]) + 
        inprod(X[i,1:ncolsX],alphaP[1:ncolsX]) + wP[i,t] + thetaP*wP[i,t-1]
      
      # Antibody levels (A)
      muA[i,t] <- inprod(M[i,1:ncolsM,t-1],betaA[1:ncolsM]) + 
        inprod(X[i,1:ncolsX],alphaA[1:ncolsX]) + wA[i,t] + thetaA*wA[i,t-1]
      
      # Disease status (D)
      muD1[i,t] <- inprod(M[i,1:ncolsM,t-1],betaD1[1:ncolsM]) + 
        inprod(X[i,1:ncolsX],alphaD1[1:ncolsX]) + wD[i,t] + thetaD*wD[i,t-1]
      muD2[i,t] <- inprod(M[i,1:ncolsM,t-1],betaD2[1:ncolsM]) + 
        inprod(X[i,1:ncolsX],alphaD2[1:ncolsX]) + wD[i,t] + thetaD*wD[i,t-1]
      muD3[i,t] <- inprod(M[i,1:ncolsM,t-1],betaD3[1:ncolsM]) + 
        inprod(X[i,1:ncolsX],alphaD3[1:ncolsX]) + wD[i,t] + thetaD*wD[i,t-1]
      
      # Inflammatory-Regulatory Responses (IR)
      IRmean[i,1:dimIR,t-1] <- muIR[i,1:dimIR,t-1] + wIR[i,1:dimIR,t] + 
        (thetaIR[1:dimIR,1:dimIR] %*% wIR[i,1:dimIR,t-1])
      muIR[i,1:dimIR,t-1] <- c(inprod(M[i,1:ncolsM,t-1],betaI1[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaI1[1:ncolsX]),
                               inprod(M[i,1:ncolsM,t-1],betaI2[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaI2[1:ncolsX]),
                               inprod(M[i,1:ncolsM,t-1],betaI3[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaI3[1:ncolsX]),
                               inprod(M[i,1:ncolsM,t-1],betaR1[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaR1[1:ncolsX]),
                               inprod(M[i,1:ncolsM,t-1],betaR2[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaR2[1:ncolsX]),
                               inprod(M[i,1:ncolsM,t-1],betaR3[1:ncolsM]) + 
                                 inprod(X[i,1:ncolsX],alphaR3[1:ncolsX]))[1:dimIR]
    }
  }
  
  # Autoregressive Matrix of Longitudinal Outcomes
  for(i in 1:N){
    for(t in 1:Time){
      D1[i,t] <- 1*(D[i,t]==1)
      D2[i,t] <- 1*(D[i,t]==2)
      D3[i,t] <- 1*(D[i,t]==3)
      refD[i,t] <- (1-D1[i,t]-D2[i,t]-D3[i,t])
      for(k in 1:ncolsM){
        M[i,k,t] <- c(refD[i,t]*P[i,t],
                      refD[i,t]*A[i,t],
                      refD[i,t]*IR[i,1,t],
                      refD[i,t]*IR[i,2,t],
                      refD[i,t]*IR[i,3,t],
                      refD[i,t]*IR[i,4,t],
                      refD[i,t]*IR[i,5,t],
                      refD[i,t]*IR[i,6,t],
                      D1[i,t]*P[i,t],
                      D2[i,t]*P[i,t],
                      D3[i,t]*P[i,t],
                      D1[i,t]*A[i,t],
                      D2[i,t]*A[i,t],
                      D3[i,t]*A[i,t],
                      D1[i,t]*IR[i,1,t],
                      D2[i,t]*IR[i,1,t],
                      D3[i,t]*IR[i,1,t],
                      D1[i,t]*IR[i,2,t],
                      D2[i,t]*IR[i,2,t],
                      D3[i,t]*IR[i,2,t],
                      D1[i,t]*IR[i,3,t],
                      D2[i,t]*IR[i,3,t],
                      D3[i,t]*IR[i,3,t],
                      D1[i,t]*IR[i,4,t],
                      D2[i,t]*IR[i,4,t],
                      D3[i,t]*IR[i,4,t],
                      D1[i,t]*IR[i,5,t],
                      D2[i,t]*IR[i,5,t],
                      D3[i,t]*IR[i,5,t],
                      D1[i,t]*IR[i,6,t],
                      D2[i,t]*IR[i,6,t],
                      D3[i,t]*IR[i,6,t])[k]
      }
    }
  }
  
})
BayesianModel

parameters <- c("betaP", "alphaP", "sigmaP", "betaA", "alphaA", "sigmaA", 
                "betaD1", "betaD2", "betaD3", "alphaD1", "alphaD2", "alphaD3",
                "betaI1", "betaI2", "betaI3", "betaR1", "betaR2", "betaR3", 
                "alphaI1", "alphaI2", "alphaI3", "alphaR1", "alphaR2", 
                "alphaR3", "SigmaIR", "P", "D", "A", "IR", "wP", "wA", "wD", 
                "wIR", "thetaP", "thetaA", "thetaD", "thetaIR", "sigma.wP", 
                "sigma.wA", "Sigma.wIR", "sigma.wD")

nimbledata <- list(P = dt_Pc, 
                   A = dt_A,
                   D = dt_D,
                   IR = dt_IR,
                   Censored = Censored)

nimbleconstants <- list(N = N,
                        Time = Time,
                        first_timepoint = first_timepoint,
                        last_timepoint = last_timepoint,
                        ncolsM = ncolsM,
                        ncolsX = ncolsX,
                        dimIR = dimIR,
                        X = dt_X,
                        covSampleIR = covSampleIR,
                        LOD = logit_DLt)

nimbleinits <- function() {
  
  # Initial Values for Missing Values
  init.A <- matrix(NA,N,Time)
  init.P <- matrix(NA,N,Time)
  init.D <- matrix(NA,N,Time)
  init.IR <- nimArray(value = NA, dim = c(N,dimIR,Time))
  init.Censored <- matrix(NA,N,Time)
  
  init.A[missing_ind_A] <- rnorm(length(missing_ind_A))
  init.D[missing_ind_D] <-  sample(c(1,2,3,4),length(missing_ind_D),replace = TRUE)
  init.IR[missing_ind_IR] <- rnorm(length(missing_ind_IR))
  init.P[missing_ind_P] <- runif(length(missing_ind_P),logit_DLt+0.01,-1*logit_DLt/4)
  init.P[censored_idxP] <- runif(length(censored_idxP),2*logit_DLt,logit_DLt)
  init.Censored[missing_ind_P] <- rep(1,length(missing_ind_P))
  
  # All initial values
  list(betaP = rep(0,ncolsM),
       betaA = rep(0,ncolsM), 
       betaD1 = rep(0,ncolsM), 
       betaD2 = rep(0,ncolsM), 
       betaD3 = rep(0,ncolsM),
       betaI1 = rep(0,ncolsM), 
       betaI2 = rep(0,ncolsM),
       betaI3 = rep(0,ncolsM),
       betaR1 = rep(0,ncolsM),
       betaR2 = rep(0,ncolsM),
       betaR3 = rep(0,ncolsM),
       alphaP = runif(ncolsX),
       alphaA = runif(ncolsX),
       alphaD1 = runif(ncolsX),
       alphaD2 = runif(ncolsX),
       alphaD3 = runif(ncolsX),
       alphaI1 = runif(ncolsX),
       alphaI2 = runif(ncolsX),
       alphaI3 = runif(ncolsX),
       alphaR1 = runif(ncolsX),
       alphaR2 = runif(ncolsX),
       alphaR3 = runif(ncolsX),
       sigma2P = runif(1),
       sigma2A = runif(1),
       sigma2.wA = runif(1),
       sigma2.wP = runif(1),
       sigma2.wD = runif(1),
       sigma2.betaP = 1,
       sigma2.betaA = 1,
       sigma2.betaD1 = 1,
       sigma2.betaD2 = 1,
       sigma2.betaD3 = 1,
       sigma2.betaI1 = 1,
       sigma2.betaI2 = 1,
       sigma2.betaI3 = 1,
       sigma2.betaR1 = 1,
       sigma2.betaR2 = 1,
       sigma2.betaR3 = 1,
       sigma2.alphaP = 1,
       sigma2.alphaA = 1,
       sigma2.alphaD1 = 1,
       sigma2.alphaD2 = 1,
       sigma2.alphaD3 = 1,
       sigma2.alphaI1 = 1,
       sigma2.alphaI2 = 1,
       sigma2.alphaI3 = 1,
       sigma2.alphaR1 = 1,
       sigma2.alphaR2 = 1,
       sigma2.alphaR3 = 1,
       thetaP  = runif(1),
       thetaA  = runif(1),
       thetaD  = runif(1),
       thIR = runif(dimIR*dimIR),
       wP = matrix(runif(N*Time),N,Time),
       wA = matrix(runif(N*Time),N,Time),
       wD = matrix(runif(N*Time),N,Time),
       wIR = nimArray(value = runif(N*dimIR*Time), dim = c(N,dimIR,Time)),
       PrecIR = diag(dimIR),
       Prec.wIR = diag(dimIR),
       R.wIR = diag(dimIR),
       P = init.P,
       A = init.A,
       D = init.D,
       IR = init.IR,
       Censored = init.Censored)
}


########################################################################
# Running Code
########################################################################

time1 <- Sys.time() 
results_fit <- foreach(x = as.numeric(paste0("62022",1:ncore)), 
                       .packages = "nimble", .verbose = TRUE) %dopar% {
  nimbleMCMC(code = BayesianModel,
             constants = nimbleconstants,
             data = nimbledata,
             inits = nimbleinits,
             monitors = parameters,
             niter = 100000,
             nburnin = 50000,
             nchains = 1,
             thin = 2,
             progressBar = TRUE,
             samplesAsCodaMCMC = TRUE,
             setSeed = x)
}
time2 <- Sys.time() 
(runtime <- time2-time1)
saveRDS(results_fit, file = "results_fit.rds", version = 2)
stopImplicitCluster()
stopCluster(cl)
 
########################################################################
# MCMC Results
########################################################################
 
# Results
results_file <- readRDS(file = "results_fit.rds")
results_mcmc <- as.mcmc.list(lapply(1:3, function(x){as.mcmc(results_file[[x]])}))

parm.interest <- c("betaP", "alphaP", "sigmaP", "betaA", "alphaA", "sigmaA", 
                   "betaD1", "betaD2", "betaD3", "alphaD1", "alphaD2", "alphaD3",
                   "betaI1", "betaI2", "betaI3", "betaR1", "betaR2", "betaR3", 
                   "alphaI1", "alphaI2", "alphaI3", "alphaR1", "alphaR2", 
                   "alphaR3", "SigmaIR", "wP", "wA", "wD", 
                   "wIR", "thetaP", "thetaA", "thetaD", "thetaIR", "sigma.wP", 
                   "sigma.wA", "Sigma.wIR", "sigma.wD")


## PLOTS ====================

# Parameters
MCMCtrace(results_mcmc, 
          params = parm.interest,
          ISB = TRUE,
          pdf = TRUE,
          Rhat = TRUE, 
          iter = 20000)

# Latent quantities
MCMCtrace(results_mcmc, 
          params = c(vec_idx_miss_A,
                     vec_idx_miss_P,
                     vec_idx_miss_D,
                     vec_idx_miss_IR),
          ISB = FALSE,
          pdf = TRUE,
          Rhat = TRUE,
          iter = 20000)


