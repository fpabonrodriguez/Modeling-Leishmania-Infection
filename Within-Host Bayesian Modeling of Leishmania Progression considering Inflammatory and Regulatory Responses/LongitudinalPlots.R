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
                 "lubridate", "ggpubr", "stringr", "nimble", "gridExtra",
                 "igraph", "parallel", "doParallel", "MCMCvis")
invisible(install_load(my_packages))
nimbleOptions(verbose = FALSE)

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


ggbg <- function(color) {
  points(0, 0, pch=16, cex=1e6, col=color)
  grid(col="white", lty=1)
}

pcr_load <- as.data.frame(with(leish,cbind(ID,PCR_Load_TP1, PCR_Load_TP2, PCR_Load_TP3,
                                           PCR_Load_TP4, PCR_Load_TP5, PCR_Load_TP6,
                                           PCR_Load_TP7)))
antib <- as.data.frame(with(leish, cbind(ID, SLA.ELISA.OD.ratio_TP1,SLA.ELISA.OD.ratio_TP2,
                                         SLA.ELISA.OD.ratio_TP3,SLA.ELISA.OD.ratio_TP4,
                                         SLA.ELISA.OD.ratio_TP5,SLA.ELISA.OD.ratio_TP6,
                                         SLA.ELISA.OD.ratio_TP7)))

infl1 <- cbind(leish$ID,dt_I1)
infl2 <- cbind(leish$ID,dt_I2)
infl3 <- cbind(leish$ID,dt_I3)
reg1 <- cbind(leish$ID,dt_R1)
reg2 <- cbind(leish$ID,dt_R2)
reg3 <- cbind(leish$ID,dt_R3)
reg2$X.PD1..of.CD4.CD49dhi_4 <- rep(NA, nrow(reg2)) 
reg3$X.PD1..of.CD8.CD49dhi_4 <- rep(NA, nrow(reg3))

long_pcrload <- reshape(pcr_load, 
                        direction = "long",
                        varying = list(names(pcr_load)[2:8]),
                        v.names = "Value",
                        idvar = c("ID"),
                        timevar = "TimePoint",
                        times = 1:7)
long_pcrload <- long_pcrload[order(long_pcrload$ID),] 


long_antib <- reshape(antib, 
                      direction = "long",
                      varying = list(names(antib)[2:8]),
                      v.names = "Value",
                      idvar = c("ID"),
                      timevar = "TimePoint",
                      times = 1:7)
long_antib <- long_antib[order(long_antib$ID),] 


long_I1 <- reshape(infl1, 
                   direction = "long",
                   varying = list(names(infl1)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_I1 <- long_I1[order(long_I1$ID),] 


long_I2 <- reshape(infl2, 
                   direction = "long",
                   varying = list(names(infl2)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_I2 <- long_I2[order(long_I2$ID),] 


long_I3 <- reshape(infl3, 
                   direction = "long",
                   varying = list(names(infl3)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_I3 <- long_I3[order(long_I3$ID),] 


long_R1 <- reshape(reg1, 
                   direction = "long",
                   varying = list(names(reg1)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_R1 <- long_R1[order(long_R1$ID),] 


long_R2 <- reshape(reg2, 
                   direction = "long",
                   varying = list(names(reg2)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_R2 <- long_R2[order(long_R2$ID),] 


long_R3 <- reshape(reg3, 
                   direction = "long",
                   varying = list(names(reg3)[2:8]),
                   v.names = "Value",
                   idvar = c("ID"),
                   timevar = "TimePoint",
                   times = 1:7)
long_R3 <- long_R3[order(long_R3$ID),] 


g1 <- ggplot(data = long_pcrload, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(255,0,0,max=255,alpha=150,names = "red"), size=1) +
  labs(title = "Pathogen Load (P) \n(Number of parasites per mL of blood)",
       x="Time Point", y="Value")

g2 <- ggplot(data = long_antib, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0,102,0,max=255,alpha=150,names="darkgreen"), size=1) +
  labs(title = "Antibody Levels (A) \n (ELISA SLA OD ratio)",
       x="Time Point", y="Value")


long_I1$Value <- as.numeric(as.character(long_I1$Value))
long_R1$Value <- as.numeric(as.character(long_R1$Value))
long_I2$Value <- as.numeric(as.character(long_I2$Value))
long_R2$Value <- as.numeric(as.character(long_R2$Value))
long_I3$Value <- as.numeric(as.character(long_I3$Value))
long_R3$Value <- as.numeric(as.character(long_R3$Value))

longR2_reduced <- long_R2[long_R2$TimePoint %in% c(3,5),]
longR3_reduced <- long_R3[long_R3$TimePoint %in% c(3,5),]

df2 <- data.frame(
  id = unique(long_R2$`leish$ID`),
  x = rep(3,length(unique(long_R2$`leish$ID`))),
  xend = rep(5,length(unique(long_R2$`leish$ID`))), 
  y = longR2_reduced[longR2_reduced$TimePoint==3,]$Value,
  yend = longR2_reduced[longR2_reduced$TimePoint==5,]$Value  
)

df3 <- data.frame(
  id = unique(long_R3$`leish$ID`),
  x = rep(3,length(unique(long_R3$`leish$ID`))),
  xend = rep(5,length(unique(long_R3$`leish$ID`))), 
  y = longR3_reduced[longR3_reduced$TimePoint==3,]$Value,
  yend = longR3_reduced[longR3_reduced$TimePoint==5,]$Value  
)

g3 <- ggplot(data = long_I1, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0,0,1, 0.5), size=1) +
  labs(title = "Inflammatory Response (I1) \n (%IFN-gamma+ of CD4+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1)

g4 <- ggplot(data = long_R1, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0.6,0.2,0.9, 0.7), size=1) +
  labs(title = "Regulatory Response (R1) \n (%IL-10+ of CD4+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1)

g5 <- ggplot(data = long_I2, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0,0,1, 0.5), size=1) +
  labs(title = "Inflammatory Response (I2) \n (%CSFElo of CD4+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1)

g6 <- ggplot(data = long_R2, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0.6,0.2,0.9, 0.7), size=1) +
  labs(title = "Regulatory Response (R2) \n (%PD1+ of CD4+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1) + 
  geom_segment(data = df2, aes(x, y, xend = xend, yend = yend, group = id),
               alpha = 0.35, colour=rgb(0.6,0.2,0.9, 0.7))

g7 <- ggplot(data = long_I3, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0,0,1, 0.5), size=1) +
  labs(title = "Inflammatory Response (I3) \n (%CSFElo of CD8+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1) 

g8 <- ggplot(data = long_R3, aes(x = factor(TimePoint), y = Value, group = ID)) +
  geom_line(aes(color=ID), colour=rgb(0.6,0.2,0.9, 0.7), size=1) +
  labs(title = "Regulatory Response (R3) \n (%PD1+ of CD8+ CD49dhi)",
       x="Time Point", y="Percentage") + ylim(0,1) + 
  geom_segment(data = df3, aes(x, y, xend = xend, yend = yend, group = id),
               alpha = 0.35, colour=rgb(0.6,0.2,0.9, 0.7))


grid.arrange(g1,g3,g5,g7,
             g2,g4,g6,g8, ncol=4)

 