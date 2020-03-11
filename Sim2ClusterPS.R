##############################################################################
#                                                                            #
#   Simulation 2 - Propensity Score Analysis with Cluster-Level Treatment    #
#                                                                            #
##############################################################################

################################
#### Packages and Functions ####

## Packages
library(simstudy)
library(tictoc)
library(dplyr)
library(tidyr)
library(MatchIt)
library(WeightIt)
library(cobalt)
library(lme4)

#### Converting logits to probabilities and vice versa ####
LogitToProb <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


##########################################


##########################################
#### Defining Conditions and Outcomes ####

## Study Conditions
DataGenConds <- crossing(nsub = c(20, 60, 100),           # number of subjects per cluster
                         nclus = c(20, 60, 100),          # number of clusters
                         icc = c(.05, .1, .2),            # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                         aggvar = c(.1,.3,.5,1),          # error/1-% sampled in observed L2 covariates
                         psmod = c("Cluster","Ignore","Fixed","Random"),   # PS model specification
                         mw = c("Matching", "Weighting")) # PS conditioning method

nreps <- 10  # number of replications

## Simulation dependent variables
DVnames <- c("PS_Converged", "Y_Converged","Abs_Std_Diff", "Variance_Ratio", "Balanced_ASD", "Balanced_VR",
             "PS_Bias", "PS_MAE", "PS_RMSE", "Matched", "TE_Bias", "TE_MAE", "TE_RMSE")
## Creating blank matrix for summary statistics for each data generation condition
DVSummaryStats <- matrix(-999, ncol = length(DVnames) + ncol(DataGenConds), nrow = nrow(DataGenConds)) # Summary stats for generated datasets


############################


#############################
#### Running Simulation ####

set.seed(48226) # Seed for reproducing results; 1213611 if using diff seed

## Variables consistent across all conditions
delta <- .5                                         # Average Treatment Effect ($\delta$)
marg <- -1.0986                                     # marginal probability of treatment;intercept of -1.0986 produces a marginal probability of .25; LogitToProb(-1.0986)
L1names <- paste0("x",1:20)                         # names of Level 1 covariates
trueL2names <- paste(L1names[1:10], "c", sep = "_") # names of true level 2 covariates
obsL2names <- paste(trueL2names, "o", sep = "_")    # names of observed level 2 covariates


#### Starting Simulation ####

con <- 1

for(con in 1:nrow(DataGenConds)){
  
  tic(as.character(con))          # Record time to run all replications for each condition
  
  nsub <- DataGenConds[con,][[1]]   # Number of Level 1 units (i.e., subjects)
  nclus <- DataGenConds[con,][[2]]  # Number of Level 2 units (i.e., clusters)
  ntot <- nclus*nsub                # Total sample size
  tau00 <- DataGenConds[con,][[3]]    # between cluster (L2) variance ($\tau_{00}$)
  sigma2 <- 1 - tau00               # within cluster (L1) variance ($\sigma^2$)
  aggvar <- DataGenConds[con,][[4]]  # error for creating observed L2 covariates
  psmod <- DataGenConds[con,][[5]]  # PS model specification
  mw <- DataGenConds[con,][[6]]     # PS conditioning method

  ## Outcome (DV) values for each replication of a condition
  DVRepStats <- matrix(-999, ncol = length(DVnames), nrow = nreps)
  
  for(r in 1:nreps){
  
  ##########################
  #### Generating data  ####
  
    #### V1: Generate L1, aggregate true L2, then add random error for observed L2 ####
    
    ## Level 1 covariates - 20 continuous covariates from a standard normal distribution with correlation of .2
    GenData <- genCorData(n = ntot, mu = rep(0, 20), sigma = sqrt(sigma2), rho = 0.2, corstr = "cs", 
                          cnames = L1names, idname = "sid") %>%
      mutate(rij = rnorm(n = ntot, mean = 0, sd = sqrt(sigma2)),                   # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\sigma^2$ = 1)
             cid = as.character(rep(1:nclus, each = nsub))) %>%                    # creating cluster ID
      group_by(cid) %>%
      mutate_at(vars(x1:x20), ~. + rnorm(n = 1, mean = 0, sd = sqrt(tau00))) %>%    # Adding cluster variation to maintain ICC(1)
      mutate_at(vars(x1:x10), list(c = ~mean(.))) %>%                               # aggregating L1 to true L2 value
      mutate_at(vars(x1_c:x10_c), list(o = ~. + rnorm(n = 1, mean = 0, sd = aggvar))) %>%  # Adding random variation to true cluster mean; selects a random number for each cluster for each covariate
      mutate(zuj = rlogis(n = 1, location = 0, scale = 1),                         # random error for PS model (probability of treatment) from logistic distribution with mean = 0  and variance of $\pi^2 / 3$
             uj = rnorm(n = 1, mean = 0, sd = sqrt(tau00))) %>%                    # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\tau_{00} depends on condition)
      ungroup()
    
  
  #### Calculating true PS (and logit), treatment assignment (z), and outcome (Y) ####
  SampData <- GenData %>%
    mutate(TrueLogit = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +                              # True logit of the probability of treatment exposure (i.e. propensity score)
             .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + zuj,
           TruePS = LogitToProb(TrueLogit),                                                         # true PS
           z = ifelse(TrueLogit > 0, 1, 0),                                                                  # Determining treatment exposure
           Yij = .5*x1 + .5*x2 + .5*x3 + .5*x4 + .5*x5 + .5*x6 + .5*x7 + .5*x8 + .5*x9 + .5*x10 +            # Generating outcome values
             .5*x11 + .5*x12 + .5*x13 + .5*x14 + .5*x15 + .5*x16 + .5*x17 + .5*x18 + .5*x19 + .5*x20 + rij + # Level 1 covariates
             .5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c + 
             .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + delta*z + uj,                                # Level 2 covariates
           cid = factor(cid))                                                                                # converting cluster id from character to factor
           
    #### Estimating observed PS ####
  
  if(psmod == "Cluster"){ 
    
    ## Model appraises treatment at the cluster-level and only uses L2 covariates
    PS.mod <- glm(as.formula(paste0("z ~ 1 + ", paste(obsL2names, collapse = " + "))),
                  data = SampData, family = binomial("logit"))
    
    # Did the PS model converge?
    PSModConverge <- ifelse(PS.mod$converged == TRUE, 1, 0)
    
  } else if(psmod == "Ignore"){
    
    ## Model ignores clusters and uses L1 and L2 covariates, the latter conceptually treated as L1 covariates
    PS.mod <- glm(as.formula(paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "))),
                  data = SampData, family = binomial("logit"))
    
    # Did the PS model converge?
    PSModConverge <- ifelse(PS.mod$converged == TRUE, 1, 0)
  
  } else if(psmod == "Fixed"){
    
    ## Model includes fixed effects for clusters and L1 covariates
    PS.mod <- glm(as.formula(paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + cid")),
                  data = SampData, family = binomial("logit"))
    
    # Did the PS model converge?
    PSModConverge <- ifelse(PS.mod$converged == TRUE, 1, 0)
    
  } else if(psmod == "Random"){
    
    ## Model includes random effects for clusters and covariates at both L1 and L2
    # default iterations is 10,000 I believe
    # use optimizers? options are: http://svmiller.com/blog/2018/06/mixed-effects-models-optimizer-checks/
    # check for false positives?
    # set nAGQ=0 instead of the default nAGQ = 1; makes the laplace approximation less precise; see here for further explanation: https://stats.stackexchange.com/questions/77313/why-cant-i-match-glmer-family-binomial-output-with-manual-implementation-of-g
    PS.mod <- glmer(as.formula(paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "), " + (1|cid)")),
                    data = SampData, family = binomial("logit"), control=glmerControl(optimizer="bobyqa"), nAGQ = .5) 
    
    # Did the PS model converge?
    PSModConverge <- PS.mod@optinfo$conv$opt
    
  } else {stop("something went wrong")}
  
  #### Adding estimated PS to dataframe ####
  SampData$PS <- fitted(PS.mod)                             # observed PS 
  SampData$Logit <- predict(PS.mod)                         # observed logit of the PS
  SampData$LogitDiff <- SampData$Logit - SampData$TrueLogit # Error in observed (estimated) and true logit of the PS
  SampData$PSDiff <- SampData$PS - SampData$TruePS          # Error in observed and true PS
  
  #########################################
  
  ######################################
  #### Analyzing Generated Datasets ####
  
  
  #### Conditioning on the PS ####
  
  #### Assessing Balance ####
  
  #### Estimating Treatment Effect ####
  
  #######################################
  
  #### Dependent Variable Values for each replication ####
  
  }

  ## logging time to run replication
  toc(quiet = TRUE, log = TRUE)
  
  # #### DVs averaged across replications ####
  #
  # colnames(DataGenRepStats) <- DGSnames
  # DataGenStats <- data.frame(DataGenRepStats) %>%
  #   summarize_all(mean, na.rm = TRUE) %>%
  #   bind_cols(DataGenConds[con, ])
  # 
  # DGS[con, ] <- as.matrix(DataGenStats[1, ])
  # 
  # if(con == 1){
  #   Con1SampDat <- SampData
  #   Con1DataGenRepStats <- DataGenRepStats
  # }
  
}