#######################################################################
#                                                                     #
#   Simulation 1 - Procedures for Simulating Aggregated Covariates    #
#                                                                     #
#######################################################################

## Loading packages
library(simstudy)
library(tictoc)
library(dplyr)
library(tidyr)

#### Converting logits to probabilities and vice versa ####
# can use plogis() instead
LogitToProb <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#############################
#### Defining Conditions ####

## Study Conditions
DataGenConds <- crossing(nsub = c(20, 60, 100),           # number of subjects per cluster
                         nclus = c(20, 60, 100),          # number of clusters
                         delta = .5, #c(.2, .5, .8),      # treatment effect size
                         icc = c(.05, .1, .2),            # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                         aggvar = c(.1,.3,.5,1),          # error/1-% sampled in observed L2 covariates
                         method = c("Formative.RE","Formative.Sample","Reflective.RE","Reflective.Sample"))   # data generation method

nreps <- 1000  # number of replications

############################

##########################
#### Generating data  ####


## Creating blank matrix for summary statistics for each data generation condition
DGSnames <- c("Prop_Z1","True_Logit","True_PS","Obs_Logit","Obs_PS",
              "Logit_Bias","Logit_MAE","Logit_RMSE","PS_Bias","PS_MAE","PS_RMSE",
              "Y_ICC1","Y_Variance","L1_Variance","L1_Cor",
              "True_L2_Variance","True_L2_Cor","Obs_L2_Variance","Obs_L2_Cor",
              "True.Obs_L2_Cor", "X_ICC1","X_ICC2","Convergence")
DGS <- matrix(-999, ncol = length(DGSnames) + ncol(DataGenConds), nrow = nrow(DataGenConds)) # Summary stats for generated datasets

## Variables consistent across all conditions
marg <- -1.0986                  # marginal probability of treatment;intercept of -1.0986 produces a marginal probability of .25; LogitToProb(-1.0986)
L1names <- paste0("x", 1:20)                             # names of L1 covariates
trueL2names <- paste(L1names[1:10], "c", sep = "_")      # names of true L2 covariates
obsL2names <- paste(trueL2names, "o", sep = "_")         # names of aggregated L2 covariates

#### Starting Generation ####
set.seed(48226) # Seed for reproducing results

for(con in 1:nrow(DataGenConds)){
  
  tic(as.character(con))          # Record time to run all replications for each condition

  nsub <- DataGenConds[con,][[1]]  # Number of Level 1 units (i.e., subjects)
  nclus <- DataGenConds[con,][[2]]  # Number of Level 2 units (i.e., clusters)
  ntot <- nclus*nsub                # Total sample size
  icc <- DataGenConds[con,][[4]]   # Intraclass correlation coefficient
  aggvar <- DataGenConds[con,][[5]]  # error for creating observed L2 covariates
  method <- DataGenConds[con,][[6]] # method for generating L1, true L2, and obs L2 covariates
  tau00 <- icc                     # between cluster (L2) variance ($\tau_{00}$)
  sigma2 <- 1 - tau00              # within cluster (L1) variance ($\sigma^2$)

  
  ## Stats for each data generation replication
  DataGenRepStats <- matrix(-999, ncol = length(DGSnames), nrow = nreps)
  
  ###################################################################
  #### Testing various ways of creating aggregated L2 covariates ####
  
  for(r in 1:nreps){
    
    if(method == "Formative.RE"){
      #### V1: Generate L1, aggregate true L2, then add random error for observed L2 ####
      
      ## Level 1 covariates - 20 continuous covariates from a standard normal distribution with correlation of .2
      GenData <- genCorData(n = ntot, mu = rep(0, 20), sigma = sqrt(sigma2), rho = 0.2, corstr = "cs",  # compound symmetry correlation structure based on sigma and rho
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
      
    } else if(method == "Formative.Sample"){
      
      #### V2: Generate L1, aggregate true L2, sample from L1 and aggregate observed L2 ####
      
      ## Level 1 covariates - 20 continuous covariates from standard normal distribution with correlation of .2
      GenL1 <- genCorData(n = ntot, mu = rep(0, 20), sigma = sqrt(sigma2), rho = 0.2, corstr = "cs", 
                         cnames = L1names, idname = "sid") %>%
        mutate(rij = rnorm(n = ntot, mean = 0, sd = sqrt(sigma2)),    # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\sigma^2$ = 1)
               cid = as.character(rep(1:nclus, each = nsub))) %>%     # creating cluster ID
        group_by(cid) %>%
        mutate_at(vars(x1:x20), ~. + rnorm(n = 1, mean = 0, sd = sqrt(tau00))) %>%    # Adding cluster variation to maintain ICC(1)
        mutate_at(vars(x1:x10), list(c = ~mean(.))) %>%               # aggregating L1 to true L2 value
        mutate(zuj = rlogis(n = 1, location = 0, scale = 1),         # random error for PS model (probability of treatment) from logistic distribution with mean = 0  and variance of $\pi^2 / 3$
               uj = rnorm(n = 1, mean = 0, sd = sqrt(tau00)))        # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\tau_{00} depends on condition)
      
      ## Random sample within each cluster and aggregate for "observed" L2 value  
      GenL2 <- GenL1 %>%
        group_by(cid) %>%
        sample_frac(size = ifelse(aggvar == 1, .3, 1-aggvar), replace = FALSE) %>%
        mutate_at(vars(x1:x10), list(c_o = ~mean(.))) %>%
        select(cid, x1_c_o:x10_c_o) %>%
        unique()
      
      ## Join observed covariates to initial data and calculate probability of treatment and outcome
      GenData <- GenL1 %>%
        left_join(GenL2, by = "cid") %>%
        ungroup()
      
    } else if(method == "Reflective.RE"){
      
      #### V3: Generate true L2, generate L1 from L2, then aggregate observed L2 plus random error ####
      
      ## Level 2 covariates - continuous covariates from normal distribution with correlation of .2
      GenL2 <- genCorData(n = nclus, mu = rep(0,20), sigma = sqrt(tau00), rho = .2, corstr = "cs",
                         cnames = paste0(L1names, "_c"), idname = "cid") %>%    # trueL2names; Generate 20 L2 covariates to maintain correlation amongst all L1 instead of half being uncorrelated with the other
        mutate(zuj = rlogis(n = nclus, location = 0, scale = 1),                # random error for PS model (probability of treatment) from logistic distribution with mean = 0  and variance of $\pi^2 / 3$
               uj = rnorm(nclus, mean = 0, sd = sqrt(tau00)))                   # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\tau_{00} depends on condition)
      
      ## Generate L1 and observed L2 covariates and calculate probability of treatment and outcome
      GenData <- GenL2[rep(seq_len(nrow(GenL2)), each = nsub), ] %>%        # repeating true L2 covariates nsub times
        mutate(sid = c(1:ntot)) %>%                                         # adding subject id
        group_by(cid) %>%
        mutate_at(vars(x1_c:x20_c), list(t = ~rnorm(nsub, mean = mean(.), sd = sqrt(sigma2)))) %>%  # Generate 20 L1 variables based on mean for each cluster
        rename_at(vars(x1_c_t:x20_c_t), ~stringr::str_remove(., "_c_t")) %>%         
        mutate_at(vars(x1_c:x10_c), list(o = ~. + rnorm(1,mean = 0,sd = aggvar))) %>%          # Adding random variation to true cluster mean; selects a random number for each cluster for each covariate
        ungroup() %>%
        mutate(rij = rnorm(n = ntot, mean = 0, sd = sqrt(sigma2)))              # random L1 error for outcome model
        # inner_join(GenL1, by = "sid")
      
    } else if(method == "Reflective.Sample"){
      
      #### V4: Generate true L2, generate L1 from L2, sample from L1 then aggregate observed L2 ####
      
      ## Level 2 covariates - continuous covariates from standard normal distribution with correlation of .2
      GenL2 <- genCorData(n = nclus, mu = rep(0,20), sigma = sqrt(tau00), rho = .2, corstr = "cs",
                         cnames = paste0(L1names, "_c"), idname = "cid") %>%    # Generate 20 L2 covariates to maintain correlation amongst all L1 instead of half being uncorrelated with the other
        mutate(zuj = rlogis(n = nclus, location = 0, scale = 1),                # random error for PS model (probability of treatment) from logistic distribution with mean = 0  and variance of $\pi^2 / 3$
               uj = rnorm(nclus, mean = 0, sd = sqrt(tau00)))                   # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\tau_{00} depends on condition)
      
      ## Generate L1 covariates
      AggL1 <- GenL2[rep(seq_len(nrow(GenL2)),each = nsub),] %>%
        mutate(sid = c(1:ntot)) %>%
        group_by(cid) %>%
        mutate_at(vars(x1_c:x20_c), list(t = ~rnorm(nsub, mean = mean(.), sd = sqrt(sigma2)))) %>%  # Generate 20 L1 variables based on mean for each cluster
        rename_at(vars(x1_c_t:x20_c_t), ~stringr::str_remove(., "_c_t")) %>%         
        ungroup() %>%
        mutate(rij = rnorm(n = ntot, mean = 0, sd = sqrt(sigma2)))              # random L1 error for outcome model
      
      ## Random sample within each cluster and aggregate for "observed" value  
      AggL2 <- AggL1 %>%
        group_by(cid) %>%
        sample_frac(size = ifelse(aggvar == 1, .3, 1-aggvar), replace = FALSE) %>%
        mutate_at(vars(x1:x10), list(c_o = ~mean(.))) %>%
        select(cid,x1_c_o:x10_c_o) %>%
        unique()
      
      ## Join datasets and calculate probability of treatment and outcome value
      GenData <- AggL1 %>%
        # inner_join(GenL1, by = "sid") %>%
        left_join(AggL2, by = "cid") %>%
        ungroup()
      
      
    } else {stop("something went wrong")}
    
    
    #### Calculating Y, treatment assignment, true and observed log PS, and ICCs ####
    SampData <- GenData %>%
      mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +              # True logit of the probability of treatment exposure (i.e. propensity score)
               .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + zuj,     
             z = ifelse(logitz1 > 0, 1, 0),                                                 # Determining treatment exposure
             Yij = .5*x1 + .5*x2 + .5*x3 + .5*x4 + .5*x5 + .5*x6 + .5*x7 + .5*x8 + .5*x9 + .5*x10 +            # Generating outcome values
               .5*x11 + .5*x12 + .5*x13 + .5*x14 + .5*x15 + .5*x16 + .5*x17 + .5*x18 + .5*x19 + .5*x20 + rij + # Level 1 covariates
               .5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c + 
               .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + .5*z + uj)                                # Level 2 covariates
    
    ## Estimating observed PS
    # need to investigate why some models do not converge, which then affects the bias and RMSE calculations; Not sure if this is actually a problem         
    PS.mod <- glm(as.formula(paste0("z ~ 1 + ",paste(obsL2names,collapse = " + "))), data = SampData, family = binomial("logit"))
    
    ## Adding estimated PS to dataframe
    SampData$PS <- fitted(PS.mod)                           # observed PS 
    SampData$Logit <- predict(PS.mod)                       # observed logit of the PS
    SampData$TruePS <- LogitToProb(SampData$logitz1)        # true PS
    SampData$LogitDiff <- SampData$Logit - SampData$logitz1 # Error in observed (estimated) and true logit of the PS
    SampData$PSDiff <- SampData$PS - SampData$TruePS        # Error in observed and true PS
    
    ## Calculating ICC(1) and (2)
    ICC1 <- purrr::map_dbl(paste0("x",1:10), ~ICC::ICCbareF(factor(cid), quo_name(.x), SampData)) # ICC(1) of aggregated covariates
    ICC2 <- (nsub*ICC1) / (1 + (nsub - 1)*ICC1)

    #### Description of characteristics for each replication #### 
    DataGenRepStats[r,1] <- mean(SampData$z)                             # proportion with treatment exposure
    DataGenRepStats[r,2] <- mean(SampData$logitz1)                       # mean true logit of treatment exposure
    DataGenRepStats[r,3] <- mean(SampData$TruePS)                        # mean true probability of treatment exposure
    DataGenRepStats[r,4] <- mean(SampData$Logit)                         # mean observed logit of treatment exposure
    DataGenRepStats[r,5] <- mean(SampData$PS)                            # mean observed probability of treatment exposure
    DataGenRepStats[r,6] <- mean(SampData$LogitDiff)                            # bias of logit
    DataGenRepStats[r,7] <- mean(abs(SampData$LogitDiff))                       # mae of logit
    DataGenRepStats[r,8] <- sqrt(mean(SampData$LogitDiff^2))                    # rmse of logit
    DataGenRepStats[r,9] <- mean(SampData$PSDiff)                               # bias of PS
    DataGenRepStats[r,10] <- mean(abs(SampData$PSDiff))                         # mae of PS
    DataGenRepStats[r,11] <- sqrt(mean(SampData$PSDiff^2))                      # rmse of PS
    DataGenRepStats[r,12] <- ICC::ICCbareF(factor(cid), Yij, SampData)          # ICC(1) of outcome Yij
    DataGenRepStats[r,13] <- var(SampData$Yij)                                  # variance of the outcome
    DataGenRepStats[r,14] <- cov(SampData[,L1names]) %>% diag() %>% mean()              # mean variance of L1 covariates
    DataGenRepStats[r,15] <- cor(SampData[,L1names]) %>% .[lower.tri(.)] %>%    # mean correlation of L1 covariates
      psych::fisherz() %>% mean()
    DataGenRepStats[r,16] <- cov(SampData[,trueL2names]) %>% diag() %>% mean()         # mean variance of true L2 covariates
    DataGenRepStats[r,17] <- cor(SampData[,trueL2names]) %>% .[lower.tri(.)] %>%       # mean correlation of true L2 covariates
      psych::fisherz() %>% mean()
    DataGenRepStats[r,18] <- cov(SampData[,obsL2names]) %>% diag() %>% mean()           # mean variance of observed L2 covariates
    DataGenRepStats[r,19] <- cor(SampData[,obsL2names]) %>% .[lower.tri(.)] %>%         # mean correlation of observed L2 covariates
      psych::fisherz() %>% mean()
    DataGenRepStats[r,20] <- purrr::map2_dbl(.x = trueL2names,.y = obsL2names,~cor(SampData[,.x],SampData[,.y])) %>% # correlation b/t true and observed L2 covariates
      psych::fisherz() %>% mean()
    DataGenRepStats[r,21] <- mean(ICC1)                                   # mean ICC(1) of L1 X covariates
    DataGenRepStats[r,22] <- mean(ICC2)                                   # mean ICC(2) of L2 X covariates
    DataGenRepStats[r,23] <- ifelse(PS.mod$converged == TRUE, 1, 0)       # Did the PS model converge?
    
  }
  
  ## logging time to run condition
  toc(quiet = TRUE, log = TRUE)
  
  #### Sample characteristics averaged across replications ####
  colnames(DataGenRepStats) <- DGSnames
  DataGenStats <- data.frame(DataGenRepStats) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    bind_cols(DataGenConds[con, ])
  
  DGS[con, ] <- as.matrix(DataGenStats[1, ])
  
  if(con == 1){
    Con1SampDat <- SampData
    Con1DataGenRepStats <- DataGenRepStats
  }
}

#######################################################

#################################
#### Saving Initial Results  ####

colnames(DGS) <- colnames(DataGenStats)
DGSbyCond <- data.frame(DGS, stringsAsFactors = FALSE) %>%
  filter(nsub != -999) %>%                                     # Quality control check; should have all 432 rows after filtering, indicating all conditions were completed
  mutate_at(vars(one_of(DGSnames)), as.numeric) %>%            # Converting from character to numeric
  mutate_at(vars(contains("_Cor")), ~psych::fisherz2r(.))      # Converting z scores back to correlations

#### Extracting timing ####
# Character vector
DataGenTimingLog <- tic.log(format = TRUE) %>% unlist()

# Converted to dataframe and cleaned for analysis
DataGenTLogDF <- data.frame(temp = DataGenTimingLog, stringsAsFactors = FALSE) %>%
  separate(temp, c("Row", "Time"), sep = ": ") %>%
  mutate(Time = as.numeric(gsub(" sec elapsed", "", Time, fixed = TRUE))) %>%
  bind_cols(DataGenConds)

# Total time (in min) to run the simulation
sum(DataGenTLogDF$Time) / 60    

# Clearing time log
tic.clearlog()

#### Saving Simulation Results ####

## Intermediate information from condition 432 (aggvar = 1; icc = .2; nclus = 100; nsub = 100; method = "Reflective.Sample")
Con432SampDat <- SampData # saved sample data from last rep of condition 432 
Con432DataGenRepStats <- data.frame(DataGenRepStats, stringsAsFactors = FALSE) # saved summary statistics for each rep of condition 432

save(Con432SampDat, Con432DataGenRepStats,
     Con1SampDat, Con1DataGenRepStats,
     DGSbyCond, DataGenTLogDF, DGSnames,
     file = "Sim1_ThousandReps.RData")

##################################################

