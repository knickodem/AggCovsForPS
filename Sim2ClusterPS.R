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

## Is Get_BalanceTableForSim actually necessary or just helps me visualize prior to running the code? Not necessary

#### Tailoring bal.tab output for my purposes ####
# UnAdj refers to unadjusted differences vs adjusted (matched, weighted)
# Get_BalanceTableForSim <- function(btab){
#   
#   TheBT <- btab[[1]] %>%
#     tibble::rownames_to_column("Covariate") %>%
#     select(Covariate, Type, starts_with("Diff"), starts_with("V.Ratio"))
#   names(TheBT) <- stringr::str_replace_all(names(TheBT), "\\.", "_")
#   
#   return(TheBT)
# }


##########################################


##########################################
#### Defining Conditions and Outcomes ####

## Study Conditions
DataGenConds <- crossing(nsub = c(20, 60, 100),            # number of subjects per cluster
                         nclus = c(60, 100, 140),           # number of clusters
                         icc = c(.05, .1, .2),             # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                         aggvar = c(.1, .3, .5, 1),        # error/1-% sampled in observed L2 covariates
                         psmod = c("Cluster","Subject"),   # PS estimation and conditioning at cluster or subject level
                         mw = c("Matching", "Weighting"))  # PS conditioning method

nreps <- 10  # number of replications

## Simulation dependent variables
DVnames <- c("Prop_Treated", "PS_Converged", # Proportion in treatment group and Yes/No for PS model convergence
             "PS_Bias", "PS_MAE", "PS_RMSE", "Logit_Bias", "Logit_MAE", "Logit_RMSE",          # Bias and Error in estimation of PS and logit of PS
             "Analytic_Subjects", "Analytic_Clusters", "True.Agg_L2_Cor", "Y_ICC1", "X_ICC1", "X_ICC2",    # Characteristics of sample used in outcome analysisNum subjects & clusters and correlation b/t true and obs L2 in sample used in outcome model
             "Un.ASD", "Un.VR", "Un.ASD_Balanced_Count", "Un.VR_Balanced_Count",               # For unadjusted sample - mean Absolute standardized difference & Variance ratio; Count of the 30 covariates balanced (on ASD or VR) at subject level
             "Adj.ASD", "Adj.VR", "Adj.ASD_Balanced_Count", "Adj.VR_Balanced_Count",           # For adjusted sample - mean ASD & VR; Count of the 30 covariates balanced (on ASD or VR) at subject level
             "Y_Converged", "TE_Bias", "TE_MAE", "TE_RMSE")                                    # Yes/No for Outcome model convergence; If yes, Bias, MAE and RMSE for estimation of the treatment effect

## Creating blank matrix for summary statistics for each data generation condition
DVSummaryStats <- matrix(-999, ncol = length(DVnames) + ncol(DataGenConds), nrow = nrow(DataGenConds)) # Summary stats for generated datasets


############################


#############################
#### Running Simulation ####

set.seed(1213611) # Seed for reproducing results; 1213611 if using diff seed; 48226 if same

## Variables consistent across all conditions
delta <- .50                                        # Average Treatment Effect ($\delta$)
marg <- -1.386294                                   # marginal probability of treatment;intercept of -1.0986 produces a marginal probability of .25; stats::plogis(-1.0986); -1.386294 = .2
L1names <- paste0("x", 1:20)                         # names of Level 1 covariates
trueL2names <- paste(L1names[1:10], "c", sep = "_") # names of true level 2 covariates
obsL2names <- paste(trueL2names, "o", sep = "_")    # names of observed level 2 covariates



#### Starting Simulation ####

con <- 1
# r <- 1

for(con in 1:nrow(DataGenConds)){
  
  tic(as.character(con))          # Record time to run all replications for each condition
  
  ## Defining characteristics of the simulation condition
  nsub <- DataGenConds[con,][[1]]     # Number of Level 1 units (i.e., subjects)
  nclus <- DataGenConds[con,][[2]]    # Number of Level 2 units (i.e., clusters)
  ntot <- nclus*nsub                  # Total sample size
  tau00 <- DataGenConds[con,][[3]]    # between cluster (L2) variance ($\tau_{00}$)
  sigma2 <- 1 - tau00                 # within cluster (L1) variance ($\sigma^2$)
  aggvar <- DataGenConds[con,][[4]]   # error for creating observed L2 covariates
  psmod <- DataGenConds[con,][[5]]    # PS estimation and equating at cluster or subject level
  mw <- DataGenConds[con,][[6]]       # PS equating method
  
  ## Dependent Variable (DV) values for each replication of a condition
  DVReplicationStats <- matrix(-999, ncol = length(DVnames), nrow = nreps)
  
  for(r in 1:nreps){
    
    ##########################
    #### Generating data  ####
    
    #### Formative-RE: Generate L1, aggregate true L2, then add random error for observed L2 ####
    
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
             TruePS = stats::plogis(TrueLogit),                                                                  # true PS
             z = ifelse(TrueLogit > 0, 1, 0),                                                                  # Determining treatment exposure
             Yij = .5*x1 + .5*x2 + .5*x3 + .5*x4 + .5*x5 + .5*x6 + .5*x7 + .5*x8 + .5*x9 + .5*x10 +            # Generating outcome values
               .5*x11 + .5*x12 + .5*x13 + .5*x14 + .5*x15 + .5*x16 + .5*x17 + .5*x18 + .5*x19 + .5*x20 + rij + # Level 1 covariates
               .5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c + 
               .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + delta*z + uj,                                # Level 2 covariates
             cid = factor(cid))                                                                                # converting cluster id from character to factor
    
    #### Estimating observed PS ####
    
    if(psmod == "Cluster"){ 
      
      ## Model appraises treatment at the cluster-level and only uses L2 covariates
      ThePSModel <- paste0("z ~ 1 + ", paste(obsL2names, collapse = " + "))
      PS.mod <- glm(as.formula(ThePSModel), data = SampData, family = binomial("logit"))
      
    } else if(psmod == "Subject"){
      
      ## Model ignores clusters and uses L1 and L2 covariates, the latter conceptually treated as L1 covariates
      ThePSModel <- paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "))
      PS.mod <- glm(as.formula(ThePSModel), data = SampData, family = binomial("logit"))
      
    } else {stop("something went wrong")}
    
    ## Did the PS model converge?
    PSModConverge <- ifelse(PS.mod$converged == TRUE, 1, 0)
    
    #### Adding estimated PS to dataframe ####
    SampData$PS <- fitted(PS.mod)                             # observed PS 
    SampData$Logit <- predict(PS.mod)                         # observed logit of the PS
    SampData$LogitDiff <- SampData$Logit - SampData$TrueLogit # Error in observed (estimated) and true logit of the PS
    SampData$PSDiff <- SampData$PS - SampData$TruePS          # Error in observed and true PS
    
    #### Calculating the pre-conditioning DVs that utilize the full sample ####
    DVReplicationStats[r,1] <- mean(SampData$z)                   # proportion with treatment exposure
    DVReplicationStats[r,2] <- PSModConverge                      # PS model convergence
    DVReplicationStats[r,3] <- mean(SampData$PSDiff)              # bias of PS
    DVReplicationStats[r,4] <- mean(abs(SampData$PSDiff))         # mae of PS
    DVReplicationStats[r,5] <- sqrt(mean(SampData$PSDiff^2))      # rmse of PS
    DVReplicationStats[r,6] <- mean(SampData$LogitDiff)           # bias of logit
    DVReplicationStats[r,7] <- mean(abs(SampData$LogitDiff))      # mae of logit
    DVReplicationStats[r,8] <- sqrt(mean(SampData$LogitDiff^2))   # rmse of logit
    
    #################################
    
    #### Overlap assumption ####
    ## Removing cases where PS is equal to 1 or 0
    AnalyticSample <- SampData %>%
      filter(PS > .005 & PS < .995)
    
    if(PSModConverge == 0 | nrow(AnalyticSample) == 0){
      
      ## Don't run conditioning and outcome model. Make remaining DV values NA
      DVReplicationStats[r,9:26] <- NA
      
    } else {
      
      ######################################
      #### Analyzing Generated Datasets ####
      
      #### Conditioning on the PS ####
      
      if(mw == "Matching"){
        if(psmod == "Cluster"){
          
          ## Extracting cluster-level data for matching; i.e., if nclus = 100, nrow(ClusterData) = 100 unless a cluster had a PS == 0|1
          ClusterData <- AnalyticSample %>%
            select(cid, one_of(obsL2names), z, PS, Logit) %>%
            unique()
          
          ## Running matching algorithm - 1:1 nearest neighbor w/o replacement with caliper of .2 SD of the logit of the PS
          TheMatches <- matchit(as.formula(ThePSModel), data = ClusterData,
                                method = "nearest", replace = FALSE, caliper = .2,
                                distance = ClusterData$Logit)
          
          ## Extracting weight
          ClusterData$weight <- TheMatches$weights
          
          ## Adding weight to subject-level data
          AnalyticSample <- inner_join(AnalyticSample, ClusterData %>% select(cid,weight), by = "cid")
          
        } else if(psmod == "Subject"){
          
          ## Running matching algorithm - 1:1 nearest neighbor w/o replacement with caliper of .2 SD of the logit of the PS
          TheMatches <- matchit(as.formula(ThePSModel), data = AnalyticSample,
                                method = "nearest", replace = FALSE, caliper = .2,
                                distance = AnalyticSample$Logit)
          
          ## Extracting weight
          AnalyticSample$weight <- TheMatches$weights
          
          
        } else {stop("something went wrong")}
        
      } else if(mw == "Weighting"){
        
        ## I only need the weights, not the entire matchit or weightit object
        # TheWeights <- weightit(as.formula(TheModel), data = SampData,
        #                        method = "ps", estimand = "ATE", ps = SampData$PS)
        AnalyticSample$weight <- get_w_from_ps(ps = AnalyticSample$PS, estimand = "ATE",
                                               treat = AnalyticSample$z, treated = 1)
        
      } else {stop("something went wrong")}
      
      
      #### Assessing Balance ####
      ## Balance is only evaluated at the subject-level in accordance with WWC guidelines, but includes all covariates
      TheBalance <- bal.tab(formula = as.formula(paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "))),
                            data = AnalyticSample,
                            continuous = "std", binary = "std", s.d.denom = "all",
                            abs = TRUE, un = TRUE, quick = TRUE,
                            disp.means = FALSE, disp.sd = FALSE, disp.v.ratio = TRUE,
                            weights = AnalyticSample$weight, method = tolower(mw))
      
      ## Calculating ICC(1) and (2) of 10 aggregated covariates
      ICC1 <- purrr::map_dbl(paste0("x",1:10), ~ICC::ICCbareF(factor(cid), quo_name(.x), AnalyticSample)) # ICC(1) of aggregated covariates
      ICC2 <- (nsub*ICC1) / (1 + (nsub - 1)*ICC1)
      
      #### Calculating the post-conditioning DVs for each replication ####
      DVReplicationStats[r,9] <- nrow(AnalyticSample[AnalyticSample$weight != 0,])                      # Number of subjects in analytic sample
      DVReplicationStats[r,10] <- unique(AnalyticSample[AnalyticSample$weight != 0,]$cid) %>% length()  # Number of clusters in analytic sample
      DVReplicationStats[r,11] <-  purrr::map2_dbl(.x = trueL2names, .y = obsL2names,
                                                   ~cor(AnalyticSample[, .x], AnalyticSample[, .y])) %>% mean() # correlation b/t true and observed L2 covariates
      DVReplicationStats[r,12] <- ICC::ICCbareF(factor(cid), Yij, AnalyticSample)  # ICC(1) of outcome Yij
      DVReplicationStats[r,13] <- mean(ICC1)                                       # mean ICC(1) of L1 X covariates
      DVReplicationStats[r,14] <- mean(ICC2)                                       # mean ICC(2) of L2 X covariates
      DVReplicationStats[r,15] <- mean(TheBalance[[1]]$Diff.Un)                    # mean Absolute Standardized Difference (ASD) before conditioning (Unadjusted)
      DVReplicationStats[r,16] <- mean(TheBalance[[1]]$V.Ratio.Un)                 # mean Variance Ration (VR) before conditioning (Unadjusted)
      DVReplicationStats[r,17] <- sum(TheBalance[[1]]$Diff.Un < .10)               # Unadjusted count of covariates (out of 30) that were balanced at L1 based on ASD
      DVReplicationStats[r,18] <- sum(TheBalance[[1]]$V.Ratio.Un < 2)              # Unadjusted count of covariates (out of 30) that were balanced at L1 based on VR
      DVReplicationStats[r,19] <- mean(TheBalance[[1]]$Diff.Adj)                   # mean ASD after conditioning (Adjusted)
      DVReplicationStats[r,20] <- mean(TheBalance[[1]]$V.Ratio.Adj)                # mean VR after conditioning (Adjusted)
      DVReplicationStats[r,21] <- sum(TheBalance[[1]]$Diff.Adj < .10)              # Adjusted count of covariates (out of 30) that were balanced at L1 based on ASD
      DVReplicationStats[r,22] <- sum(TheBalance[[1]]$V.Ratio.Adj < 2)             # Adjusted count of covariates (out of 30) that were balanced at L1 based on VR
      
      #### Estimating Treatment Effect ####
      TheOutcomeModel <- paste0("Yij ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "), " + z + (1|cid)")
      Out.mod <- lmer(as.formula(TheOutcomeModel), data = AnalyticSample[AnalyticSample$weight != 0, ], # removes observations where weight = 0 (these were the unmatched observations)
                      weights = AnalyticSample[AnalyticSample$weight != 0, ]$weight) 
      
      ## Did the Outcome model converge?
      OutModConverge <- ifelse(is.null(Out.mod@optinfo$conv$lme4$code), 1, 0) # Need to double check this
      DVReplicationStats[r,23] <- OutModConverge  # Outcome model convergence
      
      if(OutModConverge == 0 | !("z" %in% attr(fixef(Out.mod), "names"))){ 
        
        ## Don't calculate TE and make remaining DV values NA
        DVReplicationStats[r,24:26] <- NA
        
      } else {
        
        ## Difference between treatment effect estimate and true delta (i.e. bias)
        TEDiff <- fixef(Out.mod)[["z"]] - delta
        
        #### Treatment Effect DVs for each replication ####
        DVReplicationStats[r,24] <- TEDiff          # bias of TE
        DVReplicationStats[r,25] <- abs(TEDiff)     # mae of TE
        DVReplicationStats[r,26] <- TEDiff^2        # rmse of TE
      }
      
    }
  }
  
  ## logging time to run condition
  toc(quiet = TRUE, log = TRUE)
  
  #### DVs averaged across replications ####
  colnames(DVReplicationStats) <- DVnames
  DVConditionStats <- data.frame(DVReplicationStats) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(TE_RMSE = sqrt(TE_RMSE)) %>%
    bind_cols(DataGenConds[con, ])
  
  DVSummaryStats[con, ] <- as.matrix(DVConditionStats[1, ])
  
  if(con == 1){
    
    Con1SampDat <- SampData
    Con1AnalyticSample <- AnalyticSample
    Con1DVReplicationStats <- DVReplicationStats
  }
  
}

#################################
#### Saving Initial Results  ####

colnames(DVSummaryStats) <- colnames(DVConditionStats)
DVsbyCond <- data.frame(DVSummaryStats, stringsAsFactors = FALSE) %>%
  filter(nsub != -999) %>%                                     # Quality control check; should have all 432 rows after filtering, indicating all conditions were completed
  mutate_at(vars(one_of(DVnames)), as.numeric)                # Converting from character to numeric

#### Extracting timing ####
# Character vector
Sim2TimingLog <- tic.log(format = TRUE) %>% unlist()

# Converted to dataframe and cleaned for analysis
Sim2TLogDF <- data.frame(temp = Sim2TimingLog, stringsAsFactors = FALSE) %>%
  separate(temp, c("Row", "Time"), sep = ": ") %>%
  mutate(Time = as.numeric(gsub(" sec elapsed", "", Time, fixed = TRUE))) %>%
  bind_cols(DataGenConds)

# Total time (in min) to run the simulation
sum(Sim2TLogDF$Time) / 60    # for nrep = 1, 8.3 min so for 1000 I project 138.5 hrs or almost 6 days
                             # nrep = 10 took 72.2 min which projects to 120.3 hrs or 5 days

# Clearing time log
# tic.clearlog()

#### Saving Simulation Results ####

## Intermediate information from condition 432
Con432SampDat <- SampData   # saved sample data from last rep of condition 432 
Con432AnalyticSample <- AnalyticSample   # saved summary statistics for each rep of condition 432
Con432DVReplicationStats <- data.frame(DVReplicationStats, stringsAsFactors = FALSE)

save(Con432SampDat, Con432AnalyticSample, Con432DVReplicationStats,
     Con1SampDat, Con1AnalyticSample, Con1DVReplicationStats,
     DVsbyCond, Sim2TLogDF, DVnames,
     file = "Sim2_TenRep.RData")

##################################################

# load("Sim2_con1and432sample.RData")
# load("Sim2_OneRep.RData")
