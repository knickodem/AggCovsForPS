##################################################################################
#                                                                                #
#   Simulation 2 - Cluster-Level Treatment Exposure and Subject-Level Outcome    #
#                                                                                #
##################################################################################

## Loading Packages
library(simstudy)
library(tictoc)
library(dplyr)
library(tidyr)
library(MatchIt)
library(cobalt)
library(lme4)


##########################################
#### Defining Conditions and Outcomes ####

## Study Conditions
DataGenConds <- crossing(nsub = c(20, 60, 100),            # number of subjects per cluster
                         nclus = c(60, 100, 140),           # number of clusters
                         icc = c(.05, .1, .2),             # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                         aggvar = c(.1, .3, .5, 1),        # error/1-% sampled in observed L2 covariates
                         psmod = c("Cluster","Subject"),   # PS estimation and conditioning at cluster or subject level
                         mw = c("Matching", "Weighting"))  # PS conditioning method

nreps <- 1000  # number of replications

## Simulation dependent variables
DVnames <- c("Prop_Treated", "PS_Converged", "Baseline_True.Agg_L2_Cor",               # Proportion in treatment group and Yes/No for PS model convergence
             "PS_Bias", "PS_MAE", "PS_RMSE", "Logit_Bias", "Logit_MAE", "Logit_RMSE",  # Bias and Error in estimation of PS and logit of PS
             "Analytic_Subjects", "Analytic_Clusters", "True.Agg_L2_Cor",              # Num subjects & clusters and correlation b/t true and obs L2 in sample used in outcome model
             "Y_ICC1", "X_ICC1", "X_ICC2",                                             # Characteristics of sample used in outcome analysis
             "Un.ASD", "Un.VR", "Un.ASD_Balanced_Count", "Un.VR_Balanced_Count",       # For unadjusted sample - mean Absolute standardized difference & Variance ratio; Count of the 30 covariates balanced (on ASD or VR) at subject level
             "Adj.ASD", "Adj.VR", "Adj.ASD_Balanced_Count", "Adj.VR_Balanced_Count",   # For adjusted sample - mean ASD & VR; Count of the 30 covariates balanced (on ASD or VR) at subject level
             "Y_Converged", "TE_Bias", "TE_MAE", "TE_RMSE",                            # Yes/No for Outcome model convergence; If yes, Bias, MAE and RMSE for estimation of the treatment effect
             "Baseline_Converged", "Baseline_Bias", "Baseline_MAE", "Baseline_RMSE")   # TE estimation from full (unadjusted) sample  


## Creating blank matrix for summary statistics for each data generation condition
DVSummaryStats <- matrix(-999, ncol = length(DVnames), nrow = nrow(DataGenConds)) # Summary stats for generated datasets


## Variables consistent across all conditions
delta <- .50                                        # Average Treatment Effect ($\delta$)
marg <- -1.386294                                   # marginal probability of treatment;intercept of -1.0986 produces a marginal probability of .25; stats::plogis(-1.0986); -1.386294 = .2
L1names <- paste0("x", 1:20)                         # names of Level 1 covariates
trueL2names <- paste(L1names[1:10], "c", sep = "_") # names of true level 2 covariates
obsL2names <- paste(trueL2names, "o", sep = "_")    # names of aggregated level 2 covariates

############################


#############################
#### Running Simulation ####

set.seed(1213611) # Seed for reproducing results

tic("Total")
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
      ungroup() %>%
      mutate(TrueLogit = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +                              # True logit of the probability of treatment exposure (i.e. propensity score)
               .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + zuj,
             TruePS = stats::plogis(TrueLogit),                                                                  # true PS
             z = ifelse(TrueLogit > 0, 1, 0),                                                                  # Determining treatment exposure
             Yij = .5*x1 + .5*x2 + .5*x3 + .5*x4 + .5*x5 + .5*x6 + .5*x7 + .5*x8 + .5*x9 + .5*x10 +            # Generating outcome values
               .5*x11 + .5*x12 + .5*x13 + .5*x14 + .5*x15 + .5*x16 + .5*x17 + .5*x18 + .5*x19 + .5*x20 + rij + # Level 1 covariates
               .5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c + 
               .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + delta*z + uj,                                # Level 2 covariates
             cid = factor(cid))                                                                                # converting cluster id from character to factor
    
    #### Estimating PS ####
    
    if(psmod == "Cluster"){ 
      
      ## Model appraises treatment at the cluster-level and only uses L2 covariates
      ThePSModel <- paste0("z ~ 1 + ", paste(obsL2names, collapse = " + "))
      PS.mod <- glm(as.formula(ThePSModel), data = GenData, family = binomial("logit"))
      
    } else if(psmod == "Subject"){
      
      ## Model ignores clusters and uses L1 and L2 covariates, the latter conceptually treated as L1 covariates
      ThePSModel <- paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "))
      PS.mod <- glm(as.formula(ThePSModel), data = GenData, family = binomial("logit"))
      
    } else {stop("something went wrong")}
    
    ## Did the PS model converge?
    PSModConverge <- ifelse(PS.mod$converged == TRUE, 1, 0)
    
    ## Adding estimated PS to dataframe
    GenData$PS <- fitted(PS.mod)
    GenData$Logit <- predict(PS.mod)                       # observed logit of the PS
    GenData$LogitDiff <- GenData$Logit - GenData$TrueLogit # Error in observed (estimated) and true logit of the PS
    GenData$PSDiff <- GenData$PS - GenData$TruePS          # Error in observed and true PS
    
    ## Saving DVs
    DVReplicationStats[r,1] <- mean(GenData$z)            # proportion with treatment exposure
    DVReplicationStats[r,2] <- PSModConverge              # PS model convergence
    
    #### Checking Overlap assumption ####
    ## Removing cases where PS is equal to 1 or 0
    AnalyticSample <- GenData %>%
      filter(PS > .001 & PS < .999)
    
    if(PSModConverge == 0 | nrow(AnalyticSample) == 0){
      
      ## Don't run conditioning and outcome model. Make remaining DV values NA
      DVReplicationStats[r,3:31] <- NA
      
    } else {
      
      #### Calculating the pre-conditioning DVs that utilize the full sample ####
      DVReplicationStats[r,3] <-  purrr::map2_dbl(.x = trueL2names, .y = obsL2names,
                                                   ~cor(GenData[,.x],GenData[,.y])) %>%        # correlation b/t true and observed L2 covariates in full sample
        psych::fisherz() %>% ifelse(. > 7.3, 7.3, .) %>% ifelse(. < -7.3, -7.3, .) %>% mean()  # z = +- 7.3 is equivalent to r = .9999991; when r = 1, z = Inf which throws of calculation of mean correlation
      DVReplicationStats[r,4] <- mean(GenData$PSDiff)              # bias of PS
      DVReplicationStats[r,5] <- mean(abs(GenData$PSDiff))         # mae of PS
      DVReplicationStats[r,6] <- sqrt(mean(GenData$PSDiff^2))      # rmse of PS
      DVReplicationStats[r,7] <- mean(GenData$LogitDiff)           # bias of logit
      DVReplicationStats[r,8] <- mean(abs(GenData$LogitDiff))      # mae of logit
      DVReplicationStats[r,9] <- sqrt(mean(GenData$LogitDiff^2))   # rmse of logit
      
      #################################
      
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
          tryCatch(expr = {TheMatches <- matchit(as.formula(ThePSModel), data = ClusterData,
                                                 method = "nearest", replace = FALSE, caliper = .2,
                                                 distance = ClusterData$Logit)},
                   error = function(e){TheMatches <<- TRUE; return(TheMatches)})
          
          ## Extracting the weights or assigning weights to 0 if no matches
          if(class(TheMatches) == "matchit"){
            
            ClusterData$weight <- TheMatches$weights
            
          } else {
            
            ClusterData$weight <- 0
          }
          
          ## Adding weight to subject-level data
          AnalyticSample <- inner_join(AnalyticSample, ClusterData %>% select(cid,weight), by = "cid")
          
        } else if(psmod == "Subject"){
          
          ## Running matching algorithm - 1:1 nearest neighbor w/o replacement with caliper of .2 SD of the logit of the PS
          # Catch error if one should occur
          tryCatch(expr = {TheMatches <- matchit(as.formula(ThePSModel), data = AnalyticSample,
                                                 method = "nearest", replace = FALSE, caliper = .2,
                                                 distance = AnalyticSample$Logit)},
                   error = function(e){TheMatches <<- TRUE; return(TheMatches)})
          
          ## Extracting the weights or assigning weights to 0 if no matches
          if(class(TheMatches) == "matchit"){
            
            AnalyticSample$weight <- TheMatches$weights
            
          } else {
            
            AnalyticSample$weight <- 0
          }
          
        } else {stop("something went wrong")}
        
      } else if(mw == "Weighting"){
        
        AnalyticSample$weight <- WeightIt::get_w_from_ps(ps = AnalyticSample$PS, estimand = "ATE",
                                                         treat = AnalyticSample$z, treated = 1)
        
      } else {stop("something went wrong")}
      
      ## If no matches occured
      if(sum(AnalyticSample$weight) == 0){
        
        ## Don't run conditioning and outcome model. Make remaining DV values NA
        DVReplicationStats[r,2] <- 0      # Change PS convergence code to 0
        DVReplicationStats[r,3:31] <- NA  # make remaining DV values NA
        
      } else {
        ## With weighting or when matching was successful
        
        #### Assessing Balance ####
        ## Balance is only evaluated at the subject-level in accordance with WWC guidelines, but includes all covariates
        TheBalance <- bal.tab(formula = as.formula(paste0("z ~ 1 + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "))),
                              data = AnalyticSample,
                              continuous = "std", binary = "std", s.d.denom = "pooled",
                              abs = TRUE, un = TRUE, quick = TRUE,
                              disp.means = FALSE, disp.sd = FALSE, disp.v.ratio = TRUE,
                              weights = AnalyticSample$weight, method = tolower(mw))
        
        ## Calculating ICC(1) and (2) of 10 aggregated covariates
        ICC1 <- purrr::map_dbl(paste0("x",1:10), ~ICC::ICCbare(factor(cid), quo_name(.x), AnalyticSample[AnalyticSample$weight != 0,])) %>% # ICC(1) of aggregated covariates
          ifelse(. < 0, 0, .)                       # negatives constrained to 0
        ICC2 <- (nsub*ICC1) / (1 + (nsub - 1)*ICC1)
        
        #### Calculating the post-conditioning DVs for each replication ####
        DVReplicationStats[r,10] <- nrow(AnalyticSample[AnalyticSample$weight != 0,])                      # Number of subjects in analytic sample
        DVReplicationStats[r,11] <- unique(AnalyticSample[AnalyticSample$weight != 0,]$cid) %>% length()  # Number of clusters in analytic sample
        DVReplicationStats[r,12] <-  purrr::map2_dbl(.x = trueL2names, .y = obsL2names,
                                                     ~cor(AnalyticSample[AnalyticSample$weight != 0, .x],
                                                          AnalyticSample[AnalyticSample$weight != 0, .y])) %>% # correlation b/t true and observed L2 covariates
          psych::fisherz() %>% ifelse(. > 7.3, 7.3, .) %>% ifelse(. < -7.3, -7.3, .) %>% mean()  # z = +- 7.3 is equivalent to r = .9999991; when r = 1, z = Inf which throws of calculation of mean correlation
        DVReplicationStats[r,13] <- ICC::ICCbare(factor(cid), Yij, AnalyticSample[AnalyticSample$weight != 0,]) %>%
          ifelse(. < 0, 0, .)                                                        # ICC(1) of outcome Yij (negatives constrained to 0)
        DVReplicationStats[r,14] <- mean(ICC1)                                       # mean ICC(1) of L1 X covariates
        DVReplicationStats[r,15] <- mean(ICC2)                                       # mean ICC(2) of L2 X covariates
        DVReplicationStats[r,16] <- mean(TheBalance[[1]]$Diff.Un)                    # mean Absolute Standardized Difference (ASD) before conditioning (Unadjusted)
        DVReplicationStats[r,17] <- mean(TheBalance[[1]]$V.Ratio.Un)                 # mean Variance Ration (VR) before conditioning (Unadjusted)
        DVReplicationStats[r,18] <- sum(TheBalance[[1]]$Diff.Un < .10)               # Unadjusted count of covariates (out of 30) that were balanced at L1 based on ASD
        DVReplicationStats[r,19] <- sum(TheBalance[[1]]$V.Ratio.Un < 2)              # Unadjusted count of covariates (out of 30) that were balanced at L1 based on VR
        DVReplicationStats[r,20] <- mean(TheBalance[[1]]$Diff.Adj)                   # mean ASD after conditioning (Adjusted)
        DVReplicationStats[r,21] <- mean(TheBalance[[1]]$V.Ratio.Adj)                # mean VR after conditioning (Adjusted)
        DVReplicationStats[r,22] <- sum(TheBalance[[1]]$Diff.Adj < .10)              # Adjusted count of covariates (out of 30) that were balanced at L1 based on ASD
        DVReplicationStats[r,23] <- sum(TheBalance[[1]]$V.Ratio.Adj < 2)             # Adjusted count of covariates (out of 30) that were balanced at L1 based on VR
        
        #### Estimating Treatment Effect ####
        ## The Outcome Model
        TheOutcomeModel <- paste0("Yij ~ 1 + z + ", paste(L1names, collapse = " + "), " + ", paste(obsL2names, collapse = " + "), " + (1|cid)")
        # Running model and catching errors if necessary
        tryCatch(expr = {Out.mod <- lmer(as.formula(TheOutcomeModel), data = AnalyticSample[AnalyticSample$weight != 0, ], # removes observations where weight = 0 (these were the unmatched observations)
                        weights = AnalyticSample[AnalyticSample$weight != 0, ]$weight)},
                 error = function(e){Out.mod <<- TRUE; return(Out.mod)})
        
        if(class(Out.mod) == "lmerMod"){
          ## Did the Outcome model converge?
          OutModConverge <- ifelse(is.null(Out.mod@optinfo$conv$lme4$code), 1, 0)
        } else{
          ## Did the Outcome model converge?
          OutModConverge <- 0
        }
        
        ## Recording Outcome model convergence
        DVReplicationStats[r,24] <- OutModConverge
        
        if(is.logical(Out.mod)){
          
          ## Don't calculate TE and make remaining DV values NA
          DVReplicationStats[r,25:31] <- NA
          
        } else if(OutModConverge == 0 | !("z" %in% attr(fixef(Out.mod), "names"))){ 
          
          ## Don't calculate TE and make remaining DV values NA
          DVReplicationStats[r,25:31] <- NA
          
        } else {
          
          ## Difference between treatment effect estimate and true delta (i.e. bias)
          TEDiff <- fixef(Out.mod)[["z"]] - delta
          
          #### Treatment Effect DVs for each replication from outcome model ####
          DVReplicationStats[r,25] <- TEDiff          # bias of TE
          DVReplicationStats[r,26] <- abs(TEDiff)     # mae of TE
          DVReplicationStats[r,27] <- TEDiff^2        # rmse of TE
          
          #### Estimating Treatment Effect with full sample ####
          ## Establishes baseline to which the PS methods can be compared
          # (i.e, is the PS even worth it)
          Base.mod <- lmer(as.formula(TheOutcomeModel), data = GenData)
          
          ## Did the Baseline model converge?
          BaseModConverge <- ifelse(is.null(Base.mod@optinfo$conv$lme4$code), 1, 0)
          DVReplicationStats[r,28] <- BaseModConverge  # Baseline model convergence
          
          if(BaseModConverge == 0 | !("z" %in% attr(fixef(Base.mod), "names"))){ 
            
            ## Don't calculate TE
            DVReplicationStats[r,29:31] <- NA
            
          } else {
            
            ## Difference between treatment effect estimate from baseline model and the true delta (i.e. bias)
            BaseTEDiff <- fixef(Base.mod)[["z"]] - delta
            
            #### Treatment Effect DVs for each replication ####
            DVReplicationStats[r,29] <- BaseTEDiff          # bias of TE from baseline model
            DVReplicationStats[r,30] <- abs(BaseTEDiff)     # mae of TE from baseline model
            DVReplicationStats[r,31] <- BaseTEDiff^2        # rmse of TE from baseline model
            
          }  # ends BaseModConverge evaluation
          
        } # ends OutModConverge evaluation
        
      } # ends TheMatches evaluation (i.e., when no matches, skip remaining analysis)
    }  # ends PSModConverge evaluation
  }  # ends replication
  
  ## logging time to run condition
  toc(quiet = TRUE, log = TRUE)
  
  #### DVs averaged across replications ####
  DVSummaryStats[con, ] <- colMeans(DVReplicationStats, na.rm = TRUE)
  
  if(con == 1){
    
    Con1GenDatP <- GenData                     # saves dataset from the 1000th rep  
    Con1AnalyticSampleP <- AnalyticSample      # saves dataset from the 1000th rep
    Con1DVReplicationStatsP <- DVReplicationStats
  }
  
}
toc()

#################################
#### Saving Initial Results  ####

save(DVSummaryStats, file = "Sim2_ThousandRep_DVSummaryOnly_Pooled.RData")

colnames(DVSummaryStats) <- DVnames
Sim2DVsbyCondPooled <- data.frame(DVSummaryStats, stringsAsFactors = FALSE) %>%
  mutate_all(as.numeric) %>%                # Converting from character to numeric
  mutate(Logit_RMSE = sqrt(Logit_RMSE),
         TE_RMSE = sqrt(TE_RMSE),
         Baseline_RMSE = sqrt(Baseline_RMSE),
         True.Agg_L2_Cor = psych::fisherz2r(True.Agg_L2_Cor),                       # Converting z-scores to correlations
         Baseline_True.Agg_L2_Cor = psych::fisherz2r(Baseline_True.Agg_L2_Cor)) %>% # Converting z-scores to correlations
  tibble::rownames_to_column("Con") %>%
  left_join(DataGenConds %>% tibble::rownames_to_column("Con"), by = "Con") %>%
  select(-Con)

#### Extracting timing ####
# Character vector
Sim2TimingLog <- tic.log(format = TRUE) %>% unlist()

# Converted to dataframe and cleaned for analysis
Sim2TLogDFPooled <- data.frame(temp = Sim2TimingLog, stringsAsFactors = FALSE) %>%
  separate(temp, c("Row", "Time"), sep = ": ") %>%
  mutate(Time = as.numeric(gsub(" sec elapsed", "", Time, fixed = TRUE)))

# Total time (in min) to run the simulation
sum(Sim2TLogDFPooled$Time) / 60

# Clearing time log
tic.clearlog()

#### Saving Simulation Results ####

## Intermediate information from condition 432
Con432GenDatP <- GenData   # saved sample data from last rep of condition 432 
Con432AnalyticSampleP <- AnalyticSample   # saved summary statistics for each rep of condition 432
Con432DVReplicationStatsP <- data.frame(DVReplicationStats, stringsAsFactors = FALSE)


save(Con432GenDatP, Con432AnalyticSampleP, Con432DVReplicationStatsP,
     Con1GenDatP, Con1AnalyticSampleP, Con1DVReplicationStatsP,
     DVSummaryStats, Sim2DVsbyCondPooled, Sim2TLogDFPooled, DVnames,
     file = "Sim2_ThousandRep_Pooled.RData")

##################################################
