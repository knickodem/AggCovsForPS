#### Examining different ways of generating true PS ####

# aam = Arpino and Mealli (2011), also Leite et al (2015)
# aus = Austin's favored approach
# irt = method commonly used in irt simulations for generating item responses

library(simstudy)
library(tictoc)
library(dplyr)
library(tidyr)

#### Converting logits to probabilities and vice versa ####
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
                         icc = c(.05, .1, .2),            # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                         aggvar = c(.1,.3,.5,1),          # error/1-% sampled in observed L2 covariates
                         method = c("aam", "aus", "irt"),
                         model = c("main","interact"))
# method = c("Formative.RE","Formative.Sample","Reflective.RE","Reflective.Sample"))   # data generation method

nreps <- 100  # number of replications

############################

## Creating blank matrix for summary statistics for each data generation condition
DGSnames <- c("Prop_Z1","True_Logit","True_PS","Obs_Logit","Obs_PS",
              "Logit_Bias","Logit_MAE","Logit_RMSE","PS_Bias","PS_MAE","PS_RMSE",
              "L1_Variance","L1_Cor",
              "True_L2_Variance","True_L2_Cor","Obs_L2_Variance","Obs_L2_Cor",
              "True.Obs_L2_Cor", "X_ICC1","X_ICC2","Converged")
DGS <- matrix(-999, ncol = length(DGSnames) + ncol(DataGenConds), nrow = nrow(DataGenConds)) # Summary stats for generated datasets
SampDatabyCond <- list()  # save sample data from first rep of each condition

## Variables consistent across all conditions
delta <- .5 #c(.2, .5, .8),      # treatment effect size
marg <- -1.0986                  # marginal probability of treatment;intercept of -1.0986 produces a marginal probability of .25; LogitToProb(-1.0986)
L1names <- paste0("x",1:20)
trueL2names <- paste(L1names[1:10], "c", sep = "_")
obsL2names <- paste(trueL2names, "o", sep = "_")

#### Starting Generation ####
set.seed(48226) # Seed for reproducing results

for(con in 1:nrow(DataGenConds)){
  
  tic(as.character(con))          # Record time to run all replications for each condition
  
  nsub <- DataGenConds[con,][[1]]  # Number of Level 1 units (i.e., subjects)
  nclus <- DataGenConds[con,][[2]]  # Number of Level 2 units (i.e., clusters)
  ntot <- nclus*nsub                # Total sample size
  icc <- DataGenConds[con,][[3]]   # Intraclass correlation coefficient
  aggvar <- DataGenConds[con,][[4]]  # error for creating observed L2 covariates
  method <- DataGenConds[con,][[5]] # method for generating L1, true L2, and obs L2 covariates
  model <- DataGenConds[con,][[6]]
  tau00 <- icc                     # between cluster (L2) variance ($\tau_{00}$)
  sigma2 <- 1 - tau00              # within cluster (L1) variance ($\sigma^2$)
  
  
  ## Stats for each data generation replication
  DataGenRepStats <- matrix(-999,ncol = length(DGSnames), nrow = nreps)
  
  ##################################################
  #### Testing various ways of creating true PS ####
  
  for(r in 1:nreps){
    
    ## Level 1 covariates - 20 continuous covariates from a standard normal distribution with correlation of .2
    GenData <- genCorData(n = ntot, mu = rep(0, 20), sigma = sqrt(sigma2), rho = 0.2, corstr = "cs", 
                          cnames = L1names, idname = "sid") %>%
      mutate(rij = rnorm(n = ntot, mean = 0, sd = sqrt(sigma2)),                   # random error for outcome model from normal distribution (i.e, $\mu$ = 0, $\sigma^2$ = 1)
             cid = as.character(rep(1:nclus, each = nsub))) %>%                    # creating cluster ID
      group_by(cid) %>%
      mutate_at(vars(x1:x20), ~. + rnorm(n = 1, mean = 0, sd = sqrt(tau00))) %>%    # Adding cluster variation to maintain ICC(1)
      mutate_at(vars(x1:x10), list(c = ~mean(.))) %>%                               # aggregating L1 to true L2 value
      mutate_at(vars(x1_c:x10_c), list(o = ~. + rnorm(n = 1, mean = 0, sd = aggvar))) %>% # Adding random variation to true cluster mean; selects a random number for each cluster for each covariate
      ungroup()
    
    if(method == "aam"){
      
      #### Calculating Y, treatment assignment, true and observed log PS, and ICCs ####
      SampData <- GenData %>%
        group_by(cid) %>%
        mutate(zuj = rlogis(n = 1, location = 0, scale = 1)) %>%                    # random error for true PS model (probability of treatment) from logistic distribution with mean = 0  and variance of $\pi^2 / 3$
        ungroup()
        
        if(model == "main"){
          
          SampData <- SampData %>%
            mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                     .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c + zuj)
          
        } else if(model == "interact"){
          
          SampData <- SampData %>%
            mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                     .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c +
                     .5*(x1_c*x2_c) + .5*(x1_c*x3_c) + .5*(x1_c*x4_c) + .5*(x2_c*x3_c) +
                     .5*(x3_c*x4_c) + zuj)
          
        } else {stop("wrong")}
      
      SampData <- SampData %>%
        mutate(TruePS = LogitToProb(logitz1),                              # true PS
               z = ifelse(logitz1 > 0, 1, 0))                                       # Determining treatment exposure
      
    } else if(method == "aus"){
      
      if(model == "main"){
        
        SampData <- GenData %>%
          mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                   .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c)
        
      } else if(model == "interact"){
        
        SampData <- GenData %>%
          mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                   .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c +
                   .5*(x1_c*x2_c) + .5*(x1_c*x3_c) + .5*(x1_c*x4_c) + .5*(x2_c*x3_c) +
                   .5*(x3_c*x4_c))
        
      } else {stop("wrong")}
      
      SampData <- SampData %>%
        mutate(TruePS = LogitToProb(logitz1)) %>%                           # true PS
        group_by(cid) %>%
        mutate(z = rbinom(n = 1, size = 1, prob = LogitToProb(logitz1))) %>%         # treatment exposure selected from Bernoulli distribution where probability = true PS
        ungroup()
      
    } else if(method == "irt"){
      
      SampData <- GenData %>%
        group_by(cid) %>%
        mutate(compz = runif(n = 1, min = 0, max = 1)) %>%                       # Random value from a uniform distribution for each cluster
        ungroup()
      
      if(model == "main"){
        
        SampData <- SampData %>%
          mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                   .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c)
        
      } else if(model == "interact"){
        
        SampData <- SampData %>%
          mutate(logitz1 = marg + 5*x1_c + .5*x2_c + .5*x3_c + .5*x4_c + .5*x5_c +    # True logit of the probability of treatment exposure (i.e. propensity score)
                   .5*x6_c + .5*x7_c + .5*x8_c + .5*x9_c + .5*x10_c +
                   .5*(x1_c*x2_c) + .5*(x1_c*x3_c) + .5*(x1_c*x4_c) + .5*(x2_c*x3_c) +
                   .5*(x3_c*x4_c))
        
      } else {stop("wrong")}
      
      SampData <- SampData %>%
        mutate(TruePS = LogitToProb(logitz1),                           # true PS
               z = ifelse(logitz1 > compz, 1, 0))
      
      
    } else {stop("something went wrong")}
    
    ## Estimating observed PS
    # need to investigate why some models do not converge, which then affects the bias and RMSE calculations; Not sure if this is actually a problem         
    PS.mod <- glm(as.formula(paste0("z ~ 1 + ",paste(obsL2names,collapse = " + "))), data = SampData, family = binomial("logit"))
    
    ## Adding estimated PS to dataframe
    SampData$PS <- fitted(PS.mod)                           # observed PS 
    SampData$Logit <- predict(PS.mod)                       # observed logit of the PS
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
    DataGenRepStats[r,8] <- sqrt(mean(SampData$LogitDiff^2))                    # rmsea of logit
    DataGenRepStats[r,9] <- mean(SampData$PSDiff)                               # bias of PS
    DataGenRepStats[r,10] <- mean(abs(SampData$PSDiff))                         # mae of PS
    DataGenRepStats[r,11] <- sqrt(mean(SampData$PSDiff^2))                      # rmsea of PS
    DataGenRepStats[r,12] <- cov(SampData[,L1names]) %>% diag() %>% mean()              # mean variance of L1 covariates
    DataGenRepStats[r,13] <- cor(SampData[,L1names]) %>% .[lower.tri(.)] %>% mean()    # mean correlation of L1 covariates
    DataGenRepStats[r,14] <- cov(SampData[,trueL2names]) %>% diag() %>% mean()         # mean variance of true L2 covariates
    DataGenRepStats[r,15] <- cor(SampData[,trueL2names]) %>% .[lower.tri(.)] %>% mean() # mean correlation of true L2 covariates
    DataGenRepStats[r,16] <- cov(SampData[,obsL2names]) %>% diag() %>% mean()           # mean variance of observed L2 covariates
    DataGenRepStats[r,17] <- cor(SampData[,obsL2names]) %>% .[lower.tri(.)] %>% mean()  # mean correlation of observed L2 covariates
    DataGenRepStats[r,18] <- purrr::map2_dbl(.x = trueL2names,.y = obsL2names,~cor(SampData[,.x],SampData[,.y])) %>% mean() # correlation b/t true and observed L2 covariates
    DataGenRepStats[r,19] <- mean(ICC1)                                   # mean ICC(1) of L1 X covariates
    DataGenRepStats[r,20] <- mean(ICC2)                                   # mean ICC(2) of L2 X covariates
    DataGenRepStats[r,21] <- ifelse(PS.mod$converged == TRUE, 1, 0)       # Did the PS model converge?
    
    
    if(r == 1){
      SampDatabyCond[[con]] <- SampData # saving sample data from first repetition from each condition
    }
  }
  
  ## logging time to run replication
  toc(quiet = TRUE, log = TRUE)
  
  #### Sample characteristics averaged across replications ####
  colnames(DataGenRepStats) <- DGSnames
  DataGenStats <- data.frame(DataGenRepStats) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    bind_cols(DataGenConds[con, ])
  
  DGS[con, ] <- as.matrix(DataGenStats[1, ])
  
  if(con == 1){
    Con1DataGenRepStats <- DataGenRepStats
  }
  
}

#################################
#### Saving Initial Results  ####

colnames(DGS) <- colnames(DataGenStats)
DGSbyCond <- data.frame(DGS, stringsAsFactors = FALSE) %>%
  filter(nsub != -999) %>%                                     # Quality control check; should have all 432 rows after filtering, indicating all conditions were completed
  mutate_at(vars(one_of(DGSnames)), as.numeric)                # Converting from character to numeric
# DGSbyCond %>% mutate_if(is.numeric,~round(.,2)) %>% View()


#### Extracting timing ####
# Character vector
DataGenTimingLog <- tic.log(format = TRUE) %>% unlist()

# Converted to dataframe and cleaned for analysis
DataGenTLogDF <- data.frame(temp = DataGenTimingLog, stringsAsFactors = FALSE) %>%
  separate(temp, c("Row", "Time"), sep = ": ") %>%
  mutate(Time = as.numeric(gsub(" sec elapsed", "", Time, fixed = TRUE)))

DataGenTLogDF <- DataGenTLogDF[649:1296, ] %>%
  bind_cols(DataGenConds)

# Total time (in min) to run the simulation
sum(DataGenTLogDF$Time) / 60

# Clearing time log
tic.clearlog()

## Sample data from each condition into a data frame
AllSampData <- purrr::map_dfr(c(1:nrow(DGSbyCond)), ~bind_cols(SampDatabyCond[[.x]], as.data.frame(lapply(DataGenConds[.x, ], rep, nrow(SampDatabyCond[[.x]])))))

#### Saving Objects ####
save(Con1DataGenRepStats,
     DGSbyCond, DataGenTLogDF, DGSnames,
     SampData,
     file = "GenerateTruePS_Results.RData")



# ## testing method for replicating cluster information
# microbenchmark::microbenchmark(
#   DataGenConds[rep(seq_len(nrow(DataGenConds)), each = 20000), ],  #faster
#   as.data.frame(lapply(DataGenConds, rep, each = 20000)),
#   times = 10
#   )

###########################################################

##################
#### Analysis ####

## Converting independent variables to factors
DGSbyCond <- DGSbyCond %>%
  mutate_at(vars(nsub:model), ~forcats::as_factor(.))

## Order independent variables (and interactions) for use in regressions or long formatted data
ConditionVars <- c("nsub", "nclus", "icc", "aggvar", "method","model",
                   "nsub:nclus", "nsub:icc", "nsub:aggvar", "nsub:method","nsub:model",
                   "nclus:icc", "nclus:aggvar", "nclus:method","nclus:model",
                   "icc:aggvar", "icc:method", "icc:model",
                   "aggvar:method","aggvar:model","method:model")

## True and Observed PS distributions by z (treatment exposure) and generation method
library(ggplot2)
library(purrr)

InteractionPlot <- AllSampData %>%
  mutate(z = factor(z)) %>%
  select(cid, method, model, z, True = TruePS, Obs = PS) %>%
  filter(model == "interact") %>%
  gather(Type, PS, True, Obs) %>%
  ggplot(aes(x = PS, group = z, fill = z)) +
  geom_density(alpha = .5) +
  theme_bw() +
  facet_grid(Type ~ method)

MainOnlyPlot <- AllSampData %>%
  mutate(z = factor(z)) %>%
  select(cid, method, model, z, True = TruePS, Obs = PS) %>%
  filter(model == "main") %>%
  gather(Type, PS, True, Obs) %>%
  ggplot(aes(x = PS, group = z, fill = z)) +
  geom_density(alpha = .5) +
  theme_bw() +
  facet_grid(Type ~ method)


ConvergedPlot <- DGSbyCond %>%
ggplot(aes(x = nclus, y = Converged, fill = method)) +
  geom_boxplot(notch = FALSE) +
  theme_bw(base_size = 20) +
  facet_wrap(~aggvar)
  

table(DGSbyCond$Converged, DGSbyCond$model)
table(DGSbyCond$Converged, DGSbyCond$method)

#### In Sample Characteristics, then binding timing ####
## Regressions
Thelms <- map(DGSnames,~lm(as.formula(paste(.x, "~", paste(ConditionVars, collapse = " + "))),
                           data = DGSbyCond)) %>%
  set_names(DGSnames)
TheOmegas <- map_dfr(Thelms, ~sjstats::omega_sq(.x, partial = TRUE, ci.lvl = .95), .id = "Characteristic")

TheOmegas %>% filter(Characteristic == "Converged") %>%
  mutate_if(is.numeric, ~round(.,3)) %>%
  arrange(desc(partial.omegasq))
