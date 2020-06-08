################################################
#                                              #
#   Simulation 2 - Cluster-level Treatment    #
#           Analysis and Plotting              #
#                                              #
################################################


##################################################################
#### Loading packages, functions, and saved data for analysis ####

library(purrr)
library(stringr)
library(ggplot2)
library(forcats)
library(dplyr)
library(tidyr)

#### Simulation Results ####
load("Saved Sim2 Results/Sim2_ThousandRep_Pooled.RData")

#### Creates QQ, fitted vs residual, and Cook's d ####
ModelCheckPlots <- function(model, smooth_method = "loess", ...){
  
  ## Extracting residuals
  # Dependent variable name
  DV <- names(model$model)[[1]]
  
  # fit information
  aug <- broom::augment(model) #glm. arguments ignored for lm objects
  
  ## Adding id column
  aug <- tibble::rowid_to_column(aug, ".id")
  
  ## Q-Q Plot
  qq <- ggplot(aug, aes(sample = .resid)) +  ## standardizing doesn't change the shape of the distribution so it doesn't matter of .resid or .std.resid
    stat_qq(shape = 1) + stat_qq_line() +
    labs(x = "Theoretical Quantiles", y = "Residual") +
    theme_bw()
  
  ## Fitted vs. Residuals Plot
  frplot <- ggplot(aug, aes(x = .fitted, y = .resid)) +
    geom_hline(aes(yintercept = 0), color = "grey", size = 1, linetype = 1) +
    geom_point(shape = 1) +
    geom_smooth(method = smooth_method, ...) +
    labs(x = "Fitted", y = "Residual") +
    theme_bw()
  
  MCPlots <- cowplot::plot_grid(qq, frplot, nrow = 1, align = "h", labels = "AUTO")
  
  return(MCPlots)
}

## Shortcut for investigating factors with highest omegas
Arrange_Omegas <- function(Outcome){
  
  ot <- enquo(Outcome)
  Sim2Omegas_Wide %>% select(Factor, !!ot) %>% arrange(desc(!!ot))
  
}

#### Shortcut for line and boxplots ####
Plot_Shortcut <- function(dat, xv, yv, gpv, linebox = c("line","box", "linebox"), smooth_method = "loess"){
  
  if(linebox == "line"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv, group = gpv, color = gpv)) +
      geom_point(size = 3) + geom_smooth(method = smooth_method, se = FALSE, size = 2) +
      theme_bw(base_size = 20)
    
  } else if(linebox == "box"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv, fill = gpv)) +
      geom_boxplot(notch = FALSE) +
      theme_bw(base_size = 20)
    
  } else if(linebox == "linebox"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv)) +
      geom_smooth(aes_string(group = gpv, color = gpv), method = "loess", se = FALSE, size = 2) +
      geom_boxplot(aes_string(fill = gpv), notch = FALSE) +
      theme_bw(base_size = 20)
    
  }
  
  return(deplot)
}

## Character string of dependent variables
DVnames <- names(Sim2DVsbyCondPooled)[1:31]

########################################################

########################################################
#### Defining Factors for analysis and presentation ####

Sim2DVsbyCondPooledOrig <- Sim2DVsbyCondPooled

## Converting independent variables to factors
Sim2DVsbyCondPooled <- Sim2DVsbyCondPooledOrig %>%
  drop_na() %>%                                      # 432 observations should remain
  mutate_at(vars(nsub:mw), ~as_factor(.))


## Order independent variables (and interactions) for use in regressions or long formatted data
Sim2ConditionVars <- c("nsub", "nclus", "icc", "aggvar", "psmod", "mw",
                       "nsub:nclus", "nsub:icc", "nsub:aggvar", "nsub:psmod", "nsub:mw",
                       "nclus:icc", "nclus:aggvar", "nclus:psmod", "nclus:mw",
                       "icc:aggvar", "icc:psmod", "icc:mw",
                       "aggvar:psmod", "aggvar:mw", "psmod:mw")

## Converting characters for output
Sim2ConditionFactor <- str_replace(Sim2ConditionVars, "nsub", "Subjects") %>%
  str_replace("nclus","Clusters") %>%
  str_replace("icc","ICC(1)") %>%
  str_replace("aggvar","Error Magnitude") %>%
  str_replace("psmod","PS Model") %>%
  str_replace("mw","PS Conditioning")

## Renamed levels
Sim2ConditionRename <- Sim2ConditionVars
names(Sim2ConditionRename) <- Sim2ConditionFactor

#### Recoding factor levels ####
ICCLevels <- c(`.20` = "0.2", `.10` = "0.1", `.05` = "0.05")
SubjectLevels <- c("20", "60", "100")
ClusterLevels <- c("60", "100", "140")


###########################################################################

###########################################################################
#### Variation in outcomes explained by condition (partial-$\Omega^2$) ####

Sim2TLogDFPooled <- Sim2TLogDFPooled %>%
  bind_cols(crossing(nsub = c(20, 60, 100),            # number of subjects per cluster
                     nclus = c(60, 100, 140),           # number of clusters
                     icc = c(.05, .1, .2),             # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
                     aggvar = c(.1, .3, .5, 1),        # error/1-% sampled in observed L2 covariates
                     psmod = c("Cluster","Subject"),   # PS estimation and conditioning at cluster or subject level
                     mw = c("Matching", "Weighting")))

#### In timing ####
## Descriptives
psych::describeBy(Sim2TLogDFPooled$Time, list(Sim2TLogDFPooled$nclus), mat = TRUE, digits = 3) #Mdn = 4, 6, 8 but max is 19, 47, and 81 seconds
psych::describeBy(Sim2TLogDFPooled$Time, list(Sim2TLogDFPooled$psmod), mat = TRUE, digits = 3) #Cluster: Mdn = 4.7, max = 12; Subject: Mdn = 8.1, max = 84.6
psych::describeBy(Sim2TLogDFPooled$Time, list(Sim2TLogDFPooled$mw), mat = TRUE, digits = 3) #Maching: Mdn = 7.0, max = 84.6; Weight: Mdn = 2.5, max = 13.1

## Regression
Sim2Time <- lm(as.formula(paste0("Time ~ ", paste(Sim2ConditionVars, collapse = " + "))),
               data = Sim2TLogDFPooled)
## Effect Size
Sim2TimeOmega <- sjstats::omega_sq(Sim2Time, partial = TRUE) %>%
  mutate(Characteristic = "Run Time")

#### In Sample Characteristics, then binding timing ####
## Regressions
Sim2lms <- map(DVnames, ~lm(as.formula(paste(.x, "~", paste(Sim2ConditionVars, collapse = " + "))),
                            data = Sim2DVsbyCondPooled)) %>%
  set_names(DVnames)

## Effect size
#Note: Albers & Lakens (2017, p.9) says its okay that there are negative numbers even if it is theoretically impossible
Sim2Omegas <- map_dfr(Sim2lms, ~sjstats::omega_sq(.x, partial = TRUE), .id = "Characteristic") %>%
  bind_rows(Sim2TimeOmega) %>%
  mutate(term = factor(term, levels = Sim2ConditionVars),
         Factor = fct_recode(term, !!!Sim2ConditionRename)) %>%
  mutate_if(is.numeric, ~round(., 2))

## Wide format; for final paper, bold values > .14, which are considered large
# Table is exported at end of code
Sim2Omegas_Wide <- Sim2Omegas %>%
  select(Factor, Characteristic, partial.omegasq) %>%
  spread(Characteristic, partial.omegasq) %>%
  select(Factor, one_of(DVnames),`Run Time`) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_replace_all("\\.", "-"))

#### Model Checking ####
## QQ and FittedvResid plots for each characteristic
Sim2Model_Check_Plots <- map(Sim2lms, ~ModelCheckPlots(.x, se = FALSE, color = "black", size = 1))

## vif for the models
# outcome characteristic doesn't matter, so only really need to run this for one model since the factors are consistent across models
VIFs <- car::vif(Sim2lms[[1]]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Factor") %>%
  mutate(VIF = `GVIF^(1/(2*Df))`^2) %>%       # for approximate interpretation: https://stats.stackexchange.com/questions/70679/which-variance-inflation-factor-should-i-be-using-textgvif-or-textgvif
  mutate_if(is.numeric, ~round(., digits = 3))

######################################################

##############################################
#### Model Convergence  and PS Estimation ####

#### Convergence ####
## Omegas
Arrange_Omegas(`PS Converged`)        # aggvar, clusters, icc
Arrange_Omegas(`Y Converged`)         # All low, but psmod, clusters
Arrange_Omegas(`Baseline Converged`)  # aggvar, subjects, clusters

## Descriptives
psych::describeBy(Sim2DVsbyCondPooled$PS_Converged, list(Sim2DVsbyCondPooled$nclus, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3) %>%
  filter(group1 == "60")
psych::describeBy(Sim2DVsbyCondPooled$Y_Converged, list(Sim2DVsbyCondPooled$nclus,Sim2DVsbyCondPooled$psmod), mat = TRUE, digits = 3) %>%
  filter(group1 == "60")
psych::describeBy(Sim2DVsbyCondPooled$Baseline_Converged, list(Sim2DVsbyCondPooled$nclus), mat = TRUE, digits = 3)
# min = 90% for 100 and 140 clusters and 60% for 60 clusters

## Plot
# I don't think it makes sense to plot Baseline because the treatment appraisal has no bearing on this
Convergence_BoxPlot <- Sim2DVsbyCondPooled %>%
  select(nsub:mw, contains("_Converged")) %>%
  gather(Model, Value, PS_Converged, Y_Converged) %>%                          #contains("_Converged")
  mutate(Model = str_remove(Model, "_Converged") %>%
           fct_recode(Outcome = "Y") %>%
           fct_relevel("PS","Outcome"),  #,"Baseline"
         ICC = fct_recode(icc, !!!ICCLevels)) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "Proportion Converged", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99")) +
  facet_grid(Model ~ Clusters, labeller = labeller(Clusters = label_both)) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "top")


##############################################################

######################################################
#### Sample Characteristics after PS conditioning ####
## Omegas
Arrange_Omegas(`Analytic Subjects`) # Subjects (98), Clusters (96), mw (96), interactions
Arrange_Omegas(`Analytic Clusters`) # Clusters (98), mw (92), psmod (87), 
Arrange_Omegas(`True-Agg L2 Cor`)   # errmag (1), icc (.99), subjects (93)
Arrange_Omegas(`Y ICC1`)
Arrange_Omegas(`X ICC1`)
Arrange_Omegas(`X ICC2`)

#### Proportion Treated ####
Arrange_Omegas(`Prop Treated`) # ICC (1.0), Subject (.97)
psych::describe(Sim2DVsbyCondPooled$Prop_Treated) # .26 - .33 overall

## Ratio of PS sample to baseline sample
SamplingRatio <- Sim2DVsbyCondPooled %>%
  mutate(Subjects_Per_Cluster = Analytic_Subjects / Analytic_Clusters,
         Sampling_Ratio = Subjects_Per_Cluster / as.numeric(levels(nsub))[nsub]) %>%
  select(nsub:mw, Analytic_Subjects, Analytic_Clusters, Subjects_Per_Cluster, Sampling_Ratio)

SamplingRatioSummary <- SamplingRatio %>%
  group_by(psmod, mw) %>%
  summarize(Median = median(Sampling_Ratio),
            Min = min(Sampling_Ratio),
            Max = max(Sampling_Ratio))

#### Analytic Clusters ####
psych::describeBy(Sim2DVsbyCondPooled$Analytic_Clusters, list(Sim2DVsbyCondPooled$nclus, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$psmod),
                  mat = TRUE, digits = 3)

# Plot
Clusters_BoxPlot <- Sim2DVsbyCondPooled %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Analytic_Clusters", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "PS Clusters", limits = c(0, 140), breaks = c(20, 60, 100, 140)) +
  xlab("Baseline Clusters") +
  scale_fill_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99")) + #c("#e66101","#fdb863","#5e3c99","#b2abd2")
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### Analytic Subjects ####
psych::describeBy(Sim2DVsbyCondPooled$Analytic_Subjects, list(Sim2DVsbyCondPooled$nsub, Sim2DVsbyCondPooled$nclus,
                                                               Sim2DVsbyCondPooled$mw), mat = TRUE, digits = 3) %>%
  filter(group2 == 60 & group1 == 20)

## Subjects
BaselineSubjects <- crossing(nsub = c(20, 60, 100),            # number of subjects per cluster
                             nclus = c(60, 100, 140))  %>%         # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
  mutate(Baseline = (nsub*nclus),
         nsub = as_factor(nsub),
         nclus = as_factor(nclus))

# Not sure yet if I like this as a line or a box plot
Subjects_BoxPlot <- Sim2DVsbyCondPooled %>%
  left_join(BaselineSubjects, by = c("nsub","nclus")) %>%
  # mutate(Baseline = as_factor(Baseline)) %>%
  Plot_Shortcut(xv = "Baseline", yv = "Analytic_Subjects", gpv = "psmod", linebox = "line") +
  scale_y_continuous(name = "PS Subjects", limits = c(0, 14050), breaks = seq(0, 14000, 2000)) +
  scale_x_continuous(name = "Baseline Subjects", limits = c(0, 14050), breaks = seq(0, 14000, 2000)) +
  scale_color_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99")) + #c("#e66101","#fdb863","#5e3c99","#b2abd2")
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#### ICCs ####
psych::describeBy(Sim2DVsbyCondPooled$X_ICC1, list(Sim2DVsbyCondPooled$icc, Sim2DVsbyCondPooled$mw),
                  mat = TRUE, digits = 3)



ExpectedICC2 <- crossing(nsub = c(20, 60, 100),           # number of subjects per cluster
                         icc = c(.05, .1, .2))  %>%         # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
  mutate(Expected = (nsub*icc) / (1 + (nsub - 1)*icc),
         nsub = factor(nsub),
         icc = factor(icc))
# Final table
ICC2Table <- Sim2DVsbyCondPooled %>%
  group_by(icc, nsub, mw) %>%
  summarize(Median = median(X_ICC2),#) %>%
            Min = min(X_ICC2),
            Max = max(X_ICC2)) %>%
  ungroup() %>%
  # bind_rows(ExpectedICC2) %>%
  inner_join(ExpectedICC2, by = c("icc", "nsub")) %>%
  mutate_if(is.numeric, ~round(.,2)) %>%
  # spread(mw, Median)
  unite(col = "Range", Min, Max, sep = " - ") %>%
  select(ICC1 = icc, Subjects = nsub, Expected, mw, Median, Range)
ICC2Table <- left_join(ICC2Table %>% filter(mw == "Matching") %>% select(ICC1, Subjects, Expected, Median_m = Median, Range_m = Range),
                       ICC2Table %>% filter(mw == "Weighting") %>% select(ICC1, Subjects, Expected, Median_w = Median, Range_w = Range))

#### True-Agg Correlation ####
psych::describeBy(Sim2DVsbyCondPooled$True.Agg_L2_Cor, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$True.Agg_L2_Cor, list(Sim2DVsbyCondPooled$icc), mat = TRUE, digits = 3) 
psych::describeBy(Sim2DVsbyCondPooled$True.Agg_L2_Cor, list(Sim2DVsbyCondPooled$nsub), mat = TRUE, digits = 3) 

## Pooling the correlation in the baseline sample across the matching and weighting conditions
BaselineCorrs <- Sim2DVsbyCondPooled %>%
  group_by(nsub, nclus, icc, aggvar) %>%
  summarize(True.Agg_L2_Cor = mean(Baseline_True.Agg_L2_Cor, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(psmod = "Baseline")
BaselineCorrs <- bind_rows(BaselineCorrs %>% mutate(mw = "Matching"), BaselineCorrs %>% mutate(mw = "Weighting")) %>%
  mutate(mw = as_factor(mw))

## Plot
TrueAggCor_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineCorrs) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
         ICC = fct_recode(icc, !!!ICCLevels)) %>%
  Plot_Shortcut(xv = "aggvar", yv = "True.Agg_L2_Cor", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "Correlation", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) + # "#fdb863","#b2abd2"
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))
#############################################


###########################################
#### Covariate Balance (RQ1, 2, and 3) ####

## Omegas 
# the unadjusted values should be extremely similar between weighted and matched conditions; differences may be due to trimming
Arrange_Omegas(`Un-ASD`)                # Subjects, Clusters, ICC, (all .97,.98) then aggvar (.94), interacts, psmodel (.55)
Arrange_Omegas(`Un-VR`)                 # Clusters (1), Subjects (.86), ICC (.76), aggvar (.70)
Arrange_Omegas(`Un-ASD Balanced Count`) # Subjects, Clusters, ICC (all .98,.99), psmodel (.65), aggvar (.46)
Arrange_Omegas(`Un-VR Balanced Count`)  # Clusters (.99), aggvar (.62), ICC
Arrange_Omegas(`Adj-ASD`)               # psmod, Clusters, ICC, mw, aggvar
Arrange_Omegas(`Adj-VR`)                # Nothing? have an outlier affecting variance
Arrange_Omegas(`Adj-ASD Balanced Count`)# psmod, Clusters, ICC, mw, aggvar
Arrange_Omegas(`Adj-VR Balanced Count`) # Clusters, aggvar, Subjects, psmod, mw

## Descriptives
# ASD
psych::describeBy(Sim2DVsbyCondPooled$Adj.ASD, list(Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$mw), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$Adj.ASD_Balanced_Count, list(Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$mw),mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$Adj.ASD, list(Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$aggvar),
                  mat = TRUE, digits = 3) %>% filter(group1 == "Subject" & group2 == "Weighting")
psych::describeBy(Sim2DVsbyCondPooled$Adj.ASD, list(Sim2DVsbyCondPooled$nclus), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$Adj.ASD_Balanced_Count, list(Sim2DVsbyCondPooled$nclus),mat = TRUE, digits = 3)

# VR
Sim2DVsbyCondPooled %>% select(Adj.VR, Adj.VR_Balanced_Count, nsub:mw) %>% View() # filter(Adj.VR > 2.0)
psych::describeBy(Sim2DVsbyCondPooled$Adj.VR, list(Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$nclus), mat = TRUE, digits = 3) %>%
  filter(group3 == "60")
psych::describeBy(Sim2DVsbyCondPooled$Adj.VR_Balanced_Count, list(Sim2DVsbyCondPooled$psmod), mat = TRUE, digits = 3)

psych::describeBy(Sim2DVsbyCondPooled$Adj.VR, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$Adj.VR_Balanced_Count, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)

#### Unadjusted Values ####
## Mean
UnadjustedBalanceMean <- Sim2DVsbyCondPooled %>%
  gather(Variable, Temp, starts_with("Un.")) %>%
  filter(!(str_detect(Variable, "Balanced_"))) %>%
  mutate(Measure = case_when(str_detect(Variable, "ASD") ~ "Abs. Std. Difference",
                             str_detect(Variable, "VR") ~ "Variance Ratio",
                             TRUE ~ NA_character_),
         psmod = "Baseline") %>%
  group_by(nsub, nclus, icc, aggvar, psmod, Measure) %>%
  summarize(Value = mean(Temp, na.rm = TRUE)) %>% ungroup()
UnadjustedBalanceMean <- bind_rows(UnadjustedBalanceMean %>% mutate(mw = "Matching"), UnadjustedBalanceMean %>% mutate(mw = "Weighting")) %>%
  mutate(mw = as_factor(mw))

## Count
UnadjustedBalanceCount <- Sim2DVsbyCondPooled %>%
  gather(Variable, Temp, starts_with("Un.")) %>%
  filter(str_detect(Variable, "Balanced_")) %>%
  mutate(Measure = case_when(str_detect(Variable, "ASD") ~ "Abs. Std. Difference",
                             str_detect(Variable, "VR") ~ "Variance Ratio",
                             TRUE ~ NA_character_),
         psmod = "Baseline") %>%
  group_by(nsub, nclus, icc, aggvar, psmod, Measure) %>%
  summarize(Value = mean(Temp, na.rm = TRUE)) %>% ungroup()
UnadjustedBalanceCount <- bind_rows(UnadjustedBalanceCount %>% mutate(mw = "Matching"), UnadjustedBalanceCount %>% mutate(mw = "Weighting")) %>%
  mutate(mw = as_factor(mw))

psych::describe(UnadjustedBalanceCount[UnadjustedBalanceCount$Measure == "Abs. Std. Difference",]$Value)
psych::describe(UnadjustedBalanceCount[UnadjustedBalanceCount$Measure == "Variance Ratio",]$Value)

#### Plotting Mean ASD and Count ####
# Mean
MeanASD_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(UnadjustedBalanceMean %>% rename(Adj.ASD = Value)%>%
              filter(Measure != "Variance Ratio")) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline")) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Adj.ASD", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "Mean", limits = c(0, .22), breaks = seq(0, .20, .05)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) + #c("#e66101","#fdb863","#5e3c99","#b2abd2")
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        legend.position = "top",
        axis.title.x = element_blank())
# Count
CountASD_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(UnadjustedBalanceCount %>% rename(Adj.ASD_Balanced_Count = Value) %>%
              filter(Measure != "Variance Ratio")) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline")) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Adj.ASD_Balanced_Count", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "Covariate Count", limits = c(0, 30), breaks = seq(0, 30, 5)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE, name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) + #c("#e66101","#fdb863","#5e3c99","#b2abd2")
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))


# Juxtaposition of ASD Mean and Count
ASD_BoxPlot <- cowplot::plot_grid(MeanASD_BoxPlot, CountASD_BoxPlot, labels = "AUTO", align = "v", ncol = 1, rel_heights = c(1, .9))


#####################################################


###############################################################
#### Treatment Effect Estimation (addresses RQ1, 2, and 3) ####

Arrange_Omegas(`TE Bias`) # Error Magnitude (.97), Subjects (.62), ICC (.61), and interactions of the 3 ~.5
Arrange_Omegas(`TE MAE`) # Error Magnitude (.79), ICC (.40), Subjects (.32), clusters (.21) and interactions of with errormag
Arrange_Omegas(`TE RMSE`) # Clusters (.09) and interactions
Arrange_Omegas(`Baseline Bias`) # Error Magnitude (1), ICC (.94), Subjects (.93) and interactions of the 3 ~.9
Arrange_Omegas(`Baseline MAE`) # Error Magnitude (1), ICC (.96), Subjects (.94) and interactions of the 3 ~.9
Arrange_Omegas(`Baseline RMSE`) # Error Magnitude (1), ICC (.97), Subjects (.94) and interactions of the 3 ~.9

#### Descriptives ####
## The good conditions at 0.1
Sim2DVsbyCondPooled %>% filter(aggvar == "0.1") %>% filter(!(psmod == "Cluster" & mw == "Matching" & nclus == "60")) %>%
  gather(Measure, Value, TE_Bias, TE_MAE, TE_RMSE) %>%
  group_by(Measure) %>%
  summarize(Median = median(Value),
            SD = sd(Value))
## The bad at 0.1
psych::describeBy(Sim2DVsbyCondPooled$TE_Bias, list(Sim2DVsbyCondPooled$nclus, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3) %>%
  filter(group4 == "0.1" & group1 == "60")
psych::describeBy(Sim2DVsbyCondPooled$TE_MAE, list(Sim2DVsbyCondPooled$nclus, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3) %>%
  filter(group4 == "0.1" & group1 == "60")
psych::describeBy(Sim2DVsbyCondPooled$TE_RMSE, list(Sim2DVsbyCondPooled$nclus, Sim2DVsbyCondPooled$mw, Sim2DVsbyCondPooled$psmod, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3) %>%
  filter(group4 == "0.1" & group1 == "60")

## by aggvar
psych::describeBy(Sim2DVsbyCondPooled$TE_Bias, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$TE_MAE, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$TE_RMSE, list(Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)

## Gathering Baseline Info
BaselineTE <- Sim2DVsbyCondPooled %>%
  gather(Measure, Temp, Baseline_Bias, Baseline_MAE, Baseline_RMSE) %>%
  mutate(Measure = str_remove(Measure, "Baseline_"),
         psmod = "Baseline") %>%
  group_by(nsub, nclus, icc, aggvar, psmod, Measure) %>%
  summarize(Value = mean(Temp, na.rm = TRUE)) %>% ungroup()
BaselineTE <- bind_rows(BaselineTE %>% mutate(mw = "Matching"), BaselineTE %>% mutate(mw = "Weighting")) %>%
  mutate(mw = as_factor(mw))
# left_join(TempBaselineCorrs, by = c("nsub", "nclus", "icc", "aggvar", "psmod", "mw"))

psych::describeBy(BaselineTE[BaselineTE$mw=="Matching",]$Value, list(BaselineTE[BaselineTE$mw=="Matching",]$Measure), mat = TRUE, digits = 3)

#### Using Error Magnitude ####
# These get at all the RQs (PSmod = RQ1, correlation = RQ2, ICC or Subjects = RQ3)
## Bias
TEBias_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline")) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_Bias", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "Bias") +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) +
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), axis.title.x = element_blank(), legend.position = "top")
## MAE
TEMAE_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline")) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_MAE", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "MAE", limits = c(0, 1.06), breaks = seq(0, 1, .20)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE,name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) +
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), axis.title.x = element_blank())
## RMSE
TERMSE_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline")) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_RMSE", gpv = "psmod", linebox = "box") +
  scale_y_continuous(name = "RMSE")+#, limits = c(0, 1.06), breaks = seq(0, 1, .20)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE,name = "Treatment Appraisal", values = c("#e66101","#5e3c99","#4daf4a")) +
  facet_grid(. ~ mw) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

# Juxtaposition of Treatment Effect Plots
TE_BoxPlot <- cowplot::plot_grid(TEBias_BoxPlot, TEMAE_BoxPlot, TERMSE_BoxPlot,
                                 labels = "AUTO", align = "v", ncol = 1, rel_heights = c(1, .9, .9))

#### by ICC ####
## Bias
TEBias_ICC_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
         ICC = fct_recode(icc, !!!ICCLevels)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_Bias", gpv = "ICC", linebox = "box") +
  scale_y_continuous(name = "Bias") +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE, name = "ICC(1)", values = c("#bae4bc","#7bccc4","#2b8cbe")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), axis.title.x = element_blank())
## MAE
TEMAE_ICC_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
         ICC = fct_recode(icc, !!!ICCLevels)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_MAE", gpv = "ICC", linebox = "box") +
  scale_y_continuous(name = "MAE", limits = c(0, 1.06), breaks = seq(0, 1, .20)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE, name = "ICC(1)", values = c("#bae4bc","#7bccc4","#2b8cbe")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))
## RMSE
TERMSE_ICC_BoxPlot <- Sim2DVsbyCondPooled %>%
  bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
  mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
         ICC = fct_recode(icc, !!!ICCLevels)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "TE_RMSE", gpv = "ICC", linebox = "box") +
  scale_y_continuous(name = "RMSE")+#, limits = c(0, 1.06), breaks = seq(0, 1, .20)) +
  scale_x_discrete(name = "Error Magnitude (SD)", labels = c("0.1","0.3", "0.5", "1.0")) +
  scale_fill_manual(guide = FALSE,name = "ICC(1)", values = c("#bae4bc","#7bccc4","#2b8cbe")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

# Get Legend
TE_ICC_Legend <- cowplot::get_legend(
  TERMSE_ICC_BoxPlot + 
    guides(fill = guide_legend(title = "ICC(1)")) +
    theme(legend.position = "right")
)

# Juxtaposition of Treatment Effect Plots
TE_ICC_BoxPlot <- cowplot::plot_grid(TEBias_ICC_BoxPlot, TE_ICC_Legend, TEMAE_ICC_BoxPlot, TERMSE_ICC_BoxPlot,
                                     labels = c("A","","B","C"), axis = "l", align = "v", ncol = 2, rel_heights = c(.9, .9, 1, 1))

psych::describeBy(Sim2DVsbyCondPooled$TE_Bias, list(Sim2DVsbyCondPooled$icc, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$TE_MAE, list(Sim2DVsbyCondPooled$icc, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)
psych::describeBy(Sim2DVsbyCondPooled$TE_RMSE, list(Sim2DVsbyCondPooled$icc, Sim2DVsbyCondPooled$aggvar), mat = TRUE, digits = 3)


#### Using True-Agg Correlation ####
# # Note, can't really do box plots with a continuous x
# ## Bias
# TEBias_LinePlot <- Sim2DVsbyCondPooled %>%
#   bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
#   mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
#          ICC = fct_recode(icc, !!!ICCLevels),
#          `Error Magnitude` = fct_rev(aggvar)) %>%
#   rename(Subjects = nsub) %>%
#   Plot_Shortcut(xv = "True.Agg_L2_Cor", yv = "TE_Bias", gpv = "psmod", linebox = "line", smooth_method = "lm") +
#   scale_y_continuous(name = "Bias") +
#   scale_x_continuous(name = "True-Aggregated L2 Correlation") +
#   scale_color_manual(name = "PS Model", values = c("#e66101","#5e3c99","#4daf4a")) +
#   facet_grid(. ~ ICC, labeller = label_both) +                            # Subjects or ICC; don't need both
#   theme(strip.background = element_rect(colour = "black", fill = "white"))
# 
# ## MAE
# TEmae_LinePlot <- Sim2DVsbyCondPooled %>%
#   bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
#   mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
#          ICC = fct_recode(icc, !!!ICCLevels)) %>%
#   rename(Subjects = nsub) %>%
#   Plot_Shortcut(xv = "True.Agg_L2_Cor", yv = "TE_MAE", gpv = "psmod", linebox = "line", smooth_method = "lm") +
#   scale_y_continuous(name = "MAE") +
#   scale_x_continuous(name = "True-Aggregated L2 Correlation") +
#   scale_color_manual(name = "Model Specification", values = c("#e66101","#5e3c99","#4daf4a")) +
#   facet_grid(. ~ Subjects, labeller = label_both) +                            # Subjects or ICC; don't need both
#   theme(strip.background = element_rect(colour = "black", fill = "white"))
# 
# ## RMSE
# TErmse_LinePlot <- Sim2DVsbyCondPooled %>%
#   bind_rows(BaselineTE %>% mutate(Measure = paste("TE", Measure, sep = "_")) %>% spread(Measure, Value)) %>%
#   mutate(psmod = fct_relevel(psmod, "Cluster", "Subject", "Baseline"),
#          ICC = fct_recode(icc, !!!ICCLevels)) %>%
#   rename(Subjects = nsub) %>%
#   Plot_Shortcut(xv = "True.Agg_L2_Cor", yv = "TE_RMSE", gpv = "psmod", linebox = "line", smooth_method = "lm") +
#   scale_y_continuous(name = "RMSE") +
#   scale_x_continuous(name = "True-Aggregated L2 Correlation") +
#   scale_color_manual(name = "Model Specification", values = c("#e66101","#5e3c99","#4daf4a")) +
#   facet_grid(. ~ Subjects, labeller = label_both) +                            # Subjects or ICC; don't need both
#   theme(strip.background = element_rect(colour = "black", fill = "white"))

############################################



#############################################################
#### Saving/Exporting Summary Objects, Plots, and Tables ####

#### Saving summary objects ####
save(Sim2Time, Sim2lms, Sim2Omegas, Sim2Omegas_Wide,
     Sim2ConditionRename, ICCLevels, Sim2Model_Check_Plots, VIFs,
     Convergence_BoxPlot, Subjects_BoxPlot, Clusters_BoxPlot, # Clusters_Table, Subjects_Table, 
     TrueAggCor_BoxPlot, ICC2Table,
     UnadjustedBalanceCount, UnadjustedBalanceMean, BaselineTE, TempBaselineCorrs,
     MeanASD_BoxPlot, CountASD_BoxPlot, ASD_BoxPlot,
     TEBias_BoxPlot, TEMAE_BoxPlot, TERMSE_BoxPlot, TE_BoxPlot,
     TEBias_ICC_BoxPlot, TEMAE_ICC_BoxPlot, TERMSE_ICC_BoxPlot, TE_ICC_BoxPlot,
     file = "Sim2_Thousand_SummaryResults.RData")

#### Exporting Model Check Plots ####
map2(.x = Sim2Model_Check_Plots, .y = names(Sim2Model_Check_Plots),
     ~ggsave(.x, file = paste0("Saved Sim2 Results/Model_Check_Plot_Pooled_", .y, ".png"),
             height = 5, width = 10))

#### Exporting Research Question Plots ####
## Convergence
ggsave(Convergence_BoxPlot, file = "Sim2_Pooled_Convergence_BoxPlot.png", width = 10, height = 8)

## Analytic Sample
ggsave(Clusters_BoxPlot, file = "Sim2_Pooled_Clusters_BoxPlot.png", width = 10, height = 5)
ggsave(Subjects_BoxPlot, file = "Sim2_Pooled_Subjects_BoxPlot.png", width = 10, height = 5)
ggsave(TrueAggCor_BoxPlot, file = "Sim2_Pooled_TrueAggCor_BoxPlot.png", width = 10, height = 5)

## ASD
ggsave(ASD_BoxPlot, file = "Sim2_Pooled_ASD_BoxPlot.png", width = 9, height = 9)

## Treatment Effect
ggsave(TE_BoxPlot, file = "Sim2_Pooled_TE_BoxPlot.png", width = 9, height = 12)
ggsave(TE_ICC_BoxPlot, file = "Sim2_Pooled_TE_ICC_BoxPlot.png", width = 10, height = 8)


## Export Omega Table
openxlsx::write.xlsx(Sim2Omegas_Wide, file = "Sim2_Pooled_OmegaTable.xlsx")
