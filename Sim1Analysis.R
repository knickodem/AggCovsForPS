################################################
#                                              #
#   Simulation 1 - Generate L2 Aggregations    #
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
# load("Sim1_10Reps.RData")
load("Saved Sim1 Results/Sim1_ThousandReps.RData")


#### Converting logits to probabilities and vice versa ####
LogitToProb <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#### Creates QQ, fitted vs residual, and fitted vs observed plots ####
## works for lm, lmer, glm, glmer (I think)
ModelCheckPlots <- function(model, smoother = "loess"){
  
  FR <- data.frame(Fitted = fitted(model),
                   Residuals = resid(model)) %>%
    tibble::rownames_to_column("ID")
  
  qq <- ggplot(FR, aes(sample = Residuals)) +
    stat_qq(shape = 1) + stat_qq_line() +
    labs(x = "Theoretical", y = "Observed") +
    theme_bw(base_size = 18)
  
  ## Fitted vs. Residuals Plot
  frplot <- FR %>%
    ggplot(aes(x = Fitted, y = Residuals)) +
    geom_point(shape = 1) +
    geom_hline(aes(yintercept = 0), color = "grey", size = 1, linetype = 1) +
    geom_smooth(method = smoother, se = FALSE, color = "black", size = 2, linetype = 1) +
    theme_bw(base_size = 18)
  
  MCPlots <- cowplot::plot_grid(qq, frplot, nrow = 1, align = "h", labels = "AUTO")
  
  return(MCPlots)
}

## Shortcut for investigating factors with highest omegas
Arrange_Omegas <- function(Outcome){
  
  ot <- enquo(Outcome)
  TheOmegas_Wide %>% select(Factor, !!ot) %>% arrange(desc(!!ot))
  
}

#### Shortcut for line and boxplots ####
Plot_Shortcut <- function(dat, xv, yv, gpv, linebox = c("line","box")){
  
  if(linebox == "line"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv, group = gpv, color = gpv)) +
      geom_point(size = 3) + geom_smooth(method = "loess", se = FALSE, size = 2) +
      theme_bw(base_size = 20)
    
  } else if(linebox == "box"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv, fill = gpv)) +
      geom_boxplot(notch = FALSE) +
      theme_bw(base_size = 20)
    
  }
  
  return(deplot)
}

#### Getting summary descriptives on a continuous variable ####
## unquoted grouping variables should be specified in ...
Get_Descriptives <- function(data, ContinuousVariable, ..., digits=5, AllContinuous=TRUE){
  
  groups <- quos(...)
  CV <- enquo(ContinuousVariable)
  
  data_descrip <- data %>% group_by(!!!groups) %>%
    summarize(n = sum(!is.na(!!CV)),
              Median = median(!!CV, na.rm=TRUE),
              Mean = mean(!!CV, na.rm=TRUE),
              SD = sd(!!CV, na.rm=TRUE),
              Min = min(!!CV, na.rm=TRUE),
              Max = max(!!CV, na.rm=TRUE)) %>% ungroup()
  
  data_descrip <- data_descrip %>%
    mutate_if(is.numeric, ~round(., digits = digits))
  
  return(data_descrip)
}

########################################################



########################################################
#### Defining Factors for analysis and presentation ####


## Converting independent variables to factors
DGSbyCond <- DGSbyCond %>%
  mutate_at(vars(nsub:method), ~as_factor(.))

## Order independent variables (and interactions) for use in regressions or long formatted data
ConditionVars <- c("nsub", "nclus", "icc", "aggvar", "method",
                   "nsub:nclus", "nsub:icc", "nsub:aggvar", "nsub:method",
                   "nclus:icc", "nclus:aggvar", "nclus:method",
                   "icc:aggvar", "icc:method", "aggvar:method")

## Converting characters for output
ConditionFactor <- str_replace(ConditionVars, "nsub", "Subjects") %>%
  str_replace("nclus","Clusters") %>%
  str_replace("icc","ICC(1)") %>%
  str_replace("aggvar","Error Magnitude") %>%
  str_replace("method","Method")

## Renamed levels
ConditionRename <- ConditionVars
names(ConditionRename) <- ConditionFactor

#### Recoding factor levels ####

## ICC
ICCfac <- c(`.20` = "0.2", `.10` = "0.1", `.05` = "0.05")

## Correlations
Corrfac <- c(`Level 1` = "L1_Cor", `True Level 2` = "True_L2_Cor", `Observed Level 2` = "Obs_L2_Cor")


##################################################################


###########################################################################
#### Variation in outcomes explained by condition (partial-$\Omega^2$) ####

#### In timing ####
TheTimes <- lm(as.formula(paste0("Time ~ ", paste(ConditionVars, collapse = " + "))),
               data = DataGenTLogDF)
TimesAnova <- anova(TheTimes) # I don't actually use this
TimesOmega <- sjstats::omega_sq(TheTimes, partial = TRUE, ci.lvl = .95) %>%
  mutate(Characteristic = "Run Time")

#### In Sample Characteristics, then binding timing ####
## Regressions
Thelms <- map(DGSnames,~lm(as.formula(paste(.x, "~", paste(ConditionVars, collapse = " + "))),
                           data = DGSbyCond)) %>%
  set_names(DGSnames)
TheAnovas <- map(Thelms,~anova(.x)) # I don't actually use this

## Effect size
#Note: Albers & Lakens (2017, p.9) says its okay that there are negative numbers even if it is theoretically impossible
TheOmegas <- map_dfr(Thelms, ~sjstats::omega_sq(.x, partial = TRUE, ci.lvl = .95), .id = "Characteristic") %>%
  bind_rows(TimesOmega) %>%
  mutate(term = factor(term, levels = ConditionVars),
         Factor = fct_recode(term, !!!ConditionRename)) %>%
  mutate_if(is.numeric, ~round(.,2))

# Both packages give similar results (CIs aren't expected to be exactly the same given the random draws from the bootstrap)
# Based on my hand calculations, both  appear to use the formulas from Albers & Lakens (2017, p. 35). Not sure about the bootstrapped CIs
# partial = FALSE constrains lower bound to 0, but still gives odd CIs
# In summary, I don't think the CIs can be trusted from either package despite giving similar ranges
# sjstats::omega_sq(Thelms$Y_Variance, partial = TRUE, ci.lvl = .95, method = "quantile")
# groupedstats::lm_effsize_ci(Thelms$Y_Variance, effsize = "omega",partial = TRUE, conf.level = 0.95,nboot = 1000, method = "quantile")

## Wide format; for final paper, bold values > .14, which are considered large
# Table is exported at end of code
TheOmegas_Wide <- TheOmegas %>%
  select(Factor,Characteristic,partial.omegasq) %>%
  spread(Characteristic, partial.omegasq) %>%
  select(Factor,one_of(DGSnames),`Run Time`) %>%
  rename_all(~str_replace_all(.,"_"," ") %>% str_replace_all("\\.","-"))

#### Omega Plot ####
# Note: does not help organize and interpret information
# # set.seed(3854)
# OmegaPlot <- TheOmegas %>%
#   mutate(Characteristic = str_replace_all(Characteristic,"_"," ") %>% str_replace_all("\\.","-"),
#          Characteristic = as_factor(Characteristic),
#          Factor = fct_rev(Factor)) %>%
#   ggplot(aes(x = Factor, y = partial.omegasq)) + #, shape = Factor
#   geom_hline(yintercept = c(.01, .06, .14), linetype = "solid", size = 1, color = "grey" ) +
#   geom_point(size = 2) +
#   # geom_jitter(width = .1, size = 4) +
#   # geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
#   scale_y_continuous(name = expression(omega^2*italic(p)), breaks = seq(0, 1, .1)) + #, limits = c(-0.01,.51)
#   # scale_color_grey() +
#   # scale_shape_manual(values = c(16, 1, 8, 17, 2, 18, 4, 9)) +
#   theme_classic(base_size = 16) +
#   theme(legend.position = "none") +
#   # theme(legend.title = element_blank(), axis.title.y = element_blank(), #, axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#   #       legend.justification=c(1,0), legend.position=c(1,.03)) +
#   coord_flip() +
#   facet_wrap(~Characteristic)
# 
# ggsave(OmegaPlot, file = "Omega2p.png", width = 12, height = 6)

#### Model Checking ####

## QQ and FittedvResid plots for each characteristic
Model_Check_Plots <- map(Thelms, ~ModelCheckPlots(.x))

## vif for the models
# outcome characteristic doesn't matter, so only really need to run this for one model since the factors are consistent across models
VIFs <- map_dfr(Thelms, ~car::vif(.x) %>% as.data.frame() %>% tibble::rownames_to_column("Factor"),.id = "Characteristic") %>%
  mutate(VIF = `GVIF^(1/(2*Df))`^2) %>%       # for approximate interpretation: https://stats.stackexchange.com/questions/70679/which-variance-inflation-factor-should-i-be-using-textgvif-or-textgvif
  mutate_if(is.numeric,~round(., digits = 3))

## Note: Create tables only if needed to report specific values in text, but not needed for display
#### Summary descriptives of all depenendent variables ####
DVSummary <- DGSbyCond %>%
  gather(Characteristic, Value, one_of(DGSnames)) %>%
  filter(!str_detect(Characteristic, "Logit")) %>%
  Get_Descriptives(Value, Characteristic, digits = 3)

######################################################


#############################################
#### RQ1: True and Obs L2 Corr by Method ####

## TheOmegas shows that the True-Obs L2 Corr is associated with everything but nclusters
# summary(Thelms$True.Obs_L2_Cor) # doesn't tell me much
Arrange_Omegas(`True-Obs L2 Cor`)

## Box Plot
TrueObsCorbyMethodError_BoxPlot <- DGSbyCond %>%
  Plot_Shortcut(xv = "aggvar", yv = "True.Obs_L2_Cor", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "True and Observed L2 Correlation", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  theme(legend.justification = c(0, 0), legend.position = c(.1, .1))

#########################################################

#########################################################
#### RQ2: True and Obs L2 Corr by Contextual Factors ####

## Loess plot including all manipulated factors except number of clusters
# Note: doing a boxplot isn't helpful because each box would only summarize 4 data points
TrueObsCor_LoessPlot <- DGSbyCond %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "True.Obs_L2_Cor", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "True and Observed L2 Correlation", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(ICC ~ Subjects, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "top")

#########################################################

################################################
#### RQ3: Sample characteristic differences ####

# variances of and correlations between L1 covariates, true L2 covariates, and observed L2 covariates;
# proportion of subjects in treatment and control groups; and variance and ICC(1) of the subject-level outcome

## Sample Characteristic Values in Long Format
DGSbyCond_Long <- DGSbyCond %>%
  gather(Characteristic, Value,one_of(DGSnames))

#### Summary Tables by IVs ####
# Note: does not include logit characteristics
## Summarized by Method
CharacteristicsByMethod <- DGSbyCond_Long %>%
  filter(!str_detect(Characteristic, "Logit")) %>%
  Get_Descriptives(Value, method, Characteristic, digits = 3)

## Summarized by ICC
CharacteristicsByICC <- DGSbyCond_Long %>%
  filter(!str_detect(Characteristic, "Logit")) %>%
  Get_Descriptives(Value, icc, Characteristic, digits = 3)

## Summarized by Method and ICC
CharacteristicsByMethodICC <- DGSbyCond_Long %>%
  filter(!str_detect(Characteristic, "Logit")) %>%
  Get_Descriptives(Value, method, icc, Characteristic, digits = 3)

#### Variances ####
Arrange_Omegas(`Y Variance`)       # 100% Method and ICC, also subjects and clusters ; Expected to vary based on ICC (and clusters?)
Arrange_Omegas(`L1 Variance`)      # Clusters and ICC                                ; Expected to vary by ICC, subjects and clusters?
Arrange_Omegas(`True L2 Variance`) # ~100% ICC, Method, and Subjects, also clusters  ; Expected to vary by ICC and clusters
Arrange_Omegas(`Obs L2 Variance`)  # ~100% aggvar, Method, also ICC and subjects     ; Expected to vary by ICC and clusters


# ## Observed L2
# L2Variance_BoxPlot <- DGSbyCond_Long %>%
#   filter(str_detect(Characteristic, "Obs_L2_Variance")) %>%
#   mutate(ICC = fct_recode(icc, !!!ICCfac),
#          Subjects = fct_rev(nsub)) %>%
#   rename(Clusters = nclus) %>%
#   Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "method", linebox = "box") +
#   scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
#   scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
#   scale_y_continuous(name = "Observed L2 Variance", limits = c(0, 1.23), breaks = seq(0, 1.2, .2)) +
#   theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "none")
# 
# ## Observed L2 juxtaposed with true-observed L2 correlation
# CorrVar_BoxPlot <- cowplot::plot_grid(TrueObsCorbyMethodError_BoxPlot, L2Variance_BoxPlot, labels = "AUTO", align = "hv", ncol = 1)
# ggsave(CorrVar_BoxPlot, file = "Juxtapose_L2Corr_and_Var_BoxPlot.png", width = 10, height = 14)


## Loess plot including all manipulated factors except number of clusters
ObsL2Var_LoessPlot <- DGSbyCond %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Obs_L2_Variance", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "Observed L2 Variance", limits = c(0, 1.23), breaks = seq(0, 1.2, .2)) +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(ICC ~ Subjects, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "none")

## Observed L2 juxtaposed with true-observed L2 correlation
CorrVar_LoessPlot <- cowplot::plot_grid(TrueObsCor_LoessPlot, ObsL2Var_LoessPlot,
                                        labels = "AUTO", align = "hv", ncol = 1, rel_heights = c(1, .9))
ggsave(CorrVar_LoessPlot, file = "Juxtapose_L2Corr_and_Var_LoessPlot.png", width = 10, height = 14)

## Y
YVariance_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Y_Variance")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Subjects = fct_rev(nsub)) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  scale_y_continuous(name = "Variance of Y") +
  xlab("ICC(1) Condition") +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#### ICCs ####
Arrange_Omegas(`Y ICC1`)           # 100% ICC and Method, also subjects and clusters ; Expected to be equal to ICC condition
Arrange_Omegas(`X ICC1`)           # 100% ICC, then clusters                         ; Expected to be equal to ICC condition
Arrange_Omegas(`X ICC2`)           # 100% ICC and Subjects, also clusters            ; Expected to vary based on ICC and nsubjects conditions

ExpectedICC2 <- crossing(nsub = c(20, 60, 100),           # number of subjects per cluster
                         icc = c(.05, .1, .2))  %>%         # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
  mutate(Value = (nsub*icc) / (1 + (nsub - 1)*icc),
         nclus = "Expected",
         nsub = factor(nsub),
         icc = factor(icc))

## X_ICC2
XICC2_LoessPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "X_ICC2")) %>%
  bind_rows(ExpectedICC2) %>%
  mutate(Clusters = as_factor(nclus),
         ICC = fct_recode(icc, !!!ICCfac)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "Clusters", linebox = "line") +
  scale_color_manual(values = c("#cccccc","#969696", "#636363", "#525252")) + # Grey scale:
  scale_y_continuous(name = "Mean ICC(2)") +
  xlab("ICC(1) Condition") +
  facet_grid(. ~ Subjects, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        legend.justification = c(1, 0), legend.position = c(.95, .05))

## Y_ICC1
YICC1_LoessPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Y_ICC1")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Subjects = fct_rev(nsub)) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "line") +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  scale_y_continuous(name = "Observed ICC(1) of Y", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  xlab("ICC(1) Condition") +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        legend.justification = c(0, 1), legend.position = c(.05, .95))

#### Correlations ####
Arrange_Omegas(`L1 Cor`)           # 100% Methods, then interact w/ ICC and clusters ; Expected to be .2 across all conditions
Arrange_Omegas(`True L2 Cor`)      # Methods, subjects and ICC                       ; Expected to be .2 across all conditions
Arrange_Omegas(`Obs L2 Cor`)       # Method, aggvar, interact w/ ICC and subjects    ; Expected to be lower than .2 given the introduction of error

## Prep
CorrsPlotPrep <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "_Cor") & !str_detect(Characteristic, "\\.")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.),
         Characteristic = as_factor(Characteristic) %>% fct_recode(!!!Corrfac))

## Plot
Corrs_LoessPlot <- CorrsPlotPrep %>%
  Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "Correlation") +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  # scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(ICC ~ Characteristic) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))


#### Proportion in Treatment group ####
Arrange_Omegas(`Prop Z1`)          # ICC, Subjects, Method

PropZ_BoxPlot <- DGSbyCond %>%
  Plot_Shortcut(xv = "nsub", yv = "Prop_Z1", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "Proportion in Treatment Group") +
  scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(~icc)




#### RQ4: bias and precision of PS estimation ####
# scale_fill_manual(name = "ICC(1)", values = c("#b2e2e2","#66c2a4","#238b45")) + #BlueGreen ; Grey scale: c("#cccccc","#969696","#525252")

Arrange_Omegas(`Logit Bias`) # mainly clusters, some ICC and error mag
Arrange_Omegas(`Logit MAE`) # Likewise
Arrange_Omegas(`Logit RMSE`) # Likewise
Arrange_Omegas(`PS Bias`) # icc, then subjects and interactions
Arrange_Omegas(`PS MAE`) # icc, aggvar, clusters, method in that order
Arrange_Omegas(`PS RMSE`) # Clusters, aggvar, icc, method in that order


## Bias Box Plot
Bias_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_Bias")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Subjects = fct_rev(nsub)) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "Bias of PS") +
  xlab("ICC(1) Condition") +
  scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(. ~ Subjects) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))


## MAE Box Plot
MAEbyClusterICCError_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_MAE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Clusters = fct_rev(nclus)) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "MAE of PS") + #, limits = c(0, 1), breaks = seq(0,1,.2)) +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(Clusters ~ ICC) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

## RMSE Box Plot
RMSEbyClusterICCError_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_RMSE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Clusters = fct_rev(nclus)) %>%
  Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "MAE of PS") + #, limits = c(0, 1), breaks = seq(0,1,.2)) +
  scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(Clusters ~ ICC) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

###################################################################

#############################################################
#### Saving/Exporting Summary Objects, Plots, and Tables ####

## Saving summary objects
save(ConditionRename, TheTimes, Thelms, TheOmegas,
     Model_Check_Plots, VIFs, DVSummary,
     TrueObsCorbyMethodError_BoxPlot,
     file = "Sim1_ThousandSummaryResults.RData")

## Exporting Model Check Plots
map2(.x = Model_Check_Plots, .y = names(Model_Check_Plots),
     ~ggsave(.x, file = paste0("Saved Sim1 Results/Model_Check_Plot_", .y, ".png"), height = 5, width = 10))

## Exporting Research Questin Plots
ggsave(TrueObsCorbyMethodError_BoxPlot, file = "True_Obs_Cor_by_Method_Error_BoxPlot.png", height = 8, width = 8)
ggsave(TrueObsCorLoessPlot, file = "True_Obs_Cor_LoessPlot.png", height = 12, width = 12)

ggsave(TrueObsCorbyMethodErrorPlot,file = "TrueObsL2_Correlation_Plot.png",height = 7, width = 9)
ggsave(OutcomesbyMethodErrorPlot,file = "OutcomesbyMethodErrorPlot.png",height = 12, width = 12)
ggsave(OutcomesbyMethodICCPlot,file = "OutcomesbyMethodICCPlot.png",height = 12, width = 12)

## Export Omega Table
OmegaTable <- openxlsx::createWorkbook()
openxlsx::addWorksheet(OmegaTable,sheetName = "Omegas")
openxlsx::writeDataTable(OmegaTable, "Omegas", TheOmegas_Wide, rowNames = FALSE, colNames = TRUE,
                         tableStyle = "none", withFilter = FALSE, keepNA = FALSE)
openxlsx::saveWorkbook(OmegaTable, file = "Sim1_OmegaTable.xlsx", overwrite = TRUE)

# write.csv(DGSbyMethod, file = "SamplebyMethod.csv", row.names = FALSE)

