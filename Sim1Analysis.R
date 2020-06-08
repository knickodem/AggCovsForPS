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
load("Saved Sim1 Results/Sim1_ThousandRepsUpdate.RData")

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
  TheOmegas_Wide %>% select(Factor, !!ot) %>% arrange(desc(!!ot))
  
}

#### Shortcut for line and boxplots ####
Plot_Shortcut <- function(dat, xv, yv, gpv, linebox = c("line","box", "linebox")){
  
  if(linebox == "line"){
    
    deplot <- ggplot(dat, aes_string(x = xv, y = yv, group = gpv, color = gpv)) +
      geom_point(size = 3) + geom_smooth(method = "loess", se = FALSE, size = 2) +
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

#### Getting summary descriptives on a continuous variable ####
## unquoted grouping variables should be specified in ...
Get_Descriptives <- function(data, ContinuousVariable, ..., digits = 5, AllContinuous = TRUE){
  
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
DGSbyCondUpdate <- DGSbyCondUpdate %>%
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
  str_replace("method","Procedure")

## Renamed levels
ConditionRename <- ConditionVars
names(ConditionRename) <- ConditionFactor

#### Recoding factor levels ####

## ICC
ICCfac <- c(`.20` = "0.2", `.10` = "0.1", `.05` = "0.05")

## Correlations
Corrfac <- c(`Level 1` = "L1_Cor", `True Level 2` = "True_L2_Cor", `Aggregated Level 2` = "Obs_L2_Cor")


##################################################################


###########################################################################
#### Variation in outcomes explained by condition (partial-$\Omega^2$) ####

#### In timing ####
TheTimes <- lm(as.formula(paste0("Time ~ ", paste(ConditionVars, collapse = " + "))),
               data = DataGenTLogDF)
TimesAnova <- anova(TheTimes) # I don't actually use this
TimesOmega <- sjstats::omega_sq(TheTimes, partial = TRUE) %>%
  mutate(Characteristic = "Run Time")

#### In Sample Characteristics, then binding timing ####
## Regressions
Thelms <- map(DGSnameswCnvg, ~lm(as.formula(paste(.x, "~", paste(ConditionVars, collapse = " + "))),
                           data = DGSbyCondUpdate)) %>%
  set_names(DGSnameswCnvg)
TheAnovas <- map(Thelms,~anova(.x)) # I don't actually use this

## Effect size
#Note: Albers & Lakens (2017, p.9) says its okay that there are negative numbers even if it is theoretically impossible
TheOmegas <- map_dfr(Thelms, ~sjstats::omega_sq(.x, partial = TRUE), .id = "Characteristic") %>%
  bind_rows(TimesOmega) %>%
  mutate(term = factor(term, levels = ConditionVars),
         Factor = fct_recode(term, !!!ConditionRename)) %>%
  mutate_if(is.numeric, ~round(., 2))

## Wide format; for final paper, bold values > .14, which are considered large
# Table is exported at end of code
TheOmegas_Wide <- TheOmegas %>%
  select(Factor, Characteristic, partial.omegasq) %>%
  spread(Characteristic, partial.omegasq) %>%
  select(Factor, one_of(DGSnameswCnvg),`Run Time`) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_replace_all("\\.", "-"))

#### Model Checking ####
## QQ and FittedvResid plots for each characteristic
Model_Check_Plots <- map(Thelms, ~ModelCheckPlots(.x, se = FALSE, color = "black", size = 1))

## vif for the models
# outcome characteristic doesn't matter, so only really need to run this for one model since the factors are consistent across models
VIFs <- car::vif(Thelms[[1]]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Factor") %>%
  mutate(VIF = `GVIF^(1/(2*Df))`^2) %>%       # for approximate interpretation: https://stats.stackexchange.com/questions/70679/which-variance-inflation-factor-should-i-be-using-textgvif-or-textgvif
  mutate_if(is.numeric, ~round(., digits = 3))

######################################################


#############################################
#### RQ1: True and Obs L2 Corr by Method ####

Arrange_Omegas(`True-Obs L2 Cor`) #everything but nclusters

## Box Plot
TrueObsCorbyMethodError_BoxPlot <- DGSbyCondUpdate %>%
  mutate(method = str_replace(method, "\\.", "-")) %>%
  Plot_Shortcut(xv = "aggvar", yv = "True.Obs_L2_Cor", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "True and Aggregated L2 Correlation", limits = c(0, 1), breaks = seq(0, 1, .2)) +
  scale_x_discrete(name = "Error SD (or Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2"))
  # theme(legend.justification = c(0, 0), legend.position = c(.1, .1))

## Descriptives
# Note: Create tables only if needed to report specific values in text, but not needed for display
psych::describeBy(DGSbyCondUpdate$True.Obs_L2_Cor, list(DGSbyCondUpdate$method, DGSbyCondUpdate$aggvar), mat = TRUE, digits = 3)

#########################################################

#########################################################
#### RQ2: True and Obs L2 Corr by Contextual Factors ####

psych::describeBy(DGSbyCondUpdate$True.Obs_L2_Cor, list(DGSbyCondUpdate$icc), mat = TRUE, digits = 3)
psych::describeBy(DGSbyCondUpdate$True.Obs_L2_Cor, list(DGSbyCondUpdate$nsub), mat = TRUE, digits = 3)

## Loess plot including all manipulated factors except number of clusters
# Note: doing a boxplot isn't helpful because each box would only summarize 4 data points
# I don't end up using this
# TrueObsCor_LoessPlot <- DGSbyCondUpdate %>%
#   mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.)) %>%
#   rename(Subjects = nsub) %>%
#   Plot_Shortcut(xv = "aggvar", yv = "True.Obs_L2_Cor", gpv = "method", linebox = "line") +
#   scale_y_continuous(name = "True and Aggregated L2 Correlation", limits = c(0, 1), breaks = seq(0, 1, .2)) +
#   scale_x_discrete(name = "Error SD (or Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
#   scale_color_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
#   facet_grid(ICC ~ Subjects, labeller = label_both) +
#   theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "top")

#########################################################

################################################
#### RQ3: Sample characteristic differences ####

# variances of and correlations between L1 covariates, true L2 covariates, and observed L2 covariates;
# proportion of subjects in treatment and control groups; and variance and ICC(1) of the subject-level outcome

## Sample Characteristic Values in Long Format
DGSbyCond_Long <- DGSbyCondUpdate %>%
  gather(Characteristic, Value,one_of(DGSnameswCnvg))

#### ICCs ####
Arrange_Omegas(`Y ICC1`)           # 100% ICC and Method, also subjects and clusters ; Expected to be equal to ICC condition
Arrange_Omegas(`X ICC1`)           # 100% ICC, then clusters                         ; Expected to be equal to ICC condition
Arrange_Omegas(`X ICC2`)           # 100% ICC and Subjects, also clusters            ; Expected to vary based on ICC and nsubjects conditions

## Table for ICC2
# Expected values
ExpectedICC2 <- crossing(nsub = c(20, 60, 100),           # number of subjects per cluster
                         icc = c(.05, .1, .2))  %>%         # intraclass correlation [ICC(1)] for outcome, Y, and reflective L2 aggregations
  mutate(Value = (nsub*icc) / (1 + (nsub - 1)*icc),
         nclus = "Expected",
         nsub = factor(nsub),
         icc = factor(icc))
# Final table
ICC2Table <- DGSbyCondUpdate %>%
  group_by(icc, nsub) %>%
  summarize(Median = median(X_ICC2),
            Min = min(X_ICC2),
            Max = max(X_ICC2)) %>% ungroup() %>%
  inner_join(ExpectedICC2, by = c("icc", "nsub")) %>%
  mutate_if(is.numeric, ~round(.,2)) %>%
  unite(col = "Range",Min, Max, sep = " - ") %>%
  select(ICC1 = icc, Subjects = nsub, Expected = Value, Median, Range)

## Y_ICC1
YICC1_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Y_ICC1")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Subjects = fct_rev(nsub),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2"), guide = "none") +
  scale_y_continuous(name = "ICC(1) of Y", limits = c(0, .8), breaks = seq(0, .8, .1)) +
  xlab("ICC(1) Condition") +
  theme(strip.background = element_rect(colour = "black", fill = "white")) #legend.justification = c(0, 1), legend.position = c(.05, .95)

psych::describeBy(DGSbyCondUpdate$Y_ICC1, list(DGSbyCondUpdate$method), mat = TRUE, digits = 3)

#### Variances ####
Arrange_Omegas(`Y Variance`)       # 100% Method and ICC, also subjects and clusters ; Expected to vary based on ICC (and clusters?)
Arrange_Omegas(`L1 Variance`)      # Clusters and ICC                                ; Expected to vary by ICC, subjects and clusters?
Arrange_Omegas(`True L2 Variance`) # ~100% ICC, Method, and Subjects, also clusters  ; Expected to vary by ICC and clusters
Arrange_Omegas(`Obs L2 Variance`)  # ~100% aggvar, Method, also ICC and subjects     ; Expected to vary by ICC and clusters

## Descriptives
psych::describe(DGSbyCondUpdate$Y_Variance)
psych::describeBy(DGSbyCondUpdate$True_L2_Variance, list(DGSbyCondUpdate$icc), mat = TRUE, digits = 3) %>% View()

## Y
YVariance_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Y_Variance")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         Subjects = fct_rev(nsub),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  scale_y_continuous(name = "Variance of Y", limits = c(5, 30), breaks = seq(5, 30, 5)) +
  xlab("ICC(1) Condition") +
  theme(strip.background = element_rect(colour = "black", fill = "white"))


# ## True L2 variance
TrueL2Var_Box <- DGSbyCondUpdate %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac)) %>%
  rename(Subjects = nsub) %>%
  Plot_Shortcut(xv = "Subjects", yv = "True_L2_Variance", gpv = "method", linebox = "box") +
  scale_color_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(. ~ ICC, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "none")


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

# ## Loess plot including all manipulated factors except number of clusters
# ObsL2Var_LoessPlot <- DGSbyCondUpdate %>%
#   mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.)) %>%
#   rename(Subjects = nsub) %>%
#   Plot_Shortcut(xv = "aggvar", yv = "Obs_L2_Variance", gpv = "method", linebox = "line") +
#   scale_y_continuous(name = "Observed L2 Variance", limits = c(0, 1.23), breaks = seq(0, 1.2, .2)) +
#   scale_x_discrete(name = "Error SD or (Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
#   scale_color_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
#   facet_grid(ICC ~ Subjects, labeller = label_both) +
#   theme(strip.background = element_rect(colour = "black", fill = "white"), legend.position = "none")



#### Correlations ####
Arrange_Omegas(`L1 Cor`)           # 100% Methods, then interact w/ ICC and clusters ; Expected to be .2 across all conditions
Arrange_Omegas(`True L2 Cor`)      # Methods, subjects and ICC                       ; Expected to be .2 across all conditions
Arrange_Omegas(`Obs L2 Cor`)       # Method, aggvar, interact w/ ICC and subjects    ; Expected to be lower than .2 given the introduction of error

## Descriptives
psych::describeBy(DGSbyCondUpdate$L1_Cor, list(DGSbyCondUpdate$method), mat = TRUE, digits = 2)
psych::describeBy(DGSbyCondUpdate$True_L2_Cor, list(DGSbyCondUpdate$method), mat = TRUE, digits = 2)


## Prep
CorrsPlotPrep <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "_Cor") & !str_detect(Characteristic, "\\.")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac) %>% fct_rev(.),
         Characteristic = as_factor(Characteristic) %>% fct_recode(!!!Corrfac),
         method = str_replace(method, "\\.", "-"))

## Plot
Corrs_LoessPlot <- CorrsPlotPrep %>%
  Plot_Shortcut(xv = "aggvar", yv = "Value", gpv = "method", linebox = "line") +
  scale_y_continuous(name = "Correlation") +
  scale_x_discrete(name = "Error SD (or Sampling Ratio)", labels = c(".1 (.9)",".3 (.7)", ".5 (.5)", "1 (.3)")) +
  # scale_fill_manual(name = "Method", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  scale_color_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(ICC ~ Characteristic, labeller = labeller(ICC = label_both, Characteristic = label_value)) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))


#### Proportion in Treatment group ####
Arrange_Omegas(`Prop Z1`)          # ICC, Subjects, Method

## Descriptives
psych::describeBy(DGSbyCondUpdate$Prop_Z1, list(DGSbyCondUpdate$method), mat = TRUE, digits = 4) %>% View()  
median(DGSbyCondUpdate$Prop_Z1)

PropTreatBoxPlot  <- DGSbyCondUpdate %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Subjects = nsub) %>%
  ggplot(aes(x = Subjects, y = Prop_Z1, fill = method)) +
  geom_boxplot(notch = FALSE) +
  scale_y_continuous(name = "Proportion in Treatment Group", limits = c(.29,.38), breaks = seq(.30, .38, .02)) +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  theme_bw(base_size = 18) +
  facet_grid(. ~ ICC, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

#######################################################

##################################################
#### RQ4: bias and precision of PS estimation ####
# scale_fill_manual(name = "ICC(1)", values = c("#b2e2e2","#66c2a4","#238b45")) + #BlueGreen ; Grey scale: c("#cccccc","#969696","#525252")

Arrange_Omegas(`Logit Bias`) # mainly clusters, some ICC and error mag
Arrange_Omegas(`Logit MAE`) # Likewise
Arrange_Omegas(`Logit RMSE`) # Likewise
Arrange_Omegas(`PS Bias`) # icc, then subjects and interactions
Arrange_Omegas(`PS MAE`) # icc, aggvar, clusters, method in that order
Arrange_Omegas(`PS RMSE`) # Clusters, aggvar, icc, method in that order

#### Convergence ####
Arrange_Omegas(Convergence)

## Descriptives
psych::describeBy(DGSbyCondUpdate$Convergence, list(DGSbyCondUpdate$nclus), mat = TRUE, digits = 3) #, DGSbyCondUpdate$method

ConvergedPlot <- DGSbyCondUpdate %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Convergence", gpv = "method", linebox = "box") +
  ylab("Convergence Rate") +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2"))
  # facet_grid(. ~ Clusters, labeller = label_both) +
  # theme(strip.background = element_rect(colour = "black", fill = "white"))


#### Bias ####
## Descriptives
psych::describeBy(DGSbyCondUpdate$Logit_Bias, list(DGSbyCondUpdate$nclus, DGSbyCondUpdate$method), mat = TRUE, digits = 3) %>% View() #, DGSbyCondUpdate$method
  
## Bias Box Plot
# PS
BiasPS_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_Bias")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "Bias of PS", limits = c(-.045, 0), breaks = seq(0,-.04, -.01)) +
  scale_fill_manual(name = "Procedure", values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(. ~ Clusters, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_blank(), legend.position = "top")
# Logit
BiasLogit_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Logit_Bias")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "ICC", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "Bias of Logit", breaks = c(0, -4.0e+12, -8.0e+12, -1.2e+13), labels = c("0", "-4 trillion", "-8 trillion", "-12 trillion")) +
  xlab("ICC(1) Condition") +
  scale_fill_manual(guide = FALSE, values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  facet_grid(. ~ Clusters, labeller = label_both) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

# Juxtaposition PS and Logit
Bias_BoxPlot <- cowplot::plot_grid(BiasPS_BoxPlot, BiasLogit_BoxPlot,
                                   labels = "AUTO", align = "v", ncol = 1, rel_widths = c(.9, 1))



#### Bias, MAE, and RMSE Table ####
## MAE and RMSE Median (Range) by ICC and Cluster
# psych::describeBy(DGSbyCondUpdate$PS_RMSE, list(DGSbyCondUpdate$icc, DGSbyCondUpdate$nclus), mat = TRUE, digits = 3) %>% View()
# PSLogitDescrips <- DGSbyCondUpdate %>%
#   select(ICC = icc, Clusters = nclus, contains("RMSE"), contains("MAE")) %>%
#   gather(Temp, Value, -ICC, -Clusters) %>%
#   separate(Temp, into = c("Measure", "Error"), sep = "_") %>%
#   mutate(Clusters = fct_recode(Clusters, '8100' = "100")) %>%
#   unite(ICCclus, ICC, Clusters, sep = "_") %>%
#   group_by(ICCclus, Error, Measure, Error) %>%
#   summarize(Median = median(Value) %>% round(., 2) %>% format(., scientific = FALSE),
#             Range = paste0("(", format(round(min(Value), 2), scientific = FALSE), " - ", format(round(max(Value), 2), scientific = FALSE), ")")) %>%
#   gather(MdnRng, Value, Median, Range) %>%
#   spread(ICCclus, Value) %>%
#   arrange(Error, Measure)

## Logit Bias, MAE, and RMSE Median (Range) by generation procedure for 100 clusters
Logit100ClusterDescrips <- DGSbyCondUpdate %>%
  filter(nclus == "100") %>%
  select(Procedure = method, Logit_Bias, Logit_MAE, Logit_RMSE) %>%
  gather(Temp, Value, -Procedure) %>%
  mutate(Procedure = str_replace(Procedure, "\\.", "-"),
         Temp = str_remove(Temp, "Logit_")) %>%
  group_by(Procedure, Temp) %>%
  summarize(Median = paste0(format(round(median(Value), 2), scientific = FALSE), " (", format(round(min(Value), 2), scientific = FALSE), " - ", format(round(max(Value), 2), scientific = FALSE), ")")) %>%
  spread(Temp, Median)

#### MAE and RMSE Plots ####
## MAE Box Plots
# PS
MAEPS_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_MAE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "MAE of PS", limits = c(0, .33), breaks = seq(0, .30, .05)) +
  xlab("Clusters") +
  scale_fill_manual(guide = FALSE, values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_blank())
# Logit
MAELogit_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Logit_MAE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "MAE of Logit", limits = c(0, 3.0e+13), breaks = c(0, 1.0e+13, 2.0e+13, 3.0e+13), labels = c("0", "10 trillion", "20 trillion", "30 trillion")) +
    scale_fill_manual(guide = FALSE, values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  # facet_grid(. ~ Clusters) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

## RMSE Box Plot
# PS
RMSEPS_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "PS_RMSE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "RMSE of PS", limits = c(0, .33), breaks = seq(0, .30, .05)) +
  scale_fill_manual(guide = FALSE, values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        axis.title.x = element_blank())
# Logit
RMSELogit_BoxPlot <- DGSbyCond_Long %>%
  filter(str_detect(Characteristic, "Logit_RMSE")) %>%
  mutate(ICC = fct_recode(icc, !!!ICCfac),
         method = str_replace(method, "\\.", "-")) %>%
  rename(Clusters = nclus) %>%
  Plot_Shortcut(xv = "Clusters", yv = "Value", gpv = "method", linebox = "box") +
  scale_y_continuous(name = "RMSE of Logit", limits = c(0, 3.0e+13), breaks = c(0, 1.0e+13, 2.0e+13, 3.0e+13), labels = c("0", "10 trillion", "20 trillion", "30 trillion")) +
  scale_fill_manual(guide = FALSE, values = c("#e66101","#fdb863","#5e3c99","#b2abd2")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

# Get Legend
MAERMSE_Legend <- cowplot::get_legend(
  RMSELogit_BoxPlot + 
    guides(fill = guide_legend(nrow = 1, title = "Procedure")) +
    theme(legend.position = "top")
)
# Juxtaposition PS and Logit
MAERMSE_BoxPlot <- cowplot::plot_grid(MAEPS_BoxPlot, RMSEPS_BoxPlot, MAELogit_BoxPlot, RMSELogit_BoxPlot,
                                   labels = "AUTO", align = "hv", ncol = 2)
MAERMSE_BoxPlot <- cowplot::plot_grid(MAERMSE_Legend, MAERMSE_BoxPlot, ncol = 1, rel_heights = c(.1, 1))

###################################################################

#############################################################
#### Saving/Exporting Summary Objects, Plots, and Tables ####

#### Saving summary objects ####
save(TheTimes, Thelms, TheOmegas, TheOmegas_Wide,
     ConditionRename, ICCfac, Model_Check_Plots, VIFs, 
     ICC2Table, TrueL2Var_Box, 
     TrueObsCorbyMethodError_BoxPlot, Corrs_LoessPlot,
     YICC1_BoxPlot, YVariance_BoxPlot, PropTreatBoxPlot,
     BiasPS_BoxPlot, BiasLogit_BoxPlot, Bias_BoxPlot,
     MAEPS_BoxPlot, MAELogit_BoxPlot, RMSEPS_BoxPlot, RMSELogit_BoxPlot, MAERMSE_BoxPlot,
     file = "Sim1_Thousand_SummaryResultsUpdate.RData")

#### Exporting Model Check Plots ####
map2(.x = Model_Check_Plots, .y = names(Model_Check_Plots),
     ~ggsave(.x, file = paste0("Saved Sim1 Results/Model_Check_Plot_", .y, ".png"), height = 5, width = 10))


#### Exporting Research Question Plots ####
## RQ1 and 2
ggsave(TrueObsCorbyMethodError_BoxPlot, file = "True_Obs_Cor_by_Method_Error_BoxPlot.png", height = 6, width = 10)

## RQ3
# Juxtaposition of ICC(1) and Variance of Y
ICCVarY_BoxPlot <- cowplot::plot_grid(YICC1_BoxPlot, YVariance_BoxPlot,
                                      labels = "AUTO", align = "h", ncol = 2, rel_widths = c(.7, 1))
ggsave(ICCVarY_BoxPlot, file = "Juxtapose_ICC_Var_Y_BoxPlot.png", width = 12, height = 6)

# Correlations
ggsave(Corrs_LoessPlot, file = "L1L2Corrs_LoessPlot.png", height = 8, width = 14)

# Pr(Z = 1)
ggsave(PropTreatBoxPlot, file = "PropTreat_BoxPlot.png", height = 6, width = 12)

## RQ4
# Bias
ggsave(Bias_BoxPlot, file = "Bias_BoxPlot.png", width = 12, height = 10)

# MAE and RMSE
ggsave(MAERMSE_BoxPlot, file = "MAERMSE_BoxPlot.png", width = 10, height = 8)

## Export Omega Table
openxlsx::write.xlsx(TheSim2Omegas_Wide, file = "Sim1_OmegaTable.xlsx")


