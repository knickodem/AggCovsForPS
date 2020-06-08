#####################################
#                                   #
#     SROs and Student Outcomes     #
#   Outcome HLM Model and Results   #
#                                   #
#####################################

## Loading packages, functions, and PS results
source("Applied_PF.R")
load("Saved Applied Results/PSbalance.RData")
rm(OutcomeBalance, CovBalAllPlot, CovBalanceObjects)

#### Prepping Data ####

## Identifying variables to use in the outcome regression
HLMvars <- c(StudentCovNames, SchoolCovNames)

## Defining models
# UNCmod <- "~ 1 + (1|level2)"
# HLMmod <- paste0(" ~ 1 + ", paste(HLMvars[!(HLMvars %in% ExcludeCovs)], collapse = " + "), " + (1|level2)")
SROmod <- paste0(" ~ 1 + ", paste(HLMvars[!(HLMvars %in% ExcludeCovs)], collapse = " + "), " + SRO + (1|level2)")

## Dataset updates
FinalSamp16wPS$level2 <- as.character(FinalSamp16wPS$level2) # Converting school id from numeric to character
FinalSamp16wPS$Baseline_match <- 1                           # Indicator that all observations should be used in baseline models


################################
#### Running Outcome Models ####

#### Baseline ####
tic()
BaselineOutcomeMods <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod)), data =  FinalSamp16wPS[FinalSamp16wPS[["Baseline_match"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 51 - 97 sec

map_dbl(BaselineOutcomeMods, ~getME(.x,"n")) # 98503

#### Subject Level ####
tic()
SLOutcomeMods <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod)), data =  FinalSamp16wPS[FinalSamp16wPS[["SL_match"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 15 - 23 sec

map_dbl(SLOutcomeMods, ~getME(.x,"n")) # 30464

#### Cluster Level ####
tic()
CLOutcomeMods <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod)), data =  FinalSamp16wPS[FinalSamp16wPS[["CL_match"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 15 - 24 sec

map_dbl(CLOutcomeMods, ~getME(.x,"n")) #29592

#######################################


############################
#### Extracting Results ####

#### Fixed Effects - for full SRO model only ####
## Baseline
BaselineFixedEff <- map_dfr(BaselineOutcomeMods, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
  select(term, estimate, conf.low, conf.high), .id = "Outcome")

## Subject Level
SLFixedEff <- map_dfr(SLOutcomeMods, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
                              select(term, estimate, conf.low, conf.high), .id = "Outcome")

## Cluster Level
CLFixedEff <- map_dfr(CLOutcomeMods, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
                        select(term, estimate, conf.low, conf.high), .id = "Outcome")


#### R2 / f2 ####
tic(); BaselineR2 <- map2_dfr(.x = BaselineOutcomeMods, .y = OutcomeVars, ~Get_f2(.x, .y, "Baseline_match")); toc() # ~56 sec
tic(); SLR2 <- map2_dfr(.x = SLOutcomeMods, .y = OutcomeVars, ~Get_f2(.x, .y, "SL_match")); toc() # ~ 17 sec
tic(); CLR2 <- map2_dfr(.x = CLOutcomeMods, .y = OutcomeVars, ~Get_f2(.x, .y, "CL_match")); toc() # ~ 16 sec

#### Final Tables ####
## SRO Effect
SROEffectTable <-  bind_rows(BaselineR2, SLR2, CLR2) %>%
  select(Sample, Outcome, f2) %>%
  mutate(Sample = str_remove(Sample, "_match"),
         Outcome = str_replace_all(Outcome, "\\.", "/") %>% str_replace_all("_", " "),
         f2 = round(f2, 2)) %>%
  left_join(bind_rows(mutate(BaselineFixedEff, Sample = "Baseline"),
                      mutate(SLFixedEff, Sample = "SL"),
                      mutate(CLFixedEff, Sample = "CL")) %>%
              filter(term == "SRO") %>%
              mutate_if(is.numeric, ~round(., 2)) %>%
              mutate(CI95 = paste0("[", conf.low, ", ", conf.high, "]")) %>%
              select(Sample, Outcome, estimate, CI95),
            by = c("Sample", "Outcome")) %>%
  gather(Measure, Value, estimate, CI95, f2) %>%
  unite("OM", Outcome, Measure, sep = "_") %>%
  mutate(OM = as_factor(OM),
         Sample = as_factor(Sample)) %>%
  spread(OM, Value) %>%
  select(Sample, starts_with("GPA"), starts_with("Commitment"), starts_with("Empowerment"), starts_with("Teacher"))

## All Fixed Effects
AllFE <- map2(.x = list(BaselineFixedEff, SLFixedEff, CLFixedEff), .y = c("Baseline", "SL", "CL"),
              ~mutate_if(.x, is.numeric, ~round(., 2)) %>%
                mutate(CI95 = paste0("[", conf.low, ", ", conf.high, "]")) %>%
                select(Outcome, Covariate = term, estimate, CI95) %>%
                rename(!!paste0(.y, "_Est") := estimate, !!paste0(.y, "_CI95") := CI95)) %>%
  reduce(left_join, by = c("Outcome", "Covariate")) %>%
  mutate(Covariate = recode(Covariate, !!!SchoolCovsRecode) %>%
           recode(., !!!StudentCovsRecode) %>%
           PrintableCovariatesShortcut() %>%
           factor(levels = c("(Intercept)", StudentFacOrder, SchoolFacOrder, "SRO")))
                

##################################################


########################
#### Model Checking ####

#### VIFs ####
VIFs <- map2(.x = list(BaselineOutcomeMods[["GPA"]], SLOutcomeMods[["GPA"]], CLOutcomeMods[["GPA"]]),
             .y = c("Baseline", "SL", "CL"),
             ~vif.table(.x, .y, multilevel = TRUE)) %>%
  reduce(left_join, by = "Predictor")


#### Residual Plots ####
## Baseline
tic()
BaselineModelCheck <- map(BaselineOutcomeMods, .y = OutcomeVars,
                    ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                      cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = BaselineModelCheck, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_Baseline_", .y, ".png"),
             height = 10, width = 8))
toc() # ~773 sec

## Subject Level
# Generating Plots
tic()
SLModelCheck <- map(SLOutcomeMods, .y = OutcomeVars,
                    ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                      cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = SLModelCheck, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_SL_", .y, ".png"),
             height = 10, width = 8))
toc() # ~94 sec

## Cluster Level
# Generating Plots
tic()
CLModelCheck <- map(CLOutcomeMods, .y = OutcomeVars,
                    ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                      cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = CLModelCheck, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_CL_", .y, ".png"),
             height = 10, width = 8))
toc() # ~82 sec

##################################################

#### Saving Models and Results ####
save(BaselineOutcomeMods, SLOutcomeMods, CLOutcomeMods,
     BaselineR2, SLR2, CLR2,
     BaselineFixedEff, SLFixedEff, CLFixedEff,
     SROEffectTable, AllFE,
     VIFs, BaselineModelCheck, SLModelCheck, CLModelCheck,
     file = "Saved Applied Results/HLM_Models_and_Output.RData")

