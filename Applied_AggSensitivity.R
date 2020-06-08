####################################################
#                                                  #
#          SROs and Student Outcomes               #
#  Use of Aggregated Covariates Sensitivity Check  #
#                                                  #
####################################################

## Loading packages, functions, and cleaned dataset
source("Applied_PF.R")
load("Saved Applied Results/CompleteData.RData")

## Redefining School Covariates
SchoolCovNames_Sen <- SchoolCovNames[!(str_detect(SchoolCovNames, "_Perc"))]


## Recoding covariate names for figures and tables
StudentCovsRecode <- c(FRL = "Free/Reduced Priced Lunch", College = "College Aspiration",
                       Gambling = "Gambled", Vandalism = "Vandalized", Theft = "Stole", Diet = "Eat Produce",
                       Family.Community_Support_Center = "Family/Community Support", Trauma1 = "Experienced Trauma")
SchoolCovsRecode_Sen <- c(Total_Students10_Scaled = "Total Students", Fifth_L2 = "Elementary", Eighth_L2 = "Middle School",
                      HS_L2= "High School", Charter_Magnet = "Charter/Magnet", TC = "Urban",
                      p16.FRL = "% Free/Reduced Priced Lunch", p16.SPED = "% Special Education", p16.ELL = "% English Language Learner",
                      Proficient_Read = "% Proficient - Reading", Proficient_Math = "% Proficient - Math",
                      Prop_FTE_Absent = "Teacher Absences", Student_Teacher_Ratio_Scaled = "Student/Teacher Ratio", Guards = "Security Guard")

OutcomeFacOrder <- PrintableCovariatesShortcut(OutcomeVars)
StudentFacOrder <- recode(StudentCovNames, !!!StudentCovsRecode) %>%
  PrintableCovariatesShortcut()
SchoolFacOrder_Sen <- recode(SchoolCovNames_Sen, !!!SchoolCovsRecode_Sen) %>%
  PrintableCovariatesShortcut()
SampleFacOrder <- c("Full","SL","CL")

###################################
#### Propensity Score Modeling ####

## Variables excluded from the analysis due to redundancy
ExcludeCovs <- c("White","p16.White",    # American_Indian, #p16.American_Indian?   # Redundant with other race/ethnic variables
                 "Eleven","Ninth","Fifth_L2","HS_L2")  # Makes interpretation middle school vs high school even though L2 is not mutually exclusive

#### Cluster Level Treatment ####
## SRO
tic()
CLPSModel_Sen <- f.build("SRO",SchoolCovNames_Sen[!(SchoolCovNames_Sen %in% ExcludeCovs)])
CLSRO.Mod_Sen <- glm(CLPSModel_Sen, family = binomial("logit"),
                 data = CompleteSamp16L2Only)
toc()

tic()
SLPSModel_Sen <- f.build("SRO",c(StudentCovNames[!(StudentCovNames %in% ExcludeCovs)],
                             SchoolCovNames_Sen[!(SchoolCovNames_Sen %in% ExcludeCovs)]))
SLSRO.Mod_Sen <- glm(as.formula(SLPSModel_Sen), family = binomial("logit"),
                 data = FinalSample)
toc() #13.2 sec
# summary(SLSRO.Mod)
# car::vif(SLSRO.Mod) # some high numbers for L2 vars (same as CL model)

#### Adding Propensity Scores to dataset ####
CompleteSamp16L2OnlywPS_Sen <- CompleteSamp16L2Only %>%
  mutate(Logit_CL_Sen = predict(CLSRO.Mod_Sen),
         PS_CL_Sen = fitted(CLSRO.Mod_Sen))
FinalSamp16wPS_Sen <- FinalSample %>%
  mutate(Logit_SL_Sen = predict(SLSRO.Mod_Sen),
         PS_SL_Sen = fitted(SLSRO.Mod_Sen)) %>%
  left_join(CompleteSamp16L2OnlywPS_Sen %>% select(level2, Logit_CL_Sen, PS_CL_Sen), by = "level2")

############################
#### Conditioning on PS ####

#### Running Matching Algorithm ####
## Cluster Level - matching schools
Match.CL_Sen <- matchit(as.formula(CLPSModel_Sen), data = CompleteSamp16L2OnlywPS_Sen,
                    method = "nearest", replace = FALSE, caliper = .2, distance = CompleteSamp16L2OnlywPS_Sen$Logit_CL_Sen)

## Subject Level - matching students
tic()
Match.SL_Sen <- matchit(as.formula(SLPSModel_Sen), data = FinalSamp16wPS_Sen,
                    method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS_Sen$Logit_SL_Sen)
toc() # 756 - 1020 sec

#### Adding Indicator to dataset ####
# Cluster level
CompleteSamp16L2OnlywPS_Sen <- CompleteSamp16L2OnlywPS_Sen %>%
  mutate(CL_match_Sen = Match.CL_Sen$weights,)
# Subject level, then adding cluster info
FinalSamp16wPS_Sen <- FinalSamp16wPS_Sen %>%
  mutate(SL_match_Sen = Match.SL_Sen$weights) %>%
  left_join(CompleteSamp16L2OnlywPS_Sen %>% select(level2, CL_match_Sen), by = "level2")

#### Saving Matches ####
save(StudentCovNames, StudentCovsRecode, StudentFacOrder,            # Covariate Names
     SchoolCovNames_Sen, SchoolCovsRecode_Sen, SchoolFacOrder_Sen,
     ExcludeCovs, OutcomeVars, OutcomeFacOrder,                      # Outcome names
     CLPSModel_Sen, SLPSModel_Sen, CLSRO.Mod_Sen, SLSRO.Mod_Sen, PSDistributionPlot_Sen, # PS model
     Match.CL_Sen, Match.SL_Sen,                       # Conditioning methods
     FinalSamp16wPS_Sen, CompleteSamp16L2OnlywPS_Sen,                        # Datasets
     file = "Saved Applied Results/PSModels_and_Matches_Sen.RData")

##########################################

###########################
#### Assessing Balance ####

#### Calculating Covariate Descriptives and Balance ####
# Calculating standardized mean difference for each of the conditioning methods using cobalt
CovBalanceObjects_Sen <- map(.x = c("CL_match_Sen","SL_match_Sen"),
                         ~bal.tab(formula = f.build("SRO",c(StudentCovNames, SchoolCovNames_Sen)),
                                  data = FinalSamp16wPS_Sen,
                                  continuous = "std", binary = "std", s.d.denom = "pooled",
                                  abs = TRUE, un = TRUE,
                                  disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = TRUE,
                                  weights = FinalSamp16wPS_Sen[[.x]], method = "matching")) %>%
  set_names(c("CL", "SL"))

## Cleaning up table for presentation
CovBalanceTable_Sen <- map2(.x = CovBalanceObjects_Sen, .y = names(CovBalanceObjects_Sen),
                        ~Get_BalanceTable(.x, UnAdj = "Adj", renameUnAdj = .y) %>%
                          select(Covariate, starts_with("Diff"), starts_with("V.Ratio"))) %>%
  reduce(left_join, by = "Covariate") %>%
  left_join(Get_BalanceTable(CovBalanceObjects_Sen$SL, UnAdj = "Un", renameUnAdj = "Baseline") %>%
              select(Covariate, starts_with("Diff"), starts_with("V.Ratio")), by = "Covariate") %>%
  mutate(Covariate = recode(Covariate, !!!SchoolCovsRecode_Sen) %>%
           recode(., !!!StudentCovsRecode) %>%
           PrintableCovariatesShortcut()) %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  select(Covariate,Diff.Baseline, Diff.SL, Diff.CL, V.Ratio.Baseline, V.Ratio.SL, V.Ratio.CL) %>%
  mutate(Covariate = factor(Covariate, levels = c(StudentFacOrder, SchoolFacOrder_Sen)))

## Checking that unadjusted values were equivalent
map_lgl(c("M.0.Un","SD.0.Un","M.1.Un","SD.1.Un","Diff.Un","V.Ratio.Un"),
        ~identical(CovBalanceObjects_Sen$SL$Balance[[.x]], CovBalanceObjects_Sen$CL$Balance[[.x]]))   # identical

CovBalanceTableLong_Sen <- CovBalanceTable_Sen %>%
  select(Covariate, starts_with("Diff"), starts_with("V.Ratio")) %>%
  gather(Temp, Value, -Covariate) %>%
  mutate(Temp = str_replace(Temp, "Diff", "ASD") %>%
           str_replace(., "V.Ratio", "VR")) %>%
  separate(Temp, into = c("Measure", "Sample"), sep = "\\.") %>%
  spread(Measure, Value)

#### Number of Matches for each Method ####
# Students
StudentSampleSizes_Sen <- FinalSamp16wPS_Sen %>%
  group_by(SRO) %>%
  summarize(SL_match_Sen = sum(SL_match_Sen),
            CL_match_Sen = sum(CL_match_Sen),
            Baseline = n()) %>%
  gather(Sample, Students, -SRO) %>%
  mutate(SRO = as.character(SRO))
# Schools
SchoolSampleSize_Sen <- map_dfr(c("CL_match_Sen", "SL_match_Sen"),
                            ~FinalSamp16wPS_Sen[FinalSamp16wPS_Sen[[.x]] != 0, ] %>%
                              group_by(SRO) %>%
                              summarize(Schools = length(unique(level2)),
                                        Sample = .x)) %>%
  mutate(SRO = as.character(SRO)) %>%
  bind_rows(as.data.frame(table(CompleteSamp16L2OnlywPS_Sen$SRO)) %>%
              mutate(Sample = "Baseline") %>% rename(SRO = Var1, Schools = Freq))
# Combined Table
SampleSizes_Sen <- left_join(StudentSampleSizes_Sen, SchoolSampleSize_Sen, by = c("SRO", "Sample")) %>%
  mutate(Sample = str_remove(Sample, "_match_Sen"),     #str_replace(Sample, "SL_weight", "Subject") %>% str_replace(., "CL_weight", "Cluster"),
         SRO = recode(SRO, `0` = "TNo", `1` = "SRO")) %>%
  gather(Temp, Value, Students, Schools) %>%
  group_by(Sample, Temp) %>%
  mutate(Total = sum(Value),
         Percent = round((Value / Total) * 100, 1),
         Keep = paste0(Value, " (", Percent, ")")) %>% ungroup() %>%
  unite(col = "SRO_Type", SRO, Temp) %>%
  select(Sample, SRO_Type, Keep) %>%
  spread(SRO_Type, Keep)
#### Balance Summaries ####
UnderThresholdTable_Sen <- CovBalanceTableLong_Sen %>%
  group_by(Sample) %>%
  summarize(ASD_Count = sum(ASD < 0.10),
            ASD_Percent = round((ASD_Count / nrow(CovBalanceTable_Sen))*100, 1),
            VR_Count = sum(VR < 1.50, na.rm = TRUE),
            VR_Percent = round((VR_Count / sum(!is.na(CovBalanceTable_Sen$V.Ratio.Baseline)))*100, 1)) %>%
  ungroup() %>%
  mutate(ASD = paste0(ASD_Count, " (", ASD_Percent, ")"),
         VR = paste0(VR_Count, " (", VR_Percent, ")"))

SSUTTable_Sen <- left_join(SampleSizes_Sen, UnderThresholdTable_Sen, by = "Sample") %>%
  select(Sample, ends_with("Students"), ends_with("Schools"), ASD, VR) %>%
  mutate(Sample = factor(Sample, levels = SampleFacOrder)) %>%
  arrange(Sample)

## Ordering variables based on largest difference in the original sample for plotting purposes
CovBalVarOrder_Sen <- CovBalanceTable_Sen %>%
  mutate(Covariate = fct_reorder(as_factor(Covariate), Diff.Baseline , .desc = FALSE))
CovBalVarOrder_Sen <- levels(CovBalVarOrder_Sen$Covariate)


## Balance Plot for all Samples
# Need to clean up covariate names
CovBalAllPlot_Sen <- CovBalanceTableLong_Sen %>%
  mutate(Covariate = factor(Covariate, levels = CovBalVarOrder_Sen)) %>%
  ggplot(aes(x = Covariate, y = ASD, group = Sample)) +
  geom_hline(yintercept = .1, linetype = 1, color = "black", size = 2) +
  geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = 2) +
  geom_point(aes(shape = Sample, color = Sample), size = 5, position = position_dodge(1)) +
  ylab("Standardized Difference") + xlab("") +
  scale_shape_manual(values = c(19, 18, 17)) + # 15,7,3
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = rev(c("#e66101","#5e3c99","#4daf4a"))) +
  coord_flip() +
  theme_bw(base_size = 22) +
  theme(panel.grid.major = element_line(color = "gray87"),
        panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.justification = c(1, 0), legend.position = c(.9, .1))

ggsave(CovBalAllPlot_Sen, file = "Covariate_Balance_Plot_Sen.png", width = 15, height = 15)

### Outcome Descriptives and Effect Size Difference ####
## Provides the standardized difference without adjusting for covariates
OutcomeBalance_Sen <- map(.x = c("CL_match_Sen", "SL_match_Sen"),
                      ~bal.tab(formula = f.build("SRO",c(OutcomeVars)),
                               data = FinalSamp16wPS_Sen,
                               continuous = "std", binary = "std", s.d.denom = "pooled",
                               abs = FALSE, un = TRUE,
                               disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = FALSE,
                               weights = FinalSamp16wPS_Sen[[.x]], method = "matching")) %>%
  set_names(c("CL", "SL"))

## Cleaning up table for presentation
OutcomeBalanceTable_Sen <- Get_BalanceTable(OutcomeBalance_Sen$SL, UnAdj = "Un", renameUnAdj = "X") %>% mutate(Sample = "Baseline") %>%
  bind_rows(Get_BalanceTable(OutcomeBalance_Sen$SL, UnAdj = "Adj", renameUnAdj = "X") %>% mutate(Sample = "SL")) %>%
  bind_rows(Get_BalanceTable(OutcomeBalance_Sen$CL, UnAdj = "Adj", renameUnAdj = "X") %>% mutate(Sample = "CL")) %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  mutate(SRO.X = paste0(M.1.X, " (", SD.1.X, ")"),
         No_SRO.X = paste0(M.0.X, " (", SD.0.X, ")"),
         Covariate = factor(Covariate, levels = OutcomeVars)) %>%
  select(Covariate, Sample, SRO.X, No_SRO.X, ASD.X = Diff.X) %>%
  gather(Temp, Value, ends_with("X")) %>%
  arrange(Covariate) %>%
  mutate(Temp = str_replace(Temp, "X", as.character(Covariate)) %>% as_factor()) %>%
  select(-Covariate) %>%
  spread(Temp, Value) %>%
  mutate(Sample = factor(Sample, levels = SampleFacOrder))
#####################################################################


#### Saving data, models, and balance results ####
save(StudentCovNames, StudentCovsRecode, StudentFacOrder,            # Covariate Names
     SchoolCovNames_Sen, SchoolCovsRecode_Sen, SchoolFacOrder_Sen,
     ExcludeCovs, OutcomeVars, OutcomeFacOrder,                      # Outcome names
     FinalSamp16wPS_Sen, CompleteSamp16L2OnlywPS_Sen,
     CovBalanceObjects_Sen, CovBalanceTable_Sen, CovBalanceTableLong_Sen, OutcomeBalanceTable_Sen,
     SampleSizes_Sen, UnderThresholdTable_Sen, SSUTTable_Sen, CovBalVarOrder_Sen, CovBalAllPlot_Sen,
     file = "Saved Applied Results/PSbalance_Sen.RData")

############################################################
#####   Outcome HLM Model and Results for Applied Study ####

load("Saved Applied Results/PSbalance_Sen.RData")
rm(CovBalAllPlot_Sen, CovBalanceObjects_Sen)

#### Prepping Data ####

## Identifying variables to use in the outcome regression
HLMvars_Sen <- c(StudentCovNames, SchoolCovNames_Sen)

## Defining models
# UNCmod <- "~ 1 + (1|level2)"
# HLMmod <- paste0(" ~ 1 + ", paste(HLMvars[!(HLMvars %in% ExcludeCovs)], collapse = " + "), " + (1|level2)")
SROmod_Sen <- paste0(" ~ 1 + ", paste(HLMvars_Sen[!(HLMvars_Sen %in% ExcludeCovs)], collapse = " + "), " + SRO + (1|level2)")

## Dataset updates
FinalSamp16wPS_Sen$level2 <- as.character(FinalSamp16wPS_Sen$level2) # Converting school id from numeric to character
FinalSamp16wPS_Sen$Baseline_match_Sen <- 1                           # Indicator that all observations should be used in baseline models


################################
#### Running Outcome Models ####

#### Baseline ####
tic()
BaselineOutcomeMods_Sen <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod_Sen)),
                                                       data =  FinalSamp16wPS_Sen[FinalSamp16wPS_Sen[["Baseline_match_Sen"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 45 sec

map_dbl(BaselineOutcomeMods_Sen, ~getME(.x,"n")) # 98503

#### Subject Level ####
tic()
SLOutcomeMods_Sen <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod_Sen)),
                                                 data =  FinalSamp16wPS_Sen[FinalSamp16wPS_Sen[["SL_match_Sen"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 13 sec

map_dbl(SLOutcomeMods_Sen, ~getME(.x,"n")) # 32372; orig = 30464

#### Cluster Level ####
tic()
CLOutcomeMods_Sen <- map(.x = OutcomeVars, ~lmer(as.formula(paste(.x, SROmod_Sen)),
                                                 data =  FinalSamp16wPS_Sen[FinalSamp16wPS_Sen[["CL_match_Sen"]] != 0, ])) %>%
  set_names(OutcomeFacOrder)
toc() # ~ 12 sec

map_dbl(CLOutcomeMods_Sen, ~getME(.x,"n")) # 32238; orig = 29592

#######################################


############################
#### Extracting Results ####

#### Fixed Effects - for full SRO model only ####
## Baseline
BaselineFixedEff_Sen <- map_dfr(BaselineOutcomeMods_Sen, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
                              select(term, estimate, conf.low, conf.high), .id = "Outcome")

## Subject Level
SLFixedEff_Sen <- map_dfr(SLOutcomeMods_Sen, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
                        select(term, estimate, conf.low, conf.high), .id = "Outcome")

## Cluster Level
CLFixedEff_Sen <- map_dfr(CLOutcomeMods_Sen, ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE, conf.level = .95, conf.method = "Wald") %>%
                        select(term, estimate, conf.low, conf.high), .id = "Outcome")


#### R2 / f2 ####
tic(); BaselineR2_Sen <- map2_dfr(.x = BaselineOutcomeMods_Sen, .y = OutcomeVars, ~Get_f2(.x, .y, "Baseline_match_Sen")); toc() # ~ 42 sec
tic(); SLR2_Sen <- map2_dfr(.x = SLOutcomeMods_Sen, .y = OutcomeVars, ~Get_f2(.x, .y, "SL_match_Sen")); toc() # ~ 15 sec
tic(); CLR2_Sen <- map2_dfr(.x = CLOutcomeMods_Sen, .y = OutcomeVars, ~Get_f2(.x, .y, "CL_match_Sen")); toc() # ~ 13 sec

#### Final Tables ####
## SRO Effect
SROEffectTable_Sen <-  bind_rows(BaselineR2_Sen, SLR2_Sen, CLR2_Sen) %>%
  select(Sample, Outcome, f2) %>%
  mutate(Sample = str_remove(Sample, "_match_Sen"),
         Outcome = str_replace_all(Outcome, "\\.", "/") %>% str_replace_all("_", " "),
         f2 = round(f2, 2)) %>%
  left_join(bind_rows(mutate(BaselineFixedEff_Sen, Sample = "Baseline"),
                      mutate(SLFixedEff_Sen, Sample = "SL"),
                      mutate(CLFixedEff_Sen, Sample = "CL")) %>%
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

########################
#### Model Checking ####

#### VIFs ####
VIFs_Sen <- map2(.x = list(BaselineOutcomeMods_Sen[["GPA"]], SLOutcomeMods_Sen[["GPA"]], CLOutcomeMods_Sen[["GPA"]]),
             .y = c("Baseline", "SL", "CL"),
             ~vif.table(.x, .y, multilevel = TRUE)) %>%
  reduce(left_join, by = "Predictor")

#### Residual Plots ####
## Baseline
tic()
BaselineModelCheck_Sen <- map(BaselineOutcomeMods_Sen, .y = OutcomeVars,
                          ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                            cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = BaselineModelCheck_Sen, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_Sen_Baseline_", .y, ".png"),
             height = 10, width = 8))
toc() # ~1575 sec

## Subject Level
# Generating Plots
tic()
SLModelCheck_Sen <- map(SLOutcomeMods_Sen, .y = OutcomeVars,
                    ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                      cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = SLModelCheck_Sen, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_Sen_SL_", .y, ".png"),
             height = 10, width = 8))
toc() # ~217 sec

## Cluster Level
# Generating Plots
tic()
CLModelCheck_Sen <- map(CLOutcomeMods_Sen, .y = OutcomeVars,
                    ~ModelCheckPlots(.x, se = FALSE, color = "blue", size = 1, smooth_method = "loess") %>%
                      cowplot::plot_grid(plotlist = ., labels = "AUTO", ncol = 2))
# Exporting Plots
map2(.x = CLModelCheck_Sen, .y = OutcomeVars,
     ~ggsave(.x, file = paste0("Saved Applied Results/Model_Check_Plot_Sen_CL_", .y, ".png"),
             height = 10, width = 8))
toc() # ~214 sec

##################################################

#### Saving Models and Results ####
save(BaselineOutcomeMods_Sen, SLOutcomeMods_Sen, CLOutcomeMods_Sen,
     BaselineR2_Sen, SLR2_Sen, CLR2_Sen,
     BaselineFixedEff_Sen, SLFixedEff_Sen, CLFixedEff_Sen,
     SROEffectTable_Sen, OutcomeBalanceTable_Sen,
     VIFs_Sen, BaselineModelCheck_Sen, SLModelCheck_Sen, CLModelCheck_Sen,
     file = "Saved Applied Results/HLM_Models_and_Output_Sen.RData")

load("Saved Applied Results/HLM_Models_and_Output_Sen.RData")

