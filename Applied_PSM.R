#########################################
#                                       #
#       SROs and Student Outcomes       #
#       Propensity Score Matching       #
#                                       #
#########################################

## Loading packages, functions, and cleaned dataset
source("Applied_PF.R")
load("Saved Applied Results/CompleteData.RData")

## Recoding covariate names for figures and tables
StudentCovsRecode <- c(FRL = "Free/Reduced Priced Lunch", College = "College Aspiration",
                       Gambling = "Gambled", Vandalism = "Vandalized", Theft = "Stole", Diet = "Eat Produce",
                       Family.Community_Support_Center = "Family/Community Support", Trauma1 = "Experienced Trauma")
SchoolCovsRecode <- c(Total_Students10_Scaled = "Total Students", Fifth_L2 = "Elementary", Eighth_L2 = "Middle School",
                      HS_L2= "High School", Charter_Magnet = "Charter/Magnet", TC = "Urban",
                      p16.FRL = "% Free/Reduced Priced Lunch", p16.SPED = "% Special Education", p16.ELL = "% English Language Learner",
                      College_Perc = "% College Aspiration", Proficient_Read = "% Proficient - Reading", Proficient_Math = "% Proficient - Math",
                      Gambling_Perc = "% Gambled", Vandalism_Perc = "% Vandalized", Theft_Perc = "% Stole",
                      Trauma1_Perc = "% Experienced Trauma", Diet_Perc = "% Eat Produce",
                      Prop_FTE_Absent = "Teacher Absences", Student_Teacher_Ratio_Scaled = "Student/Teacher Ratio", Guards = "Security Guard",
                      Positive_Identity_Perc_Center = "School Positive Identity", Social_Competence_Perc_Center = "School Social Competence",
                      Family.Community_Support_Perc_Center = "School Family/Community Support")

OutcomeFacOrder <- PrintableCovariatesShortcut(OutcomeVars)
SampleFacOrder <- c("Full","SL","CL")
StudentFacOrder <- recode(StudentCovNames, !!!StudentCovsRecode) %>%
  PrintableCovariatesShortcut()
SchoolFacOrder <- recode(SchoolCovNames, !!!SchoolCovsRecode) %>%
  PrintableCovariatesShortcut()


###############################
#### Baseline Descriptives ####

#### Student Level ####
# Note: Getting the baseline descriptives in the original scale rather than the centered and standardized scales used latter
# Check if unadjusted Diff is the same for centered/scaled. If different, then use the latter to be comparable with the PS models
# Slightly, so I'll only keep the descriptives here
## All Covariates by SRO exposure
BaselineStudentDescrips <- bal.tab(formula = f.build("SRO", c(StudentPSVarNames, SchoolPSVarNames)),
                                   data = FinalSample,
                                   continuous = "std", binary = "std", s.d.denom = "pooled",
                                   abs = TRUE, un = TRUE,
                                   disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = FALSE) %>%
  Get_BalanceTable("Un", "Full") %>%
  mutate(Covariate = recode(Covariate, Positive_Identity_Perc = "School Positive Identity", Social_Competence_Perc = "School Social Competence",
                            Family.Community_Support_Perc = "School Family/Community Support", Student_Teacher_Ratio = "Student/Teacher Ratio") %>%
           recode(., !!!SchoolCovsRecode) %>%
           recode(., !!!StudentCovsRecode) %>%
           PrintableCovariatesShortcut() %>%
           str_remove("10") %>% as_factor()) %>%
  select(-Diff.Full, -V.Ratio.Full) %>%
  mutate_if(is.numeric, ~ifelse(Covariate == "Total Students", .*10, .) %>%
              round(., 2))
## Total + by SRO
BaselineStudentDescripsTotal <- FinalSample %>%
  select(studentnumber, level2, SRO, one_of(StudentPSVarNames), one_of(SchoolPSVarNames)) %>%
  gather(Covariate, Value, -studentnumber, -level2) %>%
  Get_Descriptives(Value, Covariate, digits = 6, AllContinuous = TRUE) %>%
  mutate(Covariate = recode(Covariate, Positive_Identity_Perc = "School Positive Identity", Social_Competence_Perc = "School Social Competence",
                            Family.Community_Support_Perc = "School Family/Community Support", Student_Teacher_Ratio = "Student/Teacher Ratio") %>%
           recode(., !!!SchoolCovsRecode) %>%
           recode(., !!!StudentCovsRecode) %>%
           PrintableCovariatesShortcut() %>%
           str_remove("10")) %>%
  mutate_if(is.numeric, ~ifelse(Covariate == "Total Students", .*10, .) %>%
              round(., 2)) %>%
  left_join(BaselineStudentDescrips, by = "Covariate") %>%
  mutate(Total = ifelse(Type == "Binary" | is.na(Type), as.character(Mean * 100), paste0(Mean, " (", SD, ")")),
         SRO = ifelse(Type == "Binary" | is.na(Type), as.character(M.1.Full * 100), paste0(M.1.Full, " (", SD.1.Full, ")")),
         No_SRO = ifelse(Type == "Binary" | is.na(Type), as.character(M.0.Full * 100), paste0(M.0.Full, " (", SD.0.Full, ")")),
         Covariate = factor(Covariate, levels = c("SRO", StudentFacOrder, SchoolFacOrder)))
  

#### School Level #### 
## Only school covariates by SRO exposure
BaselineSchoolDescrips <- bal.tab(formula = f.build("SRO", c(SchoolPSVarNames)),
                                        data = CompleteSamp16L2Only,
                                        continuous = "std", binary = "std", s.d.denom = "pooled",
                                        abs = TRUE, un = TRUE,
                                        disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = TRUE) %>%
  Get_BalanceTable("Un", "Baseline") %>%
  mutate(Covariate = recode(Covariate, Positive_Identity_Perc = "School Positive Identity", Social_Competence_Perc = "School Social Competence",
                            Family.Community_Support_Perc = "School Family/Community Support", Student_Teacher_Ratio = "Student/Teacher Ratio") %>%
           recode(., !!!SchoolCovsRecode) %>%
           PrintableCovariatesShortcut() %>%
           str_remove("10") %>% as_factor()) %>%
  mutate_at(vars(M.0.Baseline:SD.1.Baseline), ~ifelse(Covariate == "Total Students", .*10, .) %>%
              round(., 2))

################################################
  
###########################
#### Aggregate Quality ####

# Aggregate values included responses from students with missing data and from grade 5 students

#### ICCs and Sampling Ratio ####
## Intercept only multilevel model for each agg covariate
ICCMods <- map(str_remove_all(StudentCovNames, "_Center"), ~lmer(as.formula(paste0(.x," ~ 1 + (1|level2)")), data = FinalSample)) %>%
  set_names(str_remove_all(StudentCovNames, "_Center"))

## Calculate ICC(1) and ICC(2) for each agg covariate
ICC12 <- map_df(ICCMods, ~Get_ICC(.x, type = "1")) %>% mutate(Type = "ICC(1)") %>%     # Extracting ICC(1) from each model
  bind_rows(map_df(ICCMods, ~Get_ICC(.x, type = "2")) %>% mutate(Type = "ICC(2)")) %>% # Extracting ICC(2) from each model
  gather(Variable, Value, -Type)

## Calculating sampling ratio by school on each covariate, then averaging
SamplingRatio <- SL16no5 %>%                                         # Using dataset with missing data because
  select(level2, TotalEnrollment16, Total_Enrollment = mss_sample,   # mde reported size & mss sample size for each school
         one_of(paste0(str_remove_all(StudentCovNames, "_Center"),"_n"))) %>%                   # n's for each student covariate
  filter(level2 %in% FinalSample$level2) %>%                         # only keeping schools in FinalSample
  gather(Variable, n, Total_Enrollment, ends_with("_n")) %>%
  mutate(SR = n / TotalEnrollment16) %>%                             # calculates sampling ratio for each school
  group_by(Variable) %>%
  summarize(Value = mean(SR, na.rm = TRUE)) %>% ungroup() %>%        # average sampling ratio for each covariate
  mutate(Variable = str_remove(Variable, "_n"),
         Type = "SR")

#### Correlation between aggregate and true L2 ####

# School-level covariates do not include the replaced values and instead have original missingness
# The schools included, however, are only those in the Complete dataset

## Identifying covariates to compare
CompareOCRMDEMSS <- SL16no5 %>%
  select(level2,Grade_Levels,TotalEnrollment16,Total_Enrollment,mss_sample,
         p16.Female:p16.Homeless,Female_Perc:Homeless_Perc,-LGB_Perc,
         Skipped_Class_Perc:Out.School_Suspension_Perc, Bullying_Perc, Bullied_Perc,
         starts_with("Prop"), ends_with("Per_Student")) %>%
  filter(level2 %in% FinalSample$level2)                        # only keeping schools in FinalSample
  

AggTrueCorrs <- CompareOCRMDEMSS %>% select(Total_Enrollment = mss_sample, MDE = TotalEnrollment16) %>%           ## Total Students
  rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs") %>%
  bind_rows(select(CompareOCRMDEMSS, American_Indian_Perc, MDE = p16.American_Indian) %>%                         ## Race/Ethnic Proportions (MDE and MSS)
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%  # before combining Hmong and Somali into Asian and Black: .925,.546,.952,.794,.967,.251
  bind_rows(select(CompareOCRMDEMSS, Asian_Perc, MDE = p16.Asian) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Black_Perc, MDE = p16.Black) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, White_Perc, MDE = p16.White) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Latinx_Perc, MDE = p16.Latino) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Multiracial_Perc, MDE = p16.Multiracial) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Female_Perc, MDE = p16.Female) %>%                                           ## School Composition
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Special_Education_Perc, MDE = p16.SPED) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, FRL_Perc, MDE = p16.FRL) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Homeless_Perc, MDE = p16.Homeless) %>%
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Skipped_School_Perc, CRDC = Prop_Truant) %>%                                 ## Misbehavior and discipline
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%   # interpreted as % of students with [disciplinary action]
  bind_rows(select(CompareOCRMDEMSS, Skipped_Class_Perc, CRDC = Prop_Truant) %>%                                  # Truant is absent 15+ days per year whereas skipped class/school is once or more in last 30 days
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, In.School_Suspension_Perc, CRDC = Prop_InSchool_Suspensions) %>%       
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  bind_rows(select(CompareOCRMDEMSS, Out.School_Suspension_Perc, CRDC = Prop_OutSchool_Suspensions) %>% 
              rquery.cormat(.,type = "flatten", graph = FALSE, digits = 3, usena = "pairwise.complete.obs")) %>%
  mutate(row = str_remove(row,"_Perc")) %>%
  select(Variable = row,column,cor)

#### Combined Table ####
AggQualityTable <- ICC12 %>%
  bind_rows(SamplingRatio) %>%
  mutate(Value = round(Value, 3)) %>%
  spread(Type, Value) %>%
  left_join(AggTrueCorrs, by = "Variable") %>%
  mutate(Variable = str_replace_all(Variable,"_"," ") %>%
           str_replace(.,"\\.","-"))

## Correlations between MDE and CRDC
MDECRDCcor <- c(SPED = cor(SL16no5$Prop_Enrollment_IDEA, SL16no5$p16.SPED, use = "pairwise.complete.obs"),
                ELL = cor(SL16no5$Prop_Enrollment_LEP, SL16no5$p16.ELL, use = "pairwise.complete.obs"))

#### Initial Descriptives and Quality of Aggregated Covariates ####
save(BaselineStudentDescripsTotal, BaselineSchoolDescrips,
     CompareOCRMDEMSS, ICCMods, AggQualityTable, MDECRDCcor,
     file = "Saved Applied Results/Descrips_and_AggQuality.RData")
  
##############################################################################


###################################
#### Propensity Score Modeling ####

#### Prep ####
## Inversing Treatment and Control Designation
# Since SRO (treatment) > Not_SRO, this allows for more SRO units to be retained, thus increasing overall sample size (only matters if doing 1:many matching)
FinalSample <- FinalSample %>%
  mutate(Not_SRO = ifelse(SRO == 1, 0, 1))
CompleteSamp16L2Only <- CompleteSamp16L2Only %>%
  mutate(Not_SRO = ifelse(SRO == 1, 0, 1))

## Variables excluded from the analysis due to redundancy
ExcludeCovs <- c("White","p16.White",    # American_Indian, #p16.American_Indian?   # Redundant with other race/ethnic variables
                   "Eleven","Ninth","Fifth_L2","HS_L2")  # Makes interpretation middle school vs high school even though L2 is not mutually exclusive

#### Cluster Level Treatment ####
## SRO
tic()
CLPSModel <- f.build("SRO",SchoolCovNames[!(SchoolCovNames %in% ExcludeCovs)])
CLSRO.Mod <- glm(CLPSModel, family = binomial("logit"),
                 data = CompleteSamp16L2Only)
toc()


#### Subject Level Treatment ####
## SRO
tic()
SLPSModel <- f.build("SRO",c(StudentCovNames[!(StudentCovNames %in% ExcludeCovs)],
                                 SchoolCovNames[!(SchoolCovNames %in% ExcludeCovs)]))
SLSRO.Mod <- glm(as.formula(SLPSModel), family = binomial("logit"),
                 data = FinalSample)
toc() #13.2 sec

#### Adding Propensity Scores to dataset and Save Results ####
CompleteSamp16L2OnlywPS <- CompleteSamp16L2Only %>%
  mutate(Logit_CL = predict(CLSRO.Mod),
         PS_CL = fitted(CLSRO.Mod))
FinalSamp16wPS <- FinalSample %>%
  mutate(Logit_SL = predict(SLSRO.Mod),
         PS_SL = fitted(SLSRO.Mod)) %>%
  left_join(CompleteSamp16L2OnlywPS %>% select(level2, Logit_CL, PS_CL), by = "level2")

## Distribution of PS
psych::describeBy(FinalSamp16wPS$PS_SL, list(FinalSamp16wPS$SRO), mat = TRUE, digits = 4)
psych::describeBy(FinalSamp16wPS$PS_CL, list(FinalSamp16wPS$SRO), mat = TRUE, digits = 4)
max(FinalSamp16wPS$Logit_SL)
max(FinalSamp16wPS$Logit_CL)

# Used for my own sake, not necessarily for the manuscript
PSDistributionPlot <- FinalSamp16wPS %>%
  select(studentnumber, level2, SRO, starts_with("PS_")) %>%
  gather(Sample, PS, starts_with("PS")) %>%
  mutate(Sample = str_remove(Sample, "PS_"),
         Condition = case_when(SRO == 1 ~ "Yes",
                               SRO == 0 ~ "No") %>%
           as_factor()) %>%
  ggplot(aes(x = PS, fill = Condition)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(name = "SRO", values = c("#e66101", "#5e3c99")) +  #c("#3C3C3C","#5E0000")
  theme_bw(base_size = 20) +
  facet_wrap(~Sample) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))

## Exporting plot
ggsave(PSDistributionPlot, file = "PSDistributionPlot.png", width = 10, height = 5)

###################################

############################
#### Conditioning on PS ####

#### Running Matching Algorithm ####
## Cluster Level - matching schools
Match.CL <- matchit(as.formula(CLPSModel), data = CompleteSamp16L2OnlywPS,
                    method = "nearest", replace = FALSE, caliper = .2, distance = CompleteSamp16L2OnlywPS$Logit_CL)

## Subject Level - matching students
tic()
Match.SL <- matchit(as.formula(SLPSModel), data = FinalSamp16wPS,
                    method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS$Logit_SL)
toc() # 756 - 1020 sec

#### Weighting ####
## Cluster Level
Weight.CL <- WeightIt::get_w_from_ps(ps = CompleteSamp16L2OnlywPS$PS_CL, estimand = "ATE",
                        treat = CompleteSamp16L2OnlywPS$SRO, treated = 1)
## Subject Level
Weight.SL <- WeightIt::get_w_from_ps(ps = FinalSamp16wPS$PS_SL, estimand = "ATE",
                                     treat = FinalSamp16wPS$SRO, treated = 1)

#### Adding Indicator to dataset ####
# Cluster level
CompleteSamp16L2OnlywPS <- CompleteSamp16L2OnlywPS %>%
  mutate(CL_match = Match.CL$weights,
         CL_weight = Weight.CL)
# Subject level, then adding cluster info
FinalSamp16wPS <- FinalSamp16wPS %>%
  mutate(SL_match = Match.SL$weights,
         SL_weight = Weight.SL) %>%
  left_join(CompleteSamp16L2OnlywPS %>% select(level2, CL_match, CL_weight), by = "level2")

#### Saving Matches ####
save(StudentCovNames, StudentCovsRecode, StudentFacOrder,            # Covariate Names
     SchoolCovNames, SchoolCovsRecode, SchoolFacOrder,
     ExcludeCovs, OutcomeVars, OutcomeFacOrder, SampleFacOrder,      # Outcome names
     CLPSModel, SLPSModel, CLSRO.Mod, SLSRO.Mod, PSDistributionPlot, # PS model
     Match.CL, Match.SL, Weight.CL, Weight.SL,                       # Conditioning methods
     FinalSamp16wPS, CompleteSamp16L2OnlywPS,                        # Datasets
     file = "Saved Applied Results/PSModels_and_Matches.RData")

##########################################

###########################
#### Assessing Balance ####

load("Saved Applied Results/PSModels_and_Matches.RData")
load("Saved Applied Results/Descrips_and_AggQuality.RData")

### Outcome Descriptives and Effect Size Difference ####
## Provides the standardized difference without adjusting for covariates
OutcomeBalance <- map(.x = c("CL_match", "SL_match"),
                      ~bal.tab(formula = f.build("SRO",c(OutcomeVars)),
                               data = FinalSamp16wPS,
                               continuous = "std", binary = "std", s.d.denom = "pooled",
                               abs = FALSE, un = TRUE,
                               disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = FALSE,
                               weights = FinalSamp16wPS[[.x]], method = "matching")) %>%
  set_names(c("CL", "SL"))

## Cleaning up table for presentation
OutcomeBalanceTable <- Get_BalanceTable(OutcomeBalance$SL, UnAdj = "Un", renameUnAdj = "X") %>% mutate(Sample = "Full") %>%
  bind_rows(Get_BalanceTable(OutcomeBalance$SL, UnAdj = "Adj", renameUnAdj = "X") %>% mutate(Sample = "SL")) %>%
  bind_rows(Get_BalanceTable(OutcomeBalance$CL, UnAdj = "Adj", renameUnAdj = "X") %>% mutate(Sample = "CL")) %>%
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

#### Calculating Covariate Descriptives and Balance ####
# Calculating standardized mean difference for each of the conditioning methods using cobalt
CovBalanceObjects <- map(.x = c("CL_match","SL_match"),
                      ~bal.tab(formula = f.build("SRO",c(StudentCovNames, SchoolCovNames)),
        data = FinalSamp16wPS,
        continuous = "std", binary = "std", s.d.denom = "pooled",
        abs = TRUE, un = TRUE,
        disp.means = TRUE, disp.sd = TRUE, disp.v.ratio = TRUE,
        weights = FinalSamp16wPS[[.x]], method = "matching")) %>%
  set_names(c("CL", "SL"))

## Cleaning up table for presentation
CovBalanceTablePrep <- map2(.x = CovBalanceObjects, .y = names(CovBalanceObjects),
                        ~Get_BalanceTable(.x, UnAdj = "Adj", renameUnAdj = .y) %>%
                          select(Covariate, starts_with("Diff"), starts_with("V.Ratio"))) %>%
  reduce(left_join, by = "Covariate") %>%
  left_join(Get_BalanceTable(CovBalanceObjects$SL, UnAdj = "Un", renameUnAdj = "Full") %>%
              select(Covariate, starts_with("Diff"), starts_with("V.Ratio")), by = "Covariate") %>%
  mutate(Covariate = recode(Covariate, !!!SchoolCovsRecode) %>%
           recode(., !!!StudentCovsRecode) %>%
           PrintableCovariatesShortcut())

CovBalanceTable <- BaselineStudentDescripsTotal %>%
  left_join(CovBalanceTablePrep %>% mutate_if(is.numeric, ~round(., 2)), by = "Covariate") %>%
  select(Covariate, Total, SRO, No_SRO,Diff.Full, Diff.SL, Diff.CL,
         V.Ratio.Full, V.Ratio.SL, V.Ratio.CL) %>%
  mutate(Covariate = factor(Covariate, levels = c("SRO", StudentFacOrder, SchoolFacOrder)))
  
# ## Checking that unadjusted values were equivalent
map_lgl(c("M.0.Un","SD.0.Un","M.1.Un","SD.1.Un","Diff.Un","V.Ratio.Un"),
        ~identical(CovBalanceObjects$SL$Balance[[.x]], CovBalanceObjects$CL$Balance[[.x]]))   # identical

CovBalanceTableLong <- CovBalanceTablePrep %>%
  select(Covariate, starts_with("Diff"), starts_with("V.Ratio")) %>%
  gather(Temp, Value, -Covariate) %>%
  mutate(Temp = str_replace(Temp, "Diff", "ASD") %>%
           str_replace(., "V.Ratio", "VR")) %>%
  separate(Temp, into = c("Measure", "Sample"), sep = "\\.") %>%
  spread(Measure, Value)

#### Number of Matches for each Method ####
# Students
StudentSampleSizes <- FinalSamp16wPS %>%
  group_by(SRO) %>%
  summarize(SL_match = sum(SL_match),
            CL_match = sum(CL_match),
            # SL_weight = sum(SL_weight != 0),
            # CL_weight = sum(CL_weight != 0),
            Full = n()) %>%
  gather(Sample, Students, -SRO) %>%
  mutate(SRO = as.character(SRO))
# Schools
SchoolSampleSize <- map_dfr(c("CL_match", "SL_match"),
                            ~FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ] %>%
                              group_by(SRO) %>%
                              summarize(Schools = length(unique(level2)),
                                        Sample = .x)) %>%
  mutate(SRO = as.character(SRO)) %>%
  bind_rows(as.data.frame(table(CompleteSamp16L2OnlywPS$SRO)) %>%
              mutate(Sample = "Full") %>% rename(SRO = Var1, Schools = Freq))
# Combined Table
SampleSizes <- left_join(StudentSampleSizes, SchoolSampleSize, by = c("SRO", "Sample")) %>%
  mutate(Sample = str_remove(Sample, "_match"),     #str_replace(Sample, "SL_weight", "Subject") %>% str_replace(., "CL_weight", "Cluster"),
         SRO = recode(SRO, `0` = "TNo", `1` = "SRO")) %>%
  gather(Temp, Value, Students, Schools) %>%
  group_by(Sample, Temp) %>%
  mutate(Total = sum(Value),
         Percent = round((Value / Total) * 100, 1),
         Keep = paste0(Value, " (", Percent, ")")) %>% ungroup() %>%
  unite(col = "SRO_Type", SRO, Temp) %>%
  select(Sample, SRO_Type, Keep) %>%
  spread(SRO_Type, Keep)
# SampleSizes <- bind_rows(SampleSizes,  data.frame(SRO = "Total", as.list(colSums(SampleSizes[-1])))) # Adding a Total row

#### Balance Summaries ####
UnderThresholdTable <- CovBalanceTableLong %>%
  group_by(Sample) %>%
  summarize(ASD_Count = sum(ASD < 0.10),
            ASD_Percent = round((ASD_Count / nrow(CovBalanceTablePrep))*100, 1),
            VR_Count = sum(VR < 2.0, na.rm = TRUE),
            VR_Percent = round((VR_Count / sum(!is.na(CovBalanceTablePrep$V.Ratio.Full)))*100, 1)) %>%
  ungroup() %>%
  mutate(ASD = paste0(ASD_Count, " (", ASD_Percent, ")"),
         VR = paste0(VR_Count, " (", VR_Percent, ")"))

SSUTTable <- left_join(SampleSizes, UnderThresholdTable, by = "Sample") %>%
  select(Sample, ends_with("Students"), ends_with("Schools"), ASD, VR) %>%
  mutate(Sample = factor(Sample, levels = SampleFacOrder)) %>%
  arrange(Sample)
  

## Ordering variables based on largest difference in the original sample for plotting purposes
CovBalVarOrder <- CovBalanceTable %>%
  mutate(Covariate = fct_reorder(as_factor(Covariate), Diff.Full , .desc = FALSE))
CovBalVarOrder <- levels(CovBalVarOrder$Covariate)


## Balance Plot for all Samples
# Need to clean up covariate names
CovBalAllPlot <- CovBalanceTableLong %>%
  mutate(Covariate = factor(Covariate, levels = CovBalVarOrder),
         Sample = factor(Sample, levels = SampleFacOrder)) %>%
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

ggsave(CovBalAllPlot, file = "Covariate_Balance_Plot.png", width = 15, height = 21)


#### Saving data, models, and balance results ####
save(StudentCovNames, StudentCovsRecode, StudentFacOrder,            # Covariate Names
     SchoolCovNames, SchoolCovsRecode, SchoolFacOrder,
     ExcludeCovs, OutcomeVars, OutcomeFacOrder, SampleFacOrder,     # Outcome names
     FinalSamp16wPS, CompleteSamp16L2OnlywPS,
     CovBalanceObjects, CovBalanceTablePrep, CovBalanceTable, CovBalanceTableLong,
     OutcomeBalance, OutcomeBalanceTable,
     SampleSizes, UnderThresholdTable, SSUTTable,
     CovBalVarOrder, CovBalAllPlot,
     file = "Saved Applied Results/PSbalance.RData")
