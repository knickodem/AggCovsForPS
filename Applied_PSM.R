###################################################
#                                                 #
#       PS Matching for Applied Study             #
#                                                 #
###################################################

## Loading packages, functions, and cleaned dataset
source("Applied_PF.R")
load("Saved Applied Results/CompleteData.RData")


#######################
####  Descriptives ####

#### Descriptives at the Student-level ####
## Includes school characteristics - Aligns with WWC standards (p. 25) - Actually I'm not entirely sure this is what WWC implies
## Across all schools
DescripsAtL1 <- FinalSample %>%
  select(studentnumber, level2, one_of(StudentCovNames), SRO, one_of(SchoolCovNames)) %>%
  gather(Variable, Value, -studentnumber, -level2) %>%
  Get_Descriptives(Value, Variable, AllContinuous = FALSE)

## By SRO status
balance.CompL1 <- bal.tab(formula = f.build("SRO", c(StudentCovNames, SchoolCovNames)),
                          data = FinalSample,
                          continuous = "std", binary = "std", s.d.denom = "all",
                          abs = TRUE, disp.v.ratio = TRUE)
balCompL1 <- Get_BalanceTable(balance.CompL1, "Un", replace01 = c("_No_SRO", "_SRO"), ModName = "Original")
#mutate(Variable = fct_reorder(as_factor(Variable),Std_Diff,.desc=TRUE))


#### Descriptives at the School-level ####
## Only includes school characteristics
# Values differ from those calculated at the student level because the number of students per school is unequal.
DescripsAtL2 <- CompleteSamp16L2Only %>%
  select(level2, SRO, one_of(SchoolCovNames)) %>%
  gather(Variable, Value, -level2) %>%
  Get_Descriptives(Value, Variable, AllContinuous = TRUE)

## By SRO status
balance.CompL2 <- bal.tab(formula = f.build("SRO",SchoolCovNames),
                           data = CompleteSamp16L2Only,
                          continuous = "std", binary = "std", s.d.denom="all",
                          abs = TRUE, disp.v.ratio = TRUE)
balCompL2 <- Get_BalanceTable(balance.CompL2, "Un", replace01 = c("_No_SRO", "_SRO"), ModName = "Original")
         #mutate(Variable = fct_reorder(as_factor(Variable),Std_Diff,.desc=TRUE))

#################################################

###########################
#### Aggregate Quality ####

# Aggregate values included responses from students with missing data and from grade 5 students

#### ICCs and Sampling Ratio ####
## Calculate ICC(1) and ICC(2) for each agg covariate
ICCMods <- map(str_remove_all(StudentCovNames, "_Center"), ~lmer(as.formula(paste0(.x," ~ 1 + (1|level2)")), data = FinalSample)) %>%
  set_names(str_remove_all(StudentCovNames, "_Center"))
ICC12 <- map_df(ICCMods, ~Get_ICC(.x, type = "1")) %>% mutate(Type = "ICC(1)") %>%     # Extracting ICC(1) from each model
  bind_rows(map_df(ICCMods, ~Get_ICC(.x, type = "2")) %>% mutate(Type = "ICC(2)")) %>% # Extracting ICC(2) from each model
  gather(Variable, Value, -Type)

## Calculating sampling ratio by school on each covariate, then averaging
SamplingRatio <- SL16no5 %>%
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

#### Initial Descriptives and Quality of Aggregated Covariates ####
save(DescripsAtL1, DescripsAtL2, balCompL1, ICCMods, AggQualityTable,
     file = "Saved Applied Results/Descrips_and_AggQuality.RData")
  
##############################################################################


###################################
#### Propensity Score Modeling ####

#### Prep ####
## Inversing Treatment and Control Designation
# Since SRO (treatment) > Not_SRO, this allows for more SRO units to be retained, thus increasing overall sample size
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
toc() #7.75 sec when using FinalSample; .5 sec when using CompleteSamp16L2Only?
# summary(CLSRO.Mod)
# car::vif(CLSRO.Mod) # High for PIO, SC, and FCS; does this matter though? Does only balance matter?

## Not SRO
NotCLPSModel <- f.build("Not_SRO",SchoolCovNames[!(SchoolCovNames %in% ExcludeCovs)])
NotCLSRO.Mod <- glm(NotCLPSModel, family = binomial("logit"),
                 data = CompleteSamp16L2Only)

#### Ignore Cluster with Subject Level Treatment ####
## Not SRO
tic()
ICPSModel <- f.build("SRO",c(StudentCovNames[!(StudentCovNames %in% ExcludeCovs)],
                                 SchoolCovNames[!(SchoolCovNames %in% ExcludeCovs)]))
ICSRO.Mod <- glm(as.formula(ICPSModel), family = binomial("logit"),
                 data = FinalSample)
toc() #13.2 sec
# summary(ICSRO.Mod)
# car::vif(ICSRO.Mod) # some high numbers for L2 vars (same as CL model)

## Not SRO
NotICPSModel <- f.build("Not_SRO",c(StudentCovNames[!(StudentCovNames %in% ExcludeCovs)],
                                    SchoolCovNames[!(SchoolCovNames %in% ExcludeCovs)]))
NotICSRO.Mod <- glm(as.formula(NotICPSModel), family = binomial("logit"),
                 data = FinalSample)

# #### I'm starting suspect these will never work
# ## Fixed Effects with Subject Level Treatment
# tic()
# FESRO.Mod <- glm(as.formula(paste0(StudentPSModel," + factor(level2)")), family = binomial("logit"),
#                  data = FinalSample)
# toc() # 1231.9 sec and model did not converge
# # summary(FESRO.Mod)
# # car::vif(FESRO.Mod) # 
# 
# 
# ## Random Effects with Subject Level Treatment
# tic()
# RESRO.Mod <- glmer(as.formula(paste(ICPSModel, "+ (1|level2)")),
#                    family = binomial("logit"), control = glmerControl(optimizer = "bobyqa"), 
#                    data = FinalSample)
# toc() # sec
# summary(RESRO.Mod)
# car::vif(RESRO.Mod) # 
# 
# tt <- getME(FESRO.Mod,"theta")
# ll <- getME(FESRO.Mod,"lower")
# min(tt[ll==0])


#### Adding Propensity Scores to dataset and Save Results ####
CompleteSamp16L2OnlywPS <- CompleteSamp16L2Only %>%
  mutate(Logit_CL = predict(CLSRO.Mod),
         PS_CL = fitted(CLSRO.Mod),
         Logit_NotCL = predict(NotCLSRO.Mod),
         PS_NotCL = fitted(NotCLSRO.Mod))
FinalSamp16wPS <- FinalSample %>%
  mutate(Logit_IC = predict(ICSRO.Mod),
         PS_IC = fitted(ICSRO.Mod),
         Logit_NotIC = predict(NotICSRO.Mod),
         PS_NotIC = fitted(NotICSRO.Mod),
         # Logit_FE = predict(FESRO.Mod),
         # PS_FE = fitted(FESRO.Mod),
         # Logit_RE = predict(RESRO.Mod),
         # PS_RE = fitted(RESRO.Mod)
         ) %>%
  left_join(CompleteSamp16L2OnlywPS %>% select(level2, Logit_CL, PS_CL, Logit_NotCL, PS_NotCL), by = "level2")

## Distribution of PS
PSDistributionPlot <- FinalSamp16wPS %>%
  select(studentnumber, level2, SRO, starts_with("PS_")) %>%
  gather(Sample, PS, starts_with("PS")) %>%
  mutate(Sample = str_remove(Sample, "PS_"),
         Condition = case_when(SRO == 1 ~ "Has SRO",
                               SRO == 0 ~ "No SRO") %>%
           as_factor()) %>%
  ggplot(aes(x = PS, fill = Condition)) +
  geom_density(alpha = 0.4) +
  # geom_vline(data = PlotMeans, aes(xintercept = Mean, color = Retained), size = 2) +
  # scale_x_continuous(breaks = seq(0,21,3)) +
  scale_fill_manual(values = c("#3C3C3C","#5E0000")) +
  # scale_color_manual(values = c("#3C3C3C","#5E0000")) +
  theme_classic(base_size = 20) +
  # theme(legend.justification = c(1,1), legend.position = c(.95,.95)) +
  facet_wrap(~Sample)

## Exporting plot
ggsave(PSDistributionPlot, file = "PSDistributionPlot.png", width = 10, height = 5)

## Saving results
save(ExcludeCovs, CLPSModel, ICPSModel, NotCLPSModel, NotICPSModel,
     CLSRO.Mod, ICSRO.Mod, NotCLSRO.Mod, NotICSRO.Mod,#FESRO.Mod, RESRO.Mod, StudentPSModel,
     FinalSamp16wPS, CompleteSamp16L2OnlywPS, PSDistributionPlot,
     file = "Saved Applied Results/PSModels.RData")


###################################

##################
#### Matching ####

#### Running Matching Algorithm ####

## Cluster Level - matching schools
# SRO
Match.CL <- matchit(as.formula(CLPSModel), data = CompleteSamp16L2OnlywPS,
                    method = "nearest", replace = FALSE, caliper = .2, distance = CompleteSamp16L2OnlywPS$Logit_CL)
# Not SRO
NotMatch.CL <- matchit(as.formula(NotCLPSModel), data = CompleteSamp16L2OnlywPS,
                    ratio = 2, method = "nearest", replace = FALSE, caliper = .2, distance = CompleteSamp16L2OnlywPS$Logit_NotCL)

## Subject Level - matching students (ignoring clusters)
# SRO
tic()
Match.IC <- matchit(as.formula(ICPSModel), data = FinalSamp16wPS,
                    method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS$Logit_IC)
toc() # 1020 sec
# Not SRO
tic()
NotMatch.IC <- matchit(as.formula(NotICPSModel), data = FinalSamp16wPS,
                       ratio = 2, method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS$Logit_NotIC)
toc() # 394 sec


# # Fixed effects model
# Match.FE <- matchit(as.formula(NotICPSModel), data = FinalSamp16wPS,
#                     method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS$Logit_FE)
# 
# # Random effects model
# Match.RE <- matchit(as.formula(ICPSModel), data = FinalSamp16wPS,   # Model is ignored since distance is provided
#                     method = "nearest", replace = FALSE, caliper = .2, distance = FinalSamp16wPS$Logit_RE)


#### Adding Indicator to dataset ####
# Cluster level
CompleteSamp16L2OnlywPS <- CompleteSamp16L2OnlywPS %>%
  mutate(CL_1to1 = Match.CL$weights,
         CL_1to2 = NotMatch.CL$weights)
# Subject level, then adding cluster info
FinalSamp16wPS <- FinalSamp16wPS %>%
  mutate(IC_1to1 = Match.IC$weights,
         IC_1to2 = NotMatch.IC$weights
         # Matched_FE = Match.FE$weights,
         # Matched_RE = Match.RE$weights
         ) %>%
  left_join(CompleteSamp16L2OnlywPS %>% select(level2, CL_1to1, CL_1to2), by = "level2")

## Student
table(FinalSamp16wPS$SRO,FinalSamp16wPS$IC_1to1) # Retains an extra (69398 - 15209) = 54189 treatment students
table(FinalSamp16wPS$IC_1to1,FinalSamp16wPS$IC_1to2 != 0) # Retains an extra (69398 - 15209) = 54189 treatment students

## School
table(CompleteSamp16L2OnlywPS$Not_SRO, CompleteSamp16L2OnlywPS$CL_1to2) # Retains an extra (166 - 114) = 52 treatment schools
table(CompleteSamp16L2OnlywPS$CL_1to1, CompleteSamp16L2OnlywPS$CL_1to2 != 0) # although same number of control schools retained (114), they werent the same ones


#### Saving Matches ####
save(FinalSamp16wPS, CompleteSamp16L2OnlywPS,
     Match.CL, Match.IC, NotMatch.CL, NotMatch.IC, #Match.FE,Match.RE,
     file = "Saved Applied Results/Matches.RData")

##########################################

###########################
#### Assessing Balance ####

# load("Saved Applied Results/PSModels.RData")
load("Saved Applied Results/Matches.RData")


#### Standardized Mean Differences ####
# Calculating standardized mean difference for each of the matching methods
# "pooled" uses the standard deviation from the two groups from the full, unadjusted sample
# baltab uses weights from MatchIt object

## At the School Level
# Equating
SchoolBalanceObjects <- map(c("CL_1to1","CL_1to2"),
                            ~bal.tab(formula = f.build("SRO",c(SchoolCovNames)), data = CompleteSamp16L2OnlywPS,
                                     continuous = "std", binary = "std", s.d.denom = "all",
                                     disp.v.ratio = TRUE, abs = TRUE,
                                     weights = CompleteSamp16L2OnlywPS[[.x]], method = "matching"))

# Combining into one table
SchoolBalAll <- map2_dfr(.x = SchoolBalanceObjects, .y = c("CL_1to1","CL_1to2"),
                          ~Get_BalanceTable(.x, UnAdj = "Adj", replace01 = c("_No_SRO","_SRO"), ModName = .y)) %>%
  bind_rows(balCompL2) %>%
  rename(Method = Model)


## At the Student Level
# Equating
StudentBalanceObjects <- map(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
                             ~bal.tab(formula = f.build("SRO",c(StudentCovNames,SchoolCovNames)), data = FinalSamp16wPS,
                                      continuous = "std", binary = "std", s.d.denom = "all",
                                      disp.v.ratio = TRUE, abs = TRUE,
                                      weights = FinalSamp16wPS[[.x]], method = "matching")) %>%
  set_names(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"))

# Combining into one table
StudentBalAll <- map2_dfr(.x = StudentBalanceObjects, .y = names(StudentBalanceObjects),
                           ~Get_BalanceTable(.x, UnAdj = "Adj", replace01 = c("_No_SRO","_SRO"), ModName = .y)) %>%
  bind_rows(balCompL1) %>%
  rename(Method = Model)


#### Number of Matches for each Method ####
## Students
StudentSampleSizes <- map_dfr(StudentBalanceObjects,
                             ~tibble::rownames_to_column(.x[["Observations"]],"Observations"),.id = "Method") %>%
  bind_rows(as.data.frame(balance.CompL1$Observations) %>% mutate(Method = "Original")) %>%
  mutate(Total = Control + Treated)

# StudentSampleSizes %>% filter(Observations == "Matched (Unweighted)")

## Schools - only provides unweighted
SchoolSampleSizes <- map_dfr(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
                             ~FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ] %>% group_by(SRO) %>%
                               summarize(schools = length(unique(level2)),
                                         Method = .x) %>%
                               spread(SRO, schools)) %>%
  rename(Control = `0`, Treated = `1`) %>%
  bind_rows(as.data.frame(balance.CompL2$Observations) %>% mutate(Method = "Original")) %>%
  mutate(Total = Control + Treated)


## Combined for output table
ComboSampleSizes <- StudentSampleSizes %>%
  filter(Observations == "Matched (ESS)" | is.na(Observations)) %>%     # Could use (Unweighted) too
  mutate_at(vars(Treated,Total),~round(.,1)) %>%   # Only need this if using (ESS)
  inner_join(SchoolSampleSizes, by = "Method", suffix = c(".Stu",".Sch")) %>%
  mutate(Control = paste0(Control.Stu, " (", Control.Sch,")"),
         Treated = paste0(Treated.Stu, " (", Treated.Sch,")"),
         Total = paste0(Total.Stu, " (", Total.Sch,")"))


#### Balance Summaries ####
## Percent of variables under .1 threshold
UnderThreshold <- StudentBalAll %>%
  filter(!(Variable %in% c("distance"))) %>% # ExcludeCovs
  group_by(Method) %>% summarize(Count = sum(Diff < .10),
                                Variables = n()) %>%
  mutate(Percent = round(Count / Variables*100, 1))

## Percent of continuous variables with variance ratios < 1.2
# Note: bal.tab automatically reverses all ratios to be > 1
VarRatio <- StudentBalAll %>%
  filter(!(Variable %in% c("distance")) & !is.na(V.Ratio)) %>%
  group_by(Method) %>% summarize(Count = sum(V.Ratio < 1.20), # A Cohen's d of .10 is equivalent ot an Odds Ratio of 1.20
                                Variables = n()) %>%
  mutate(Percent = round(Count / Variables*100, 1))

## Model with best balance on the most variables
MostBalanced <- StudentBalAll %>% filter(!(Variable %in% c("distance"))) %>%
  group_by(Variable) %>%
  filter(Diff == min(Diff)) %>%
  group_by(Method) %>% summarize(Count = n()) %>%
  mutate(Variables = sum(Count),
         Percent = round(Count / Variables*100, 1))


## Combining Under Threshold, Most Balanced, and number of matches
# Need to report NumVarRatioVariables in table notes
NumVarRatioVariables <- unique(VarRatio$Variables)
UTMBVRSS <- left_join(UnderThreshold, MostBalanced, by = c("Method","Variables"), suffix=c(".UT",".MB")) %>%
  left_join(VarRatio %>% select(-Variables), by = "Method") %>%
  left_join(ComboSampleSizes %>% select(Method, Control, Treated, Total), by = "Method")
  # select(Method = Model, Matched,`Num Variables` = Variables,
  #        Count.UT, Percent.UT, Count.VR = Count, Percent.VR = Percent, Count.MB, Percent.MB)


## Ordering variables based on largest difference in the original sample for plotting purposes
BalVarOrder <- StudentBalAll %>% filter(Method == "Original") %>%
  mutate(Variable = fct_reorder(as_factor(Variable), Diff , .desc=TRUE))
BalVarOrder <- levels(BalVarOrder$Variable)


## Balance Plot for all Models
BalAllPlotGrey <- StudentBalAll %>% filter(!(Variable %in% c("distance"))) %>%
  mutate(Variable = factor(Variable, levels = rev(BalVarOrder))) %>%
  ggplot(aes(x = Variable, y = Diff, group = Method)) +
  geom_hline(yintercept = .1, linetype = 1, color = "gray5") +
  geom_point(aes(shape = Method, color = Method),size = 4 ,position = position_dodge(1)) +
  ylab("Standardized Difference") + xlab("") +
  scale_shape_manual(values=c(19,18,17,15,7)) + #,3
  scale_color_grey() +
  #scale_color_brewer(palette="Dark2") +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_line(color = "gray87"),
        panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.justification=c(1,0),legend.position=c(.9,.1))

#### Student-Level Demographics from Original and Best fitting sample ####
## Might move this to HLM script


#### Saving data, models, and balance results ####
save(FinalSamp16wPS, CompleteSamp16L2OnlywPS,
     StudentCovNames, SchoolCovNames, ExcludeCovs, OutcomeVars,
     SchoolBalanceObjects, SchoolBalAll, StudentBalanceObjects, StudentBalAll,
     ComboSampleSizes, VarRatio, UnderThreshold, MostBalanced, UTMBVRSS,
     BalVarOrder, BalAllPlotGrey,
     file = "Saved Applied Results/PSbalance.RData")
