###################################################
#                                                 #
#       Data Prep for Applied Study               #
#                                                 #
###################################################

## Loading packages and functions
source("Applied_PF.R")

###############################################################
#### Loading Data and transforming student-level variables ####

#### MDE School Enrollment Data ####
MDEenrollOrig <- readxl::read_excel("enrollment_public_file_2016.xlsx", sheet = "School")
names(MDEenrollOrig) <- str_replace_all(names(MDEenrollOrig)," ","_")
MDEenroll <- MDEenrollOrig %>%
  filter(Grade == "All Grades") %>%
  mutate(District_Number = ifelse(District_Name=="MINNEAPOLIS PUBLIC SCHOOL DIST.","999",District_Number),
         Temp = paste0(District_Number,District_Type,School_Number) %>% as.numeric(),
         level2 = Temp*3+1,
         p16.Asian = Total_Asian_Percent + Total_Native_Hawaiian_or_Pacific_Islander_Percent,
         p16.FRL = str_remove(`Total_Students_Eligible_for_Free_or_Reduced-Priced_Meals_Percent`,"%") %>% as.numeric(),
         p16.SPED = str_remove(Total_Students_Receiving_Special_Education_Services_Percent,"%") %>% as.numeric(),
         p16.ELL = str_remove(Total_English_learner_Identified_Percent,"%") %>% as.numeric(),
         p16.Homeless = str_remove(Total_Students_Experiencing_Homelessness_Percent,"%") %>% as.numeric()) %>%
  select(level2,TotalEnrollment16 = Total_Enrollment,p16.Female = Total_Female_Percent,
         p16.American_Indian = Total_American_Indian_or_Alaska_Native_Percent,
         p16.Asian, p16.Latino = Total_Hispanic_or_Latino_Percent,p16.Black = Total_Black_or_African_American_Percent,
         p16.White = Total_White_Percent,p16.Multiracial = Total_Two_or_More_Races_Percent,
         p16.FRL:p16.Homeless) %>%
  mutate_at(vars(starts_with("p16")),~./100)

#### Minnesota Student Survey ####
## Original data
mssOrig <- read_sav("G:/My Drive/MyDrG/Data/MSS2013-2016-SECURE-20180701.sav")
# names(mssOrig)

## Checking alpha of SRO_Discomfort scale
# mssOrig %>% filter(Year==2016) %>% select(starts_with("W26")) %>% psych::alpha()

## Keeping only 2016, dichtomizing, and renaming student-level variables
mss16 <- mssOrig %>%
  filter(Year == 2016) %>%
  mutate(Fifth = ifelse(Y2==5, 1, 0),
         Eighth = ifelse(Y2==8, 1, 0),
         Ninth = ifelse(Y2==9, 1, 0),
         Eleven = ifelse(Y2==11, 1, 0),
         Female = ifelse(Y1==2, 1, 0),
         # American_Indian = ifelse(racegraphs==0,NA,ifelse(racegraphs==1,1,0)), # using raceethnic designation that aligns more with MDE,
         # Asian = ifelse(racegraphs==0,NA,ifelse(racegraphs==2,1,0)),           # but also accounting for students who said Hmong or Somali only
         # Black = ifelse(racegraphs==0,NA,ifelse(racegraphs==3,1,0)),
         # White = ifelse(racegraphs==0,NA,ifelse(racegraphs==5,1,0)),
         # Multiracial = ifelse(racegraphs==0,NA,ifelse(racegraphs==6,1,0)),
         # Latinx = ifelse(racegraphs==0,NA,ifelse(racegraphs==7,1,0)),
         # Somali = ifelse(racegraphs==0,NA,ifelse(racegraphs==8,1,0)),
         # Hmong = ifelse(racegraphs==0,NA,ifelse(racegraphs==9,1,0)),
         American_Indian = ifelse(raceethnic == 1, 1, 0),
         Asian = ifelse(raceethnic == 9 & racegraphs == 9, 1, ifelse(raceethnic %in% c(2, 4), 1, 0)),
         Black = ifelse(raceethnic == 9 & racegraphs == 8, 1, ifelse(raceethnic == 3, 1, 0)),
         White = ifelse(raceethnic == 5, 1, 0),
         Multiracial = ifelse(raceethnic == 6, 1, 0),
         Latinx = ifelse(raceethnic == 7, 1, 0),
         Safe_Travel = ifelse(is.na(Y22a),NA,ifelse(Y22a %in% c(3,4),1,0)), # to and from school; binary agree/disagree
         Unsupervised = Y30 - 1,                                            # only asked of Grade5
         Sports_Team = ifelse(W33a==1,0,1),                                  # Participated at least 1 day per week
         School_Club = ifelse(Y34c==1,0,1),
         Community_Club = ifelse(Y34g==1,0,1),
         Tutoring = ifelse(Y34d==1,0,1),
         Leadership = ifelse(Y34e==1,0,1),
         Artistic_Lessons = ifelse(W33e==1,0,1),
         Physical_Lessons = ifelse(W33f==1,0,1),
         Religious_Activity = ifelse(Y34h==1,0,1),
         Has_SRO = case_when(W25==1 ~ 1,W25==2 ~ 0, TRUE ~ NA_real_),
         SRO_Discomfort = W26a + W26b + W26c,                            # Not planning to use this b/c of missingness and lack of psychometric evidence
         Poor_Health = Y36 - 1,                                          # Treating as continuous; Excellent (0) to Poor (4)
         Vandalism = ifelse(Y77b == 1, 0, 1),                            # At least once in last 12 months
         Theft = ifelse(Y77d == 1, 0, 1)) %>%
  mutate_at(vars(American_Indian:Latinx),~ifelse(racegraphs == 0, NA, .)) %>%
  rename(GPA = Grades, Grade = Y2, Special_Education = Y11, FRL = Y12,
         Skipped_Class = Y15d, Skipped_School = Y16d, Sick_3_Days = stayedhome3, # All variables were dichotomized as at least once in last 30 days (even iss and oss),
         Disciplined = disc, In.School_Suspension = iss, Out.School_Suspension = oss,  # expect sick is 3 days and skipped meals was asked as Yes/No; Suspended = suspended
         Skipped_Meal = Y48, Sleep_8_Hrs = Sleep8, Smoked_Cigarette = dCig,
         Smoked_eCig = dEcig, Smoked_Marijuana = dMar, Drank_Alcohol = dAlc, Took_Rx = dRx,
         Commitment_to_Learning = CtL, Positive_Identity = PI, Social_Competence = SC,
         Empowerment = EM, Family.Community_Support = FCS, Teacher.School_Support = TSS, OST_Experiences = PosEx,
         Bullied = BD, Bullying = BLY, Mental_Distress = MD, Family_Violence = FV)

# # items about school safety/criminal activity Y23 and 24? OST Y31 and Y32? Only asked in 2013
#
# prop.table(table(mssOrig$Y31,mssOrig$Y2,useNA = "always"),2)
# table(mss16$Y77b,mss16$Grade,useNA = "always")
# table(mssOrig$Y32,mssOrig$Year,useNA = "always")


## Extracting student-level information
StudentVars <- mss16 %>%
  select(studentnumber,level2,Grade,Has_SRO,                                     # IDs and SRO; SRO_Discomfort,
         Fifth:Latinx,Special_Education,FRL,LGB,Homeless,Parentjail,Trauma1,     # Background Characteristics ;Hmong
         Unsupervised:Religious_Activity,OST_Experiences,Safe_Travel,            # Out of School Time
         Poor_Health,Sick_3_Days,PhyEd,Skipped_Meal,Diet,Sleep_8_Hrs,            # Health and Diet
         Skipped_Class,Skipped_School,Gambling,Vandalism,Theft,                  # Truancy/Crime  
         Disciplined,In.School_Suspension,Out.School_Suspension,                 # Discipline
         Smoked_Cigarette, Smoked_eCig,Smoked_Marijuana,Drank_Alcohol,Took_Rx,   # Drug Use
         Bullied, Bullying, Mental_Distress, Family_Violence,                    # Developmental Challenges
         GPA,College,Commitment_to_Learning,Positive_Identity,Social_Competence, # Academic and SE Skills Outcomes
         Empowerment,Family.Community_Support,Teacher.School_Support             # Support Outcomes
  )

# ## Proportion missing by Grade
# StuVarPropMiss <- StudentVars %>%
#   group_by(Grade) %>%
#   summarize_at(vars(-group_cols(),-studentnumber,-level2),~round(sum(is.na(.))/n(),3))
# NotG5 <- StuVarPropMiss %>%
#   gather(Variable,Percent,-Grade) %>% filter(Grade==5 & Percent > .99) # Items not asked of grade 5 students


#########################################################

###########################################
#### Extracting School-Level Variables ####

## Identifying variables to aggregate
AggVars <- StudentVars %>%
  select(-studentnumber,-level2,-Grade) %>%
  names()

## Calculates school mean and n for each covariate
mssAgg <- mss16 %>%
  select(level2,one_of(AggVars)) %>%
  group_by(level2) %>%
  summarise_at(vars(-group_cols(),one_of(AggVars)),
               list(Perc = ~mean(.,na.rm=TRUE),    # For most covariates the mean in the Percent of students with characteristic
                    n = ~sum(!is.na(.)))) %>%
  ungroup() %>%
  mutate_all(~ifelse(is.nan(.),NA,.))

## Determining grade levels present for each school (when grades 9 and 11 are combined into HS)
GradeLevels16 <- mss16 %>%
  mutate(Grade = ifelse(Grade %in% c(9,11), "HS",as.character(Grade))) %>%
  select(Grade,level2) %>%
  unique() %>%                    # 1081 unique schools
  mutate(G2 = Grade) %>%
  spread(Grade, G2) %>%
  unite("Grade_Levels", `5`:HS,sep = ",", remove = TRUE) %>%
  mutate(Grade_Levels = str_replace_all(Grade_Levels, "NA,?", "") %>% str_remove(",$"))

## Extracting school-level information and joining grade levels and aggregated variables
SchoolLevel16 <- mss16 %>% 
  select(Year,level2:mss_sample,Diversity.16,medianincome,TC,                   ## Variables from MSS, MDE (non-MCA), or ACS; TotalSOC16:
         LEA_HBPOLICY_IND:Expenses_Per_Student,                                 ## Variables from OCR
         starts_with("pLevel"),mean.read05:sd.read10,mean.math05,sd.math11) %>% ## MCA variables from MDE
  unique() %>%                                                                  # 1082 schools (school 1815483031 is listed twice, but has NA on all OCR vars for one case)
  filter(!(level2 == 1815483031 & is.na(SCH_FTECOUNSELORS))) %>%                # a random OCR variable was used here
  left_join(GradeLevels16, by = "level2") %>%
  left_join(mssAgg, by = "level2") %>%
  left_join(MDEenroll, by = "level2") %>%
  select(Year,Grade_Levels,everything())

# SchoolLevel16 %>% select(level2, Grade_Levels, TotalStudents16, TotalEnrollment16,
#                          starts_with("mde"),starts_with('p16')) %>% View()

###############################################################

#######################################################
#### Determining school-level academic performance ####

# Need to use % proficient given that the scale scores are not comparable across grades
# Interpretation becomes the % of students who performed higher than the cut score designated for their grade
# Although % proficent reduces variability, I think it is at least appropriate given that all schools use the same standards and cut scores

## At each MCA grade
GradeAP <- SchoolLevel16 %>% 
  select(level2, Grade_Levels, starts_with("pLevel")) %>%
  gather(Temp, Percent, -level2, -Grade_Levels) %>%
  mutate(Performance_Level = str_sub(Temp, start = 7, end = 7),
         MCA_Grade = str_sub(Temp, start = -2),
         MCA_Grade = ifelse(MCA_Grade %in% c("10","11"), "HS",str_remove(MCA_Grade, "0")),
         Subject = str_sub(Temp, start = -6, end = -3) %>%
           str_to_title()) %>%
  filter(!is.na(Percent))   # 10092
#filter(str_detect(Grade_Levels,MCA_Grade)) # 10028 removes 64 records (1 grade each from 8 schools) see below for details

# ## When MCA and MSS disagree on grade level represented
# GradeAP_MCAnoMSS <- GradeAP %>%  filter(!str_detect(Grade_Levels,MCA_Grade)) # 8 schools
# Investigate <- GradeAP %>% filter(level2 %in% unique(GradeAP_MCAnoMSS$level2))
# Investigate2 <- SL16 %>% filter(level2 %in% unique(GradeAP_MCAnoMSS$level2))
#  # l2        ; GL   ; MCAGL ; According to MDE's 2016 Student Enrollment file (x - 1 / 3)
#  # 303004    ; HS   ; 8, HS ; grades 7 - 12, so did not give MSS to grade 8?
#  # 86119321  ; 8    ; 5     ; grades 1 - 8, so did not give MSS to grade 5? Enrollment too low to report MCA for grade 8; should probably just eliminate this school
#  # 187503634 ; HS   ; 8, HS ; grades 6 - 12, so did not give MSS to grade 8?
#  # 299709319 ; 5    ; 5, 8  ; grades K - 8, so did not give MSS to grade 8?
#  # 299709676 ; 5    ; 5, 8  ; grades K - 8, so did not give MSS to grade 8?
#  # 866403010 ; HS   ; 8, HS ; grades 7 - 12, so did not give MSS to grade 8?
#  # 869103010 ; 8    ; 5, 8  ; grades 5 - 8, so did not give MSS to grade 5?
#  # 870903061 ; HS   ; 8, HS ; grades 7 - 12, so did not give MSS to grade 8?

## School-level averages (if > 1 grade is represented)
# Contains all proficiency levels
SLAP <- GradeAP %>%
  group_by(level2,Grade_Levels,Subject,Performance_Level) %>%
  summarize(n_MCA_Grades = n(),
            Mean  = mean(Percent)) %>%
  mutate(Sum = sum(Mean)) %>%
  group_by(level2) %>%
  mutate(check = n())       # 12 schools are missing either reading or math scores while 27 schools have neither; check <- SchoolLevel16 %>% filter(level2 %in% unique(SLAP[SLAP$check==4,]$level2))

## Calculating % proficient for math and for reading
SLAP_Proficiency <- SLAP %>%
  filter(Performance_Level %in% c("M","E")) %>%
  unite("Temp",Subject, Performance_Level) %>%
  select(level2,Grade_Levels,Temp,Mean) %>%
  spread(Temp,Mean) %>%
  mutate(Proficient_Math = (Math_E + Math_M)/100,
         Proficient_Read = (Read_E + Read_M)/100)

#######################################################


#### Identifying Possible Treatments ####
# ## Distributions of potential 'treatment' variables
# # For most of these I need to examine the extant literature on what is a meaningful cut point
# table(SchoolLevel16$SCH_FTECOUNSELORS >= 1, useNA = "always")
# table(SchoolLevel16$SCH_FTESERVICES_NUR >= 1, useNA = "always")
# table(SchoolLevel16$SCH_FTESERVICES_PSY >= 1, useNA = "always")
# table(SchoolLevel16$SCH_FTESERVICES_SOC >= 1, useNA = "always")
# table(SchoolLevel16$SCH_FTESECURITY_GUA, useNA = "always")
# table(SchoolLevel16$SCH_FTESECURITY_LEO, useNA = "always")
# table(SchoolLevel16$Prop_FTE_Under_2yrs > .1, useNA = "always")

# ## Need to determine what the cutoffs should be for designating has SRO
# ggplot(SchoolLevel16,aes(x=SCH_FTESECURITY_LEO,y = Has_SRO_Perc)) +
#   geom_point() + #geom_density() +
#   theme_bw()
# table(SchoolLevel16$Has_SRO_Perc>.5,SchoolLevel16$SCH_FTESECURITY_LEO <.2, useNA="always")
# table(SchoolLevel16$SCH_FTESECURITY_LEO >= .2, useNA = "always") %>% sum()
# SchoolLevel16 %>% select(level2,Grade_Levels,mss_sample, TotalStudents16,Has_SRO_Perc,SCH_FTESECURITY_LEO) %>% 
#   left_join(mss16 %>% group_by(level2) %>% summarize(Count = sum(!is.na(Has_SRO))), by="level2") %>%
#   filter(Grade_Levels!="5" & SCH_FTESECURITY_LEO < .2) %>% View()
# # For non grade 5 only schools where OCR data said no SRO present, students largely agreed or disagreed completely
# SchoolLevel16 %>% filter(SCH_FTESECURITY_LEO==0) %>%
#   ggplot(aes(Has_SRO_Perc)) + geom_histogram()

#######################################################
#### Finalizing school-level and combined datasets ####

#### Transforming school-level variables ####
## also joining academic performance and MDE enrollment covariates
SL16 <- SchoolLevel16 %>%
  mutate(Fifth_L2 = ifelse(str_detect(Grade_Levels, "5"),1,0),
         Eighth_L2 = ifelse(str_detect(Grade_Levels, "8"),1,0),
         HS_L2 = ifelse(str_detect(Grade_Levels, "HS"),1,0),
         Charter_Magnet = ifelse(SCH_STATUS_MAGNET==1 | SType == 7, 1, 0),
         Total_Students10 = TotalEnrollment16/10,
         ExpensesPerStudent = scale(Expenses_Per_Student/1000, center = TRUE, scale = FALSE)[,],   # Don't end up using this
         Median_Income = scale(medianincome/1000,center = TRUE, scale = FALSE)[,],                 # Don't end up using this
         SRO = ifelse(Has_SRO_Perc > .5 | SCH_FTESECURITY_LEO >= .2, 1, 0) %>% replace_na(., 0),     # Identifies TRUE cases, but leaves everything else NA, so the NAs are converted to 0 which are dealt with later
         Guards = ifelse(SCH_FTESECURITY_GUA < .2, 0, 1),                                         # .2 used as cutoff as this indicates present at least 1 day per week.
         Counselor = ifelse(SCH_FTECOUNSELORS < .2, 0, 1),
         Nurse = ifelse(SCH_FTESERVICES_NUR < .2, 0, 1),
         Psychologist = ifelse(SCH_FTESERVICES_PSY < .2, 0, 1),
         Social_Worker = ifelse(SCH_FTESERVICES_SOC < .2, 0, 1),
         Attacks_Per_Student = Incidents_Attack/Total_Enrollment,        # Not going to use these Per_Student based on unreliabiliy of reporting in CRDC
         Threats_Per_Student = Incidents_Threats/Total_Enrollment,
         Bullying_Per_Student = Incidents_Bullying/Total_Enrollment) %>%
  mutate(SRO = ifelse(rowSums(is.na(select(.,Has_SRO_Perc, SCH_FTESECURITY_LEO))) == 2, NA, SRO)) %>%   # Changes 0 to NA for rows were both items had NA)
  left_join(SLAP_Proficiency %>% select(level2,Grade_Levels,Proficient_Math,Proficient_Read),
            by = c("level2","Grade_Levels"))


# table(SL16$Grade_Levels,SL16$Counselor,useNA = "always") #%>% margin.table(margin=2)
# table(SL16$Grade_Levels, useNA = "always")
# table(SL16[SL16$Grade_Levels!="5",]$SRO, useNA = "always")
# table(SL16[SL16$Grade_Levels!="5",]$SRO,SL16[SL16$Grade_Levels!="5",]$SCH_FTESECURITY_LEO < .2,useNA = "always")





#### Joining School Info to Student data ####
## and dropping unnecessary MCA variables
StuSch16 <- StudentVars %>%
  left_join(SL16, by = "level2") %>%
  select(-starts_with("pLevel"),-starts_with("mean."),-starts_with("sd."))
# check <- StudentVars %>% anti_join(SL16, by = "level2") # 0, thus all students were joined to a school

#### Saving Preliminary Data ####
save(SL16, StuSch16, AggVars,
     file = "Saved Applied Results/AppliedPrelimData.RData")

#######################################################

#######################################################
#### Covariate Selection and Missing Data Analysis ####

## Loading data if necessary
load("Saved Applied Results/AppliedPrelimData.RData")

#### Covariate Selection ####
## Student-level covariates in the PS model
# Although some covariates might be useful, they were only asked of certain grades (Unsupervised = 5; LGB = HS)
# Others have high missingness rates (OST_Experiences)
# Or have no correlation with SRO and Outcomes (PhyEd)
# Or too highly correlated with each other (Parentjail, Homeless used to calculate Trauma; skip school & class)
OutcomeVars <- c("GPA","Commitment_to_Learning","Empowerment","Teacher.School_Support")
StudentPSVarNames <- AggVars[!(AggVars %in% c(OutcomeVars,"Has_SRO","Unsupervised","LGB","Fifth",
                                              "OST_Experiences","PhyEd","Skipped_School","Bullying",
                                              "Parentjail","Homeless", "Family_Violence"))]

## School-level covariates used in the PS model
# Although some covariates might be useful, they have high rates of missingness (listed after #)[OSTExperience was dropped at L1, but retained at L2]
# Or are too highly correlated with others (Parenjail,Homeless, and FV used to calculate Trauma; Diversity and % White;
#    Skipped school w/ class[more variability], Bullied [not necessarily at school so may not be affected by SRO]w/ Bullying )
# Or are uncorrelated with SRO and outcomes (Social Worker, ExpensesPerStudent, Prop_FTE_Under2, PhyEd)
### Need to decide based on missingnes whether to use MSS, MDE or CRDC in some cases
SchoolPSVarNames <- SL16 %>%
  select(Total_Students10,Fifth_L2,Eighth_L2,HS_L2,Charter_Magnet,TC,                 # Characteristics; SRO_Discomfort
         p16.Female:p16.ELL,                                                         # Composition; LGB_Perc,Parentjail_Perc,p16.Homeless,Diversity.16
         Sports_Team_Perc:Sleep_8_Hrs_Perc,-PhyEd_Perc,                              # OST and Health/Diet ; Unsupervised_Perc, PhyEd_Perc
         Prop_Truant,Skipped_Class_Perc,Gambling_Perc:Theft_Perc,                    # Truancy and Crime; Skipped_School_Perc, Prop_Truant is defined as absent 15 or more times per school year
         Disciplined_Perc,Prop_InSchool_Suspensions,Prop_OutSchool_Suspensions,      # Discipline; Prop_Arrest (missing SPPS schools), Prop_Expulsions & Prop_Law_Referral (low variation)
         Smoked_Cigarette_Perc:Took_Rx_Perc,                                         # Drug Use
         Trauma1_Perc,Bullied_Perc,Mental_Distress_Perc,                             # Bullying & Challenges; Family_Violence_Perc, Bullying_Perc, Attacks_Per_Student,Threats_Per_Student,Bullying_Per_Student (I don't trust the reporting quality of per student measures)
         Positive_Identity_Perc,Social_Competence_Perc,Family.Community_Support_Perc, # Assets
         College_Perc,Proficient_Math,Proficient_Read,                                # Academics
         Prop_FTE_Absent,Student_Teacher_Ratio,Guards,Counselor,Nurse,Psychologist    # Resources ; Median_Income, Social_Worker,ExpensesPerStudent,Prop_FTE_Under_2yrs,
  ) %>% names()

#### Correlation of covariates with outcomes ####
## Evaluating whether covariates are related to outcome, treatment, and/or redundant with each other
# Kainz (and others) recommend selecting covariates on a theoretical basis rather than sample dependent empirical methods
# However, I have a lot of potential variables, some of which are closely related, and models take forever to run and/or won't converge
CovOutCorrs <- StuSch16 %>%
  filter(Grade_Levels!="5" & Grade!="5") %>%  # Removing grade 5 students and schools who won't be included in the analysis
  select(SRO, one_of(OutcomeVars), one_of(StudentPSVarNames), one_of(SchoolPSVarNames)) %>%
  cor(use = "pairwise.complete.obs") %>% round(.,3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")
# write.csv(CovOutCorrs, "CovOutCorrs.csv")

## Uncorrelated with SRO and Outcomes
Uncorr <- CovOutCorrs %>% select(Variable:Teacher.School_Support) %>%
  filter_at(vars(SRO:Teacher.School_Support), all_vars(abs(.) < .10))
# Mostly L2 covariates. Some with really low % (and probably little variability)

## Investigating Multicollinearity
MulticollinearityCheck <- CovOutCorrs %>% select(Variable, one_of(StudentPSVarNames,SchoolPSVarNames)) %>%
  filter(Variable %in% c(StudentPSVarNames,SchoolPSVarNames)) %>%
  filter_at(vars(-Variable), any_vars(abs(.) > .75 & abs(.) < 1)) %>%
  mutate_if(is.numeric, list(~ifelse(abs(.) < .75, NA, .))) %>%
  select_if(function(x){!all(is.na(x))})

## Correlation Notes:
# SRO not correlated with any L1 covariates; largest was PhyEd with -.09
# Some L2 OST covariates uncorrelated with both SRO and outcomes, but others both. Theory would say keep all, empirically says drop em
# The overall rate and variability for some L2 covariates might be too low even if they are theoretically relevant (e.g. % Homeless); Exclude, but make a note in text
# Outcomes tended to be correlated with same covariates
# Outcomes also tended to have stronger correlations with L1 covariates than L2
# ExpensesPerStudent uncorrelated with all; a bit surprising, but I also am not overly confident in the reliability of the measure
# Multicollinearity not a problem with L1 covariates
# Some high correlations at L2 (> .80)


#### Level 2 Missing data investigation ####
## Shortcut dataset
# Removing schools that only serve grade 5 and those missing MCA scores 
# not many consistencies as to why the 28 schools are missing data on one or both MCAs
SL16no5 <- SL16 %>% filter(Grade_Levels!="5") %>%
  filter(!(is.na(Proficient_Read)|is.na(Proficient_Math))) %>%
  filter(!(is.na(Nurse))) # Any CRDC covariate could be used instead of Nurse

## Missing summary table
L2missing <- SL16no5 %>%
  select(one_of(SchoolPSVarNames),SRO,-starts_with("Proficient")) %>% #,SRO_Discomfort_Perc
  md.pattern(plot=FALSE) %>% as.data.frame() %>% tibble::rownames_to_column(.,var="Case_Count") %>%
  mutate(Case_Count = str_remove(Case_Count,"X") %>% str_remove("\\..*"),
         Case_Count = ifelse(Case_Count=="","Total Missing",Case_Count))
names(L2missing)[[length(L2missing)]] <- "Total_Missing"
# summarize_all(~sum(is.na(.))) %>%
# gather(Variable,Count_Missing)

## Descriptives across all students and by Missing status
# I don't actually reference this table
L2MissingDescrips <- SL16 %>% filter(Grade_Levels!="5") %>%
  mutate(p16.FRL = ifelse(is.na(p16.FRL), FRL_Perc, p16.FRL),                       # Replacing missing MDE values with CRDC or MSS values
         p16.Homeless = ifelse(is.na(p16.Homeless), Homeless_Perc, p16.Homeless),
         p16.SPED = ifelse(is.na(p16.SPED), Prop_Enrollment_IDEA, p16.SPED),
         p16.ELL = ifelse(is.na(p16.ELL), Prop_Enrollment_LEP, p16.ELL)) %>%
  select(level2,SRO,one_of(SchoolPSVarNames)) %>%
  mutate(NumMissing = rowSums(is.na(.)),
         Missing = ifelse(NumMissing > 0, 1, 0)) %>%
  Compare_Missing_Shortcut(VarsToCompare = c("SRO",SchoolPSVarNames), AllCont = TRUE) %>%
  mutate(Std_Diff = abs(Std_Diff))


## L2 Notes: 1081 total schools, but 528 are 5th grade only; 553 serve target grades
# drops to 525 after removing schools without MCA scores
# 461 schools have complete data
# 31 are missing p16.FRL (MDE), 16 of which are missing only that; All 31 have values for FRL_Perc (MSS)
# 22 missing p16.ELL (MDE), 10 missing only that; All 22 have values for Prop_Enrollment_LEP (CRDC)
# All missing p16 variables have alternatives
# The remaining 2 missing CRDC variables (46 and 40 students)
# After replacing values when possible, 30 of the 553 (5.4%) of schools had missing data
# Schools with missing data tended to have fewer SROs but more Guards, have HS or 5th grade, have more experienced and less absent teachers,
# be in the TC (63 to 35%), be more diverse (e.g. %Wh: 58 to 75%), be a charter/mag (32 to 10%), smaller, higher expenses/student,
# generally higher misbehavior and drug use, higher FRL and SPED populations
# Only 9 schools listed as a traditional school, the others are Charter schools or part of Special, Intermediate or Tribal districts
# mice requires separate imputation routines at each level; imputation at L2 doesn't seem worth it.


#### Level 1 Missing data investigation ####
## Shortcut dataset
ShortcutL1dat <- StuSch16 %>%
  filter(Grade_Levels!="5" & Grade!="5") %>%
  mutate(p16.FRL = ifelse(is.na(p16.FRL), FRL_Perc, p16.FRL),                       # Replacing missing MDE values with CRDC or MSS values
         p16.Homeless = ifelse(is.na(p16.Homeless), Homeless_Perc, p16.Homeless),
         p16.SPED = ifelse(is.na(p16.SPED), Prop_Enrollment_IDEA, p16.SPED),
         p16.ELL = ifelse(is.na(p16.ELL), Prop_Enrollment_LEP, p16.ELL)) %>%       #126868 in 553; length(unique(ShortcutL1dat$level2))
  filter_at(vars(one_of(SchoolPSVarNames)), all_vars(!is.na(.)))                   #124977 in 523 


## Missing summary table
L1missing <- ShortcutL1dat %>% 
  select(one_of(OutcomeVars,StudentPSVarNames)) %>%
  md.pattern(plot = FALSE) %>% as.data.frame() %>%
  tibble::rownames_to_column(., var = "Case_Count") %>%
  mutate(Case_Count = str_remove(Case_Count,"X") %>% str_remove("\\..*"),
         Case_Count = ifelse(Case_Count=="","Total Missing",Case_Count))
names(L1missing)[[length(L1missing)]] <- "Total_Missing"

## % missing by variable
L1VarMissingness <- L1missing %>%
  filter(Case_Count=="Total Missing") %>%
  select(-Case_Count, Total_Missing) %>%
  gather(Variable, Count_Missing) %>%
  mutate(Percent_Missing = Count_Missing / nrow(ShortcutL1dat))

# ## Further investigation of missing cases
# L1missing %>% select(Case_Count,Skipped_School,Latinx:Total_Missing) %>% 
#   filter(Case_Count!="Total Missing") %>%
#   mutate(Case_Count = as.numeric(Case_Count)) %>%
#   filter(Case_Count > 20 & Total_Missing==2) %>% View()

## Descriptives across all students and by Missing status
L1MissingDescrips <- ShortcutL1dat %>%
  select(studentnumber,level2,SRO,one_of(OutcomeVars,StudentPSVarNames,SchoolPSVarNames)) %>%
  mutate(NumMissing = rowSums(is.na(.)),
         Missing = ifelse(NumMissing > 0, 1, 0)) %>%
  Compare_Missing_Shortcut(VarsToCompare = c("SRO","TC",OutcomeVars,StudentPSVarNames), AllCont = FALSE) %>%
  mutate(Std_Diff = abs(Std_Diff))

## Summarizing the criterion by which schools met SRO designation
Temp <- ShortcutL1dat %>% select(level2,SRO,Has_SRO_Perc,SCH_FTESECURITY_LEO,one_of(SchoolPSVarNames)) %>%
  unique() %>% drop_na()
table(Temp$Has_SRO_Perc > .5,Temp$SCH_FTESECURITY_LEO >= .2, useNA = "always") %>% margin.table(1)

## Investigating MCAR - Predicting Outcome missingness from covariate missingness
BinaryMissing <- ShortcutL1dat %>%
  select(one_of(OutcomeVars,StudentPSVarNames)) %>%
  mutate_all(~ifelse(is.na(.), 1, 0))              # missing value = 1; nonmissing = 0
# Does missingness in covariate predict missingness in outcome
mcar.mods <- map(OutcomeVars,~glm(formula(paste0(.x," ~ 1 + ",paste(StudentPSVarNames,collapse =" + "))), data = BinaryMissing)) %>%
  set_names(nm = OutcomeVars)
mcarModSummary <- map_dfr(mcar.mods,~broom::tidy(.x) %>% select(term,estimate,p.value),.id="Outcome") %>%
  mutate_if(is.numeric,~round(.,3))
mcarModSummary %>% select(-estimate) %>% spread(Outcome,p.value) %>%# View() 
  summarize_at(vars(-term),~sum(. < .05))
# summary(mcar.mods$Commitment_to_Learning)


## L1 Notes: 97805 (77.2%) have complete data (24% with missing)
# Bullying had the most missing with 10611 (8.5% of cases)
# 240017 of 7127223 (3.4%) total cells were missing
# missingness on outcomes: GPA (1.6%), CtL (1.9%), Emp (6.7%), TSS (5.8%); mean(BinaryMissing$Teacher.School_Support)
# Out of 43 predictors (51 total, but some had no missingness), num of predictors sig related to outcome based on missingness (1/0):  GPA (30), CtL (36), Emp (24), TSS (34)
# Conclusions: Data not MCAR. Not sure how to determine if its MAR or NMAR

##############################################

####################################
#### Centering and Final Sample ####

#### Complete case sample ####
CompleteSamp16 <- ShortcutL1dat %>%
  select(studentnumber,level2,one_of(OutcomeVars),one_of(StudentPSVarNames),SRO,one_of(SchoolPSVarNames)) %>%
  drop_na() # 98251 in 523  

CompleteSamp16L2Only <- CompleteSamp16 %>%
  select(level2,SRO,one_of(SchoolPSVarNames)) %>% unique()


#### Centering and Scaling ####

## Group mean (Within cluster) centering SEL measures
CompleteSamp16 <- CompleteSamp16 %>%
  group_by(level2) %>%
  mutate_at(vars(Bullied, Mental_Distress, Positive_Identity, Social_Competence, Family.Community_Support),
            list(Center = ~scale(., center = TRUE, scale = FALSE)[,])) %>%
  ungroup()

## Grand mean centering SEL measures, standardizing total students and student-teacher ratio
CompleteSamp16L2Only <- CompleteSamp16L2Only %>%
  mutate_at(vars(Total_Students10,Student_Teacher_Ratio), list(Scaled = ~scale(., center = TRUE, scale = TRUE)[,])) %>%
  mutate_at(vars(OST_Experiences_Perc,Bullied_Perc,Mental_Distress_Perc,
                 Positive_Identity_Perc,Social_Competence_Perc,Family.Community_Support_Perc),
            list(Center = ~scale(., center = TRUE, scale = FALSE)[,]))


## Joining scaled and centered L2 variables
FinalSample <- CompleteSamp16 %>%
  left_join(CompleteSamp16L2Only %>% select(level2, ends_with("_Scaled"), ends_with("_Center")), by = "level2")


#### Defining covariates for use in PS and outcome models ####
## Student level
StudentCovNames <- FinalSample %>%
  select(Eighth:Took_Rx,Bullied_Center:Family.Community_Support_Center,College) %>%
  names()

## School level
SchoolCovNames <- FinalSample %>%
  select(Total_Students10_Scaled,Fifth_L2:Religious_Activity_Perc,
         OST_Experiences_Perc_Center,Safe_Travel_Perc:Trauma1_Perc,
         Bullied_Perc_Center:Family.Community_Support_Perc_Center,
         College_Perc:Prop_FTE_Absent,Student_Teacher_Ratio_Scaled,
         Guards,Counselor,Nurse,Psychologist) %>%
  names()

#### Saving Missing Data Investigations and Final Dataset ####
save(SL16no5, ShortcutL1dat, CompleteSamp16, CompleteSamp16L2Only, FinalSample,          # Intermediate and final datasets
     OutcomeVars, StudentPSVarNames, SchoolPSVarNames, StudentCovNames, SchoolCovNames,  # Vectors of variable names
     L1missing, L2missing, L1MissingDescrips, L2MissingDescrips, mcarModSummary,         # missing data investigations
     CovOutCorrs, file = "Saved Applied Results/CompleteData.RData")

#######################################################################


#############################
#### Multiple Imputation ####

# load("Saved Applied Results/FinalData_and_Descrips.RData")
# 
# ## Removing Missing colums b/c I don't need them anymore
# OrigSamp16 <- OrigSamp16 %>% select(-NumMissing,-Missing)
# 
# ## Identifying groups of variables with similar needs for imputation
# L1ContVars <- c(OutcomeVars,"Bullied","Bullying","Mental_Distress","Family_Violence",  # Continuous L1 variables
#               "Positive_Identity","Social_Competence","Family.Community_Support")      # L1 variables w/ no missing values
# L1VarsNoMissing <- c("Eighth", "Ninth", "Eleven", "College")
# RaceEthVars <- OrigSamp16 %>% select(American_Indian:Latinx) %>% names()    # Vector of race variables
# 
# 
# ## Defining imputation method for each variable
# MIMethods <- data.frame(Variable = names(OrigSamp16)) %>%
#   mutate(MI_Method = ifelse(Variable %in% c(OutcomeVars,StudentPSVarNames[!(StudentPSVarNames %in% L1VarsNoMissing)]), "2l.jomo",""))
# 
#   # mutate(MI_Method = ifelse(Variable %in% L1ContVars, "2l.lmer",
#   #                           ifelse(Variable %in% StudentPSVarNames[!(StudentPSVarNames %in% c(L1ContVars,L1VarsNoMissing))], "2l.bin","")))
# 
# 
# set.seed(845972)
# ## Preliminary imputation just to get the prediction matrix, then making changes as necessary 
# predmat0 <- make.predictorMatrix(OrigSamp16)
# predmat0[,c("Eleven","studentnumber")] <- 0  # Removing some redundant variables from imputation model
# predmat0[RaceEthVars,RaceEthVars] <- 0
# predmat0[,"level2"] <- -2               # signifying level2 as cluster indicator
# # predmat0[c("studentumber","level2",L1VarsNoMissing,"SRO",SchoolPSVarNames), ] <- 0 # Setting rows to 0 where no values need imputing; I don't think this is necessary, but just in case
# diag(predmat0) <- 0                     # setting diagonals to 0
# ## Running the multiple imputation
# tic()
# TheMI <- mice(data = OrigSamp16, m = 1, method = MIMethods$MI_Method, pred = predmat0, maxit=20)
# toc()
# 
# save(TheMI,"Saved Applied Results/TrialMI.RData")
# # When keeping binary variables as numeric...
# 
# # When changing binary variables to characters...
# 
#   
#   
# 
# ## Creating the dataset that will be used in the multiple imputation
# soccapMI <- soccapOrdinal %>%
#   mutate_at(vars(Boy,English_Language_Learner,Hmong,Born_in_the_US,
#                  Live_with_Both_Parents,Best_Friend),factor)
# 
# ## (Semi-)Automatic way of determining the method of imputation for each variable
# MIMethod <- data.frame(Class = map_chr(soccapMI,class),
#                        Levels= map_int(soccapMI,~length(levels(.))),
#                        Number_Missing = map_int(soccapMI,~sum(is.na(.)))) %>%
#   rownames_to_column("Variable") %>%
#   mutate(Method = ifelse(Number_Missing==0,"",
#                          ifelse(Class=="numeric","norm",
#                                 ifelse(Class=="factor"&Levels==2,"logreg",
#                                        ifelse(Class=="factor"&Levels>2,"pmm","")))))
# ## Running the multiple imputation
# # scmi <- mice(data=soccapMI,m=60,method=MIMethod$Method,maxit=30)


################################################################
