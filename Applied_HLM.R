#######################################################
#                                                     #
#   Outcome HLM Model and Results for Applied Study   #
#                                                     #
#######################################################

## Loading packages, functions, and PS results
source("Applied_PF.R")
load("Saved Applied Results/PSbalance.RData")

################################
#### Running Outcome Models ####

## Identifying variables to use in the outcome regression
HLMvars <- c(StudentCovNames, SchoolCovNames)
HLMmod <- paste0(" ~ 1 + ", paste(HLMvars[!(HLMvars %in% ExcludeCovs)], collapse = " + "), "+ SRO + (1|level2)")

#### Commitment to Learning ####
tic()
CTLmods <- map(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
               ~lmer(as.formula(paste("Commitment_to_Learning", HLMmod)), data = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ],
                    weights = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ][[.x]])) %>%
  set_names(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"))
toc() # 31 sec

summary(test)
vif.table(test, "test", multilevel = TRUE)
testcheck <- ModelCheckPlots(test, smoother = "loess")

sjstats::omega_sq(test, partial = TRUE) #, ci.lvl = .95; adding the CI takes a hell of a lot longer

#### GPA ####
tic()
GPAmods <- map(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
               ~lmer(as.formula(paste("GPA", HLMmod)), data = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ],
                     weights = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ][[.x]])) %>%
  set_names(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"))
toc() # 24 sec

#### Empowerment ####
tic()
EMPmods <- map(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
               ~lmer(as.formula(paste("Empowerment", HLMmod)), data = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ],
                     weights = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ][[.x]])) %>%
  set_names(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"))
toc() # 21 sec

#### Teacher/School Support ####
tic()
TSSmods <- map(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"),
               ~lmer(as.formula(paste("Teacher.School_Support", HLMmod)), data = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ],
                     weights = FinalSamp16wPS[FinalSamp16wPS[[.x]] != 0, ][[.x]])) %>%
  set_names(c("IC_1to1","IC_1to2","CL_1to1","CL_1to2"))
toc() # 21 sec






############################################

############################
#### Extracting Results ####

summary(TSSmods[[1]])
(-0.014134) / sd(residuals(TSSmods[[1]])) # sd(residuals) is the RMSE/ sqrt(length(residuals(TSSmods[[1]]))) # cohensD according to https://stats.stackexchange.com/questions/171941/cohens-d-from-regression-coefficient

TSSmodcheckplots <- ModelCheckPlots(TSSmods[[1]])
TSSmodcheckplots$FittedResidualPlot
