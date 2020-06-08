###############################
#                             #
#  SROs and Student Outcomes  #
#    Packages and Functions   #
#                             #
###############################

#### Packages ####
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)
library(MatchIt)
library(cobalt)
library(ggplot2)
library(tictoc)
library(lme4)
library(broom.mixed)
library(mice)
# library(micemd)
# library(MatchIt.mice)

## Shortcut for changing covariate names to version for tables and figures
PrintableCovariatesShortcut <- function(vector){
  vector %>% str_replace("p16.", "% ") %>% ifelse(str_detect(., "_Perc"), paste("%", .), .) %>%
    str_remove("_Center") %>% str_remove("_Scaled") %>% str_remove("_Perc") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\.", "/")
}


#### Computing of correlation matrix ####
## adapted from http://www.sthda.com/english/wiki/correlation-matrix-an-r-function-to-do-all-you-need
rquery.cormat<-function(x,
                        type=c('lower', 'upper', 'full', 'flatten'),
                        digits=2,
                        usena="complete.obs",
                        graph=TRUE,
                        graphType=c("correlogram", "heatmap"),
                        col=NULL,
                        alpha = .05, ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, ...) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(
      c("#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
        "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  ## Correlation matrix
  cormat<-format(round(cor(x, use = usena, ...),digits),nsmall=digits)
  pmat<-format(round(cor.pmat(x, ...),digits),nsmall=digits)
  # Reorder correlation matrix
  # ord<-corrMatOrder(cormat, order="hclust")
  # cormat<-cormat[ord, ord]
  # pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  # sym<-symnum(cormat, abbr.colnames=FALSE)
  
  ## Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  
  ## Get lower/upper triangle or flatten
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
    return(list(r=cormat, p=pmat)) #,sym=sym
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    # sym=t(sym)
    return(list(r=cormat, p=pmat)) #,sym=sym
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat) %>%
      mutate(CorSig = paste0(cor,ifelse(as.numeric(as.character(p)) < alpha, "*", " ")))
    return(cormat)
  }
  
}

####################################
#### Propensity Score Functions ####

#### Tailoring bal.tab output for my purposes ####
# UnAdj refers to unadjusted differences vs adjusted (matched, weighted)
Get_BalanceTable <- function(btab, UnAdj = c("Un","Adj"), renameUnAdj = NULL){
  TheBT <- btab[[1]] %>%
    tibble::rownames_to_column("Covariate") %>%
    # mutate(Covariate = c(StudentCovNames, SchoolCovNames)) %>%
    select(Covariate, Type ,ends_with(UnAdj), -starts_with("KS"),-contains("Threshold"))
  
  if(!is.null(renameUnAdj)){
    TheBT <- TheBT %>%
      rename_all(~str_replace(., UnAdj, renameUnAdj))
  }
  return(TheBT)
}

#### Calculate Cohen's D from two independent or paired groups ####
## SD can be pooled by providing sd1 and sd2
Get_CohensD <- function(m1,m2,sd1,sd2=NULL,n1,n2,sample="ind",proportion=FALSE){
  
  
  # raw mean difference
  md <- (m1 - m2)
  
  # Sigma for continuous variables
  if(proportion==FALSE){
    
    # Use only SD from group 1 or an overall SD
    if(is.null(sd2)){
      
      sigma <- sd1
      
    } else {
      
      # sigma for independent groups
      if(sample=="ind"){
        
        sigmanum <- (n1 - 1) * (sd1^2) + (n2 - 1) * (sd2^2)
        sigmadenom <- (n1 + n2 - 2)
        sigma <- sqrt(sigmanum/sigmadenom)
        
      } else{ 
        
        # sigma for paired groups
        sigma <- sqrt((sd1^2 + sd2^2) / 2)
        
      }
    }
    # Sigma for dichotomous variables
  } else {
    
    # Can provide overall proportion to sd2
    if(!is.null(sd2)){
      
      sigma <- sqrt(sd2 * (1 - sd2)/ n1)
      
    } else if(sample=="ind"){
      
      # for unequal independent groups 
      sigma <- sqrt((m1 * (1 - m1) / n1) + (m2 * (1 - m2) / n2))
      
    } else {
      
      # for paired groups or groups of equal size
      sigma <- sqrt((m1 * (1 - m1) + m2 * (1 - m2)) / 2)
      
    }
  }
  
  # Calculating d
  d <- md/sigma
  
  return(d)
}

#### Getting summary descriptives on a continuous variable ####
## unquoted grouping variables should be specified in ...
# requires semTools for skew and kurtosis
Get_Descriptives <- function(data,ContinuousVariable,...,digits=5,AllContinuous=TRUE){
  
  groups <- quos(...)
  CV <- enquo(ContinuousVariable)
  
  data_descrip <- data %>% group_by(!!!groups) %>%
    summarize(n = sum(!is.na(!!CV)),
              Median = median(!!CV, na.rm=TRUE),
              Mean = mean(!!CV, na.rm=TRUE),
              SE = sd(!!CV, na.rm=TRUE)/sqrt(sum(!is.na(!!CV))),
              SD = sd(!!CV, na.rm=TRUE),
              Min = min(!!CV, na.rm=TRUE),
              Max = max(!!CV, na.rm=TRUE)) %>% ungroup()
  
  if(AllContinuous==FALSE){
    
    data_descrip <- data_descrip %>%
      mutate(SE = ifelse(Min==0 & Max==1,(Mean*(1-Mean))/sqrt(n),SE),
             SD = ifelse(Min==0 & Max==1,NA,SD))
  }
  
  data_descrip <- data_descrip %>%
    mutate_if(is.numeric,~round(.,digits=digits))
  
  return(data_descrip)
}

#### cobalt can't handle missing values so need to do comparison manually ####
# Requires Get_Descriptives and Get_CohensD
# Binary column indicating missingness must be called "Missing"
# Using SRO instead of Missing I confirmed that this code functions equivalently to bal.tab
Compare_Missing_Shortcut <- function(data,VarsToCompare, AllCont = FALSE){
  
  MissDescrips <- data %>%
    select(Missing,one_of(VarsToCompare)) %>%
    gather(Variable, Value, -Missing) %>%
    Get_Descriptives(Value, Missing, Variable, AllContinuous = AllCont) # Compare Means using absolute standardized difference
  
  DescripswMissing <- MissDescrips %>% filter(Missing == 0) %>% select(-Missing) %>%
    left_join(MissDescrips %>% filter(Missing == 1) %>% select(-Missing),by = "Variable", suffix = c("_Comp","_Miss")) %>%
    left_join(data %>% select(one_of(VarsToCompare)) %>%                                      # Descriptives across all observations
                gather(Variable, Value) %>%
                Get_Descriptives(Value, Variable, AllContinuous = AllCont), by = "Variable") %>%
    mutate(Std_Diff = ifelse(is.na(SD), Get_CohensD(m1 = Mean_Comp, m2 = Mean_Miss, sample = "paired", proportion = TRUE),
                             Get_CohensD(m1 = Mean_Comp, m2 = Mean_Miss, sd1 = SD_Comp, sd2 = SD_Miss, n1 = n_Comp, n2 = n_Miss)),
           Var_Ratio = SD_Comp / SD_Miss,
           Variable = factor(Variable,levels = c(VarsToCompare))) %>%
    select(Variable,starts_with("n"), starts_with("Mean"), starts_with("SD"), Std_Diff,Var_Ratio)
  
  return(DescripswMissing)
}


####################################################


#######################################
#### Multilevel Modeling Functions ####

#### Creates QQ, fitted vs residual, and fitted vs observed plots ####
ModelCheckPlots <- function(model,
                            glm.predict = "response", glm.resid = "deviance",
                            smooth_method = "loess", ...){
  
  library(ggplot2)
  
  ## Collects output list
  MCPlots <- list(QQ = NULL,
                  L2_QQ = NULL,
                  Fitted_and_Residual_Values = NULL,
                  Residual_Correlation = NULL,
                  Cooks_Distance = NULL,
                  ROC_Curve = NULL)
  
  ## Extracting residuals
  if("lm" %in% class(model)){
    
    ## Dependent variable name
    DV <- names(model$model)[[1]]
    
    ## fit information
    aug <- broom::augment(model, type.predict = glm.predict, type.residuals = glm.resid) #glm. arguments ignored for lm objects
    
    
  } else if(grepl("merMod", class(model))){
    
    ## Dependent variable name
    DV <- names(model@frame)[[1]]
    
    ## Level 1 fit information
    aug <- broom.mixed::augment(model, type.predict = glm.predict, type.residuals = glm.resid) #glm. arguments ignored for lmer objects
    
    ## Level 2 information and plots
    clusternm <- names(model@cnms)[[1]]
    L2aug <- broom.mixed::augment(lme4::ranef(model))
    L1L2 <- dplyr::left_join(dplyr::rename(aug, level = all_of(clusternm)), L2aug, by = "level")
    
    ## L1 and L2 Residual Correlation
    L1L2corr <- ggplot(L1L2, aes(x = estimate, y = .resid)) +
      geom_hline(aes(yintercept = 0), color = "grey", size = 1, linetype = 1) +
      geom_point(shape = 1) +
      geom_smooth(method = smooth_method, ...) +
      labs(x = "Level 2 Residual", y = "Level 1 Residual") +
      theme_bw()
    
    # Storing plot in output object
    MCPlots[["Residual_Correlation"]] <- L1L2corr
    
  } else {stop("model must be of class 'lm', 'glm', or '(g)lmerMod'")}
  
  ## Adding id column
  aug <- tibble::rowid_to_column(aug, ".id")
  
  ## Cook's d
  cookplot <- ggplot(aug, aes(x = .id, y = .cooksd)) +
    geom_bar(stat = "identity", position = "identity") +
    labs(x = "Observation", y = "Cook's distance") +
    theme_bw()
  
  # Storing plot in output object
  MCPlots[["Cooks_Distance"]] <- cookplot
  
  
  
  if("glm" %in% class(model) | class(model) == "glmerMod"){
    
    ## Fitted vs. Residuals grouped by DV
    aug[[DV]] <- factor(aug[[DV]]) 
    
    binned <- ggplot(aug, aes(x = .id, y = .resid)) + 
      geom_hline(aes(yintercept = 0), color = "grey", size = 1, linetype = 1) +
      geom_point(aes_string(color = DV)) +
      geom_smooth(aes_string(color = DV), ...) +
      labs(x = "Observation", y = paste(tools::toTitleCase(glm.resid), "Residual")) +
      theme_bw()
    
    ## ROC Analysis
    rocit <- pROC::roc(response = aug[[DV]], predictor = aug[[".fitted"]]) # glm.predict muste be "response"
    
    # Extracting ROC Results
    AUC <- round(rocit$auc[[1]], 3)
    senspec <- data.frame(Sensitivity = rocit$sensitivities,
                          Specificity = rocit$specificities)
    
    rocplot <- ggplot(senspec, aes(x = Specificity, y = Sensitivity)) + 
      geom_abline(slope = 1, intercept = 1, colour = "grey", size = 1) +
      geom_line(size = 1.5, color = "black") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
      scale_x_reverse(limits = c(1, 0), breaks = seq(1, 0, -.2)) +
      annotate(geom = "text", x = .3, y = .2, label = paste("AUC =", AUC)) +
      theme_bw()
    
    
    # Storing plots in output object
    MCPlots[["Fitted_and_Residual_Values"]] <- binned
    MCPlots[["ROC_Curve"]] <- rocplot
    
  } else if(class(model) == "lm" | class(model) == "lmerMod"){
    
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
    
    # Storing plots in output object
    MCPlots[["QQ"]] <- qq
    MCPlots[["Fitted_and_Residual_Values"]] <- frplot
    
  }
  
  if(class(model) == "lmerMod"){
    
    ## Q-Q Plot of L2 residuals
    L2qq <- ggplot(L2aug, aes(sample = estimate)) +  
      stat_qq(shape = 1) + stat_qq_line() +
      labs(x = "Theoretical Quantiles", y = "Level 2 Residual") +
      theme_bw()
    
    # Storing plot in output object
    MCPlots[["L2_QQ"]] <- L2qq
    
  }
  
  ## Removing NULL list elements 
  MCPlots <- MCPlots[lengths(MCPlots) != 0]
  
  return(MCPlots)
  
}

#### Variance Inflation Factor of single or multilevel model put into a dataframe ####
vif.table <- function (fit, modelname, multilevel = FALSE, digs = 3) {
  
  if(multilevel==TRUE){
    # Taken from here: https://jonlefcheck.net/2012/12/28/dealing-with-multicollinearity-using-variance-inflation-factors/,
    # which itself was adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(fixef(fit))
    
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
      v <- v[-(1:ns), -(1:ns), drop = FALSE]
      nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  } else{
    # v <- car::vif(fit)
  }
  
  v <- v  %>% as.data.frame() %>%
    tibble::rownames_to_column("Predictor") %>%
    mutate_if(is.numeric, ~round(., digits = digs))
  names(v) <- c("Predictor", modelname)
  
  return(v)
}

#### R2 ####
# produces the same as broom.mixed::tidy(., effects = "ran_pars", scales = "vcov"); uses formulas from Lorah (2018)
HLM_R2 <- function(unc.model, model){
  
  vc.unc <- VarCorr(unc.model)
  vc.mod <- VarCorr(model)
  
  R2 <- 1 - ((attr(vc.mod,"sc")^2 + vc.mod[[1]][[1]]) / (attr(vc.unc,"sc")^2 + vc.unc[[1]][[1]]))
  
  return(R2)
  
}


#### Calculate ICC from lmer or glmer object ####
# check performance::icc
Get_ICC <- function(model, type = c("1","2","DE"), logistic = FALSE){
  
  vc <- VarCorr(model)
  
  if(logistic == FALSE){
    residvar <- attr(vc,"sc")^2
  } else {
    residvar <- pi^2/3
  }
  
  ICC <- vc[[1]][[1]] / (vc[[1]][[1]] + residvar)
  
  if(type == "1"){
    return(ICC)
  } else {
    k <- getME(model,"n") / getME(model,"l_i")
    attributes(k) <- NULL
  }
  
  if(type == "2") {
    ICC2 <- (k*ICC) / (1 + (k - 1)*ICC)
    return(ICC2)
  } else if(type == "DE"){
    DE <- 1 + (k - 1)*ICC
    return(DE)
  }
}




#### Extracting fixed effects from single or multilevel, linear or logistic models
# requires broom, MuMIn
# Look into incorporating confint() rather than my own calculations, although so far I haven't found problems with my versions
Get_Fixed <- function(model,modelname,multilevel=FALSE,logistic=c("binary","multinomial","poisson"),
                      partialsd=TRUE,alpha=.05,digits=Inf,pnames=NULL){
  
  if(length(logistic)>1){
    logistic <- "linear"
  }
  
  ## critical value for significance
  critval <- qnorm((alpha/2),lower.tail = FALSE)
  
  ## Get unstandardized coefficients and standard errors and keep only fixed effects
  coefficients <- broom::tidy(model)
  # CI <- confint_tidy(model)
  
  if(multilevel==TRUE){
    coefficients <- coefficients %>% filter(group=="fixed")
  }
  
  if(logistic=="multinomial"){
    coefficients <- coefficients %>%
      mutate(estimate = log(estimate))
  }
  
  coefficients <- coefficients %>%
    mutate(Interval = std.error*critval,
           Lower = estimate-Interval,
           Upper = estimate+Interval,
           Sig = ifelse((abs(estimate)-abs(Interval))>0,"Yes","No")) %>%
    select(Predictor=term,Estimate=estimate,Std.Error=std.error,Interval:Sig,everything())
  
  if(logistic %in% c("binary","multinomial","poisson")){
    
    ## Add odds ratio and probability
    coefficients <- coefficients %>%
      mutate_at(vars(Estimate:Upper),list(OR=exp(.)))
    
    if(logistic=="multinomial"){
      coefficients <- coefficients %>% select(Level=y.level,everything())
    }
    
  } else {
    
    ## Standardized coefficients
    stdcoefs <- MuMIn::std.coef(model,partial.sd=partialsd) %>% data.frame() %>%
      rownames_to_column("term") %>%
      mutate(Interval_std = Std..Error.*critval,
             Lower_std = Estimate.-Interval_std,
             Upper_std = Estimate.-Interval_std,
             Sig_std = ifelse((abs(Estimate.)-abs(Interval_std))>0,"Yes","No")) %>%
      select(Predictor=term,Estimate_std=Estimate.,Std.Error_std=Std..Error.,ends_with("_std"),df)
    
    ## Joining unstandardized and standardized
    coefficients <- coefficients  %>%
      left_join(stdcoefs,by="Predictor")
  }
  
  ## Joining unstandardized and standardized
  coefficients <- coefficients  %>%
    mutate(Model=modelname) %>%
    select(Model,everything()) %>%
    mutate_if(is.numeric,~round(.,digits=digits))
  
  ## Changing Predictor names from variable codes to presentable labels
  if(!is.null(pnames)){
    coefficients <- coefficients %>%
      mutate(Predictor = pnames)
  }
  
  return(coefficients)
  
}

#### Calculating f2 for SRO fixed effect ####
Get_f2 <- function(SFit, outcome, sample){
  
  UNCmod <- "~ 1 + (1|level2)"
  Covmod <- paste0(" ~ 1 + ", paste(HLMvars[!(HLMvars %in% ExcludeCovs)], collapse = " + "), " + (1|level2)")
  
  UFit <- lmer(as.formula(paste(outcome, UNCmod)), data =  FinalSamp16wPS[FinalSamp16wPS[[sample]] != 0, ])
  CFit <- lmer(as.formula(paste(outcome, Covmod)), data =  FinalSamp16wPS[FinalSamp16wPS[[sample]] != 0, ])
  
  SR2 <- HLM_R2(UFit, SFit)
  CR2 <- HLM_R2(UFit, CFit)
  f2 <- ((SR2 - CR2) / (1 - SR2))
  
  f2df <- data.frame(Sample = sample, Outcome = outcome, C_R2 = CR2, S_R2 = SR2, f2 = f2)
  return(f2df)
}

###################################################

