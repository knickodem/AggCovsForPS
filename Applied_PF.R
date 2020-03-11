################################################
#                                              #
#   Packages and Functions for Applied Study   #
#                                              #
################################################

#### Packages ####
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)
library(MatchIt)
library(WeightIt)
library(cobalt)
library(ggplot2)
library(tictoc)
library(mice)
library(lme4)
# library(micemd)
# library(MatchIt.mice)



###########################################
#### Descriptive and Generic Functions ####

#### Converting logits to probabilities and vice versa ####
LogitToProb <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

## Converting Probability to Logit
ProbToLogit <- function(prob){
  odds <- prob / (1 - prob)
  logit <- log(odds)
  return(logit)
}

#### Getting summary descriptives on a continuous variable ####
## unquoted grouping variables should be specified in ...
Get_Descriptives <- function(data,ContinuousVariable,...,digits=5,AllContinuous = FALSE){
  
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


#### Balance Plots ####
# balAllPlotGrey <- balAll %>% filter(!(Variable %in% c("distance"))) %>%
#   mutate(Variable = factor(Variable,levels=rev(levels(balOrig$Variable)))) %>%
#   rename(Sample=Model) %>%
#   ggplot(aes(x=Variable,y=Std_Diff,group=Sample)) +
#   geom_hline(yintercept = .1, linetype = 1, color = "gray5") +
#   geom_point(aes(shape=Sample,color=Sample),size=4,position=position_dodge(1)) +
#   ylab("Standardized Difference") + xlab("") +
#   scale_shape_manual(values=c(19,18,17,15,7,3)) +
#   scale_color_grey() +
#   #scale_color_brewer(palette="Dark2") +
#   coord_flip() +
#   theme_bw(base_size = 20) +
#   theme(panel.grid.major = element_line(color = "gray87"),
#         panel.grid.minor = element_line(color = "gray90"),
#         panel.background = element_rect(fill = "white", color = "black"),
#         axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black"),
#         legend.justification=c(1,0),legend.position=c(.9,.1))
# 
# # ggsave(balAllPlotGrey, file = "Scripts/Saved Objects/balAllPlotGrey.png",height = 10, width = 15)

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



#################################################################



#######################################
#### Multilevel Modeling Functions ####

#### Calculate ICC from lmer or glmer object ####
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


#### Variance Inflation Factor of single or multilevel model put into a dataframe ####
vif.table <- function (fit, modelname, multilevel = FALSE, digs = 3) {
  
  if(multilevel==TRUE){
    ## adapted from rms::vif
    
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
    v <- car::vif(fit)
  }
  
  v <- v  %>% as.data.frame() %>%
    tibble::rownames_to_column("Predictor") %>%
    mutate_if(is.numeric, ~round(., digits = digs))
  names(v) <- c("Predictor", modelname)
  
  return(v)
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


#### Creates QQ, fitted vs residual, and fitted vs observed plots ####
## works for lm, lmer, glm, glmer (I think)
## does not work yet for multivariate, but could
# Consider adding an option for using cowplot::plot_grid rather than saving as a list
ModelCheckPlots <- function(model, smoother = "loess",
                            observed = NULL, multinomial = FALSE, byImp = FALSE){
  
  
  if(multinomial == TRUE){
    
    Residual <- resid(model) %>% data.frame() %>%
      tibble::rownames_to_column("ID") %>% gather(Level, Residuals, -ID)
    
    FR <- fitted(model) %>% data.frame() %>%
      tibble::rownames_to_column("ID") %>% gather(Level, Fitted, -ID) %>%
      left_join(Residual, by = c("ID", "Level"))
    
  } else {
    
    FR <- data.frame(Fitted = fitted(model),
                     Residuals = resid(model)) %>%
      tibble::rownames_to_column("ID")
  }
  
  ## QQ Plot - appears I no longer need to do the calculation manually with the addition of stat_qq_line()
  # y <- quantile(FR$Residuals[!is.na(FR$Residuals)], c(0.25, 0.75))
  # x <- qnorm(c(0.25, 0.75))
  # slope <- diff(y)/diff(x)
  # int <- y[1L] - slope * x[1L]
  
  qq <- ggplot(FR, aes(sample = Residuals)) +
    stat_qq(shape = 1) + stat_qq_line() +
    labs(x = "Theoretical", y = "Observed") +
    theme_bw(base_size = 18) #+geom_abline(slope = slope, intercept = int, color="black")
  
  ## Fitted vs. Residuals Plot
  frplot <- FR %>%
    ggplot(aes(x = Fitted, y = Residuals)) +
    geom_point(shape = 1) +
    geom_hline(aes(yintercept = 0), color = "grey", size = 1, linetype = 1) +
    geom_smooth(method = smoother, se = FALSE, color = "black", size = 2, linetype = 1) +
    theme_bw(base_size = 18)
  
  if(multinomial==TRUE){
    qq <- qq + facet_wrap(~Level)
    frplot <- frplot + facet_wrap(~Level)
  }
  
  if(byImp==TRUE){
    qq <- qq + facet_wrap(~`.imp`,ncol=2) + theme(strip.background = element_blank())
    frplot <- frplot + facet_wrap(~`.imp`,ncol=2) + theme(strip.background = element_blank())
  }
  
  MCPlots <- list(QQPlot = qq,
                  FittedResidualPlot = frplot)
  
  
  if(!is.null(observed)){
    ## Fitted vs. Observed Plot
    afplot <- FR %>%
      ggplot(aes(x=Fitted,y=Observed)) +
      geom_point(shape=1) +
      geom_smooth(method="lm",se=FALSE,color="black",size=2,linetype=1) +
      theme_bw(base_size = 18)
    
    if(multinomial==TRUE){
      afplot <- afplot + facet_wrap(~Level)
    }
    
    if(byImp==TRUE){
      afplot <- afplot + facet_wrap(~`.imp`,ncol=2) + theme(strip.background = element_blank())
    }
    
    MCPlots <- list(QQPlot = qq,
                    FittedResidualPlot = frplot,
                    ActualFittedPlot = afplot)
  }
  
  if(class(model)=="lmerMod"){
    
    ## QQ Plot
    L2res <- data.frame(Residuals = ranef(model)[[1]][[1]])
    L2qq <- ggplot(L2res, aes(sample = Residuals)) +
      stat_qq(shape=1) + stat_qq_line() +
      labs(x = "Theoretical", y = "Observed") +
      theme_bw(base_size = 18)
    
    MCPlots <- list(QQPlot = qq,
                    L2QQPlot = L2qq,
                    FittedResidualPlot = frplot)
  }
  
  return(MCPlots)
}

###################################################


####################################
#### Propensity Score Functions ####

#### Tailoring bal.tab output for my purposes ####
# UnAdj refers to unadjusted differences vs adjusted (matched, weighted)
Get_BalanceTable <- function(btab,UnAdj=c("Un","Adj"),replace01 = NULL, ModName=NULL){
  
  TheBT <- btab[[1]] %>%
    tibble::rownames_to_column("Variable") %>%
    select(Variable,Type,ends_with(UnAdj),-starts_with("KS"),-contains("Threshold"))
  names(TheBT) <- str_remove(names(TheBT),UnAdj) %>%
    str_remove(.,"\\.$")
  
  # Replace 0 and 1 values in column names with in 2 element character vector
  if(!is.null(replace01)){
    
    names(TheBT) <- str_replace(names(TheBT),".0",replace01[[1]]) %>%
      str_replace(.,".1",replace01[[2]])
    
  }

  if(!is.null(ModName)){
  TheBT <- TheBT %>%
    mutate(Model = ModName) %>%
    select(Model,everything())
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
    mutate(Std_Diff = ifelse(is.na(SD), Mean_Comp - Mean_Miss,
                             Get_CohensD(m1 = Mean_Comp, m2 = Mean_Miss, sd1 = SD)),
           Var_Ratio = SD_Comp / SD_Miss,
           Variable = factor(Variable,levels=c(VarsToCompare))) %>%
    select(Variable,starts_with("n"),starts_with("Mean"),starts_with("SD"),Std_Diff,Var_Ratio)
  
  return(DescripswMissing)
}


############################################
