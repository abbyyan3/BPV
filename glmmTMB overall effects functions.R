oeffectsglmmtmb.zibb <- function (outcome, variables = NULL, N, zi, data, terms, case, type = "fe"){
  if (!requireNamespace("glmmTMB")) 
    install.packages("glmmTMB")
  library(glmmTMB)
  library(MASS)
  if (missing(data)) 
    stop("'data' should be specified")
  if (missing(terms)) 
    stop("'terms' should be specified")
  
  data$group = ifelse(data[, terms] == case, 1, 0)
  #data$group = data[, terms]
  #data$group = model.matrix(~group, data = data)[,2]
  #formula = paste("y ~ group + offset(log(N))", zi, "shape~1", sep = ",")
  #count = cat(paste(y, zi, "shape~1", sep = ","), "\n")
  
  if (is.null(variables)){
    formula <- as.formula(
      paste(paste0(outcome, "/", N),
            "group", 
            sep = " ~ "))
  } else {
    formula <- as.formula(
      paste(paste0(outcome, "/", N),
            paste("group",  paste(variables, collapse = " + "), sep = " + "), 
            sep = " ~ "))
  }
  
  ziformula <- as.formula(zi)
  zinbm3 = glmmTMB(formula, ziformula, family=betabinomial(link = "logit"), data=data, weights = N)
  #zinbm3 = glmmTMB(y/N~group, zi~1, family=betabinomial(link = "logit"), data=data, weights = N)
  
  newdata0 = data
  vars_df = unique(newdata0[, "group"])
  X.cond = model.matrix(lme4::nobars(formula(zinbm3)[-2]), newdata0)
  X.cond.uniq = unique(X.cond)
  

    new = data.frame(X.cond.uniq)#[X.cond.uniq[,"group"]==vars_df[i],]
    beta.cond = fixef(zinbm3)$cond
    pred.cond = X.cond.uniq %*% beta.cond
    ziformula = zinbm3$modelInfo$allForm$ziformula #$
    X.zi = model.matrix(lme4::nobars(ziformula), new)
    beta.zi = fixef(zinbm3)$zi
    pred.zi = X.zi %*% beta.zi
    #pred.ucount = exp(pred.cond)*(1-plogis(pred.zi))
    #set.seed(101)
    sigma = vcov(zinbm3)$cond
    pred.condpar.psim = mvrnorm(2000,mu=beta.cond,Sigma=sigma)
    pred.cond.psim = X.cond.uniq %*% t(pred.condpar.psim) 
    pred.zipar.psim = mvrnorm(2000,mu=beta.zi,Sigma=vcov(zinbm3)$zi)
    pred.zi.psim = X.zi %*% t(pred.zipar.psim)
    pred.ucount.psim = exp(pred.cond.psim)*(1-plogis(pred.zi.psim))
    
    #om = apply(pred.ucount.psim, 1, mean)
    #os = apply(pred.ucount.psim, 1, sum)

  d = pred.ucount.psim[2, ] - pred.ucount.psim[1, ]#om[1] - om[2]
  r = pred.ucount.psim[2, ]/pred.ucount.psim[1, ]# os[1]/os[2]
  result <- list(d, r)
  return(result)

}

oeffectsglmmtmb.zinb <- function (outcome, variables = NULL, N, zi, data, terms, case, type = "fe"){
  if (!requireNamespace("glmmTMB")) 
    install.packages("glmmTMB")
  library(glmmTMB)
  library(MASS)
  if (missing(data)) 
    stop("'data' should be specified")
  if (missing(terms)) 
    stop("'terms' should be specified")
  
  data$group = ifelse(data[, terms] == case, 1, 0)
  #data$group = data[, terms]
  #data$group = model.matrix(~group, data = data)[,2]
  #formula = paste("y ~ group + offset(log(N))", zi, "shape~1", sep = ",")
  #count = cat(paste(y, zi, "shape~1", sep = ","), "\n")

  if (is.null(variables)){
    formula <- as.formula(
      paste(outcome, 
            paste("group", "offset(log(N))", sep = " + "), 
            sep = " ~ "))
  } else {
    formula <- as.formula(
      paste(outcome, 
            paste("group",  paste(variables, collapse = " + "), "offset(log(N))", sep = " + "), 
            sep = " ~ "))
  }
  
  
  ziformula <- as.formula(zi)
  zinbm3 = glmmTMB(formula, ziformula, family=nbinom2, data=data)
  
  newdata0 = data
  vars_df = unique(newdata0[, "group"])
  X.cond = model.matrix(lme4::nobars(formula(zinbm3)[-2]), newdata0)
  X.cond.uniq = unique(X.cond)
  
  new = data.frame(X.cond.uniq)#[X.cond.uniq[,"group"]==vars_df[i],]
  beta.cond = fixef(zinbm3)$cond
  pred.cond = X.cond.uniq %*% beta.cond
  ziformula = zinbm3$modelInfo$allForm$ziformula #$
  X.zi = model.matrix(lme4::nobars(ziformula), new)
  beta.zi = fixef(zinbm3)$zi
  pred.zi = X.zi %*% beta.zi
  #pred.ucount = exp(pred.cond)*(1-plogis(pred.zi))
  #set.seed(101)
  sigma = vcov(zinbm3)$cond
  pred.condpar.psim = mvrnorm(2000,mu=beta.cond,Sigma=sigma)
  pred.cond.psim = X.cond.uniq %*% t(pred.condpar.psim) 
  pred.zipar.psim = mvrnorm(2000,mu=beta.zi,Sigma=vcov(zinbm3)$zi)
  pred.zi.psim = X.zi %*% t(pred.zipar.psim)
  pred.ucount.psim = exp(pred.cond.psim)*(1-plogis(pred.zi.psim))
  
  #om = apply(pred.ucount.psim, 1, mean)
  #os = apply(pred.ucount.psim, 1, sum)
  
  d = pred.ucount.psim[2, ] - pred.ucount.psim[1, ]#om[1] - om[2]
  r = pred.ucount.psim[2, ]/pred.ucount.psim[1, ]# os[1]/os[2]
  result <- list(d, r)
  return(result)
  
}




