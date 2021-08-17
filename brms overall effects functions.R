oeffectsbrms.zibb <- function (outcome, variables = NULL, N, zi, data, terms, case, type = "fe"){
  if (!requireNamespace("brms")) 
    install.packages("brms")
  library(brms)
  library(MASS)
  if (missing(data)) 
    stop("'data' should be specified")
  if (missing(terms)) 
    stop("'terms' should be specified")
  
  data$group = ifelse(data[, terms] == case, 1, 0)
  #data$group = model.matrix(~group, data = data)[,2]
  #formula = paste("y ~ group + offset(log(N))", zi, "shape~1", sep = ",")
  #count = cat(paste(y, zi, "shape~1", sep = ","), "\n")
  
  if (is.null(variables)){
    formula <- as.formula(
      paste(paste0(outcome, "|vint(", N, ")"),
           "group", 
            sep = " ~ "))
  } else {
    formula <- as.formula(
      paste(paste0(outcome, "|vint(", N, ")"),
            paste("group",  paste(variables, collapse = " + "), sep = " + "), 
            sep = " ~ "))
  }
  
  # set up formula and prior
  zero_inflated_beta_binomial2 <- custom_family(
    "zero_inflated_beta_binomial2", dpars = c("mu", "phi", "zi"),
    links = c("logit", "log", "logit"), lb = c(NA, 0, 0), ub = c(NA, NA, 1),
    type = "int", vars = "vint1[n]"
  )
  
  stan_funs <- "
  real zero_inflated_beta_binomial2_lpmf(int y, real mu, real phi, 
                                       real zi, int T) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                           beta_binomial_lpmf(0 | T, mu * phi, (1 - mu) * phi)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
        beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi); 
    } 
} 
"
  stanvars <- stanvar(scode = stan_funs, block = "functions")

  bform = bf(formula, zi)
  prior = set_prior("cauchy(0,2)", class="b")
  
  f = brm(bform, data=data, family=zero_inflated_beta_binomial2, 
          stanvars = stanvars,
          prior=prior, chains=2, iter=2000, warmup = 1000)
  
  newdata0 = data
  vars_df = unique(newdata0[, "group"])
  
  listom <- list()
  listos <- list()
  for (i in 1:length(vars_df)){
    new = newdata0[newdata0[,"group"]==vars_df[i],]
    X.cond = model.matrix(formula(f)[[1]][-2], new)
    betas = cbind(posterior_samples(f, pars = paste0("b_", gsub("\\(|\\)", "", colnames(X.cond))[1])),
                  posterior_samples(f, pars = paste0("b_", colnames(X.cond)[2]))) 
    lp.mu = as.matrix(X.cond %*% t(betas))         
    mu = t(exp(lp.mu))#apply(, 2, sum)
    lp.zi = posterior_linpred(f, newdata=new, dpar="zi", re.form = NA, re_formula = NA)
    zi = exp(lp.zi)/(1+exp(lp.zi)) 
    om = apply((1-zi)*mu, 1, mean)
    os = apply((1-zi)*mu, 1, sum)
    listom[[i]] <- om
    listos[[i]] <- os
  }
  
  bd = listom[[2]] - listom[[1]]
  br = listos[[2]]/listos[[1]]
  result <- list(bd, br)
  return(result)
}


oeffectsbrms.zinb <- function (outcome, variables = NULL, N, zi, data, terms, case, type = "fe"){
  if (!requireNamespace("brms")) 
    install.packages("brms")
  library(brms)
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
    
  bform = bf(formula, zi)#, "shape~1")
  prior = set_prior("cauchy(0,2)", class="b")
  
  f = brm(bform, data=data, family=zero_inflated_negbinomial(),
          prior=prior, chains=2, iter=2000, warmup = 1000)#, control = list(adapt_delta = 0.99)
  newdata0 = data
  vars_df = sort(unique(newdata0[, "group"]))
  
  listom <- list()
  listos <- list()
  for (i in 1:length(vars_df)){
    new = newdata0[newdata0[,"group"]==vars_df[i],]
    X.cond = model.matrix(formula(f)[[1]][-2], new)
    betas = cbind(posterior_samples(f, pars = paste0("b_", gsub("\\(|\\)", "", colnames(X.cond))[1])),
                  posterior_samples(f, pars = paste0("b_", colnames(X.cond)[2]))) 
    lp.mu = as.matrix(X.cond %*% t(betas))         
    mu = t(exp(lp.mu))#apply(, 2, sum)
    lp.zi = posterior_linpred(f, newdata=new, dpar="zi", re.form = NA, re_formula = NA)
    zi = exp(lp.zi)/(1+exp(lp.zi)) 
    om = apply((1-zi)*mu, 1, mean)
    os = apply((1-zi)*mu, 1, sum)
    listom[[i]] <- om
    listos[[i]] <- os
  }
  mean(listom[[1]])
  mean(listom[[2]])
  bd = listom[[2]] - listom[[1]]
  br = listos[[2]]/listos[[1]]
  result <- list(bd, br)
  return(result)
}



