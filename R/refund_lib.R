#refund library
# Line 26 and Line 62 needs to be changed for the data (-13.38757, 13.38757,)
#######################################################
## This version is for the study having exact x_i(t) ## 
#######################################################
###################
## final version ##
###################

create_predictors = function(covariates = NULL, funcs, kz = 30, kb = 30,smooth.cov = FALSE)
{
  require(mgcv)
  kb = min(kz, kb)
  p = ifelse(is.null(covariates), 0, dim(covariates)[2])
  if (is.matrix(funcs)) {
    n = nrow(funcs)
    Funcs = list(length = 1)
    Funcs[[1]] = funcs
  } else {
    n = nrow(funcs[[1]])
    Funcs = funcs
  }
  N.Pred = length(Funcs)
  t = phi = psi = CJ = list(length = N.Pred)
  for (i in 1:N.Pred) {
    t[[i]] = seq(-13.38757, 13.38757, length = dim(Funcs[[i]])[2])# 0 to 1 before
    N_obs = length(t[[i]])
    #meanFunc = apply(Funcs[[i]], 2, mean, na.rm = TRUE)
    #resd = sapply(1:length(t[[i]]), function(u) Funcs[[i]][,u] - meanFunc[u])
    #Funcs[[i]] = resd
    NormFunc = t(apply(Funcs[[i]], 1, function(u) u/sum(u)))
    Funcs[[i]] = NormFunc
    num = kb - 2
    qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots <- quantile(t[[i]], qtiles)
    phi[[i]] = cbind(1, t[[i]], sapply(knots, function(k) ((t[[i]] - k > 0) * (t[[i]] - k))))
    CJ[[i]] = Funcs[[i]] %*% phi[[i]]
  }
  X = cbind(rep(1, n), covariates)
  for (i in 1:N.Pred) {
    X = cbind(X, CJ[[i]])
  }
  return(X)
}


new_pfr = function (Y, covariates = NULL, funcs,  kb = 30, smooth.cov = FALSE, family = "gaussian", method = "REML", ...) 
{
  require(mgcv)
 
  n = length(Y)
  p = ifelse(is.null(covariates), 0, dim(covariates)[2])
  if (is.matrix(funcs)) {
    Funcs = list(length = 1)
    Funcs[[1]] = funcs
  } else {
    Funcs = funcs
  }
  N.Pred = length(Funcs)
  t = phi = psi = CJ = list(length = N.Pred)
  for (i in 1:N.Pred) {
    t[[i]] = seq(-13.38757, 13.38757, length = dim(Funcs[[i]])[2])
    # equally space is not the best way to do this
    N_obs = length(t[[i]])
    #meanFunc = apply(Funcs[[i]], 2, mean, na.rm = TRUE)
    #resd = sapply(1:length(t[[i]]), function(u) Funcs[[i]][,u] - meanFunc[u])
    #Funcs[[i]] = resd
    NormFunc = t(apply(Funcs[[i]], 1, function(u) u/sum(u)))
    Funcs[[i]] = NormFunc
    num = kb - 2
    qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots <- quantile(t[[i]], qtiles)
    phi[[i]] = cbind(1, t[[i]], sapply(knots, function(k) ((t[[i]] - k > 0) * (t[[i]] - k))))
    CJ[[i]] = Funcs[[i]] %*% phi[[i]]
  }
  X = cbind(rep(1, n), covariates)
  for (i in 1:N.Pred) {
    X = cbind(X, CJ[[i]])
  }
  D = list(length = N.Pred)
  for (i in 1:(N.Pred)) {
    D[[i]] = diag(c(rep(0, 1 + p), rep(0, kb * (i - 1)), c(rep(0, 2), rep(1, kb - 2)), rep(0, kb * (N.Pred - i))))
  }
  fit = gam(Y ~ X-1, paraPen = list(X = D), family = family, method = "REML", ...)
  coefs = fit$coef

  fitted.vals <- as.matrix(X) %*% coefs
  beta.covariates = coefs[1:(p + 1)]
  BetaHat = varBeta = varBetaHat = Bounds = list(length(N.Pred))
  for (i in 1:N.Pred) {
    BetaHat[[i]] = phi[[i]] %*% coefs[(2 + p + kb * (i - 1)):(1 + p + kb * (i))]
    varBeta[[i]] = fit$Vp[(2 + p + kb * (i - 1)):(1 + p + kb * (i)), (2 + p + kb * (i - 1)):(1 + p + kb * (i))]
    varBetaHat[[i]] = phi[[i]] %*% varBeta[[i]] %*% t(phi[[i]])
    Bounds[[i]] = cbind(BetaHat[[i]] + 1.96 * (sqrt(diag(varBetaHat[[i]]))), 
                        BetaHat[[i]] - 1.96 * (sqrt(diag(varBetaHat[[i]]))))
  }
  ret <- list(fit, fitted.vals, BetaHat, beta.covariates, X, 
              phi, psi, varBetaHat, Bounds,coefs)
  names(ret) <- c("fit", "fitted.vals", "BetaHat", "beta.covariates", 
                  "X", "phi", "psi", "varBetaHat", "Bounds","coefs")
  return(ret)
}
