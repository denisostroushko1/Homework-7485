# Use this to test the function inside ATE_80pct
# x.list.temp <- x.list[[1]]

ATE_fwrd <- function(x.list.temp) {
  X <- x.list.temp$X
  true.X <- X[,x.list.temp$true.cov]
  A <- x.list.temp$A
  Y <- x.list.temp$Y
  ps.df.org <- data.frame(Y, A, X)
  
  start_time <- Sys.time()
  ps.ml.st <- glm(A ~ 1, data = ps.df.org[,-1], family = "binomial")
  ps.ml.end <- formula(glm(A ~ ., data = ps.df.org[,-1], family = "binomial")) 
  fwrd.res <- step(ps.ml.st, direction = "forward", scope = ps.ml.end, trace = FALSE)
  fwrd.formula <- fwrd.res$formula
  
  ml.cov <- colnames(fwrd.res$model[,-1])
  ml.cov <- as.numeric(gsub("X", "", ml.cov))
  ml.ord.idx <- order(ml.cov)
  
  ps.ml <- glm(fwrd.formula, data = ps.df.org[,-1], family = "binomial")
  end_time <- Sys.time()
  execution_time <- end_time - start_time
  cat("Forward execution time:", as.numeric(execution_time, units = "secs"), "secs(s)\n")
  
  
  # Obtain PS and Stratified
  ps.df <- ps.df.org %>% 
    mutate(ps = predict(ps.ml, type = "response"),
           ps_quintile = cut(ps, 
                             breaks = c(0, quantile(ps, p = c(0.2, 0.4, 0.6, 0.8) + small_jitter), 1), 
                             labels = 1:5),
           w0 = (1-A)/(1-ps),
           w1 = A/ps)
  
  # estimated ATE
  n <- nrow(ps.df)
  nj <- table(ps.df$ps_quintile)
  te_quintile <- 
    tapply(ps.df$Y[ps.df$A == 1], ps.df$ps_quintile[ps.df$A == 1], mean) - 
    tapply(ps.df$Y[ps.df$A == 0], ps.df$ps_quintile[ps.df$A == 0], mean)
  ATE.pss <- sum(te_quintile *nj/n, na.rm = TRUE)
  
  # IPW
  ATE.ipw2 <- weighted.mean(ps.df$Y, ps.df$w1) - weighted.mean(ps.df$Y, ps.df$w0)
  
  ### Add bootstrapping in the future ###
  ATE.pss.bt <- matrix(NA, nrow=B, ncol=1)
  ATE.ipw2.bt <- matrix(NA, nrow=B, ncol=1)
  for(i in 1:B) {
    # Resampling
    ps.df.bt <- ps.df.org[sample(1:n, n, replace = TRUE), ]
    ps.ml.bt <- glm(fwrd.formula, data = ps.df.bt[,-1], family = "binomial")
    # Obtain PS and Stratified
    ps.df.bt <- ps.df.bt %>% 
      mutate(ps = predict(ps.ml.bt, type = "response"),
             ps_quintile = cut(ps, 
                               breaks = c(0, quantile(ps, p = c(0.2, 0.4, 0.6, 0.8) + small_jitter), 1), 
                               labels = 1:5),
             w0 = (1-A)/(1-ps),
             w1 = A/ps)
    
    # calculate treatment effect 
    nj <- table(ps.df.bt$ps_quintile)
    te_quintile_bt <- 
      tapply(ps.df.bt$Y[ps.df.bt$A == 1], ps.df.bt$ps_quintile[ps.df.bt$A == 1], mean) - 
      tapply(ps.df.bt$Y[ps.df.bt$A == 0], ps.df.bt$ps_quintile[ps.df.bt$A == 0], mean)
    
    ATE.pss.bt[i,] <- sum(te_quintile_bt * nj/n, na.rm = TRUE)
    
    ATE.ipw2.bt[i,] <- weighted.mean(ps.df.bt$Y, ps.df.bt$w1) - weighted.mean(ps.df.bt$Y, ps.df.bt$w0)
  }
  SE.ATE.pss <- sd(ATE.pss.bt)
  SE.ATE.ipw2 <- sd(ATE.ipw2.bt)
  # Setting up saving results
  overlap.cov <- sum(x.list.temp$cov.idx %in% ml.cov)
  miss.cov <- sum(!(x.list.temp$cov.idx %in% ml.cov))
  add.cov <- sum(!(ml.cov %in% x.list.temp$cov.idx))
  # res.df <- data.frame(true.cov = x.list.temp$cov.idx, ml.cov = which(rnd.cov),)
  # Return the ATE
  return(list(cov.info = data.frame(method = "forward selection",
                                    n = nrow(X),
                                    m = ncol(X),
                                    s = ncol(true.X),
                                    overlap.cov = overlap.cov,
                                    miss.cov = miss.cov,
                                    add.cov = add.cov,
                                    execution_time = execution_time),
              result = data.frame(method = "forward selection",
                                  n = nrow(X),
                                  m = ncol(X),
                                  s = ncol(true.X),
                                  ATE.true = x.list.temp$true.ATE,
                                  ATE.pss = ATE.pss,
                                  SE.pss = SE.ATE.pss,
                                  lwr.pss = ATE.pss - 1*qnorm(0.975)*SE.ATE.pss,
                                  upr.pss = ATE.pss + 1*qnorm(0.975)*SE.ATE.pss,
                                  ATE.ipw2 = ATE.ipw2,
                                  SE.ipw2 = SE.ATE.ipw2,
                                  lwr.ipw2 = ATE.ipw2 - 1*qnorm(0.975)*SE.ATE.ipw2,
                                  upr.ipw2 = ATE.ipw2 + 1*qnorm(0.975)*SE.ATE.ipw2),
              true.cov = x.list.temp$cov.idx,
              ml.cov = ml.cov))
}