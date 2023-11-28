# Use this to test the function inside ATE_80pct
# x.list.temp <- x.list[[1]]

ATE_80pct <- function(x.list.temp, B) {
  X <- x.list.temp$X
  true.X <- X[,x.list.temp$true.cov]
  # randomly pick 80%
  # rnd.cov <- purrr::rbernoulli(ncol(true.X), p=0.8)
  
  rnd.cov.idx <- sample(1:ncol(true.X), ncol(true.X)*0.8)
  rnd.ord.idx <- order(rnd.cov.idx)
  # Use the indices to obtain the ordered vector
  rnd.cov.idx <- rnd.cov.idx[rnd.ord.idx]
  rnd.cov <- 1:ncol(true.X) %in% rnd.cov.idx
  ml.cov <- x.list.temp$cov.idx[which(rnd.cov)]
  
  ps.cov <- true.X[,rnd.cov]
  A <- x.list.temp$A
  Y <- x.list.temp$Y
  ps.df.org <- data.frame(A, ps.cov)
  ps.ml <- glm(A ~ ., data = ps.df.org, family = "binomial")
  
  # Obtain PS and Stratified
  small_jitter = 0
  ps.df <- ps.df.org %>% 
    mutate(Y = as.numeric(Y),
           ps = predict(ps.ml, type = "response"),
           ps_quintile = cut(ps, 
                             breaks = c(0, quantile(ps, p = c(0.2, 0.4, 0.6, 0.8 )+ small_jitter), 1), 
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
    ps.ml.bt <- glm(A ~ ., data = ps.df.bt, family = "binomial")
    # Obtain PS and Stratified
    ps.df.bt <- ps.df.bt %>% 
      mutate(Y = as.numeric(Y),
             ps = predict(ps.ml.bt, type = "response"),
             ps_quintile = cut(ps, 
                               breaks = c(0, quantile(ps, p = c(0.2, 0.4, 0.6, 0.8 + small_jitter)), 1), 
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
  return(list(cov.info = data.frame(method = "select 80%",
                                    n = nrow(X),
                                    m = ncol(X),
                                    s = ncol(true.X),
                                    overlap.cov = overlap.cov,
                                    miss.cov = miss.cov,
                                    add.cov = add.cov),
              result = data.frame(method = "select 80%",
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