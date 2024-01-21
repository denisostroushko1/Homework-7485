# Use this to test the function inside ATE_80pct
# x.list.temp <- x.list[[1]]

ATE_ttest <- function(x.list.temp) {
  X <- x.list.temp$X
  true.X <- X[,x.list.temp$true.cov]
  A <- x.list.temp$A
  Y <- x.list.temp$Y
  ps.df.org <- data.frame(Y, A)
  
  start_time <- Sys.time()
  res.ttest <- t.test(Y~A, data = ps.df.org)
  end_time <- Sys.time()
  execution_time <- end_time - start_time
  cat("T test execution time:", as.numeric(execution_time, units = "secs"), "secs(s)\n")
  
  ATE.pss <- res.ttest$estimate[2]-res.ttest$estimate[1]
  ATE.pss <- unname(ATE.pss)
  SE.ATE.pss <- res.ttest$stderr
  
  ATE.ipw2 <- NA
  SE.ATE.ipw2 <- NA
  # Setting up saving results
  overlap.cov <- NA
  miss.cov <- NA
  add.cov <- NA
  ml.cov <- NA
  # res.df <- data.frame(true.cov = x.list.temp$cov.idx, ml.cov = which(rnd.cov),)
  # Return the ATE
  return(list(cov.info = data.frame(method = "ttest",
                                    n = nrow(X),
                                    m = ncol(X),
                                    s = ncol(true.X),
                                    overlap.cov = overlap.cov,
                                    miss.cov = miss.cov,
                                    add.cov = add.cov,
                                    execution_time = execution_time),
              result = data.frame(method = "ttest",
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