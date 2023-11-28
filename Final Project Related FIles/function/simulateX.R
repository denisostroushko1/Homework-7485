# library(dplyr)
# library(tidyverse)
# library(clusterGeneration) # generate random correlation matrix
# # library(simstudy)
# library(SimDesign) # rmvnorm
# library(locfit) # expit function
# library(purrr) # rberoulii

# samples
# n <- 500
# # potential covariates
# m <- 100
# # true covariates
# s <- 20
# # replications
# n.rep <- 2

sim.x <- function(n, m, s, n.rep) {
  
  sim.x_step1 <- function(n, m, s) {
    # Covariates matrix
    # Setting up Correlation matrix Sigma
    Sigma <- rcorrmatrix(m)
    # Simulating covariates X based on multi-variate normal distribution
    X_i <- rmvnorm(n, mean = rep(0, m), sigma = Sigma) 
    
    # Proportion of TRUE covariates that will be used to build the model
    s <- s
    cov.idx <- sample(1:m, s)
    ord.idx <- order(cov.idx)
    # Use the indices to obtain the ordered vector
    cov.idx <- cov.idx[ord.idx]
    true.cov <- 1:m %in% cov.idx
    
    # Treatment Assignment based on Covariates
    temp.coef <- runif(m, min = -0.5, max = 0.5) # generate coefficients
    trt.coef <- true.cov * temp.coef # indicator * coefficients = B
    trt.prob <- expit(X_i %*% trt.coef) # Pr(A_i) is expit[logit(A_i)] where logit(A_i) = xB
    temp.trt <- data.frame(true.cov, temp.coef, trt.coef) # column 3 will be the final model
    
    # distribution of the probability of getting treatment
    # hist(trt.prob)
    
    # Getting treatment or not based on bernoulli given probability trt.prob
    A <- purrr::rbernoulli(m, trt.prob) 
    A <- ifelse(A==TRUE, 1, 0)
    
    # Outcome model based on Covariates
    # range of the coefficients
    temp.coef <- runif(m, min = -1, max = 1) # random draw for coefficients
    out.coef <- true.cov * temp.coef # indicator * coefficients = B
    temp.out <- data.frame(true.cov, temp.coef, out.coef) # column 3 will be the final model
  
    # Potential outcome Y0
    Y0.underlying <-  X_i %*% out.coef # underlying Y0 without adding error term
    return(list(Y0 = Y0.underlying, 
                X = X_i,
                A = A,
                true.cov = true.cov,
                cov.idx = cov.idx))
  }
  # get Xi and Y0.underlying
  step1.list <- replicate(n.rep, sim.x_step1(n, m, s), simplify = FALSE)
  
  # Extract Y0.underlying from each replication to calculate ATE
  Y0_values <- sapply(step1.list, function(sim) sim$Y0)
  # Calculate the standard deviation across all replications
  sd_Y0 <- sd(unlist(Y0_values))
  # ATE will be ATE/sd(Y0_underlying) ~ 0.5
  true.ATE = 0.5*sd_Y0
  # error term: normal with mean 0 and variance 5% of sd(Y0_underlying)
  e = rnorm(m, 0, sd_Y0*0.05) 
  
  step1.list <- lapply(step1.list, function(sim) {
    sim$Y1 <- sim$Y0 + true.ATE
    sim$Y <- sim$A*sim$Y1 + (1-sim$A)*sim$Y0 + e
    sim$e <- e
    sim$true.ATE <- true.ATE
    return(sim)
  })
  order <- c("Y", "A", "X", "Y1", "Y0", "true.cov", "cov.idx", "e", "true.ATE")
  # Reorder the elements
  step1.list <- lapply(step1.list, function(sim) {
    sim <- sim[order]  
    return(sim)
  })
}


# Testing
# x.list <- sim.x(500, 100, 10, 2)
