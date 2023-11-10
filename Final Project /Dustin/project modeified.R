
# clean environment
rm(list=ls())

library(SimDesign) # rmvnorm
library(locfit) # expit function
library(purrr) # rberoulii
library(glmnet)
library(tidyverse)

# samples, i.e. sample size 
n <- 500
# potential covariates, m is the total number 
m <- 100

# Covariates matrix
# Setting up Correlation matrix Sigma
A <- matrix(runif(m^2)*2-1, ncol=m) 
Sigma <- t(A) %*% A
Scaled.Sigma <- cov2cor(Sigma)

# # Examining if the correlation matrix make sense
# upper.triangle <- Scaled.Sigma[upper.tri(Scaled.Sigma)]
# min(upper.triangle)
# max(upper.triangle)

# Simulating covariates X based on multi-variate normal distribution
X_i <- rmvnorm(n, mean = rep(0, m),
               sigma = Scaled.Sigma) 
# For simplicity, we could also use normal distribution for each coavariate instead.

# Proportion of TRUE covariates that will be used to build the model
covar.prop <- 0.2
# generate a vector indicator, determine whether the covariate will be used or not.
true.cov <- purrr::rbernoulli(m, p=covar.prop)

# Treatment Assignment based on Covariates
temp.coef <- runif(m, min = -0.75, max = 0.75) # generate coefficients
trt.coef <- true.cov * temp.coef # indicator * coefficients = B
trt.prob <- expit(X_i %*% trt.coef) # Pr(A_i) is expit[logit(A_i)] where logit(A_i) = xB
temp.trt <- data.frame(true.cov, temp.coef, trt.coef) # column 3 will be the final model

# distribution of the probability of getting treatment
# hist(trt.prob) 

# Getting treatment or not based on bernoulli given probability trt.prob
A <- purrr::rbernoulli(m, trt.prob)

# Outcome model based on Covariates
# range of the coefficients
temp.coef <- runif(m, min = -1, max = 1) # random draw for coefficients
out.coef <- true.cov * temp.coef # indicator * coefficients = B
temp.out <- data.frame(true.cov, temp.coef, out.coef) # column 3 will be the final model

# Potential outcome Y0
Y0.underlying <-  X_i %*% out.coef # underlying Y0 without adding error term
# error term: normal with mean 0 and 10% of the full range of Y0 as variance
e0 = rnorm(m, 0, abs(range(Y0.underlying)[1]-range(Y0.underlying)[2])*0.1) 
Y0 <- Y0.underlying + e0

# Distribution of underlying Y0 and Y0 with error term
# hist(Y0.underlying)
# hist(Y0)

# Treatment Effect (decided based on full range of Y0) , 50% of the full range of Y0
## SET TREATMENT EFFECT TO 2 

      # trt.effect <- abs(range(Y0.underlying)[1]-range(Y0.underlying)[2])*0.5
      trt.effect <- 1
# Potential outcome Y1
Y1.underlying <-  X_i %*% out.coef + trt.effect
# error term: normal with mean 0 and 10% of the full range of Y1 as variance

    ### I think that we do not need to get e1 added, we already incorporated error -- e0 -- 
    ###   to the actual Y0, so we just add treatment effect to Y0 to get Y1
    # e1 = rnorm(m, 0, abs(range(Y1.underlying)[1]-range(Y1.underlying)[2])*0.1)
    # Y1 <- Y1.underlying + e1
  Y1 <- Y0 + trt.effect

# Distribution of underlying Y1 and Y1 with error term
# hist(Y1.underlying)
# hist(Y1)

# Observed Y
Y <- A*Y1 + (1-A)*Y0

# Final simulated df
# sim.df <- data.frame(Y, A, Y0, Y1, X_i)
sim.df <- data.frame(Y, A, X_i) %>% mutate(A = as.numeric(A))

# mean(sim.df$Y1) - mean(sim.df$Y0)
# t.test(sim.df$Y~sim.df$A)

# Things need to determine:
# How to simulate correlation matrix Sigma ?
# Min and Max of coefficients?
# Scale of noise for Y0 and Y1? Scale of Treatment effect (T)?

        ####### backward variable selection 

backward_model <- MASS::stepAIC(object = lm(Y ~ ., data = sim.df), 
                                scope = list(lower = lm(Y ~ 1, data = sim.df), 
                                             upper = lm(Y ~ ., data = sim.df)),
                                direction = "backward", 
                                trace = 0)

summary(backward_model)
length(coef(backward_model))

zeros = rep(0, length(colnames(sim.df)[2:length(sim.df)]))

### store results 
  backw_selected_coef <- zeros
  backw_selected_coef[which(
        colnames(sim.df)[2:length(sim.df)] %in% 
        names(coef(backward_model))[names(coef(backward_model)) != "(Intercept)"] 
      )] <- 1
  
  backw_signif_coef <- zeros
  backw_signif_coef[which(
        colnames(sim.df)[2:length(sim.df)] %in% 
        (summary(backward_model)$coefficients[,1][summary(backward_model)$coefficients[,4] <= 0.05] %>% names())
      )] <- 1

forward_model <- MASS::stepAIC(object = lm(Y ~ 1, data = sim.df), 
                               scope = list(lower=lm(Y ~ 1, data = sim.df), 
                                            upper = lm(Y ~ ., data = sim.df)),
                                direction = "forward", 
                                trace = 0)

summary(forward_model)
length(coef(forward_model))

### store results 
  fwd_selected_coef <- zeros
  fwd_selected_coef[which(
        colnames(sim.df)[2:length(sim.df)] %in% 
        names(coef(forward_model))[names(coef(forward_model)) != "(Intercept)"] 
      )] <- 1
  
  fwd_signif_coef <- zeros
  fwd_signif_coef[which(
        colnames(sim.df)[2:length(sim.df)] %in% 
        (summary(forward_model)$coefficients[,1][summary(forward_model)$coefficients[,4] <= 0.05] %>% names())
      )] <- 1

data.frame(
  iter = 1, 
  var_names = colnames(sim.df)[2:length(sim.df)], 
  true = c(1, ifelse(sign(out.coef) != 0, 1, 0)), 
  
  backw_selected_coef = backw_selected_coef, 
  backw_signif_coef = backw_signif_coef, 
  
  fwd_selected_coef = fwd_selected_coef, 
  fwd_signif_coef = fwd_signif_coef
) -> example_of_res
            



