# clean environment
rm(list=ls())
library(dplyr)
library(tidyverse)
library(clusterGeneration) # generate random correlation matrix
library(SimDesign) # rmvnorm
library(locfit) # expit function
library(purrr) # rberoulii
library(glmnet)

set.seed(100)
source("./function/simulateX.R")
source("./function/ATE_80pct.R")
source("./function/ATE_fwrd.R")
source("./function/ATE_bkwrd.R")
source("./function/ATE_lasso.R")
source("./function/ATE_pca.R")

# Sample Size (n)
# Number of Potential Covariates (m)
# True Covariates (s)
# Replications (n.rep)
n.rep <- 10
# Bootstrapping
B <- 100
# add small jitter to avoid unique error in the cut function 
small_jitter <- 1e-10  # Adjust the magnitude of jitter as needed

n.seq <- c(4500)
m.seq <- c(500)
s.seq <- c(10,20,30)

start_time <- Sys.time()
cov.info <- NA
result <- NA
for (fct.n in seq_along(n.seq)){
  for (fct.m in seq_along(m.seq)){
    for (fct.s in seq_along(s.seq)){
      x.list <- sim.x(n = n.seq[fct.n], m = m.seq[fct.m], s = s.seq[fct.s], n.rep = n.rep)
      # Selecting 80% Correct Covariates
      ATE_80pct_res <- lapply(x.list, ATE_80pct)
      ATE_fwrd_res <- lapply(x.list, ATE_fwrd)
      ATE_bkwrd_res <- lapply(x.list, ATE_bkwrd)
      ATE_lasso_res <- lapply(x.list, ATE_lasso)
      ATE_pca_res <- lapply(x.list, ATE_pca)
      temp.cov.info <- do.call(rbind, lapply(list(ATE_80pct_res, ATE_fwrd_res, ATE_bkwrd_res, ATE_lasso_res, ATE_pca_res), function(result_list) {
        do.call(rbind, lapply(result_list, function(result) result$cov.info))
      }))
      temp.result <- do.call(rbind, lapply(list(ATE_80pct_res, ATE_fwrd_res, ATE_bkwrd_res, ATE_lasso_res, ATE_pca_res), function(result_list) {
        do.call(rbind, lapply(result_list, function(result) result$result))
      }))
      cov.info <- rbind(cov.info, temp.cov.info) %>% filter(!is.na(method))
      result <- rbind(result, temp.result) %>% filter(!is.na(method))
      print(paste0("n: ", fct.n, "; m: ", fct.m, "; n: ", fct.s))
    }
  }  
}
end_time <- Sys.time()
execution_time <- end_time - start_time
cat("Total execution time:", as.numeric(execution_time, units = "mins"), "min(s)\n")

saveRDS(cov.info, file = "cov.info9.RDS")
saveRDS(result, file = "result9.RDS")
