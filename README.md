# Homework-7485

Homework for Fall 2023 PUBH 7485 at the University of Minnesota. I found a lot of cool ways to present tables with summary statistics, as well as test and models output. I also 
found a lot of useful packages that interact with `ggplot` mainly. 

Also, textbook and lecture materials cover way more intro-level regression methods and theory, and I wish to use this file as a table of  contents. 

By default, I use `kable`, `kableExtra`, `tidyverse`, and `ggplot2` a lot, so I will only be making notes of new packages 
I find and use for the analysis. 

# HW1

* Statistical Concepts: 
  + "Table 1" introduction. Formal introduction to data description for publications and the concept of Standardized Mean Differences (SMD)
  + Unadjsuted Treatment Effect: difference in mean of the outcome variable between two (or more) groups. 
  + Average Treatment Effect (ATE) estimation using regression adjustment approach to account for confounding variables. 
  + Obtaining variance of ATE estimator using bootstrap standard error [standard deviation of bootstrap sampling distribution]
  
* Implementation via R code: 
  + Saved 'printed' version of Table 1 can be passed on to kable for further formatting of table for publication 
  + Regression adjustment is done using base R `lm` and `glm` models. In practice, any regression approach can be used to 
    obtain expected response under treatment/no treatment for each observation 
  + Bootstrapping of ATE estimator variance is implement using base R `for` loops and sampling with replacement techniques. 
  
* Packages: 
  + `TableOne`; produce Table 1 easily. Includes options to show varying levels of factor predictors, perform t-test for means, 
    perform SMD for each variable. Other less common options of styling the table are available too. 
  
# HW2

* Statistical Concepts: 
  + Propensity Score model implementation with a logistic regression model. Any regression model that can return a probability of     being in the treatment (case) group can be used. 
  + ATE estimators using Inverse Probability Weighting (IPW version 1 and 2) and Propensity Score Stratification (PSS)
  + Bootstrap resampling to get standard errors of ATE estimators 
  + Visual examination of weights for IPW estimator as a way to check for unreasonably high or low weights. 
  
* Implementation via R code: 
  + Implementation of ATE estimators from initial propensity score models
  
# HW3

* Statistical Concepts: 
  + Doubly Robust ATE estimator example: Augmented IPW ATE estimator 
  + Using regression models for the outcome variable to obtain expected response for each observation under yes/no treatment 
  + Using regression models to obtain probability of being in the treatment group for each observation 
  + Combining results from the two regression models and using steps from IPW estimation procedure to get Augmented IPW ATE estimate
  + Doubly robust estimators will have reduced variance if both the outcome and the treatment regression models are
    correctly specified, i.e. all functional relationships are correct and all confounders are included
  
* Implementation via R code: 
  + all calculation and steps are performed via base R and `dplyr` functions. 
  + Bootstrap is done at the end of the file, seed was set, results are stored once to avoid lengthy computation every time. 
  
# HW4

* Statistical Concepts: 
  + Propensity Score regression model: propensity score and treatment indicator as predictors for the outcome variables as 
    another way of regression adjustment 
  + Average treatment effect among Treated (ATT): focusing analysis on those who received treatment, and what the effect of 
    treatment is in this population. 
  + Propensity score matching for the estimation of treatment effect: this method always!! estimates ATT since we find 
    people similar to those who were treated from the pool of non-treated. Hence, the focus is on treated and like-treated 
    observations only. 
  + Evaluation of matching quality by looking at the difference in: 
    - distance for matching among the matched pairs 
    - difference in propensity scores among the matched pairs 
  + Potential source of correlation that arises from matching is not addressed, need to review opinions and existing research to 
    determine if leveraging correlation to reduce standard errors is possible/statistically sound. 
  
* Implementation via R code: 
  + same ideas as previous sections 
  + propensity score matching using a `MatchIt` package
  
* Packages: 
  + `MatchIt` functions work with vectors of fitted probabilities. Output of said function is (are) data set(s) with 
    all matched observations. Additionally, a new variable is created to index pairs and identify members that are matched to each other. 
  
# HW5

* Statistical Concepts: 
  + First ideas in mediation analysis. 
  + Fitting mediator regression model as a function of a treatment indicator and other confusers. Using fitted value of the 
    mediator variable in further analysis. 
  + Fitting outcome regression models with observed values of the mediator variable to be used in further analysis 
  + Obtaining Controlled Direct Effect (CDE) from regression adjustment method. One addition is to use a fixed value of mediator 
    for each person. We can vary value of mediator and see how it impact the value of CDE 
  + Natural Indirect Effect estimation: 
    - using mediator regression model to get expected mediator value under treatment and no treatment 
    - using outcome regression model to get expected value of outcome under no treatment using mediator under no treatment. 
      Same for expected outcome under treatment for all patients using individual expected mediator values under treatment 

* R code implementation: all the same ideas and steps as previous sections 
  
# HW6

* Statistical Concepts: 
  + G formula computation for varying treatment regimes over time 
  + Defining outcomes at time $i$ as a function of outcomes at time $i-1$ and other covariates 
  + Performing simulations and bootstrapping to obtain outcomes at the end of the observation period to compare treatment regimes, maximize desired variable, and find the best treatment strategy 
  
* Implementation via R code: 
  + `gformula` implements G-computation using a single function. It requires a lot of inputs, examples are documented in the doc. 
  + `ipwtm` allows to implement a flexible inverse probability weighting. 
  + bootstrapping standard errors for estimates from the `ipwtm` function 
    - due to longitudinal nature of data we need to resample entire sequences observed for a given patient, and not individual 
      rows from the data. We need to preserve correlated nature of outcomes to make sure that bootstrap standard error is not 
      over-inflated 
      
* Packages: 
  + `gfoRmula`: this package is great, but unintuitive to use. Several guides are available: 
    - [original git repo](https://github.com/CausalInference/gfoRmula)
    - [user guide by authors](https://www.cell.com/patterns/fulltext/S2666-3899(20)30008-8)
    - [CRAN documentation](https://cran.r-project.org/web/packages/gfoRmula/gfoRmula.pdf)
    - [Mix of Theory and R code](https://christopherbboyer.com/files/teaching/lab3_parametric_g_formula.pdf)
    - [G formula with Multiple Imputation](https://thestatsgeek.com/2023/01/31/g-formula-for-causal-inference-via-multiple-imputation/)
  + `ipw`: package for inverse probability weighting. Allows to IP weight time varying treatment indicators, i.e. when we 
    observe of sequence of treatment indicators (decisions) for a single patient 
    - 
    
# Final Project

* Description: for our final project we decided to study how the process of variable selection in the context of causal inference 
  impacts Bias and Variance of ATE estimators. Regression model use is everywhere in causal inference, but it is not 
  clear how to pick correct confounders from a set of all possible predictors that are available in the data. 
  
  We found that in large sample sizes with large number true predictors that need to be captured in the model, forward selection, LASSO, and outcome adaptive LASSO have similar performance. However, when the number of true predictors is small, Outcome Adaptive LASSO is the best performing method in terms of MSE. It guarantees more accurate capture of true confounders and higher rate of filtering out noise. 

* Statistical Concepts: 
  + Variable selection methods such as forward selection, LASSO and outcome adaptive LASSO. 
  + Full factorial simulation study design. A total of 81 combinations of permuted variables are compared to find the best 
    performing way of constructing ATE estimator. 
  + Use of linear regression model to analyze simulation result and understand marginal impact of variable selection method for 
    model building, after adjusting for sample size, variable number, and true confounder amount. 
  
  
