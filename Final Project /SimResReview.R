
###################
# Interested in two outcomes: 
##    how do permutation factors in out experiment affect SE 
##    how do same factors impact MSE 
source("Master Packages.R")

res.1 <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result1.RDS")
)

res.2a <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result2a.RDS")
)

res.2b <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result2b.RDS")
)

res.2c <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result2c.RDS")
)

res.3a <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result3a.RDS")
)

res.3b <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result3b.RDS")
)

res.3c <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result3c.RDS")
)

res.4a <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result4a.RDS")
)

res.4b <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result4b.RDS")
)

res.5 <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result5.RDS")
)

res.6 <- readRDS(
  paste0(
    getwd(), "/Final Project /simulation results/result6.RDS")
)

all_results <- 
  rbind(res.6, 
        res.5, 
        res.4b, 
        res.4a, 
        res.3c, 
        res.3b,
        res.3a, 
        res.2c,
        res.2b,
        res.2a,
        res.1) %>% unique() %>% 
  filter(!(method %in% c("experience based", "pca - 80% of variance", "ttest"))) %>% 
  mutate(method = factor(method, 
                         levels = c("best case", setdiff(unique(method), "best case")))) %>% 
  arrange(n,m,s,method)%>% 
  mutate(ipw2_mse = (ATE.ipw2 - ATE.true)^2, 
         p_important_var = s/m)

###################################################

summary(all_results$SE.ipw2)
hist(all_results$SE.ipw2)

summary(log(all_results$SE.ipw2))
hist(log(all_results$SE.ipw2))

###################################################

summary(all_results$ipw2_mse)
hist(all_results$ipw2_mse)

summary(log(all_results$ipw2_mse))
length(which(is.na(log(all_results$ipw2_mse))))
hist(log(all_results$ipw2_mse))

###################################################
# marginal plots of main effects 

              ############################
              # BOOTSTRAP STANDARD ERROR # 
              ############################

    ggplot(data = all_results, 
           aes(x = n, y = SE.ipw2)) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    ## LOG SE vs SAMPLE SIZE COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = n, y = log(SE.ipw2))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, n
      ) %>% 
      summarize(mean_log_se = mean(log(SE.ipw2))) %>% 
      ggplot(aes(x = n, y = mean_log_se, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)
    
    ## LOG SE vs N IMPORATNT PREDICTORS COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = s, y = log(SE.ipw2))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, s
      ) %>% 
      summarize(mean_log_se = mean(log(SE.ipw2))) %>% 
      ggplot(aes(x = s, y = mean_log_se, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)
    
    ## LOG SE vs % IMPORATNT PREDICTORS COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = p_important_var, y = log(SE.ipw2))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, p_important_var
      ) %>% 
      summarize(mean_log_se = mean(log(SE.ipw2))) %>% 
      ggplot(aes(x = p_important_var, y = mean_log_se, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)


              ############################
              #             MSE          # 
              ############################

    ggplot(data = all_results, 
           aes(x = n, y = ipw2_mse)) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    ## LOG MSE vs SAMPLE SIZE COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = n, y = log(ipw2_mse))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, n
      ) %>% 
      summarize(mean_log_mse = mean(log(ipw2_mse))) %>% 
      ggplot(aes(x = n, y = mean_log_mse, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)
    
    ## LOG MSE vs N IMPORATNT PREDICTORS COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = s, y = log(ipw2_mse))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, s
      ) %>% 
      summarize(mean_log_mse = mean(log(ipw2_mse))) %>% 
      ggplot(aes(x = s, y = mean_log_mse, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)
    
    ## LOG MSE vs % IMPORATNT PREDICTORS COLORED BY METHOD
    ggplot(data = all_results, 
           aes(x = p_important_var, y = log(ipw2_mse))) + 
      theme_minimal() + 
      geom_point(aes(color = method)) + 
      geom_smooth(aes(group = method, color  = method), method = "lm")
    
    all_results %>% 
      group_by(
        method, p_important_var
      ) %>% 
      summarize(mean_log_mse = mean(log(ipw2_mse))) %>% 
      ggplot(aes(x = p_important_var, y = mean_log_mse, color = method, group = method)) + 
        theme_minimal() + 
        geom_point() + 
        geom_line() + 
        geom_smooth(aes(group = 1), se = F)


############################################
# regression models to explain the outcome 

## effect on bootsrapped standard error
se_model <- lm(log(SE.ipw2) ~ 
                 method + 
                 poly(n, 2) +
                 m + poly(s,2) + 
                 m:s + 
                 method:s, data = all_results)

options(scipen = 999)

se_model %>% summary() %>% tidy() %>% 
  mutate(ci_low = estimate - 1.96 * std.error, 
         ci_high = estimate + 1.96 * std.error) %>% 
  select(term, estimate, ci_low, ci_high) %>% 
  mutate_at(c("estimate", "ci_low", "ci_high"), ~ round(exp(.), 4)) 

(se_model %>% summary())$r.squared

ggplot(data = data.frame(pred = se_model$fitted.values, 
                         resid = rstudent(se_model)), 
       aes(x = pred, y = resid)) + 
  theme_minimal() + 
  geom_point(alpha = 0.5) + 
  geom_smooth(se = T) + 
  geom_hline(yintercept = -2, colour = "blue", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept =  2, colour = "blue", linewidth = 1, linetype = "dashed")+
  geom_hline(yintercept =  0, colour = "red", linewidth = 1, linetype = "dashed")

## effect on squared error 



lm(log(ipw2_mse) ~ method + poly(n, 2) + m + s + m:s + method:s, data = all_results) %>% summary() %>% tidy() %>% 
  mutate(ci_low = estimate - 1.96 * std.error, 
         ci_high = estimate + 1.96 * std.error) %>% 
  select(term, estimate, ci_low, ci_high) %>% 
  mutate_at(c("estimate", "ci_low", "ci_high"), ~ exp(.) %>% round(., 4)) 

(lm(log(ipw2_mse) ~ method + poly(n, 2) + m + s + m:s + method:s, data = all_results) %>% summary())$r.squared


