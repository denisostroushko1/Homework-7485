
rm(list = ls())

source(paste0(getwd(),"/Master Packages.R"))

data <- read_csv('./In Class Sim Study/dataHW3.csv')[,-1]

data$person_id = seq(from = 1, to = nrow(data), by = 1)
## right away work with LOG values of CD40, and other consequent counts of CD4
data <- 
  data %>% 
  mutate_at(vars(c("CD40", "CD41", "CD42", "CD43", "CD44")), ~ log(.))

ggplot(data = data, 
       aes(x = (CD40))) + geom_histogram()

# Time period 1
m1 <- lm(CD41 ~ A0 * CD40, data = data)
summary(m1)

glm_r1 <- glm(R1 ~ A0 * CD40, data = data, family = "binomial")
summary(glm_r1)

# Time period 2
m2 <- lm(CD42 ~ A1 * CD41 * R1, data = data)
summary(m2)

glm_r2 <- glm(R2 ~ A1 * CD41, data = data, family = "binomial")
summary(glm_r2)

# Time period 3 
m3 <- lm(CD43 ~ A2 * CD42 * R2, data = data)
summary(m3)

glm_r3 <- glm(R3 ~ A2 * CD42, data = data, family = "binomial")
summary(glm_r3)

# Terminal CD4 value 
m4 <- lm(CD44 ~ A3 * CD43 * R3, data = data)
summary(m4)

glm_r4 <- glm(R4 ~ A3 * CD43, data = data, family = "binomial")
summary(glm_r4)

simulate_trial <- 
  function(DATA_FRAME, 
           TRT_VECTOR){
    
    # assign treatment based on desired Treatment 
    DATA_FRAME <- 
      DATA_FRAME %>% 
      mutate(
        A0 = TRT_VECTOR[1],
        A1 = TRT_VECTOR[2],
        A2 = TRT_VECTOR[3],
        A3 = TRT_VECTOR[4]
      )
    # predict CD41
    DATA_FRAME$CD41 = predict(m1, DATA_FRAME, type = "response")
    # predict R1 
    pred_r1 = predict(glm_r1, DATA_FRAME, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    DATA_FRAME$pred_r1 = pred_r1
    
    DATA_FRAME$R1 = sapply(pred_r1,function(z){rbinom(n = 1,size = 1,prob = z)})
    
    DATA_FRAME <- 
      DATA_FRAME %>% 
      mutate(R1 = case_when(A0 == 0 ~ 0 ,
                            R0 == 1 ~ 1, 
                            T ~ R1)
                           # A0 == 1 & R0 == 0 ~ R1)
    )
    # predict CD42
    DATA_FRAME$CD42 = predict(m2, DATA_FRAME, type = "response")
    # predict R2
    pred_r2 = predict(glm_r2, DATA_FRAME, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    DATA_FRAME$pred_r2 = pred_r2
    DATA_FRAME$R2 = sapply(pred_r2,function(z){rbinom(n = 1,size = 1,prob = z)})
    
    DATA_FRAME <- 
      DATA_FRAME %>% 
      mutate(R2 = case_when(A1 == 0 ~ 0 ,
                            R1 == 1 ~ 1, 
                            T ~ R2)
                          #  A1 == 1 & R1 == 0 ~ R2)
    )
    
    # predict CD43
    DATA_FRAME$CD43 = predict(m3, DATA_FRAME, type = "response")
    # predict R3
    pred_r3 = predict(glm_r3, DATA_FRAME, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    DATA_FRAME$pred_r3 = pred_r3
    DATA_FRAME$R3 = sapply(pred_r3,function(z){rbinom(n = 1,size = 1,prob = z)})
    
    DATA_FRAME <- 
      DATA_FRAME %>% 
      mutate(R3 = case_when(A2 == 0 ~ 0 ,
                            R2 == 1 ~ 1, 
                            T ~ R3)
                            #A2 == 1 & R2 == 0 ~ R3)
    )
    
    # predict CD43
    DATA_FRAME$CD44 = predict(m4, DATA_FRAME, type = "response")
    # predict R3
    pred_r4 = predict(glm_r4, DATA_FRAME, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    DATA_FRAME$pred_r4 = pred_r4
    DATA_FRAME$R4 = sapply(pred_r4,function(z){rbinom(n = 1,size = 1,prob = z)})
    
    DATA_FRAME <- 
      DATA_FRAME %>% 
      mutate(R4 = case_when(A3 == 0 ~ 0 ,
                            R3 == 1 ~ 1, 
                            T ~ R4)
                            # A3 == 1 & R3 == 0 ~ R4)
    )
    
    return(DATA_FRAME)
    
  }

summarise_clinical_trial = function(COMPLETE_DATA){
  COMPLETE_DATA %>% 
    summarise(
      mean_R1 = mean(R1),
      mean_R2 = mean(R2),
      mean_R3 = mean(R3),
      mean_R4 = mean(R4),
      
      mean_CD41 = mean(CD41),
      mean_CD42 = mean(CD42),
      mean_CD43 = mean(CD43),
      mean_CD44 = mean(CD44),
      
      sd_CD41 = sd(CD41),
      sd_CD42 = sd(CD42),
      sd_CD43 = sd(CD43),
      sd_CD44 = sd(CD44)
    )
}




results = data.frame(
  iter = integer(), 
  mean_R1 = numeric(), 
  mean_R2 = numeric(), 
  mean_R3 = numeric(), 
  mean_R4 = numeric(), 
  mean_CD41 = numeric(), 
  mean_CD42 = numeric(), 
  mean_CD43 = numeric(), 
  mean_CD44 = numeric(), 
  sd_CD41 = numeric(), 
  sd_CD42 = numeric(), 
  sd_CD43 = numeric(), 
  sd_CD44 = numeric(), 
  treatment = character()
)

treatments = list(
  c(0,0,0,0), 
  c(0,0,0,1), 
  c(0,0,1,1), 
  c(0,1,1,1), 
  c(1,1,1,1)
)

set.seed(17927)
REPS = 100

for(j in 1:length(treatments)){
  
  for( i in 1:REPS){
    print(paste0('j: ',j, '. i: ', i))
    
    # data frame will come with R0 = 0, some simulated value of CD40, and Pre-Determined Treatment Sequence A0, A1, A2, A3 
    # CD41, R1
    # CD42, R2, 
    # CD43, R3, 
    # CD44 will be estimated inside the loop 
    
    sim_1 = 
      data.frame(
        CD40 = rnorm(n = 100000, 
                     mean = mean(data$CD40) + rnorm(n = 1, mean = 0, sd = 1), 
                     sd = sd(data$CD40) + abs(rnorm(n = 1, mean = 0, sd = 1))
                     ), 
        R0 = 0
      ) 
    
    summary(sim_1$CD40)
    
    sim_1_complete = simulate_trial(DATA_FRAME = sim_1, 
                                    TRT_VECTOR = treatments[[j]]
                                    )
    
    results <- 
      rbind(results, 
            
            summarise_clinical_trial(sim_1_complete) %>% 
              mutate(iter = i, 
                     treatment = paste0(treatments[[j]], collapse = ", ")
                       )
      )
  }
  
}

summary(results)

write.csv(results, './In Class Sim Study/100_iter_for_each_5_treatments.csv')

beepr::beep(4)

## Create data for summary 
results %>% 
  select(treatment, iter,mean_R1, mean_R2, mean_R3, mean_R4) %>% 
  pivot_longer(cols = c( 'mean_R1', 'mean_R2', 'mean_R3', 'mean_R4'), 
               names_to = 'char_time', 
               values_to = 'mean_R'
               ) %>% 
  arrange(treatment, iter, char_time) %>% 
  group_by(treatment, iter) %>% 
  mutate(time_num = seq(from = 1, to = 4, by = 1)) %>% 
  ungroup() -> long_R

results %>% 
  select(treatment, iter,mean_CD41, mean_CD42, mean_CD43, mean_CD44) %>% 
  pivot_longer(cols = c( 'mean_CD41', 'mean_CD42', 'mean_CD43', 'mean_CD44'), 
               names_to = 'char_time', 
               values_to = 'mean_CD'
               ) %>% 
  arrange(treatment, iter, char_time) %>% 
  group_by(treatment, iter) %>% 
  mutate(time_num = seq(from = 1, to =4, by = 1)) %>% 
  ungroup() -> long_mean_CD

results %>% 
  select(treatment, iter,sd_CD41,sd_CD42, sd_CD43, sd_CD44) %>% 
  pivot_longer(cols = c( 'sd_CD41', 'sd_CD42', 'sd_CD43', 'sd_CD44'), 
               names_to = 'char_time', 
               values_to = 'sd_CD'
               ) %>% 
  arrange(treatment, iter, char_time) %>% 
  group_by(treatment, iter) %>% 
  mutate(time_num = seq(from = 1, to = 4, by = 1)) %>% 
  ungroup() -> long_sd_CD


all_long_res = 
  long_R %>% select(-char_time) %>% 
  
  left_join(long_mean_CD %>% select(-char_time), 
            by = c("treatment", "iter","time_num")) %>% 
  
  left_join(long_sd_CD %>% select(-char_time), 
            by = c("treatment", "iter","time_num")) %>% 
  
  select(treatment, iter, time_num, everything())

### Average LOG CD4 Counts
mean_plots_data_CD = 
  all_long_res %>% 
  group_by(treatment, time_num) %>% 
  summarise(mean_plot = mean(mean_CD), 
            sd_plot = sd(mean_CD)/sqrt(100)) %>% 
  
  mutate(ci_min = mean_plot - 1.96 * sd_plot, 
         ci_max = mean_plot + 1.96 * sd_plot)

ggplot(data = mean_plots_data_CD , 
       mapping = aes(x = time_num, y = mean_plot, group = treatment, color = treatment)) + 
  theme_classic() + 
  geom_point(size = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = ci_min, ymax = ci_max), width = 0.2, alpha = 0.5)


### Average Resistance 
mean_plots_data_R = 
  all_long_res %>% 
  group_by(treatment, time_num) %>% 
  summarise(mean_plot = mean(mean_R), 
            sd_plot = sd(mean_R)/sqrt(100)) %>% 
  
  mutate(ci_min = mean_plot - 1.96 * sd_plot, 
         ci_max = mean_plot + 1.96 * sd_plot)

ggplot(data = mean_plots_data_R , 
       mapping = aes(x = time_num, y = mean_plot, group = treatment, color = treatment)) + 
  theme_classic() + 
  geom_point(size = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = ci_min, ymax = ci_max), width = 0.2, alpha = 0.5)


### Average Standard deviation of CD counts from sim iterations  
mean_plots_data_CD_SD = 
  all_long_res %>% 
  group_by(treatment, time_num) %>% 
  summarise(mean_plot = mean(sd_CD), 
            sd_plot = sd(sd_CD)/sqrt(100)) %>% 
  
  mutate(ci_min = mean_plot - 1.96 * sd_plot, 
         ci_max = mean_plot + 1.96 * sd_plot)

ggplot(data = mean_plots_data_CD_SD , 
       mapping = aes(x = time_num, y = mean_plot, group = treatment, color = treatment)) + 
  theme_classic() + 
  geom_point(size = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = ci_min, ymax = ci_max), width = 0.2, alpha = 0.5)
