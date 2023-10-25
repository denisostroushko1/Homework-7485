
rm(list = ls())

source(paste0(getwd(),"/Master Packages.R"))

data <- read_csv('./In Class Sim Study/dataHW3.csv')[,-1]

data$person_id = seq(from = 1, to = nrow(data), by = 1)

ggplot(data = data, 
       aes(x = (CD40))) + geom_histogram()

ggplot(data = data, 
       aes(x = log(CD40))) + geom_histogram()

# Time period 1
m1 <- lm(log(CD41) ~ A0*log(CD40), data = data)
summary(m1)

glm_r1 <- glm(R1 ~ A0*log(CD40), data = data, family = "binomial")
summary(glm_r1)

# Time period 2
m2 <- lm(log(CD42) ~ A1*log(CD41)*R1, data = data)
summary(m2)

glm_r2 <- glm(R2 ~ A1*log(CD41), data = data, family = "binomial")
summary(glm_r2)

# Time period 3 
m3 <- lm(log(CD43) ~ A2*log(CD42)*R2, data = data)
summary(m3)

glm_r3 <- glm(R3 ~ A2*log(CD42), data = data, family = "binomial")
summary(glm_r3)

# Terminal CD4 value 
m4 <- lm(log(CD44) ~ A3*log(CD43)*R3, data = data)
summary(m4)

glm_r4 <- glm(R4 ~ A3*log(CD43), data = data, family = "binomial")
summary(glm_r4)

simulate_trial <- 
  function(data_frame){
    
    # predict CD41
    data_frame$CD41 = exp(predict(m1, data_frame, type = "response"))
    # predict R1 
    pred_r1 = predict(glm_r1, data_frame, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    data_frame$R1 = rbinom(n = nrow(data_frame), size = 1, prob = pred_r1)
    
    data_frame <- 
      data_frame %>% 
      mutate(R1 = case_when(A0 == 0 ~ 0 ,
                            R0 == 1 ~ 1, 
                            T ~ R1)
                           # A0 == 1 & R0 == 0 ~ R1)
    )
    # predict CD42
    data_frame$CD42 = exp(predict(m2, data_frame, type = "response"))
    # predict R2
    pred_r2 = predict(glm_r2, data_frame, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    data_frame$R2 = rbinom(n = nrow(data_frame), size = 1, prob = pred_r2)
    
    data_frame <- 
      data_frame %>% 
      mutate(R2 = case_when(A1 == 0 ~ 0 ,
                            R1 == 1 ~ 1, 
                            T ~ R2)
                          #  A1 == 1 & R1 == 0 ~ R2)
    )
    
    # predict CD43
    data_frame$CD43 = exp(predict(m3, data_frame, type = "response"))
    # predict R3
    pred_r3 = predict(glm_r3, data_frame, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    data_frame$R3 = rbinom(n = nrow(data_frame), size = 1, prob = pred_r3)
    
    data_frame <- 
      data_frame %>% 
      mutate(R3 = case_when(A2 == 0 ~ 0 ,
                            R2 == 1 ~ 1, 
                            T ~ R3)
                            #A2 == 1 & R2 == 0 ~ R3)
    )
    
    # predict CD43
    data_frame$CD44 = exp(predict(m4, data_frame, type = "response"))
    # predict R3
    pred_r4 = predict(glm_r4, data_frame, type = "response")
    # draw R1 realized values from predicted probabilities of R1 
    data_frame$R4 = rbinom(n = nrow(data_frame), size = 1, prob = pred_r4)
    
    data_frame <- 
      data_frame %>% 
      mutate(R4 = case_when(A3 == 0 ~ 0 ,
                            R3 == 1 ~ 1, 
                            T ~ R4)
                            # A3 == 1 & R3 == 0 ~ R4)
    )
    
    return(data_frame)
    
  }


# data frame will come with R0 = 0, some simulated value of CD40, and Pre-Determined Treatment Sequence A0, A1, A2, A3 
# CD41, R1
# CD42, R2, 
# CD43, R3, 
# CD44 will be estimated inside the loop 

sim_1 = 
  data.frame(
    A0 = 1,
    A1 = 1, 
    A2 = 1, 
    A3 = 1, 
    R0 = 0, 
    CD40 = rnorm(n = 100000, 
                 mean = mean(data$CD40), 
                 sd = sd(data$CD40)) 
  )

sim_1$CD40 <- ifelse(sim_1$CD40 <= 50, 50, sim_1$CD40)
sim_1_complete = simulate_trial(sim_1)

summary(sim_1_complete$CD40)
summary(sim_1_complete$CD41)
summary(sim_1_complete$CD42)
summary(sim_1_complete$CD43)
summary(sim_1_complete$CD44)

summary(sim_1_complete$R0)
summary(sim_1_complete$R1)
summary(sim_1_complete$R2)
summary(sim_1_complete$R3)
summary(sim_1_complete$R4)
