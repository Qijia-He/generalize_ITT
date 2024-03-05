library(mvtnorm)
library(DescTools)
library(MASS)
library(glmnet)
library(bridgedist)
library(scales)
library(dplyr)
library(randomForest)
library(mgcv)
library(dplyr)


#---------------------------
#       read data
#---------------------------
#data.hptn0 <- read.csv("NI_trial_HPTN.csv")
#data.prep <- read.csv("historical_trial_Prep.csv")

#outcome.hptn <- read.csv("HPTN084_outcome.csv")[, -1]
names(outcome.hptn) <- c("HIV", "Followup.time", "id")
outcome.hptn$Followup.time <- outcome.hptn$Followup.time/365.25
data.hptn <- merge(data.hptn0, outcome.hptn, by = "id")

data.hptn$Education[data.hptn$Education != "not complete primary school"] <- "complete primary school"
data.prep$Education[data.prep$Education != "not complete primary school"] <- "complete primary school"

data.hptn[sapply(data.hptn, is.character)] <- lapply(data.hptn[sapply(data.hptn, is.character)], 
                                                     as.factor)
data.prep[sapply(data.prep, is.character)] <- lapply(data.prep[sapply(data.prep, is.character)], 
                                                     as.factor)

data.prep$HIV <- ifelse(data.prep$HIV == "No", 0, 1)
data.hptn <- data.hptn[-which(data.hptn$Followup.time == 0), ]

data.hptn$illness3 <- ifelse(data.hptn$Chlamydia == "Pos" | data.hptn$Gonorrhea == "Pos" |
                               data.hptn$Trichomonas == "Pos", 1, 0)

data.prep$illness3 <- ifelse(data.prep$Chlamydia == "Pos" | data.prep$Gonorrhea == "Pos" |
                               data.prep$Trichomonas == "Pos", 1, 0)

data.prep$Age <- scale(data.prep$Age)
data.hptn$Age <- scale(data.hptn$Age)

posterior.S <- function(data, s){
  dt <- data
  dt$S[dt$S0 != s] <- 0
  dt$S[dt$S0 == s] <- 1
  model.f = glm(S0 ~ Age + Employment + Education + illness3 + Syphilis, data = dt, family = 'binomial')
  posterior = predict(model.f, newdata = data, type = 'response')
  return(posterior)
}

data.hptn$S0 <- 1
data.prep$S0 <- 0
dt.combine <- rbind(data.hptn, data.prep)
dt.combine$S0 <- as.factor(dt.combine$S0)
dt.combine$ps <- posterior.S(dt.combine, 1)

data.hptn$ps <- dt.combine$ps[1:nrow(data.hptn)]
data.prep$ps <- dt.combine$ps[(nrow(data.hptn) + 1):nrow(dt.combine)]

# hist(data.prep$ps)
# sum(data.prep$ps >= 0.5) # 206
index.prep.ps.high <- which(data.prep$ps >= 0.5) # 206  206/(206 + 861) = 0.193
index.prep.ps.low <- which(data.prep$ps < 0.5) # 861  861/(206 + 861) = 0.807

################################################################################
################################################################################

# Helper functions
expit <- function(x) exp(x)/(1+exp(x))

# Define D and Y models
D_md2 <- function(Z, X){
  p <- expit(Z*(2.5 - 0.1 * X[,1] - 0.4 * X[,5] + 0.3 * X[,2]) +
               (1 - Z)*(- 0.3 * X[,1] - 0.4 * X[,5] - 1.5))
  model <- rbinom(length(Z), 1, prob = p)
  return(list(model = model, p = p))
}

# D model for D_h1
D_md1 <- function(Z, X){
  p <- expit(Z*(5 - 0.2*X[,5] + 0.2*X[,1]-3) + (1-Z)*(- 0.2*X[,5] -1))
  model <- rbinom(length(Z), 1, prob = p)
  return(list(model = model, p = p) )
}

################################################################################
# Scenario I
Y_md12<- function(X, U, Z){
  p <- expit(1.4 - 1*X[,1] - 0.6*X[,2] + 0.4*X[,3] -0.6*X[,5] + 3.5*Z)
  model <- rbinom(length(U), 1, prob = p)
  return(list(model = model, p = p) )
}

# P(Y = 1 | Z = 1, X) = expit(2.6 - 0.6*X1 - 0.8*X2 + 0.4*X3)
# P(Y = 1 | Z = 0, X) = expit(1.6 - 0.7*X1 - 0.7*X2 + 0.4*X3 - 0.2*X5)
Y_md11 <- function(X, U, Z){
  p1 <- expit(2.6 - 0.6*X[,1] - 0.8*X[,2] + 0.4*X[,3])
  p0 <- expit(1.6 - 0.7*X[,1] - 0.7*X[,2] + 0.4*X[,3] - 0.2*X[,5])
  model <- rbinom(length(U), 1, prob = ifelse(Z == 1, p1, p0 ) )
  return(list(model = model, p1 = p1, p0 = p0)  )
}

# P(Y | Z, X) in the target trial
# P(Y = 1 | Z = 0, X) = 0
# P(Y = 1 | Z = 1, X) satisfies ITT = CATE(X, h_1) * CC(X, h_2)
Y_md13 <- function(X, U, Z){
  p = (expit(2.6-0.6*X[,1]-0.8*X[,2]+0.4*X[,3])-
         expit(1.6-0.7*X[,1]-0.7*X[,2]+0.4*X[,3]-0.2*X[,5]))/(expit(2-0.2*X[,5] + 0.2*X[,1]) -
                                                                expit(-0.2*X[,4]-1))*(expit(4.5 - 0.3 * X[,1] - 0.4 * X[,5] + 0.2 * X[,1]+ 0.3 * X[,2] - 2)-
                                                                                        expit(- 0.3 * X[,1] - 0.4 * X[,5] - 1.5))
  model = rbinom(length(U), 1, prob = ifelse(Z == 1, p, 0))
  return(list(model = model, p = p) )
}

Y_md22<- function(X, U, Z){
  p <- expit(1.4 - 1*X[,1] - 0.6*X[,2]^0.5 + 0.4*X[,3] -0.6*X[,5] + 3.5*Z)
  model <- rbinom(length(U), 1, prob = p)
  return(list(model = model, p = p) )
}

Y_md21 <- function(X, U, Z){
  p1 <- expit(2.6 - 0.6*X[,1] - 0.8*X[,2] + 0.4*X[,3]^0.5)
  p0 <- expit(1.6 - 0.7*X[,1] - 0.7*X[,2]^3 + 0.4*X[,3] - 0.2*X[,5])
  model <- rbinom(length(U), 1, prob = ifelse(Z == 1, p1, p0 ) )
  return(list(model = model, p1 = p1, p0 = p0)  )
}

# P(Y | Z, X) in the target trial
# P(Y = 1 | Z = 0, X) = 0
# P(Y = 1 | Z = 1, X) satisfies ITT = CATE(X, h_1) * CC(X, h_2)

Y_md23 <- function(X, U, Z){
  #p = (expit(2.6 - 0.6*X[,1] - 0.8*X[,2] + 0.4*X[,3]^0.5)-
  #      expit(1.6 - 0.7*X[,1] - 0.7*X[,2]^3 + 0.4*X[,3] - 0.2*X[,5]))/
  # (expit(2-0.2*X[,5] + 0.2*X[,1]) - expit(-0.2*X[,5]-1))*
  # (expit(4.5 - 0.3 * X[,1] - 0.4 * X[,5] + 0.2 * X[,1]+ 0.3 * X[,2] - 2)- expit(- 0.3 * X[,1] - 0.4 * X[,5] - 1.5))
  p = (D_md2(1, X)$p - D_md2(0, X)$p )*
    ( Y_md21(X, 1, 1)$p1 - Y_md21(X, 1, 0)$p0 )/
    (D_md1(1, X)$p - D_md1(0, X)$p)
  p[p > 1] <- 1
  model = rbinom(length(U), 1, prob = ifelse(Z == 1, p, 0))
  return(list(model = model, p = p) )
}
# Generate a dataset {X, Z, D, Y} given:
# (1) D(Z) model: D | Z, X
# (3) Y(D) model: Y(D) | X, U
# sample size n, and parameter c

# overlap = 1: bad overlap
# overlap = 2: original data
# overlap = 3: good overlap
generate_data_binary <- function(overlap, n, D_md1, D_md2, Y_md1, Y_md2, Y_md3){
  
  if(overlap == 2){
    dt = data.hptn[sample(dim(data.hptn)[1], n, replace = TRUE), ]
    dt_h1 = data.prep[sample(dim(data.prep)[1], n, replace = TRUE), ]
    dt_h2 = data.prep[sample(dim(data.prep)[1], n, replace = TRUE), ]
  }
  # bad overlap
  if(overlap == 1){
    n_11 = n * 0.1
    n_12 = n * 0.9
    n_21 = n * 0.15
    n_22 = n * 0.85
    dt = data.hptn[sample(dim(data.hptn)[1], n, replace = TRUE), ]
    dt_h1 = data.prep[c(sample(index.prep.ps.high, n_11, replace = TRUE), sample(index.prep.ps.low, n_12, replace = TRUE)), ]
    dt_h2 = data.prep[c(sample(index.prep.ps.high, n_21, replace = TRUE), sample(index.prep.ps.low, n_22, replace = TRUE)), ]
  }
  
  # good overlap
  if(overlap == 3){
    n_11 = n * 0.4
    n_12 = n * 0.6
    n_21 = n * 0.5
    n_22 = n * 0.5
    dt = data.hptn[sample(dim(data.hptn)[1], n, replace = TRUE), ]
    dt_h1 = data.prep[c(sample(index.prep.ps.high, n_11, replace = TRUE), sample(index.prep.ps.low, n_12, replace = TRUE)), ]
    dt_h2 = data.prep[c(sample(index.prep.ps.high, n_21, replace = TRUE), sample(index.prep.ps.low, n_22, replace = TRUE)), ]
  }
  
  X = dt %>% select("Age", "Employment", "Education", "illness3", "Syphilis")
  X1 = dt_h1 %>% select("Age", "Employment", "Education", "illness3", "Syphilis")
  X2 = dt_h2 %>% select("Age", "Employment", "Education", "illness3", "Syphilis")
  
  X$Employment <- ifelse(X$Employment == "employed", 1, 0)
  X$Education <- ifelse(X$Education == "complete primary school", 1, 0)
  X$Syphilis <- ifelse(X$Syphilis == "Neg", 0, 1)
  
  X1$Employment <- ifelse(X1$Employment == "employed", 1, 0)
  X1$Education <- ifelse(X1$Education == "complete primary school", 1, 0)
  X1$Syphilis <- ifelse(X1$Syphilis == "Neg", 0, 1)
  
  X2$Employment <- ifelse(X2$Employment == "employed", 1, 0)
  X2$Education <- ifelse(X2$Education == "complete primary school", 1, 0)
  X2$Syphilis <- ifelse(X2$Syphilis == "Neg", 0, 1)
  
  # Generate Z ~ Bern(0.5)
  Z = rbinom(n, 1, 0.5)
  Z1 = rbinom(n, 1, 0.5)
  Z2 = rbinom(n, 1, 0.5)
  
  # Generate D | Z, X
  D = D_md1(Z, X)$model
  D1 = D_md1(Z1, X1)$model
  D2 = D_md2(Z2, X2)$model
  
  # Generate U
  #U = rbinom(n, 1, 0.5)
  U = bridgedist::rbridge(n)
  
  # Generate Y(0) and Y(1)
  Y = Y_md3(X, U, Z)$model
  Y1 = Y_md1(X1, U, Z1)$model
  Y2 = Y_md2(X2, U, Z2)$model
  
  # Output X, Z, D, Y
  return(list(dt = data.frame(X, Z, D, Y),
              dt_h1 = data.frame(X1, Z1, D1, Y1),
              dt_h2 = data.frame(X2, Z2, D2, Y2)) )
}

# Estimate the conditional compliance CC(X) given dt
# CC(X) = P(D = 1 | Z = 1, X) - P(D= 1 | Z = 0, X)
# Output an estimate CC(X) stored as a function

CC <- function(dt, outcome_model = 'glm'){
  if(outcome_model == 'glm'){
    # Fit a D | Z = 1, X model in the Z == 1 stratum
    md_Z_1 = glm(D ~ X1 + X2 + X5, data = dt, family = 'binomial', subset = (Z == 1))
    
    # Fit a D | Z = 0, X model in the Z == 0 stratum
    md_Z_0 = glm(D ~ X1 + X5, data = dt, family = 'binomial', subset = (Z == 0))
    
    # Return a function estimate of CC(X)
    CC_func <- function(dt){
      return(predict(md_Z_1, newdata = dt, type = 'response') -
               predict(md_Z_0, newdata = dt, type = 'response'))
    }
  }else if (outcome_model == 'rf'){
    # Fit a D | Z = 1, X model in the Z == 1 stratum
    md_Z_1 = randomForest(x = dt[dt$Z == 1, c(1, 2, 5)], y = as.factor(dt$D[dt$Z == 1]), ntree = 500)
    
    # Fit a D | Z = 0, X model in the Z == 0 stratum
    md_Z_0 = randomForest(x = dt[dt$Z == 0, c(1, 5)], y = as.factor(dt$D[dt$Z == 0]), ntree = 500)
    
    # Return a function estimate of CC(X)
    CC_func <- function(dt){
      return(predict(md_Z_1, newdata = dt, type = 'prob')[,2] -
               predict(md_Z_0, newdata = dt, type = 'prob')[,2])
    }
  }
}

# Wald estimate the conditional average treatment effect CATE(X)
# Return a function estimate
# Method = 'gaussian' or 'binomial'
# outcome_model = 'glm', 'rf', 'gam'
Wald <- function(dt, method = 'binomial', outcome_model = 'glm'){
  # Regression-based Wald estimate
  if(outcome_model == 'glm'){
    md_D_Z_1 = glm(D ~ X1 + X5, data = dt, family = 'binomial', subset = (Z == 1))
    md_D_Z_0 = glm(D ~ X5, data = dt, family = 'binomial', subset = (Z == 0))
    md_Y_Z_1 = glm(Y ~ 1, data = dt, family = method, subset = (Z == 1))
    md_Y_Z_0 = glm(Y ~  1, data = dt, family = method, subset = (Z == 0))
    
    # Return a function estimate of CATE(X)
    CATE_func <- function(dt){
      numerator = predict(md_Y_Z_1, newdata = dt, type = 'response') -
        predict(md_Y_Z_0, newdata = dt, type = 'response')
      denominator = predict(md_D_Z_1, newdata = dt, type = 'response') -
        predict(md_D_Z_0, newdata = dt, type = 'response')
      return(numerator/denominator)
    }
  }
  else if(outcome_model == 'rf'){
    md_D_Z_1 = randomForest(x = dt[dt$Z == 1, c(1, 5)], y = as.factor(dt$D[dt$Z == 1]), ntree = 500)
    md_D_Z_0  = randomForest(as.factor(D) ~ X5, 
                             data = dt, subset = (Z == 0), ntree = 500)
    md_Y_Z_1 = randomForest(x = dt[dt$Z == 1, 1:3], y = as.factor(dt$Y[dt$Z == 1]), ntree = 500)
    md_Y_Z_0 = randomForest(x = dt[dt$Z == 0, c(1, 2, 3, 5)], y = as.factor(dt$Y[dt$Z == 0]), ntree = 500)
    
    # Return a function estimate of CATE(X)
    CATE_func <- function(dt){
      numerator = predict(md_Y_Z_1, dt, type = 'prob')[,2] -
        predict(md_Y_Z_0, dt, type = 'prob')[,2]
      denominator = predict(md_D_Z_1, newdata = dt, type = 'prob')[,2] -
        predict(md_D_Z_0, newdata = dt, type = 'prob')[,2]
      return(numerator/denominator)
    }
  }else if(outcome_model == 'gam'){
    md_D_Z_1 = gam(D ~ s(X1) + s(X5), data = dt, family = 'binomial', subset = (Z == 1))
    md_D_Z_0 = gam(D ~ s(X5), data = dt, family = 'binomial', subset = (Z == 0))
    md_Y_Z_1 = gam(Y ~ s(X1) + s(X2) + s(X3), data = dt, family = method, subset = (Z == 1))
    md_Y_Z_0 = gam(Y ~  s(X1) + s(X2) + s(X3) + s(X5), data = dt, family = method, subset = (Z == 0))
    
    # Return a function estimate of CATE(X)
    CATE_func <- function(dt){
      numerator = predict(md_Y_Z_1, newdata = dt, type = 'response') -
        predict(md_Y_Z_0, newdata = dt, type = 'response')
      denominator = predict(md_D_Z_1, newdata = dt, type = 'response') -
        predict(md_D_Z_0, newdata = dt, type = 'response')
      return(numerator/denominator)
    }
  }
}

# Construct an efficient influence function to get EIF based estimator
# Return the EIF based estimator phi.est and confidence interval
# dt1 and dt2 are from two historical trials, and dt is from the NI trial
# oucome_model = 'glm', 'rf', 'gam'

EIF <- function(dt1, dt2, dt, outcome_model = 'glm'){
  
  data <- rbind(dt, dt1, dt2)
  n_1 = n_2 = n
  data$S <- c(rep(3, n), rep(1, n_1), rep(2, n_2) ) # 3 = NI; 1 = h1; 2 = h2
  
  # estimate posteriors f(S = s|X)
  # posterior.S <- function(data, s){
  #   model.f <- lda(S ~ .- Z - D - Y, data = data)
  #   prediction <- model.f %>% predict(data)
  #   posterior1 <- prediction$posterior
  #   return(posterior[, which(colnames(posterior) == s)])
  # }
  
  posterior.S <- function(data, s){
    dt <- data
    dt$S[dt$S != s] <- 0
    dt$S[dt$S == s] <- 1
    model.f = gam(S ~ s(X1) + X2 + X3 + X5, data = dt, family = 'binomial')
    posterior = predict(model.f, newdata = data, type = 'response')
    return(posterior)
  }
  
  # f(Z | S = s, X) = 0.5
  f.Z <- function(data, s){
    return(0.5)
  }
  
  # estimate mu0_Y(X; S = s) = E[Y |S, Z=0, X]
  # estimate mu0_D(X; S = s) = E[D |S, Z=0, X]
  if(outcome_model == 'glm'){
    mu0.D <- function(data, s){
      if(s == 2) {md_D_S = glm(D ~ X1, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) }
      if(s == 1) {md_D_S = glm(D ~ X1, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) }
      value = predict(md_D_S, newdata = data, type = 'response')
      return(value)
    }
    mu1.D <- function(data, s){
      if(s == 2) {md_D_S = glm(D ~ X1, data = data, family = 'binomial', subset = (Z == 1) & (S == s)) }
      if(s == 1) {md_D_S = glm(D ~ X1, data = data, family = 'binomial', subset = (Z == 1) & (S == s)) }
      value = predict(md_D_S, newdata = data, type = 'response')
      return(value)
    }
    mu0.Y <- function(data, s){
      md_Y_S = glm(Y ~ 1, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) 
      value = predict(md_Y_S, newdata = data, type = 'response')
      return(value)
    }
    mu1.Y <- function(data, s){
      if(s == 2) {md_Y_S = glm(Y ~ 1, data = data, family = 'binomial', subset = (Z == 1) & (S == s))}
      if(s == 1) {md_Y_S = glm(Y ~ 1, data = data, family = 'binomial', subset = (Z == 1) & (S == s))}
      value = predict(md_Y_S, newdata = data, type = 'response')
      return(value)
    }
  }
  
  if(outcome_model == 'rf'){
    mu0.D <- function(data, s){
      if(s == 2) {md_D_S = randomForest(x = data[data$Z == 0 & data$S == s, c(1, 5)], 
                                        y = as.factor(data$D[data$Z == 0 & data$S == s]), 
                                        ntree = 500, nodesize = length((data$D[data$Z == 0 & data$S == s])/10)) }
      if(s == 1) { md_D_S = randomForest(as.factor(D) ~ X5, 
                                         data = data, subset = (Z == 0) & (S == s), 
                                         ntree = 500, nodesize = length((data$D[data$Z == 0 & data$S == s])/10))}
      value = predict(md_D_S, newdata = data, type = 'prob')[,2]
      return(value)
    }
    
    mu1.D <- function(data, s){
      if(s == 2) {md_D_S = randomForest(x = data[data$Z == 1 & data$S == s, c(1, 2, 5)], 
                                        y = as.factor(data$D[data$Z == 1 & data$S == s]), 
                                        ntree = 500, nodesize = length((data$D[data$Z == 1 & data$S == s])/10))}
      
      if(s == 1) {md_D_S = randomForest(x = data[data$Z == 1 & data$S == s, c(1, 5)], 
                                        y = as.factor(data$D[data$Z == 1 & data$S == s]), 
                                        ntree = 500, nodesize = length((data$D[data$Z == 0 & data$S == s])/10))}
      value = predict(md_D_S, newdata = data, type = 'prob')[,2]
      return(value)
    }
    mu0.Y <- function(data, s){
      md_Y_S = randomForest(x = data[data$Z == 0 & data$S == s, c(1, 2, 3, 5)], 
                            y = as.factor(data$Y[data$Z == 0 & data$S == s]), ntree = 500)
      value = predict(md_Y_S, newdata = data, type = 'prob')[,2]
      return(value)
    }
    
    mu1.Y <- function(data, s){
      if(s == 2) {md_Y_S = randomForest(x = data[data$Z == 1 & data$S == s, c(1, 2, 3, 5)], 
                                        y = as.factor(data$Y[data$Z == 1 & data$S == s]), ntree = 500)}
      if(s == 1) {md_Y_S = randomForest(x = data[data$Z == 1 & data$S == s, 1:3], 
                                        y = as.factor(data$Y[data$Z == 1 & data$S == s]), ntree = 500)}
      value = predict(md_Y_S, newdata = data, type = 'prob')[,2]
      return(value)
    }
  }
  
  if(outcome_model == 'gam'){
    mu0.D <- function(data, s){
      if(s == 2) {md_D_S = gam(D ~ s(X1) + X5, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) }
      if(s == 1) {md_D_S = gam(D ~ X5, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) }
      value = predict(md_D_S, newdata = data, type = 'response')
      return(value)
    }
    mu1.D <- function(data, s){
      if(s == 2) {md_D_S = gam(D ~ s(X1) + X2 + X5, data = data, family = 'binomial', subset = (Z == 1) & (S == s)) }
      if(s == 1) {md_D_S = gam(D ~ s(X1) + X5, data = data, family = 'binomial', subset = (Z == 1) & (S == s)) }
      value = predict(md_D_S, newdata = data, type = 'response')
      return(value)
    }
    mu0.Y <- function(data, s){
      md_Y_S = gam(Y ~ s(X1) + X2 + X3 + X5, data = data, family = 'binomial', subset = (Z == 0) & (S == s)) 
      value = predict(md_Y_S, newdata = data, type = 'response')
      return(value)
    }
    mu1.Y <- function(data, s){
      if(s == 2) {md_Y_S = gam(Y ~ s(X1) + X2 + X3 + X5, data = data, family = 'binomial', subset = (Z == 1) & (S == s))}
      if(s == 1) {md_Y_S = gam(Y ~ s(X1) + X2 + X3, data = data, family = 'binomial', subset = (Z == 1) & (S == s))}
      value = predict(md_Y_S, newdata = data, type = 'response')
      return(value)
    }
  }
  
  # estimate conditional compliance for S=s: delta_D(X; s)
  # estimate ITT for S=s: delta_Y(X; s)
  # estimate wald estimator for S=s: delta(X; s) = delta_D(X; s)/delta_Y(X; s)
  
  delta_D <- function(data, s){
    value <- mu1.D(data, s) - mu0.D(data, s) 
    return(value)
  }
  
  delta_Y <- function(data, s){
    value <- mu1.Y(data, s) - mu0.Y(data, s)
    return(value)
  }
  
  delta <- function(data, s){
    delta_D <- delta_D(data, s)
    delta_Y <- delta_Y(data, s)
    return(delta_Y/delta_D)
  }
  
  delta_D_1 <- delta_D(data, 1)
  delta_D_2 <- delta_D(data, 2)
  delta_Y_1 <- delta_Y(data, 1)
  delta_1 <- delta_Y_1/delta_D_1
  posterior.S3 <-  posterior.S(data, 3)
  term1 <- (2*data$Z - 1)*(data$S == 1)/f.Z(data, 1) * 
    posterior.S3/posterior.S(data, 1) *
    delta_D_2/delta_D_1 *
    (data$Y - mu0.Y(data, 1) - (data$D - mu0.D(data, 1))*delta_1)
  
  term2 <- (2*data$Z - 1)*(data$S == 2)/f.Z(data, 2) * 
    posterior.S3/posterior.S(data, 2) *
    delta_1 *
    (data$D - mu0.D(data, 2) - delta_D_2*data$Z)
  
  term3 <- (data$S == 3)*delta_1*delta_D_2
  
  k <- term1 + term2 + term3
  phi.est <- sum(k)/length(data[which(data$S == 3), 1])
  EIF <- 3*term1 + 3*term2 + 3*(term3 - (data$S == 3)*phi.est)
  var <- mean(EIF^2)/(3*n)
  lower <- phi.est - qnorm(0.975)*sqrt(var)
  upper <- phi.est + qnorm(0.975)*sqrt(var)
  #phi.win <-  sum( Winsorize(k, probs = c(0.01, 0.99)))/length(data[which(data$S == 3), 1])
  
  return(list = c(phi.est = phi.est, lower = lower, upper = upper,
                  var = var
                  #, phi.win = phi.win
  ) )
}


# Bootstrap a confidence interval for the EIF based estimator
var_boot_EIF <- function(dt1, dt2, dt, outcome_model = 'glm', winsorize = F, n_boot = 1000){
  
  est_boot_EIF = numeric(n_boot)
  for (i in 1:n_boot){
    # Resample dtt, dt2, and dt
    n_1 = dim(dt1)[1]
    dt1_resample = dt1[sample(n_1, n_1, replace = TRUE), ]
    
    n_2 = dim(dt2)[1]
    dt2_resample = dt2[sample(n_2, n_2, replace = TRUE), ]
    
    n = dim(dt)[1]
    dt_resample = dt[sample(n, n, replace = TRUE), ]
    
    est_boot_EIF[i] <- EIF(dt1_resample, dt2_resample, dt_resample, outcome_model)[1]
    #if(is.na(est_boot_EIF[i])){ dd <- data.frame(rbind(dt1_resample, dt2_resample, dt_resample))}
  }
  k.na <- sum(is.na(est_boot_EIF))
  return(c(quantile(est_boot_EIF, 0.025), quantile(est_boot_EIF, 0.975), k.na) )
}


# Bootstrap a confidence interval for the wald estimator
var_boot <- function(dt1, dt2, dt, outcome_model = 'glm', n_boot = 1000){
  
  est_boot = numeric(n_boot)
  for (i in 1:n_boot){
    # Resample dtt, dt2, and dt
    n_1 = dim(dt1)[1]
    dt1_resample = dt1[sample(n_1, n_1, replace = TRUE), ]
    
    n_2 = dim(dt2)[1]
    dt2_resample = dt2[sample(n_2, n_2, replace = TRUE), ]
    
    n = dim(dt)[1]
    dt_resample = dt[sample(n, n, replace = TRUE), ]
    
    cc_x = CC(dt2_resample)
    cate_x = Wald(dt1_resample, method = 'binomial', outcome_model)
    est_boot[i] = mean(cc_x(dt_resample) * cate_x(dt_resample))
  }
  return(est_boot)
}

# A naive estimator that estimates E[Y | Z = 1, X] and E[Y | Z = 0, X] using data from historical trial h2
# using one historical dataset, i.e., assuming conditional constancy
cond_con_est2 <- function(dt, dt2, method = 'binomial', n_boot = 1000){
  
  # Estimate E[Y | Z = 0, X] and E[Y | Z = 1, X]
  md_Y_Z_1 = glm(Y ~ X1 + X2 + X3 + X5, data = dt2, family = method, subset = (Z == 1))
  md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt2, family = method, subset = (Z == 0))
  
  # Apply estimated functions on the target trial data
  est = mean(predict(md_Y_Z_1, newdata = dt, type = 'response') -
               predict(md_Y_Z_0, newdata = dt, type = 'response'))
  
  # bootstrap a CI
  est_boot = numeric(n_boot)
  for (i in 1:n_boot){
    
    # Resample dt2 and dt
    n_1 = dim(dt2)[1]
    dt2_resample = dt2[sample(n_1, n_1, replace = TRUE), ]
    
    n = dim(dt)[1]
    dt_resample = dt[sample(n, n, replace = TRUE), ]
    
    md_Y_Z_1 = glm(Y ~ X1 + X2 + X3 + X5, data = dt2_resample, family = method, subset = (Z == 1))
    md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt2_resample, family = method, subset = (Z == 0))
    
    # Put together the estimate
    est_boot[i] = mean(predict(md_Y_Z_1, newdata = dt_resample, type = 'response') -
                         predict(md_Y_Z_0, newdata = dt_resample, type = 'response'))
  }
  
  return(c(est
           ,quantile(est_boot, 0.025), quantile(est_boot, 0.975)
  ))
}

# A naive estimator that estimates E[Y | Z = 1, X] and E[Y | Z = 0, X] using data from historical trial h1
# using one historical dataset, i.e., assuming conditional constancy
cond_con_est1 <- function(dt, dt1, method = 'binomial', n_boot = 1000){
  
  # Estimate E[Y | Z = 0, X] and E[Y | Z = 1, X]
  md_Y_Z_1 = glm(Y ~ X1 + X2 + X3, data = dt1, family = method, subset = (Z == 1))
  md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt1, family = method, subset = (Z == 0))
  
  # Apply estimated functions on the target trial data
  est = mean(predict(md_Y_Z_1, newdata = dt, type = 'response') -
               predict(md_Y_Z_0, newdata = dt, type = 'response'))
  
  # bootstrap a CI
  est_boot = numeric(n_boot)
  for (i in 1:n_boot){
    
    # Resample dt2 and dt
    n_1 = dim(dt1)[1]
    dt1_resample = dt1[sample(n_1, n_1, replace = TRUE), ]
    
    n = dim(dt)[1]
    dt_resample = dt[sample(n, n, replace = TRUE), ]
    
    md_Y_Z_1 = glm(Y ~ X1 + X2 + X3, data = dt1_resample, family = method, subset = (Z == 1))
    md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt1_resample, family = method, subset = (Z == 0))
    
    # Put together the estimate
    est_boot[i] = mean(predict(md_Y_Z_1, newdata = dt_resample, type = 'response') -
                         predict(md_Y_Z_0, newdata = dt_resample, type = 'response'))
  }
  
  return(c(est
           ,quantile(est_boot, 0.025), quantile(est_boot, 0.975)
  ))
}


################################################################################
################################################################################
# Simulation for a binary outcome
# n_1, n_2, n are sample size of historical trial h1, historical trial h2, and the NI trial 
# n_boot is the bootstrap repetitions times
# variance = TRUE means the confidence interval would be reported
run_simu_once <- function(overlap, n, D_md1, D_md2, Y_md1, Y_md2, Y_md3, 
                          n_boot = 500,
                          variance = TRUE){
  
  # data generation
  d = generate_data_binary(overlap = overlap, n = n, D_md1, D_md2, Y_md1, Y_md2, Y_md3)
  dt_h1 = d$dt_h1
  dt_h2 = d$dt_h2
  dt = d$dt
  names(dt) = names(dt_h1) = names(dt_h2) = c("X1", "X2", "X3", "X4", "X5", "Z", "D", "Y")
  
  # the ground truth estimator
  est_ground_truth = mean(dt$Y[dt$Z==1]) - mean(dt$Y[dt$Z==0])
  
  # EIF based estimator when outcome_model = 'glm'
  est_EIF <- EIF(dt1 = dt_h1, dt2 = dt_h2, dt = dt)[1]
  # EIF based estimator when outcome_model = 'gam'
  est_EIF_gam <- EIF(dt1 = dt_h1, dt2 = dt_h2, dt = dt, outcome_model = 'gam')[1]
  
  # Calculate the historical-data-driven estimator: Wald-based when outcome_model = 'glm'
  cc_x = CC(dt_h2)
  cate_x = Wald(dt_h1, method = 'binomial')
  est_wald = mean(cc_x(dt) * cate_x(dt))
  
  # Calculate the historical-data-driven estimator: Wald-based when outcome_model = 'gam'
  # cate_x_gam = Wald(dt_h2, method = 'binomial', outcome_model = 'gam')
  # est_wald_gam = mean(cc_x(dt) * cate_x_gam(dt))
  
  md_Y_Z_1 = glm(Y ~ X1 + X2 + X3 + X5, data = dt_h2, family = 'binomial', subset = (Z == 1))
  md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt_h2, family = 'binomial', subset = (Z == 0))
  
  # Apply estimated functions on the target trial data
  est_naive2 = mean(predict(md_Y_Z_1, newdata = dt, type = 'response') -
                      predict(md_Y_Z_0, newdata = dt, type = 'response'))
  
  md_Y_Z_1 = glm(Y ~ X1 + X2 + X3, data = dt_h1, family = 'binomial', subset = (Z == 1))
  md_Y_Z_0 = glm(Y ~ X1 + X2 + X3 + X5, data = dt_h1, family = 'binomial', subset = (Z == 0))
  
  # Apply estimated functions on the target trial data
  est_naive1 = mean(predict(md_Y_Z_1, newdata = dt, type = 'response') -
                      predict(md_Y_Z_0, newdata = dt, type = 'response'))
  
  if (variance){
    # t test results
    prop_test_res = prop.test(x = c(sum(dt$Y[dt$Z==1]), sum(dt$Y[dt$Z == 0])), 
                              n = c(length(dt$Y[dt$Z==1]), length(dt$Y[dt$Z==0])))
    naive_ITT = c(mean(dt$Y[dt$Z == 1]) - mean(dt$Y[dt$Z == 0]),
                  prop_test_res$conf.int)
    est_ground_truth <- naive_ITT[1]
    est_ground_truth_CI_lower <- naive_ITT[2]
    est_ground_truth_CI_upper <- naive_ITT[3]
    
    # Estimate ITT assuming conditional constancy using dt_h2
    naive_h2 = cond_con_est2(dt, dt_h2, method = 'binomial', n_boot = n_boot)
    naive_h2_est = naive_h2[1]
    naive_h2_CI_lower = naive_h2[2]
    naive_h2_CI_upper = naive_h2[3]
    
    naive_h1 = cond_con_est1(dt, dt_h1, method = 'binomial', n_boot = n_boot)
    naive_h1_est = naive_h1[1]
    naive_h1_CI_lower = naive_h1[2]
    naive_h1_CI_upper = naive_h1[3]
    
    # Calculate the historical-data-driven estimator: Wald-based
    cc_x = CC(dt_h2)
    cate_x = Wald(dt_h1, method = 'binomial')
    est_wald = mean(cc_x(dt) * cate_x(dt))
    boot_est = var_boot(dt1 = dt_h1, dt2 = dt_h2, dt = dt, n_boot = n_boot)
    wald_CI_lower = quantile(boot_est, 0.025)
    wald_CI_upper = quantile(boot_est, 0.975)
    
    # Bootstrap the variance associated with the estimate
    est_EIF <- EIF(dt1 = dt_h1, dt2 = dt_h2, dt = dt)[1]
    boot_est_EIF = var_boot_EIF(dt1 = dt_h1, dt2 = dt_h2, dt = dt, winsorize = F, n_boot = n_boot)
    EIF_CI_lower = boot_est_EIF[1]
    EIF_CI_upper = boot_est_EIF[2]
    
    EIF.md <- EIF(dt1 = dt_h1, dt2 = dt_h2, dt = dt, outcome_model = 'gam')
    est_EIF_gam <- EIF.md[1]
    lower.b <- EIF.md[2]
    upper.b <- EIF.md[3]
    
    return(c( 
      est_ground_truth, est_ground_truth_CI_lower, est_ground_truth_CI_upper,
      naive_h1_est, naive_h1_CI_lower, naive_h1_CI_upper,
      naive_h2_est, naive_h2_CI_lower, naive_h2_CI_upper,
      est_wald,  wald_CI_lower, wald_CI_upper,
      est_EIF,  EIF_CI_lower, EIF_CI_upper,
      est_EIF_gam, lower.b, upper.b
    ))
  } else {
    return(c(est_ground_truth, 
             est_naive1, est_naive2,
             est_wald,
             est_EIF, 
             est_EIF_gam
    ))
  }
  
}


n_1 = 1000
n_2 = 1000
n = 1000

n_boot = 1000

simu_res1000 = data.frame(est = c(run_simu_once(n_1, n_2, n, 0, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T),
                                  run_simu_once(n_1, n_2, n, 0.25, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T),
                                  run_simu_once(n_1, n_2, n, 0.50, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T)) )


