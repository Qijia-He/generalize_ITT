library(mvtnorm)
library(DescTools)
library(MASS)
library(glmnet)
library(bridgedist)
library(scales)
library(dplyr)
library(randomForest)
library(mgcv)
################################################################################
################################################################################

# Helper functions
expit <- function(x) exp(x)/(1+exp(x))

# Define D and Y models
D_md2 <- function(Z, X){
  p <- expit(Z*(4.5 - 0.3 * X[,1] - 0.4 * X[,5] + 0.2 * X[,1]+ 0.3 * X[,2] - 2) +
               (1-Z)*(- 0.3 * X[,1] - 0.4 * X[,5]-1.5))
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

# Generate a dataset {X, Z, D, Y} given:
# (1) D(Z) model: D | Z, X
# (3) Y(D) model: Y(D) | X, U
# sample size n, and parameter c
generate_data_binary <- function(n, c, D_model, Y_md){
  
  # Generate X ~ MVT(mu, Id)
  # mu = (c, c, c, 0, ..., 0)
  mean_vec = c(c, c, c, rep(0, 7))
  X = rmvnorm(n = n, mean = mean_vec, diag(0.5, nrow = 10))
  X[X > 4] <- 4
  X[X < 0] <- 0
  
  # Generate Z ~ Bern(0.5)
  Z = rbinom(n, 1, 0.5)
  
  # Generate D | Z, X
  D = D_model(Z, X)$model
  
  # Generate U
  #U = rbinom(n, 1, 0.5)
  U = bridgedist::rbridge(n)
  
  # Generate Y(0) and Y(1)
  Y = Y_md(X, U, Z)$model
  
  # Output X, Z, D, Y
  return(data.frame(X, Z, D, Y))
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
run_simu_once <- function(n_1, n_2, n, c, D_md1, D_md2, Y_md1, Y_md2, Y_md3, 
                          n_boot = 500,
                          variance = TRUE){
  
  # data generation
  dt_h1 = generate_data_binary(n = n_1, c = 1.2*c, D_md1, Y_md21)
  dt_h2 = generate_data_binary(n = n_2, c = c, D_md2, Y_md22)
  dt = generate_data_binary(n = n, c = 0.8*c, D_md1, Y_md23)
  
  # the ground truth estimator
  est_ground_truth = mean(dt$Y[dt$Z==1]) - mean(dt$Y[dt$Z==0])
  
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
  
    
    return(c( 
      est_ground_truth, est_ground_truth_CI_lower, est_ground_truth_CI_upper,
      naive_h1_est, naive_h1_CI_lower, naive_h1_CI_upper,
      naive_h2_est, naive_h2_CI_lower, naive_h2_CI_upper,
      est_wald,  wald_CI_lower, wald_CI_upper
    ))
  } else {
    return(c(est_ground_truth, 
             est_naive1, est_naive2,
             est_wald
    ))
  }
  
}


# Look at the sampling distribution
# For scenario I, we plug in Y_md11, Y_md12, Y_md13
# For scenario II, we plug in Y_md21, Y_md22, Y_md23
# We consider parameter c = 0, 0.25, 0.5, sample sizes n_1 = n_2 = n = 1000, 2000, 5000
n_1 = 1000
n_2 = 1000
n = 1000
D_md1 = D_md1
D_md2 = D_md2

n_boot = 1000

simu_res1000 = data.frame(est = c(run_simu_once(n_1, n_2, n, 0, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T),
                                  run_simu_once(n_1, n_2, n, 0.25, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T),
                                  run_simu_once(n_1, n_2, n, 0.50, D_md1, D_md2, Y_md11, Y_md12, Y_md13, n_boot = n_boot, variance = T)) )


simu_res1000


