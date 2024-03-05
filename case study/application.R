
# read data - patient-level data cannot be shared
# data.hptn0 <- read.csv("NI_trial_HPTN.csv")
#data.prep <- read.csv("historical_trial_Prep.csv")

#outcome.hptn <- read.csv("HPTN084_outcome.csv")[, -1]

#----------------------------------------
#             data processing
#----------------------------------------
names(outcome.hptn) <- c("HIV", "Followup.time", "id")
outcome.hptn$Followup.time <- outcome.hptn$Followup.time/365.25
data.hptn <- merge(data.hptn0, outcome.hptn, by = "id")

data.hptn["age_group"] = cut(data.hptn$Age, c(0, 25, 30, 35, Inf), c("0-25", "26-30", "31-35", ">35"), include.lowest=TRUE)

data.prep["age_group"] = cut(data.prep$Age, c(0, 25, 30, 35, Inf), c("0-25", "26-30", "31-35", ">35"), include.lowest=TRUE)

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

data.hptn <- data.hptn[, c(1:4, 12, 13, 5:11, 14:16)]

# Regression-based Wald estimate
md_D_Z_1 = glm(compliance ~  age_group + Employment + Education1 + illness3 + Syphilis, data = data.prep, family = 'binomial', subset = (Arm == "TDF/FTC"))
md_Y_Z_1 = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment + Education1 + illness3 + Syphilis  , data = data.prep, family = poisson(link=log), subset = (Arm == "TDF/FTC")) 
md_Y_Z_0 = glm(HIV ~ offset(log(Followup.time)) +  age_group + Employment + Education1 + illness3 + Syphilis  , data = data.prep, family = poisson(link=log), subset = (Arm == "Placebo")) 

#----------------------------------------
#             wald estimator
#----------------------------------------
# Return a function estimate of CATE(X), point identification
CATE_func <- function(dt){
  numerator = predict(md_Y_Z_1, newdata = dt, type = 'response')/dt$Followup.time -
    predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time
  denominator = predict(md_D_Z_1, newdata = dt, type = 'response') - 0
  return(numerator/denominator)
}

# partial identification of CATE using prep
# dt is NI trial data HPTN 084, p is the compliance rate in Prep
CATE_func_partial <- function(dt){
  numerator = predict(md_Y_Z_1, newdata = dt, type = 'response')/dt$Followup.time -
    predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time
  denominator = predict(md_D_Z_1, newdata = dt, type = 'response')
  lower.bound <- 0 - predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time
  upper.bound <- 0
  p <- denominator
  return(list(lower = p * numerator/denominator + (1 - p) * lower.bound, 
              upper = p * numerator/denominator + (1 - p) * upper.bound))
}


model.Dz.1 <- glm(compliance ~  age_group + Employment + Education1 + illness3 + Syphilis , data = data.hptn, family = 'binomial', subset = (Arm == "TDF/FTC"))

# p = E[D|Z=Placebo, X]
cc_x <- function(dt, p = 0.05){
  predict(model.Dz.1, newdata = dt, type = 'response') - p
}

#ITT, point identification
print("point identification")
mean(cc_x(data.hptn, 0) * CATE_func(data.hptn))
mean(cc_x(data.hptn, 0.05) * CATE_func(data.hptn))
mean(cc_x(data.hptn, 0.1) * CATE_func(data.hptn))

#ITT partial identification
print("partial identification")
mean(cc_x(data.hptn, 0) * CATE_func_partial(data.hptn)$lower)
mean(cc_x(data.hptn, 0) * CATE_func_partial(data.hptn)$upper)

mean(cc_x(data.hptn, 0.05) * CATE_func_partial(data.hptn)$lower)
mean(cc_x(data.hptn, 0.05) * CATE_func_partial(data.hptn)$upper)

mean(cc_x(data.hptn, 0.1) * CATE_func_partial(data.hptn)$lower)
mean(cc_x(data.hptn, 0.1) * CATE_func_partial(data.hptn)$upper)


# bootstrap CI


ITT0 <- NULL
ITT5 <- NULL
ITT10 <- NULL

ITT0.lower <- NULL
ITT5.lower <- NULL
ITT10.lower <- NULL

ITT0.upper <- NULL
ITT5.upper <- NULL
ITT10.upper <- NULL

ITT_cab_placebo <- NULL
in_placebo <- NULL

set.seed(23)

# illness3 (0/1)
# Arm (0/1)
# 2*2 = 4 strata
# sample from each strata and combine them together

data.hptn1 <- data.hptn[data.hptn$illness3 == 1 & data.hptn$Arm == "TDF/FTC", ]
data.hptn2 <- data.hptn[data.hptn$illness3 != 1 & data.hptn$Arm == "TDF/FTC", ]
data.hptn3 <- data.hptn[data.hptn$illness3 == 1 & data.hptn$Arm != "TDF/FTC", ]
data.hptn4 <- data.hptn[data.hptn$illness3 != 1 & data.hptn$Arm != "TDF/FTC", ]

data.prep1 <- data.prep[data.prep$illness3 == 1 & data.prep$Arm == "TDF/FTC", ]
data.prep2 <- data.prep[data.prep$illness3 != 1 & data.prep$Arm == "TDF/FTC", ]
data.prep3 <- data.prep[data.prep$illness3 == 1 & data.prep$Arm != "TDF/FTC", ]
data.prep4 <- data.prep[data.prep$illness3 != 1 & data.prep$Arm != "TDF/FTC", ]


data.hptn.list <- list(data.hptn1, data.hptn2, data.hptn3, data.hptn4)
data.prep.list <- list(data.prep1, data.prep2, data.prep3, data.prep4)

resample <- function(data.list){
  data_resample <- NULL
  for (k in 1:length(data.list)) {
    dt <- data.list[[k]]
    n <- length(dt[, 1])
    dt_resample = dt[sample(n, n, replace = TRUE), ]
    data_resample <- rbind(data_resample, dt_resample)
  }
  return(data_resample)
}

set.seed(223)
for (i in 1:1000) {
  
  data.hptn_resample <- resample(data.hptn.list)
  data.prep_resample <- resample(data.prep.list)
  
  # Regression-based Wald estimate
  md_D_Z_1 = glm(compliance ~ age_group + Employment + Education1  + illness3 + Syphilis, data = data.prep_resample, family = 'binomial', subset = (Arm == "TDF/FTC"))
  md_Y_Z_1 = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment + Education1  + illness3 + Syphilis, data = data.prep_resample, family = poisson(link=log), subset = (Arm == "TDF/FTC")) 
  md_Y_Z_0 = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment  + Education1 + illness3 + Syphilis, data = data.prep_resample, family = poisson(link=log), subset = (Arm == "Placebo")) 
  md_Y_Z_cab = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment + Education1 + illness3 + Syphilis, data = data.hptn_resample, family = poisson(link=log), subset = (Arm == "CAB-LA")) 
  md_Y_Z_truvada = glm(HIV ~ offset(log(Followup.time)) +  age_group + Employment + Education1 + illness3 + Syphilis, data = data.hptn_resample, family = poisson(link=log), subset = (Arm == "TDF/FTC")) 
  
  # Return a function estimate of CATE(X), point identification
  CATE_func <- function(dt){
    numerator = pmin(1, predict(md_Y_Z_1, newdata = dt, type = 'response')/dt$Followup.time) -
      pmin(1, predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time)
    denominator = predict(md_D_Z_1, newdata = dt, type = 'response') - 0
    return(numerator/denominator)
  }
  
  CATE_func_partial <- function(dt){
    numerator = pmin(1, predict(md_Y_Z_1, newdata = dt, type = 'response')/dt$Followup.time) -
      pmin(1, predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time)
    denominator = predict(md_D_Z_1, newdata = dt, type = 'response')
    lower.bound <- 0 - predict(md_Y_Z_0, newdata = dt, type = 'response')/dt$Followup.time
    upper.bound <- 0
    p <- denominator
    return(list(lower = p * numerator/denominator + (1 - p) * lower.bound, 
                upper = p * numerator/denominator + (1 - p) * upper.bound))
  }
  
  model.Dz.1 <- glm(compliance ~ age_group + Employment + Education1  + illness3 + Syphilis, data = data.hptn_resample, family = 'binomial', subset = (Arm == "TDF/FTC"))
  
  # p = E[D|Z=Placebo, X]
  cc_x <- function(dt, p = 0.05){
    predict(model.Dz.1, newdata = dt, type = 'response') - p
  }
  
  #ITT, point identification
  ITT0 = c(ITT0, mean(cc_x(data.hptn_resample, 0) * CATE_func(data.hptn_resample)) )
  ITT5 = c(ITT5, mean(cc_x(data.hptn_resample, 0.05) * CATE_func(data.hptn_resample)) )
  ITT10 = c(ITT10, mean(cc_x(data.hptn_resample, 0.1) * CATE_func(data.hptn_resample)) )
  
  #ITT, partial identification
  ITT0.lower <- c(ITT0.lower, mean(cc_x(data.hptn_resample, 0) * CATE_func_partial(data.hptn_resample)$lower) )
  ITT0.upper <- c(ITT0.upper, mean(cc_x(data.hptn_resample, 0) * CATE_func_partial(data.hptn_resample)$upper) )
  
  ITT5.lower <- c(ITT5.lower, mean(cc_x(data.hptn_resample, 0.05) * CATE_func_partial(data.hptn_resample)$lower) )
  ITT5.upper <- c(ITT5.upper, mean(cc_x(data.hptn_resample, 0.05) * CATE_func_partial(data.hptn_resample)$upper) )
  
  ITT10.lower <- c(ITT10.lower, mean(cc_x(data.hptn_resample, 0.1) * CATE_func_partial(data.hptn_resample)$lower) )
  ITT10.upper <- c(ITT10.upper, mean(cc_x(data.hptn_resample, 0.1) * CATE_func_partial(data.hptn_resample)$upper) )
  
  ITT_cab_truvada <- predict(md_Y_Z_cab, newdata = data.hptn_resample, type = 'response')/data.hptn_resample$Followup.time - 
    predict(md_Y_Z_truvada, newdata = data.hptn_resample, 
            type = 'response')/data.hptn_resample$Followup.time
  # ITT of cab-la to placebo
  ITT_cab_placebo <- c(ITT_cab_placebo, mean(ITT_cab_truvada + cc_x(data.hptn_resample, 0) * CATE_func(data.hptn_resample) ) )
  
  # hiv incidence for truvada
  in_truvada <- mean(predict(md_Y_Z_truvada, newdata = data.hptn_resample, type = 'response')/data.hptn_resample$Followup.time)
  
  # hiv incidence for placebo
  in_placebo <- c(in_placebo, in_truvada - mean(cc_x(data.hptn_resample, 0) * CATE_func(data.hptn_resample) ) )
  
  print(i)
}
# quantile(ITT_cab_placebo, c(0.025, 0.975))

quantile(ITT0, c(0.025, 0.975))
quantile(ITT5, c(0.025, 0.975))
quantile(ITT10, c(0.025, 0.975))


quantile(ITT0.lower, 0.025)
quantile(ITT0.upper, 0.975)
quantile(ITT5.lower, 0.025)
quantile(ITT5.upper, 0.975)
quantile(ITT10.lower, 0.025)
quantile(ITT10.upper, 0.975)

quantile(ITT_cab_placebo, c(0.025, 0.975))
quantile(in_placebo, c(0.025, 0.975))

#----------------------------------------
#             EIF estimator
#----------------------------------------
EIF1 <- function(dt1, dt, outcome_model = 'glm', p = 0.05){
  n_1 <- nrow(dt1)
  n <- nrow(dt)
  
  dt1$Z = ifelse(dt1$Arm == "TDF/FTC", 1, 0)
  dt$Z = ifelse(dt$Arm == "TDF/FTC", 1, 2)
  
  data <- rbind(dt, dt1)
  data$S <- c(rep(3, n), rep(1, n_1)) # 3 = NI; 1 = h1
  
  # estimate posteriors f(S = s|X)
  posterior.S <- function(data, s){
    dt <- data
    dt$S[dt$S != s] <- 0
    dt$S[dt$S == s] <- 1
    model.f = gam(S ~ age_group + Employment + Education1 + Syphilis, data = dt,
                  family = 'binomial')
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
      if(s == 1){
        md_D_S = glm(compliance ~ age_group + Employment + Education1 + Syphilis,
                     data = data, family = 'binomial', subset = (Arm == "Placebo") & (S == s)) 
        value = predict(md_D_S, newdata = data, type = 'response')
      }
      if(s == 3){
        value =  rep(p, n + n_1)}
      return(value)
    }
    mu1.D <- function(data, s){
      md_D_S = glm(compliance ~ age_group + Employment + Education1 + Syphilis, 
                   data = data, family = 'binomial', subset = (Arm == "TDF/FTC") & (S == s))
      value = predict(md_D_S, newdata = data, type = 'response')
      return(value)
    }
    mu0.Y <- function(data, s){
      md_Y_S = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment + Education1 + Syphilis, data = data, family = poisson(link=log), subset = (Arm == "Placebo") & (S == s)) 
      value = predict(md_Y_S, newdata = data, type = 'response')/data$Followup.time
      return(value)
    }
    mu1.Y <- function(data, s){
      md_Y_S = glm(HIV ~ offset(log(Followup.time)) + age_group + Employment + Education1 + Syphilis, data = data, family = poisson(link=log), subset = (Arm == "TDF/FTC") & (S == s))
      value = predict(md_Y_S, newdata = data, type = 'response')/data$Followup.time
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
  delta_D_3 <- delta_D(data, 3)
  delta_Y_1 <- delta_Y(data, 1)
  delta_1 <- delta_Y_1/delta_D_1
  posterior.S3 <-  posterior.S(data, 3)
  term1 <- (2*data$Z - 1)*(data$S == 1)/f.Z(data, 1) * 
    posterior.S3/posterior.S(data, 1) *
    delta_D_3/delta_D_1 *
    (data$HIV - mu0.Y(data, 1) - (data$compliance - mu0.D(data, 1))*delta_1)
  
  term2 <- (data$Z == 1)*(data$S == 3)/f.Z(data, 3) * 
    delta_1 *
    (data$HIV - mu0.D(data, 3) - delta_D_3*data$compliance)
  
  term3 <- (data$S == 3)*delta_1*delta_D_3
  
  k <- term1 + term2 + term3
  phi.est <- sum(k)/length(data[which(data$S == 3), 1])
  kappa <- n/(n + n_1)
  EIF <- 1/kappa*term1 + 1/kappa*term2 + 1/kappa*(term3 - (data$S == 3)*phi.est)
  EIF[which(abs(EIF) > 20)] = 0
  var <- mean(EIF^2)/(n + n_1)
  lower <- phi.est - qnorm(0.975)*sqrt(var)
  upper <- phi.est + qnorm(0.975)*sqrt(var)
  #phi.win <-  sum( Winsorize(k, probs = c(0.01, 0.99)))/length(data[which(data$S == 3), 1])
  
  return(list = c(phi.est = phi.est, lower = lower, upper = upper,
                  var = var
  ) )
}


data.prep$Education1 <- ifelse(data.prep$Education == "not complete primary school", "Yes", "No")
data.hptn$Education1 <- ifelse(data.hptn$Education == "not complete primary school", "Yes", "No")

data.prep1 <- data.prep[data.prep$Chlamydia == "Neg" & data.prep$Gonorrhea == "Neg" & data.prep$Trichomonas == "Neg", ]
data.hptn1 <- data.hptn[data.hptn$Chlamydia == "Neg" & data.hptn$Gonorrhea == "Neg" & data.hptn$Trichomonas == "Neg", ]


md_EIF <- EIF1(dt1 = data.prep1, dt = data.hptn1, outcome_model = 'glm', p = 0)
cat("p = 0, EIF estimator = ", md_EIF[1], ", CI = [", md_EIF[2], ", ", md_EIF[3], "]", "\n")

md_EIF <- EIF1(dt1 = data.prep1, dt = data.hptn1, outcome_model = 'glm', p = 0.05)
cat("p = 0.05, EIF estimator = ", md_EIF[1], ", CI = [", md_EIF[2], ", ", md_EIF[3], "]", "\n")

md_EIF <- EIF1(dt1 = data.prep1, dt = data.hptn1, outcome_model = 'glm', p = 0.1)
cat("p = 0.1, EIF estimator = ", md_EIF[1], ", CI = [", md_EIF[2], ", ", md_EIF[3], "]", "\n")





