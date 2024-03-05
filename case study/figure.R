library(survminer)
library(survival)
library(ggsurvfit)
library(cmprsk)


# Partner PrEP data
data.prep0 <- data.prep[, 4:6]
kmfit.prep <- survfit(Surv(Followup.time, HIV == 1) ~ Arm, data = data.prep0) # 2 - truvada; 1 - CAB

p.prep <-  ggsurvplot(kmfit.prep, pval = F, size = 0.1, cumevents = F, risk.table = TRUE,
                      fun = "event", palette = c("#FF3399", "#2E9FDF"), ggtheme = theme_bw(base_size = 22), 
                      ylim = c(0, 0.08),
                      xlab = "Years since Enrollment", ylab = "Cumulative Incidence \n (person-yr)",
                      title = "Partners PrEP - Female",
                      legend.title="Arm",
                      legend.labs = c("Placebo", "TDF/FTC"),
                      risk.table.fontsize = 8,
                      risk.table.height = 0.3)

# HPTN 084 data (two figures)
splots <- list()
rand <- read.csv("rand.csv")
rand <- rand[, 2:3]
names(rand) <- c("id", "Arm")
outcome.hptn <- read.csv("HPTN084_outcome.csv")[, -1]
names(outcome.hptn) <- c("HIV", "Followup.time", "id")
outcome.hptn$Followup.time <- outcome.hptn$Followup.time/365.25
data.hptn.full <- merge(rand, outcome.hptn, by = "id")
kmfit.hptn.full <- survfit(Surv(Followup.time, HIV == 1) ~ Arm, data = data.hptn.full) # 2 - truvada; 1 - CAB
# E.1 Kaplan-Meier estimates in HPTN 084
# Figure S5: Kaplan-Meier estimates of incident HIV acquisition in the HPTN 084 trial population
p.hptn <- 
  ggsurvplot(kmfit.hptn.full, pval = F, size = 0.1, cumevents = F, risk.table = TRUE,
             fun = "event", palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(base_size = 22), 
             ylim = c(0, 0.08),
             xlab = "Years since Enrollment", ylab = "Cumulative Incidence \n (person-yr)",
             title = "HPTN 084",
             legend.title="Arm",
             legend.labs = c("CAB-LA", "TDF/FTC"),
             risk.table.fontsize = 8,
             risk.table.height = 0.3)


data.hptn0 <- data.hptn[, c(4, 12, 13)]
kmfit.hptn <- survfit(Surv(Followup.time, HIV == 1) ~ Arm, data = data.hptn0) # 2 - truvada; 1 - CAB
p.hptn.target <- 
  ggsurvplot(kmfit.hptn, pval = F, size = 0.1, cumevents = F, risk.table = TRUE,
             fun = "event", palette = c("#E7B800", "#2E9FDF"), ggtheme = theme_bw(base_size = 22), 
             ylim = c(0, 0.08),
             xlab = "Years since Enrollment", ylab = "Cumulative Incidence \n (person-yr)",
             title = "HPTN 084 - Target Population",
             legend.title="Arm",
             legend.labs = c("CAB-LA", "TDF/FTC"),
             risk.table.fontsize = 8,
             risk.table.height = 0.3)

splots[[1]] <- p.prep
splots[[2]] <- p.hptn.target

# Figure 3: Left panel: Kaplan-Meier estimates of incident HIV acquisition in the Partners PrEP
# study. Right panel: Kaplan-Meier estimates of incident HIV acquisition among HPTN 084 partici-
#   pants whose partners either living with HIV or having an unknown HIV status (target population)
p = arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 1, risk.table.height = 0.3)
ggsave("cumulative_incidence.pdf", plot = p, width = 20, height = 8)





# Plot the probability of selection as in Cole and Stuart and make a comment.
# E.3: Overlap of HPTN 084 and Partners PrEP
# Figure S5: The probability of sample selection which is the probability of a participant to be in the HPTN 084 - Target Population given the covariates.
posterior.S <- function(data, s){
  dt <- data
  dt$S[dt$S != s] <- 0
  dt$S[dt$S == s] <- 1
  model.f = glm(S ~ age_group + Employment + Education + illness3 + Syphilis, data = dt, family = 'binomial')
  posterior = predict(model.f, newdata = data, type = 'response')
  return(posterior)
}

data.hptn$S <- 1
data.prep$S <- 0
dt.combine <- rbind(data.hptn, data.prep)
dt.combine$S <- as.factor(dt.combine$S)
dt.combine$ps <- posterior.S(dt.combine, 1)
dt.combine$Population <- ifelse(dt.combine$S == 1, "HPTN 084 - Target Population", "Partners PrEP - Female")

density.p <- ggdensity(dt.combine, x = "ps", 
                       fill = "Population", palette = "jco",
                       xlab = "Probability of sample selection",
                       ylab = "Density",
                       font.label = list(size = 24, face = "plain"))
density.p



