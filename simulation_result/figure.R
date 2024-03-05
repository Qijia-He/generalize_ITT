#----------------------------------------------------#
#   Take Scenario_X1 when n = 5000 as an example     #
#----------------------------------------------------#
data <- read.csv("Scenario_X1_n5000.csv")

# 1 est_ground_truth, est_ground_truth_CI_lower, est_ground_truth_CI_upper,
# 4 naive_h1_est, naive_h1_CI_lower, naive_h1_CI_upper,
# 7 naive_h2_est, naive_h2_CI_lower, naive_h2_CI_upper,
# 10 est_wald,  wald_CI_lower, wald_CI_upper,
# 13 est_EIF,  EIF_CI_lower, EIF_CI_upper,
# 16 est_EIF_gam, lower.b, upper.b

library(ggplot2)
library(latex2exp)
k <- seq(1, 54, 3)
k1 <- k + 1
k2 <- k + 2

itt2 <- 0.1164
itt1 <- 0.1167
itt3 <- 0.1164

itt.list <- c(rep(itt1, 6), rep(itt2, 6), rep(itt3, 6))

# percentage bias %
mean = colMeans(data[, k])
round (100*(mean - itt.list)/itt.list, 1)

for (i in c(5, 10, 15)) {
  m <- data[, k[i]]
  outlier = sum(data[, k[i]] >1 | data[, k[i]] < -1)
  m[ m > 1] <- 1
  m[m < -1] <- -1
  rmse <- (m- itt.list[i])^2
  print(c(outlier, round (100*(mean(m) - itt.list[i])/itt.list[i], 1), round(mean(rmse), 3)) )
}

# 95% coverage
mat_lower <- sweep(data[, k1], 2, itt.list) # lower bound - itt
mat_upper <- sweep(data[, k2], 2, itt.list) # upper bound - itt
coverage = colMeans(mat_lower < 0 & mat_upper > 0)
round (100*coverage, 1)


# D.2: Sampling distributions of 6 estimators in Scenario X1 and Scenario Y1

data1 <- data[, k]
data2 <- as.vector(t(data1))
data3 <- data.frame(est = data2,
                    tp = rep(c('1. GT', '2. Naive 1', '3. Naive 2', '4. Wald', '5. EIF par', '6. EIF gam'), length(data2)/6),
                    overlap = rep(rep(c("1. Poor overlap", "2. Limited overlap", "3. Sufficient overlap"), each = 6), length(data2)/18))


data3[sapply(data3, is.character)] <- lapply(data3[sapply(data3, is.character)], as.factor)
data3$tp <- factor(data3$tp,labels=c('1. GT'=parse(text=TeX('$\\widehat{ITT}_{hypo}$')),
                                     '2. Naive 1'=parse(text=TeX('$\\widehat{ITT}_{const, 1}$')),
                                     '3. Naive 2'=parse(text=TeX('$\\widehat{ITT}_{const, 2}$')),
                                     '4. Wald'=parse(text=TeX('$\\widehat{ITT}_{reg, par}$')),
                                     '5. EIF par'=parse(text=TeX('$\\widehat{ITT}_{EIF, par}')),
                                     '6. EIF gam'=parse(text=TeX('$\\widehat{ITT}_{EIF, gam}'))))

GT500 <- data.frame(overlap = c("1. Poor overlap", "2. Limited overlap", "3. Sufficient overlap"),
                    true_ITT = c(itt1, itt2, itt3))

GT500[sapply(GT500, is.character)] <- lapply(GT500[sapply(GT500, is.character)], as.factor)


# D.2: Sampling distributions of 6 estimators in Scenario X1 and Scenario Y1

ggplot(data = data3, aes(x = est)) + geom_histogram() + 
  facet_grid(tp ~ overlap, labeller = labeller(tp = label_parsed)) +  
  xlim(-0.1, 0.4) +
  geom_vline(data = GT500, aes(xintercept = true_ITT), color = 'red', size = 1.5, linetype = 'dashed') +
  theme_bw(base_size = 22) + xlab('')



# Web Appendix D: Additional simulation details
# D.1: Overlap in the simulated datasets

library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

overlap = 1
n = 2000
# This function can be found in simulation_Scenario_X1.R, real data is required to run the following code.
d = generate_data_binary(overlap = overlap, n = n, D_md1, D_md2, Y_md11, Y_md12, Y_md13)
dt_h1 = d$dt_h1
dt = d$dt
names(dt_h1)[7:9] <- c("Z", "D", "Y")
data <- rbind(dt_h1, dt)
data$Population <- c(rep("1. D_h1", n), rep("2. D_target", n) )

# Create the density plot
density1 <-
  ggplot(data, aes(x = ps, fill = Population)) +
  geom_density(alpha = 0.5) +
  labs(fill = "Population") +
  scale_fill_discrete(labels = c(expression(D[h1]), expression(D[target]))) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  theme(legend.title = element_text(family = "Arial", size = 27),
        legend.text = element_text(family = "Arial", size = 27),
        legend.position="top",
        text=element_text(size=30))+
  labs(title="Poor overlap",x="Probability of trial participation", y = "Density")

overlap = 2
n = 2000
d = generate_data_binary(overlap = overlap, n = n, D_md1, D_md2, Y_md11, Y_md12, Y_md13)
dt_h1 = d$dt_h1
# dt_h2 = d$dt_h2
dt = d$dt
names(dt_h1)[7:9] <- c("Z", "D", "Y")
# names(dt_h2)[7:9] <- c("Z", "D", "Y")
data <- rbind(dt_h1, dt)
data$Population <- c(rep("1. D_h1", n), rep("2. D_target", n) )
density2 <- 
  ggplot(data, aes(x = ps, fill = Population)) +
  geom_density(alpha = 0.5) +
  labs(fill = "Population") +
  scale_fill_discrete(labels = c(expression(D[h1]), expression(D[target]))) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  theme(legend.title = element_text(family = "Arial", size = 27),
        legend.text = element_text(family = "Arial", size = 27),
        legend.position="top",
        text=element_text(size=30))+
  labs(title="Limited overlap",x="Probability of trial participation", y = "Density")


overlap = 3
n = 2000
d = generate_data_binary(overlap = overlap, n = n, D_md1, D_md2, Y_md11, Y_md12, Y_md13)
dt_h1 = d$dt_h1
# dt_h2 = d$dt_h2
dt = d$dt
names(dt_h1)[7:9] <- c("Z", "D", "Y")
# names(dt_h2)[7:9] <- c("Z", "D", "Y")
data <- rbind(dt_h1, dt)
data$Population <- c(rep("1. D_h1", n), rep("2. D_target", n) )
data$Population <- c(rep("1. \\textit{D}_{h1}", n), rep("2. \\textit{D}_{target}", n))

density3 <- 
  ggplot(data, aes(x = ps, fill = Population)) +
  geom_density(alpha = 0.5) +
  labs(fill = "Population") +
  scale_fill_discrete(labels = c(expression(D[h1]), expression(D[target]))) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  theme(legend.title = element_text(family = "Arial", size = 27),
        legend.text = element_text(family = "Arial", size = 27),
        legend.position="top",
        text=element_text(size=30))+
  labs(title="Sufficient overlap",x="Probability of trial participation", y = "Density")


p = grid.arrange(density1, density2, density3, nrow = 1)

ggsave("3density.pdf", plot = p, width = 20, height = 8)
