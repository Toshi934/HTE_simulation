####### Simulation Code #######
### Last Update: 2025/05/17 ###
rm(list = ls(all.names = TRUE))

#########################################################################
### Install packages
library(tidyverse)
library(grf)
library(lmtest)
library(sandwich)
library(broom)
library(ggplot2)
library(gtsummary)
#########################################################################
### Test data generation

# Sample size
N <- 10000

# Treatment
set.seed(1)
treatment <- rbinom(N, 1, 0.5)

# Covariates
age <- rnorm(N, 40, 10) 
systolic.blood.pressure <- rnorm(N, 120, 10) 
hba1c <- rnorm(N, 5, 1) 
eGFR <- rnorm(N, 60, 20) 
eGFR <- pmin(pmax(eGFR, 10), 120)
medication <- rbinom(N, 1, 0.5)

# Parameters for HTE
b <- (0.1 + 0.05*hba1c + 0.05*eGFR + 0.3*medication) / 20

# Outcome
regY1 <- -0.1 + b*treatment + 0.0005*(systolic.blood.pressure - mean(systolic.blood.pressure)) + 0.01*hba1c + 0.001*(eGFR- mean(eGFR)) + 0.01*medication
regY1 <- pmin(pmax(regY1, 0), 1)
outcome <- rbinom(N, 1, regY1)

# Create a data frame
data <- data.frame(treatment, outcome, age, systolic.blood.pressure, hba1c,  eGFR, medication)

#########################################################################
### Set up for cross-fitting
num.rankings <- 5
num.folds <- 8

n = data %>% nrow()
folds <- sort(seq(n) %% num.folds) + 1

### Run causal forest ###
Y <-  data$outcome
W <-  data$treatment
X <- data[, c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")]
X <- as.matrix(X)

set.seed(1)
cf <- causal_forest(X = X, # covariate matrix
                           Y = Y, # outcome vector
                           W = W, # exposure vector
                           num.trees = 2000, # grow 10000 trees
                           tune.parameters = "all", # tune parameters
                           clusters = folds) 

# Obtain estimated CATE
data$tau.hat <- predict(cf, estimate.variance = F)$predictions # predict CATEs
tau.hat = predict(cf)$predictions


#########################################################################
### Test calibration performance of causal forest ###

### 1) BLP analysis ###
test_calibration(cf)



### 2) Calibration plot ###
# Obrain CATE quintile ranking
ranking <- rep(NA, n)
for (fold in seq(num.folds)){
  tau.hat.quintiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quintiles, include.lowest=TRUE,labels=seq(num.rankings))
}

# Computing AIPW scores 
tau.hat <- data$tau.hat 
e.hat <- cf$W.hat # P [W=1|X]
m.hat <- cf$Y.hat # E [Y|X]
# Estimating mu.hat(X, 1) and mu. hat (X, 0) for observations 
mu.hat.O <- m.hat - e.hat * tau.hat 
mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat
# AIPW scores
aipw.scores <- tau.hat + data[,"treatment"]/(e.hat)*(data[,"outcome"] - mu.hat.1) - (1-data[,"treatment"])/(1-e.hat)*(data[,"outcome"] - mu.hat.O)
ols <- lm(aipw.scores ~ 0 + factor(ranking))

forest.ate <- data.frame ("AIPW", paste0("Q", seq(length(unique(ranking)))), coeftest(ols, vcov = vcovHC(ols, "HC2")) [,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL #forest.ate

# Plot estimated ATEs by CATE quintile
calibration.plot <- ggplot(forest.ate) +
  aes (x = ranking, y = estimate, group = method, color = method) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = estimate -2 * std.err, ymax = estimate +2 * std.err), width = .2, position = position_dodge (0.2)) +
  ylab("ATE (percentage point)") + xlab("CATE Ranking") + 
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) 

print(calibration.plot)



### 3) Variable importance ###
covariate.list <- c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")

# Obtain variable importance 
var.imp <- variable_importance(cf) %>% 
  tidy() %>% 
  data.frame() %>% 
  mutate(varname = covariate.list) %>% 
  arrange(desc(x))

# Order the factor levels by value descending
var.imp$varname <- factor(var.imp$varname, levels = var.imp$varname[order(var.imp$x, decreasing = FALSE)])

# Barplot of variable importance
ggplot(var.imp, aes(x=varname, y=x)) + 
  geom_bar(stat = "identity") +
  coord_flip() + 
  ylab("Importance") + 
  xlab("Variables") +
  theme_bw() + 
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )



### 4) Comparison of CATE quintile ###

# Label CATE ranking
data$ranking <- ranking

# Summarize each CATE subgroup
data[, c("ranking", "age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")] %>% 
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"),
              digits = all_continuous() ~ 2,
              by = ranking) %>% 
  bold_labels()




