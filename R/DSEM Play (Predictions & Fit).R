

library(dsem)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(phylopath)
library(dplyr)
library(ggdag)
library(dynlm)
library(ggplot2)
library(reshape)
library(ggplot2)
library(dplyr)
library(brms)
library(LaplacesDemon)


### testing with practice code/data
setwd("~/Desktop")
load("Biggs/DSEM/isle_royale.rda")
data(isle_royale)
data = ts( log(isle_royale[,2:3]), start=1959)

# original model from demo of select features
sem = "
  # Link, lag, param_name
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, MtoM
  wolves -> wolves, 1, WtoW
"
fit = dsem( sem = sem,
             tsdata = data,
             estimate_delta0 = FALSE,
             control = dsem_control(
               quiet = TRUE,
               getsd = TRUE) )
summary(fit)

# making a second model to test model comparisons
sem2 = "
  # Link, lag, param_name
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
"
fit2 = dsem( sem = sem2,
            tsdata = data,
            estimate_delta0 = FALSE,
            control = dsem_control(
              quiet = TRUE,
              getsd = TRUE) )
summary(fit2)


#### Task 1: assess prediction skill/fit for an INDIVIDUAL response variable (species)


# loo residuals function seems to produce predicted values for each observation... let's start there
loodf<-loo_residuals(fit, what="loo", track_progress=FALSE)
# it's odd to me the that the SE values are the same for each row
# in any event, this gives predicted and observed values that could be compared to one another

# for each response variable (species), get correlation coefficients between predicted and observed values
cor_df <- loodf %>%
  group_by(Var2) %>%
  summarise(cor_coeff = cor(obs, est, use = "complete.obs"))

# get RMSE for each response variable (species)
metrics_df <- loodf %>%
  group_by(Var2) %>%
  summarise(
    cor_coeff = cor(obs, est, use = "complete.obs"),
    RMSE = sqrt(mean((obs - est)^2, na.rm = TRUE))
  )

# plot observed vs predicted values for each response variable (species)
plot <- ggplot(loodf, aes(x = obs, y = est)) +
  geom_point(alpha = 0.7) +  # Scatter plot
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  
  facet_wrap(~ Var2, scales = "free") + 
  theme_minimal() +
  labs(x = "Observed Values", y = "Predicted Values", title = "Predictions vs Observations")

# add correlation coefficients as annotations
plot +
  stat_cor(aes(label = paste("r = ", round(..r.., 2), sep = "")), method = "pearson", label.x.npc = "left", label.y.npc = 0.95) +
  geom_text(data = metrics_df, aes(x = min(loodf$obs), y = max(loodf$est), label = paste0("RMSE = ", round(RMSE, 2))), hjust = 0, vjust = 1, size = 5)
  

# are my simple (and/or frequentist) approaches appropriate in this context? I am not sure?


#### Task 2: assess fit for a whole model, not just one response, (or better yet, multiple whole models)

# get predictions for each model
loodf<-loo_residuals(fit, what="loo", track_progress=FALSE)
loodf2<-loo_residuals(fit2, what="loo", track_progress=FALSE)

# i am not too familiar with WAIC values, but they may be appropriate in this (Bayesian?) context?.. although idk if they make sense for SEMs...
# here is my attempt, along with the help of chatgpt, to calculate and compare WAIC values between the two models
log_lik_approx <- dnorm(loodf$obs, mean = loodf$est, sd = loodf$se, log = TRUE)
lppd <- sum(log_lik_approx)
pWAIC <- var(log_lik_approx)
WAIC <- -2 * (lppd - pWAIC)
cat("Approximate WAIC:", WAIC, "\n")

log_lik_approx <- dnorm(loodf2$obs, mean = loodf2$est, sd = loodf2$se, log = TRUE)
lppd <- sum(log_lik_approx)
pWAIC <- var(log_lik_approx)
WAIC <- -2 * (lppd - pWAIC)
cat("Approximate WAIC:", WAIC, "\n")

# the first model has the lower WAIC model, so this one is the better fit


# RMSEA (root mean square error of approximation) has also been used to assess traditional SEM fit
# RMSEA is based on chi squared stats, and therefore cant be (traditionally) calculated with dsem output
# however, chatgpt helped me figure out approximate RMSEA based on the residual variance 
squared_errors <- (loodf$obs - loodf$est)^2
df <- nrow(residuals_df) - length(coef(fit))  
N <- nrow(residuals_df)
rmsea_approx <- sqrt(sum(squared_errors) / (df * (N - 1)))
cat("Approximate RMSEA:", rmsea_approx, "\n")

squared_errors <- (loodf2$obs - loodf2$est)^2
df <- nrow(residuals_df) - length(coef(fit))  
N <- nrow(residuals_df)
rmsea_approx <- sqrt(sum(squared_errors) / (df * (N - 1)))
cat("Approximate RMSEA:", rmsea_approx, "\n")

# the first model has the lower RMSEA value, aligning with the outcome of the WAIC test
# additionally, both models produced RMSEA values lower than 0.05, which is typically the standard for RMSEA values from SEMs
