
# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(mgcv)
library(rstan)
library(DHARMa)
library(posterior)
library(logitnorm)
library(foreach)       # parallelization of moments_logit_normal
library(doMC)
theme_set(theme_bw())

# Functions ---------------------------------------------------------------

# vectorized logit_norm mean
logit_norm_vect <- function(mu, sigma) {
  if(length(sigma) == 1) sigma <- rep(sigma, length(mu))
  pmean <- numeric(length(mu))

  for(i in 1:length(mu)) {
    pmean[i] <- momentsLogitnorm(mu[i], sigma[i])["mean"]
  }
  return(pmean)
}

# compute the mean of the logit-normal distribution in parallel.
# Intended to be used with a matrix of logit-means and matrix of logit-sd,
# of the same size. the sds can also be a vector with length == ncol(mu).
# This is thought for columns to be posterior samples.
logit_norm_mean <- function(mu, sigma, cores = parallel::detectCores() - 1) {

  registerDoMC(cores = cores)

  sigma_mat <- !is.null(dim(sigma))

  # turn columns into list elements
  if(sigma_mat) {
    arg_list <- lapply(1:ncol(mu), function(j) {
      r <- cbind(mu[, j], sigma[, j])
      colnames(r) <- c("mu", "sigma")
      return(r)
    })
  } else {
    arg_list <- lapply(1:ncol(mu), function(j) {
      r <- cbind(mu[, j], sigma[j])
      colnames(r) <- c("mu", "sigma")
      return(r)
    })
  }

  # compute moments in parallel
  means_list <- foreach(cc = arg_list) %dopar% {
    logit_norm_vect(cc[, "mu"], cc[, "sigma"])
  }

  return(do.call("cbind", means_list))
}

map_hdi <- function(x, ci = 0.95) {
  # x <- rnorm(1000); ci =  0.95
  hh <- bayestestR::hdi(x, ci = ci)
  res <- c(
    "map" = bayestestR::map_estimate(x),
    "hdi_lower" = hh$CI_low,
    "hdi_upper" = hh$CI_high
  )
  return(res)
}

# Prepare data ------------------------------------------------------------

d <- read.csv("dgapes_v1.csv", sep = ";")

# amend an error. Icterus_pyrrhopterus appear as gulper in three sources and
# as masher in only one. Consider it as gulper.
d$foraging.behaviour[d$bird.sp == "Icterus_pyrrhopterus"] <- "gulper"

#  remove NAs
d <- d[complete.cases(d[, c("foraging.behaviour", "fruitstime",
                            "gape.width", "frdiam",
                            "plant.sp", "bird.sp")]), ]

d$sizediff <- d$gape.width - d$frdiam
d$masstime <- d$fruitstime / d$frmass

# nrow(d)
# sum(is.na(d$frmass)) / nrow(d) # many NA in fruitmass, 12 %

# unique plant-bird interactions
d$plant.bird.sp <- paste(d$plant.sp, d$bird.sp, sep = "__")

# make factors
d$plant.sp <- factor(d$plant.sp, levels = unique(d$plant.sp))
d$bird.sp <- factor(d$bird.sp, levels = unique(d$bird.sp))
d$plant.bird.sp <- factor(d$plant.bird.sp, levels = unique(d$plant.bird.sp))
d$source <- factor(d$source, levels = unique(d$source))
d$foraging.behaviour <- factor(d$foraging.behaviour,
                               levels = c("gulper", "masher"),
                               labels = c("(a) Gulpers", "(b) Mashers"))

# Count Ns
d$plant.sp %>% unique %>% length # 89
d$bird.sp %>% unique %>% length  # 145
unique(d$plant.bird.sp) %>% length # 637 pairs
unique(d$source) %>% length # 39 studies
nrow(d) # 719

# some pairs of species have more than one measurement of consumption rate.
# View(d[d$plant.bird.sp == "Euterpe_edulis__Ramphastos_vitellinus", ])

t1 <- table(d$source)
t1 <- t1[order(t1, decreasing = T)]
t1rel <- t1 / sum(t1)
cumsum(t1rel)

# What if we assign zero consumption to bird-plant pairs recorded in the
# same study for which there are no pairs?
sources <- levels(d$source)

# What if we create those non-interacting pairs?
dcombs <- do.call("rbind", lapply(1:length(sources), function(s) {
  # s = 2
  print(s)
  sour <- sources[s]# s <- sources[2]
  dsub <- d[d$source == sour, ]
  dall <- expand.grid(plant.sp = unique(dsub$plant.sp),
                      bird.sp = unique(dsub$bird.sp))
  dall$plant.bird.sp <- paste(dall$plant.sp, dall$bird.sp, sep = "__")

  # keep only pairs of species absent in the database
  dall <- dall[!(dall$plant.bird.sp %in% unique(d$plant.bird.sp)), ]

  if(nrow(dall) == 0) return(NULL)

  # join to obtain data from plant and bird
  dplant <- aggregate(frdiam ~ plant.sp, dsub, mean)
  dbird <- aggregate(gape.width ~ bird.sp + foraging.behaviour, dsub, mean)

  dall_plant <- left_join(dall, dplant, by = "plant.sp")
  dout <- left_join(dall_plant, dbird, by = "bird.sp")

  dout$fruitstime <- 0
  dout$source <- sour
  dout$sizediff <- dout$gape.width - dout$frdiam

  return(dout)
}))

tcombs <- dcombs$source %>% table
# plot(tcombs)
# Gondim, bledinger y Correia son casi los únicos que aportan datos extra.

# merge with observed data
cols <- c("foraging.behaviour", "fruitstime",
          "gape.width", "frdiam", "sizediff",
          "plant.sp", "bird.sp", "plant.bird.sp",
          "source")

dbin <- rbind(d[, cols], dcombs[, cols])
dbin$cons_bin <- as.numeric(dbin$fruitstime > 0)
dbin$cons_fac <- factor(as.character(dbin$cons_bin), levels = c("0", "1"),
                        labels = c("consumed", "not consumed"))

dbin$plant.sp <- factor(dbin$plant.sp, levels = levels(d$plant.sp))
dbin$bird.sp <- factor(dbin$bird.sp, levels = levels(d$bird.sp))
dbin$plant.bird.sp <- factor(dbin$plant.bird.sp,
                             levels = levels(d$plant.bird.sp))
dbin$source <- factor(dbin$source, levels = levels(d$source))
dbin$foraging.behaviour <- factor(dbin$foraging.behaviour,
                                  levels = levels(d$foraging.behaviour))

# full data set
t2 <- table(dbin$source)
t2 <- t2[order(t2, decreasing = T)]
plot(t2)
head(t2)
t2rel <- t2 / sum(t2)
barplot(t2rel)
barplot(cumsum(t2rel))

# pseudo-absences
t2 <- table(dbin$source[dbin$cons_bin == 0])
t2 <- t2[order(t2, decreasing = T)]
plot(t2)
head(t2)
t2rel <- t2 / sum(t2)
barplot(t2rel)
barplot(cumsum(t2rel))

cumsum(t2)
cumsum(t2rel)

# Bayesian model ---------------------------------------------------------

# compile the stan model
smodel <- stan_model("consumption_probability.stan", verbose = F)

dbin$behav_num <- as.numeric(dbin$foraging.behaviour)
dbin$plant_num <- as.numeric(dbin$plant.sp)
dbin$bird_num <- as.numeric(dbin$bird.sp)

temp <- aggregate(sizediff ~ bird.sp + foraging.behaviour, dbin, mean)
group_mat <- model.matrix(~ foraging.behaviour - 1, temp)

# Data for stan is passed as a list. The names of its elements must match the
# variables defined in the data {} section of the .stan file.
sdata <- list(
  N = nrow(dbin), K = 2,
  Sb = max(dbin$bird_num), Sp = max(dbin$plant_num),
  y = dbin$cons_bin,
  sizediff = dbin$sizediff,
  group = dbin$behav_num,
  group_mat = group_mat,
  plant = dbin$plant_num,
  bird = dbin$bird_num,

  # parameters to define prior distributions
  prior_a_sd = 3,
  prior_b_sd = 10 / sd(dbin$sizediff), # this means that if sd(x) = 1, the prior sd would be 10
  prior_sigma_sd = 3,
  lower_u = 0.2
)

# sample the posterior
mp1 <- sampling(
  smodel, sdata, seed = 342534543, refresh = 100,
  iter = 4000, warmup = 1000, chains = 6, cores = 6, thin = 2,
  control = list(adapt_delta = 0.95)
)

saveRDS(mp1, "exports/consumption_probability_model.rds")
mp1 <- readRDS("exports/consumption_probability_model.rds")

# perfect, 1 min.
pairs(mp1, pars = c("b", "a", "u"))
pairs(mp1, pars = c("sigma_plant", "sigma_bird", "u"))

sm1 <- summary(mp1)[[1]]
min(sm1[, "n_eff"]) # 1172.925
max(sm1[, "Rhat"])  # 1.00753

# Extract parameters --------------------------------------------------

# extract parameters
ahat <- as.matrix(mp1, "a") %>% t
bhat <- as.matrix(mp1, "b") %>% t
uhat <- as.matrix(mp1, "u") %>% t
sigma_plant <- as.matrix(mp1, "sigma_plant") %>% as.numeric
sigma_bird <- as.matrix(mp1, "sigma_bird") %>% t
sigma_raneff <- rbind(
  sqrt(sigma_plant ^ 2 + sigma_bird[1, ] ^ 2), # gulpers total sd
  sqrt(sigma_plant ^ 2 + sigma_bird[2, ] ^ 2)  # mashers total sd
)
npost <- ncol(uhat)

e_bird_hat <- as.matrix(mp1, "e_bird") %>% t
e_plant_hat <- as.matrix(mp1, "e_plant") %>% t

# randomized random effects: same as "uncertainty" mode in brms:
# resample the fitted random effects.
bird_ids_rs <- sapply(1:npost, function(i) sample(1:nrow(e_bird_hat), size = nrow(e_bird_hat)))
plant_ids_rs <- sapply(1:npost, function(i) sample(1:nrow(e_plant_hat), size = nrow(e_plant_hat)))

e_bird_rs <- sapply(1:npost, function(i) e_bird_hat[bird_ids_rs[, i], i])
e_plant_rs <- sapply(1:npost, function(i) e_plant_hat[plant_ids_rs[, i], i])

# array of samples, to summarise later and add the R2
var_names <- c("a", "b", "u", "sigma_plant", "sigma_bird")
samples <- as_draws_matrix(mp1)
samples <- subset_draws(samples, var_names)

# Residuals ---------------------------------------------------------------

ysim <- matrix(NA, nrow(dbin), npost)
group <- sdata$group
plant <- sdata$plant
bird <- sdata$bird

# make design matrices
X_group <- model.matrix(~ foraging.behaviour - 1, dbin)
X_bird <- model.matrix(~ bird.sp - 1, dbin)
X_plant <- model.matrix(~ plant.sp - 1, dbin)

# simulate new random effects, but respecting hierarchy
e_bird_sim_raw <- matrix(rnorm(npost * sdata$Sb), sdata$Sb, npost)
e_bird_sim <- e_bird_sim_raw * (group_mat %*% sigma_bird)

e_plant_sim_raw <- matrix(rnorm(npost * sdata$Sp), sdata$Sp, npost)
e_plant_sim <- e_plant_sim_raw * (outer(rep(1, sdata$Sp), sigma_plant))

# new random effects ignoring hierarchy (non-hierarchy)
e_bird_nh_raw <- matrix(rnorm(npost * nrow(dbin)), nrow(dbin), npost)
e_plant_nh_raw <- matrix(rnorm(npost * nrow(dbin)), nrow(dbin), npost)

e_bird_nh <- e_bird_nh_raw * (X_group %*% sigma_bird)
e_plant_nh <- e_plant_nh_raw * (outer(rep(1, nrow(dbin)), sigma_plant))

# compute probability
p_raw <- plogis(
  X_group %*% ahat +
  (X_group * dbin$sizediff) %*% bhat +

  ## unconditional to raneffs, with hierarchy:
  # X_bird %*% e_bird_sim + X_plant %*% e_plant_sim

  ## conditioning on raneffs:
  X_bird %*% e_bird_hat + X_plant %*% e_plant_hat

  ## conditioning on raneffs with resampling ("uncertainty" mode in brms)
  #X_bird %*% e_bird_rs + X_plant %*% e_plant_rs

  ## unconditional to raneffs, without hierarchy
  # e_bird_nh + e_plant_nh
)
upper_rep <- X_group %*% uhat
psim <- p_raw * upper_rep

ysim_vec <- as.numeric(runif(npost * nrow(dbin)) <= as.vector(psim))
ysim <- matrix(ysim_vec, ncol = npost)

res <- createDHARMa(simulatedResponse = ysim,
                    observedResponse = dbin$cons_bin,
                    fittedPredictedResponse = rowMeans(ysim),
                    integerResponse = TRUE)
plot(res)

dbin$res <- res$scaledResiduals

ggplot(dbin, aes(sizediff, res)) +
  geom_point(alpha = 0.2, size = 1.8) +
  geom_smooth(method = "gam", method.args = list(family = mgcv::betar()),
              alpha = 0.3, linewidth = 0.5,
              fill = viridis(1, option = "C", begin = 0.2),
              color = viridis(1, option = "C", begin = 0.2)) +
  facet_wrap(vars(foraging.behaviour)) +
  ylab("Residuals") +
  xlab("Size difference (bill - fruit, mm)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 10))

ggsave("figures/consumption_probability_residuals_conditional.tiff",
       width = 1961, height = round(1961 * 0.55), units = "px",
       dpi = 300)


# unconditional: largely underestimated, no matter whether hierarchy is considered or not
# conditional: perfect.
# uncertainty: (as brms does), quite bad also.

# This is because the random effects distributions are poorly approximated by a
# Gaussian.

par(mfrow = c(1, 2))
plot(density(rowMeans(e_bird_hat), adjust = 1.5), main = "birds random effects")
plot(density(rowMeans(e_plant_hat), adjust = 1.5), main = "plants random effects")
par(mfrow = c(1, 1))

# Predictions -------------------------------------------------------------

pdata <- rbind(
  data.frame(
    sizediff = seq(min(dbin$sizediff[dbin$foraging.behaviour == "(a) Gulpers"]),
                   max(dbin$sizediff[dbin$foraging.behaviour == "(a) Gulpers"]),
                   length.out = 200),
    behav_num = 1
  ),
  data.frame(
    sizediff = seq(min(dbin$sizediff[dbin$foraging.behaviour == "(b) Mashers"]),
                   max(dbin$sizediff[dbin$foraging.behaviour == "(b) Mashers"]),
                   length.out = 200),
    behav_num = 2
  )
)

# add factor
pdata$foraging.behaviour <- factor(as.character(pdata$behav_num),
                                   levels = c("1", "2"),
                                   labels = c("(a) Gulpers", "(b) Mashers"))

X_group_pred <- model.matrix(~ foraging.behaviour - 1, pdata)
sigmas_mat_pred <- X_group_pred %*% sigma_raneff

mu_raw_logit_pred <- X_group_pred %*% ahat +
                     (X_group_pred * pdata$sizediff) %*% bhat
# mu_raw_logitnorm_pred <- logit_norm_mean(mu_raw_logit_pred,
#                                          sigmas_mat_pred)
# saveRDS(mu_raw_logitnorm_pred, "exports/consumption_probability_pred-mu-samples.rds")
mu_raw_logitnorm_pred <- readRDS("exports/consumption_probability_pred-mu-samples.rds")

mu_pred <- mu_raw_logitnorm_pred * (X_group_pred %*% uhat)

# summarize the posterior distribution of predictions
pdata$mu <- rowMeans(mu_pred)
pdata$mu_lower <- apply(mu_pred, 1, quantile, probs = 0.025, method = 8)
pdata$mu_upper <- apply(mu_pred, 1, quantile, probs = 0.975, method = 8)

# averaging data by bird species
dg <- dbin[dbin$foraging.behaviour == "(a) Gulpers", ]
dm <- dbin[dbin$foraging.behaviour == "(b) Mashers", ]

dg1 <- aggregate(cbind(cons_bin, sizediff) ~ bird.sp + foraging.behaviour, dg, mean)
dm1 <- aggregate(cbind(cons_bin, sizediff) ~ bird.sp + foraging.behaviour, dm, mean)

dagg <- rbind(dg1, dm1)
dagg$foraging.behaviour <- factor(dagg$foraging.behaviour,
                                  levels = c("(a) Gulpers", "(b) Mashers"))


# predictor densities
dens_g <- density(dbin$sizediff[dbin$foraging.behaviour == "(a) Gulpers"],
                  n = 2 ^ 10)
dens_m <- density(dbin$sizediff[dbin$foraging.behaviour == "(b) Mashers"],
                  n = 2 ^ 10)

dens_data <- rbind(
  data.frame(sizediff = dens_g$x, dens = dens_g$y,
             foraging.behaviour = "(a) Gulpers"),
  data.frame(sizediff = dens_m$x, dens = dens_m$y,
             foraging.behaviour = "(b) Mashers")
)

dens_data$foraging.behaviour <- factor(dens_data$foraging.behaviour,
                                       levels = c("(a) Gulpers", "(b) Mashers"))
dens_data$dens_scaled <- dens_data$dens * 0.9 / max(dens_data$dens)

# plot
viri_op <- "C"
ggplot(pdata, aes(sizediff, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(aes(sizediff, ymax = dens_scaled, ymin = 0), data = dens_data,
              color = NA, alpha = 0.2, inherit.aes = F,
              fill = viridis(1, option = viri_op, begin = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(color = NA, alpha = 0.3, fill = viridis(1, option = viri_op)) +
  geom_line(color = viridis(1, option = viri_op), linewidth = 0.5) +
  geom_point(aes(sizediff, cons_bin), data = dagg, inherit.aes = F,
             alpha = 0.4, size = 1.8) +
  facet_wrap(vars(foraging.behaviour)) +
  # ylim(0, 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 10)) +
  ylab("Consumption probability") +
  xlab("Size difference (bill - fruit, mm)") +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0.01, 0.01))

ggsave("figures/consumption_probability_predictions.tiff",
       width = 1961, height = round(1961 * 0.55), units = "px",
       dpi = 300)

# Width: 945 (single column), 1476 (1.5 column) or 1961 (double column) pixels
# (at 300 dpi). Resolution: 300-600 dpi.
# Size: <50 MB (for exceptions on file size, see below).

# R2 ----------------------------------------------------------------------

# for gulpers (1) and mashers (2)
r2mat <- matrix(NA, npost, 2)
colnames(r2mat) <- c("R2[1]", "R2[2]")
attr(r2mat, "nchains") <- attr(samples, "nchains")

# make design matrices
X_group <- model.matrix(~ foraging.behaviour - 1, dbin)

upper_mat <- X_group %*% uhat
sigma_mat <- X_group %*% sigma_raneff

# mu without variation among raneffs
mu_raw_logit <- X_group %*% ahat + (X_group * dbin$sizediff) %*% bhat
mu_raw_logitnorm <- logit_norm_mean(mu_raw_logit, sigma_mat)
saveRDS(mu_raw_logitnorm, "exports/consumption_probability_fitted-mu-samples.rds")
mu_raw_logitnorm <- readRDS("exports/consumption_probability_fitted-mu-samples.rds")
mu_sim <- mu_raw_logitnorm * upper_mat

var_y <- mu_sim * (1 - mu_sim)
# get var(mu) and mean(var) by group
mu_var <- aggregate(mu_sim ~ foraging.behaviour, dbin, var)[, -1] %>% as.matrix
var_mean <- aggregate(var_y ~ foraging.behaviour, dbin, mean)[, -1] %>% as.matrix

r2mat[, 1] <- mu_var[1, ] / (mu_var[1, ] + var_mean[1, ]) * 100
r2mat[, 2] <- mu_var[2, ] / (mu_var[2, ] + var_mean[2, ]) * 100

samples_full <- bind_draws(samples, r2mat)

# Coefficients table ------------------------------------------------------

summ1 <- summarise_draws(samples_full)
summ2 <- apply(samples_full, 2, map_hdi) %>% t
pgt0 <- apply(samples_full, 2, function(x) sum(x > 0) / length(x))
coef_table <- cbind(summ1[, 1], summ2, pgt0 = pgt0, summ1[, -1])
write.csv(coef_table, "exports/consumption_probability_table_summary.csv",
          row.names = F)

# Exports for plotting ---------------------------------------------------

plot_list <- list(
  pdata = pdata,
  dagg = dagg,
  dens_data = dens_data,
  data = dbin
)
saveRDS(plot_list, "exports/consumption_probability_plot_list.rds")
