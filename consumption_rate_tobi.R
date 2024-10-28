# Packages ----------------------------------------------------------------
library(devtools)
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
setwd("C:/Users/tobias/Documents/POSDOC CONICET/diferencia ingesta/manuscrito/analsisis")
d <- read.csv("dgapes_v1.csv", sep = ";")

# ammend one error. Icterus_pyrrhopterus appear as gulper in three sources and
# as masher in only one. Consider it as gulper.
d$foraging.behaviour[d$bird.sp == "Icterus_pyrrhopterus"] <- "gulper"

summary(d$fruitstime[d$foraging.behaviour == "gulper"])

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
d$plant.sp %>% unique %>% length # 95
d$bird.sp %>% unique %>% length  # 158
unique(d$plant.bird.sp) %>% length # 735 pairs
unique(d$source) %>% length # 41 studies
nrow(d) # 819

a1 <- aggregate(fruitstime ~ bird.sp + foraging.behaviour, d, mean)
table(a1$foraging.behaviour)
# (a) Gulpers (b) Mashers
# 107          51

aggregate(sizediff ~ foraging.behaviour, d, length)
# foraging.behaviour sizediff
# 1        (a) Gulpers      549
# 2        (b) Mashers      270

sum(d$sizediff[d$foraging.behaviour == "(a) Gulpers"] < 0) /
  sum(d$foraging.behaviour == "(a) Gulpers") * 100
# 7.84 % (40 / 510)
sum(d$sizediff[d$foraging.behaviour == "(b) Mashers"] < 0) /
  sum(d$foraging.behaviour == "(b) Mashers") * 100
# 18.85 % (46 / 244)

summary(d$fruitstime[d$foraging.behaviour == "(a) Gulpers"])
quantile(d$fruitstime[d$foraging.behaviour == "(a) Gulpers"],
         probs = seq(0, 1, by = 0.05))
sort(d$fruitstime[d$foraging.behaviour == "(a) Gulpers"])
plot(ecdf(d$fruitstime[d$foraging.behaviour == "(a) Gulpers"]))


summary(d$fruitstime[d$foraging.behaviour == "(b) Mashers"])
quantile(d$fruitstime[d$foraging.behaviour == "(b) Mashers"],
         probs = seq(0, 1, by = 0.05))

# some pairs of species have more than one measurement of consumption rate.
# View(d[d$plant.bird.sp == "Euterpe_edulis__Ramphastos_vitellinus", ])

tbird <- table(d$bird.sp) %>% sort
tplant <- table(d$plant.sp) %>% sort
tsource <- table(d$source) %>% sort

barplot(sort(tbird)); abline(h = 1, col = "red")
barplot(sort(tplant)); abline(h = 1, col = "red")
barplot(sort(tsource)); abline(h = 1, col = "red")


# REDOTO EDIT HERE --------------------------------------------------------

# REDOTO, acá podés meter los histogramas.
# Si querés exportar los datos, usá el data.frame "d":
# write.csv(d, "datos_limpios.csv")




#############################################
##phylogenetic distances
#birds
library(ape)

treebird<- read.nexus("output.nex") #load the phylogenetic tree

plot.phylo(treebird[[55]], type= "fan") #let's see how horrible it is
plot.phylo
abird<- multi2di(treebird)#solve the mutichotomies with a random procedure

##take the phylogenetic correlation matrices for all the hypothetic trees
treelist = vector('list', 3)
for (i in 1:length(abird)){

  treelist[[i]] = as.matrix(vcv.phylo(abird[[i]], cor=TRUE))
}
## average the phylogenetic correlation into one matrix
bircor = apply(simplify2array(treelist), 1:2, mean)


##plants


plantas<-read.table("plantas.txt", sep=";", header = TRUE)
plant<- as.data.frame(unique(plantas$plant.sp))
colnames(plant)<- "plant.sp"
genus<- matrix(nrow= nrow(plant))
for( i in 1:nrow(plant)){
  tmp = strsplit(as.character(plant$plant.sp[i]), " ") ##take the genus
  genus[i,] = sort(substring(tmp[[1]][1], 1, ))

}


devtools::install_github("ropensci/taxize")##need this to get the families
library(taxize)
plant$genus<- genus


plant$family <- tax_name(plant$genus, get = 'family', db = 'ncbi')$family  ##get the families


devtools::install_github("jinyizju/V.PhyloMaker2") ###get the megatree of plants and the tools to prune it

library(V.PhyloMaker2)

plantree<- phylo.maker(plant, scenarios = "S1")##use the plant data frame to get the trees
plot.phylo(plantree[[1]], type="fan")#see how horrible it is
plantree<- multi2di(plantree[[1]])#solve multichotomies (luckilly they are not so much)

placor<- as.matrix(vcv.phylo(plantree, cor=TRUE)) #get the phylogenetic correlation matrix

save(bircor, placor, file="filogenias.R") #save both birds and plant matrices
load("filogenias.R") #load them to happyly analyse the data

# Prior predictive check -------clr# Prior predictive check --------------------------------------------------

# check prior for phi, the dispersion parameter of the gamma distribution.

# The gamma distribution is usually parameterized in terms of shape and rate
# (alpha and beta), but it is easiear to define it in terms of mu = mean and
# phi = dispersion parameter.

# x ~ gamma(shape = 1 / phi,
#           rate = 1 / (mu * phi))

# mean(x) = mu
# var(x) = mu ^ 2 * phi
# phi = var(x) / mu ^ 2

# define a half-normal prior for phi, but choose its sigma:
prior_phi_sd <- 10

# # plot gamma density sampling from that prior, at at the observed mean
# mu <- mean(d$fruitstime)
# phi1 <- abs(rnorm(1, 0, prior_phi_sd))
# curve(dgamma(x, shape = 1 / phi1, rate = 1 / (mu * phi1)), to = 50,
#       ylim = c(0, 0.5))
#
# # repeat many times to see variability
# for(i in 1:500) {
#   phi1 <- abs(rnorm(1, 0, prior_phi_sd))
#   curve(dgamma(x, shape = 1 / phi1, rate = 1 / (mu * phi1)),
#         add = T, col = rgb(0, 0, 0, 0.05))
# }

# densities peaking at the right are those concentrated on the mean (low phi);
# they are somewhat symmetric, with mean = mode.
# densities peaking at the left are those with high variance (high phi),
# with mode << mean.

# Bayesian model -------------------------------------------------------------

# compile the stan model
smodel <- stan_model("consumption_rate.stan", verbose = F) # use this

d$behav_num <- as.numeric(d$foraging.behaviour)
d$plant_num <- as.numeric(d$plant.sp)
d$bird_num <- as.numeric(d$bird.sp)
d$source_num <- as.numeric(d$source)

temp <- aggregate(sizediff ~ bird.sp + foraging.behaviour, d, mean)
group_mat <- model.matrix(~ foraging.behaviour - 1, temp)

# Data for stan is passed as a list. The names of its elements must match the
# variables defined in the data {} section of the .stan file.
sdata <- list(
  N = nrow(d), K = 2,
  Sp = max(d$plant_num),
  Sb = max(d$bird_num),
  Ns = max(d$source_num),

  y = d$fruitstime,
  sizediff = d$sizediff,
  group = d$behav_num,
  group_mat = group_mat,
  plant = d$plant_num,
  bird = d$bird_num,
  source = d$source_num,

  # parameters to define prior distributions
  prior_a_sd = 3,
  prior_b_sd = 10 / sd(d$sizediff), # this means that if sd(x) = 1, the prior sd would be 10
  prior_phi_sd = 3,
  prior_sigma_sd = 3,
  lower_u = 1,
  upper_u = 30#quantile(d$fruitstime, 0.95)
)

# sample the posterior
mr1 <- sampling(
  smodel, sdata, seed = 342534543, refresh = 100,
  #iter = 5, chains = 1, cores = 1, # for testing
  iter = 3400, warmup = 1000, chains = 8, cores = 8, thin = 3,
  control = list(adapt_delta = 0.98)
)
saveRDS(mr1, "exports/consumption_rate_model.rds")
mr1 <- readRDS("exports/consumption_rate_model.rds")

# pairs(mr1, pars = c("u", "a", "b")) # hard to identify alpha and u
# pairs(mr1, pars = c("phi", "sigma_bird", "sigma_plant", "sigma_source"))

sm1 <- summary(mr1)[[1]]
min(sm1[, "n_eff"]) # 3045.712
max(sm1[, "Rhat"])  # 1.003287

# Extract parameters --------------------------------------------------

ahat <- as.matrix(mr1, "a") %>% t
bhat <- as.matrix(mr1, "b") %>% t
uhat <- as.matrix(mr1, "u") %>% t

phihat <- as.matrix(mr1, "phi") %>% t

sigma_plant <- as.matrix(mr1, "sigma_plant") %>% as.numeric
sigma_bird <- as.matrix(mr1, "sigma_bird") %>% t
sigma_source <- as.matrix(mr1, "sigma_source") %>% as.numeric
sigma_raneff <- rbind(
  sqrt(sigma_plant ^ 2 + sigma_source ^ 2 + sigma_bird[1, ] ^ 2), # gulpers total sd
  sqrt(sigma_plant ^ 2 + sigma_source ^ 2 + sigma_bird[2, ] ^ 2)  # mashers total sd
)

e_bird_hat <- as.matrix(mr1, "e_bird") %>% t
e_plant_hat <- as.matrix(mr1, "e_plant") %>% t
e_source_hat <- as.matrix(mr1, "e_source") %>% t

npost <- ncol(uhat)

# array of samples, to summarise later and add the R2
var_names <- c("a", "b", "u", "phi", "sigma_plant", "sigma_bird", "sigma_source")
samples <- as_draws_matrix(mr1)
samples <- subset_draws(samples, var_names)

# Residuals ---------------------------------------------------------------

ysim <- matrix(NA, nrow(d), npost)

# make design matrices
X_group <- model.matrix(~ foraging.behaviour - 1, d)
X_bird <- model.matrix(~ bird.sp - 1, d)
X_plant <- model.matrix(~ plant.sp - 1, d)
X_source <- model.matrix(~ source - 1, d)

# simulate new random effects, but respecting hierarchy
e_bird_sim_raw <- matrix(rnorm(npost * sdata$Sb), sdata$Sb, npost)
e_bird_sim <- e_bird_sim_raw * (group_mat %*% sigma_bird)

e_plant_sim_raw <- matrix(rnorm(npost * sdata$Sp), sdata$Sp, npost)
e_plant_sim <- e_plant_sim_raw * (outer(rep(1, sdata$Sp), sigma_plant))

e_source_sim_raw <- matrix(rnorm(npost * sdata$Ns), sdata$Ns, npost)
e_source_sim <- e_source_sim_raw * (outer(rep(1, sdata$Ns), sigma_source))

# compute probability
mu_raw <- plogis(
  X_group %*% ahat +
  (X_group * d$sizediff) %*% bhat +

  # # unconditional to raneffs, with hierarchy:
  # X_bird %*% e_bird_sim + X_plant %*% e_plant_sim + X_source %*% e_source_sim

  # conditioning on raneffs:
  X_bird %*% e_bird_hat + X_plant %*% e_plant_hat + X_source %*% e_source_hat
)

upper_rep <- X_group %*% uhat
mu_sim <- mu_raw * upper_rep

phi_mat <- X_group %*% phihat

for(i in 1:npost) {
  alpha <- 1 / phi_mat[, i]
  beta <- 1 / (phi_mat[, i] * mu_sim[, i])
  ysim[, i] <- rgamma(nrow(d), alpha, beta)
}

res <- createDHARMa(simulatedResponse = ysim,
                    observedResponse = d$fruitstime,
                    fittedPredictedResponse = rowMeans(ysim),
                    integerResponse = TRUE)
plot(res)

d$res <- res$scaledResiduals

ggplot(d, aes(sizediff, res)) +
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

# unconditional: largely underestimated, no matter whether hierarchy is considered or not
# conditional: perfect.

# This is because the random effects distributions are poorly approximated by a
# Gaussian.

par(mfrow = c(2, 2))
plot(density(rowMeans(e_bird_hat), adjust = 1.5), main = "birds random effects")
plot(density(rowMeans(e_plant_hat), adjust = 1.5), main = "plants random effects")
plot(density(rowMeans(e_source_hat), adjust = 1.5), main = "source random effects")
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
hist(rowMeans(e_bird_hat), breaks = 20, main = "birds random effects")
hist(rowMeans(e_plant_hat), breaks = 15, main = "plants random effects")
hist(rowMeans(e_source_hat), breaks = 15, main = "source random effects")
par(mfrow = c(1, 1))

# Predictions -------------------------------------------------------------

pdata <- rbind(
  data.frame(
    sizediff = seq(min(d$sizediff[d$foraging.behaviour == "(a) Gulpers"]),
                   max(d$sizediff[d$foraging.behaviour == "(a) Gulpers"]),
                   length.out = 200),
    behav_num = 1
  ),
  data.frame(
    sizediff = seq(min(d$sizediff[d$foraging.behaviour == "(b) Mashers"]),
                   max(d$sizediff[d$foraging.behaviour == "(b) Mashers"]),
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
# saveRDS(mu_raw_logitnorm_pred, "exports/consumption_rate_pred-mu-samples.rds")
mu_raw_logitnorm_pred <- readRDS("exports/consumption_rate_pred-mu-samples.rds")

mu_pred <- mu_raw_logitnorm_pred * (X_group_pred %*% uhat)

# summarize the posterior distribution of predictions
pdata$mu <- rowMeans(mu_pred)
pdata$mu_lower <- apply(mu_pred, 1, quantile, probs = 0.025, method = 8)
pdata$mu_upper <- apply(mu_pred, 1, quantile, probs = 0.975, method = 8)

# averaging by bird
dg <- d[d$foraging.behaviour == "(a) Gulpers", ]
dm <- d[d$foraging.behaviour == "(b) Mashers", ]

dg1 <- aggregate(cbind(fruitstime, sizediff) ~ bird.sp + foraging.behaviour, dg, mean)
dm1 <- aggregate(cbind(fruitstime, sizediff) ~ bird.sp + foraging.behaviour, dm, mean)

dagg <- rbind(dg1, dm1)
dagg$foraging.behaviour <- factor(dagg$foraging.behaviour,
                                  levels = c("(a) Gulpers", "(b) Mashers"))

# predictor densities
dens_g <- density(d$sizediff[d$foraging.behaviour == "(a) Gulpers"],
                  n = 2 ^ 10)
dens_m <- density(d$sizediff[d$foraging.behaviour == "(b) Mashers"],
                  n = 2 ^ 10)

dens_data <- rbind(
  data.frame(sizediff = dens_g$x, dens = dens_g$y,
             foraging.behaviour = "(a) Gulpers"),
  data.frame(sizediff = dens_m$x, dens = dens_m$y,
             foraging.behaviour = "(b) Mashers")
)

dens_data$foraging.behaviour <- factor(dens_data$foraging.behaviour,
                                       levels = c("(a) Gulpers", "(b) Mashers"))

ytop <- 22
dens_data$dens_scaled <- dens_data$dens * ytop * 0.9 / max(dens_data$dens)
viri_op <- "C"

ggplot(pdata, aes(sizediff, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(aes(sizediff, ymax = dens_scaled, ymin = 0), data = dens_data,
              color = NA, alpha = 0.2, inherit.aes = F,
              fill = viridis(1, option = viri_op, begin = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(color = NA, alpha = 0.3, fill = viridis(1, option = viri_op)) +
  geom_line(color = viridis(1, option = viri_op), linewidth = 0.5) +
  geom_point(aes(sizediff, fruitstime), data = dagg, inherit.aes = F,
             alpha = 0.4, size = 1.8) +
  facet_wrap(vars(foraging.behaviour)) +
  # ylim(0, 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 10)) +
  ylab("Consumption rate (fruits / min)") +
  xlab("Size difference (bill - fruit, mm)") +
  scale_y_continuous(limits = c(0, ytop),
                     expand = c(0.01, 0.01))

# R2 ----------------------------------------------------------------------

# for gulpers (1) and mashers (2)
r2mat <- matrix(NA, npost, 2)
colnames(r2mat) <- c("R2[1]", "R2[2]")
attr(r2mat, "nchains") <- attr(samples, "nchains")

# make design matrices
X_group <- model.matrix(~ foraging.behaviour - 1, d)

upper_mat <- X_group %*% uhat
phi_mat <- X_group %*% phihat
sigma_mat <- X_group %*% sigma_raneff

# mu without variation among raneffs
mu_raw_logit <- X_group %*% ahat + (X_group * d$sizediff) %*% bhat
# mu_raw_logitnorm <- logit_norm_mean(mu_raw_logit, sigma_mat)
# saveRDS(mu_raw_logitnorm, "exports/consumption_rate_pred-mu-samples.rds")
mu_raw_logitnorm <- readRDS("exports/consumption_rate_pred-mu-samples.rds")
mu_sim <- mu_raw_logitnorm * upper_mat

var_y <- mu_sim ^ 2 * phi_mat
# get var(mu) and mean(var) by group
mu_var <- aggregate(mu_sim ~ foraging.behaviour, d, var)[, -1] %>% as.matrix
var_mean <- aggregate(var_y ~ foraging.behaviour, d, mean)[, -1] %>% as.matrix

r2mat[, 1] <- mu_var[1, ] / (mu_var[1, ] + var_mean[1, ]) * 100
r2mat[, 2] <- mu_var[2, ] / (mu_var[2, ] + var_mean[2, ]) * 100

samples_full <- bind_draws(samples, r2mat)

# Coefficients table ------------------------------------------------------

summ1 <- summarise_draws(samples_full)
summ2 <- apply(samples_full, 2, map_hdi) %>% t
pgt0 <- apply(samples_full, 2, function(x) sum(x > 0) / length(x))
coef_table <- cbind(summ1[, 1], summ2, pgt0 = pgt0, summ1[, -1])
write.csv(coef_table, "exports/consumption_rate_table_summary.csv",
          row.names = F)

# Exports for plotting ---------------------------------------------------

plot_list <- list(
  pdata = pdata,
  dagg = dagg,
  dens_data = dens_data,
  data = d
)
saveRDS(plot_list, "exports/consumption_rate_plot_list.rds")
