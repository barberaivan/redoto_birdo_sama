library(tidyverse)
library(viridis)
library(ggh4x)
library(deeptime) # ggarrange2, aligns everything
library(mgcv)
theme_set(theme_classic())

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.35),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}


l1 <- readRDS("exports/consumption_rate_plot_list.rds")
l2 <- readRDS("exports/consumption_probability_plot_list.rds")

viri_op <- "C"

for(i in 1:4) {
  l1[[i]]$foraging.behaviour <- factor(as.character(l1[[i]]$foraging.behaviour),
                                       levels = c("(a) Gulpers", "(b) Mashers"),
                                       labels = c("Gulpers", "Mashers"))
  l2[[i]]$foraging.behaviour <- factor(as.character(l2[[i]]$foraging.behaviour),
                                       levels = c("(a) Gulpers", "(b) Mashers"),
                                       labels = c("Gulpers", "Mashers"))
}

# Rate
pdata1 = l1$pdata
dagg1 = l1$dagg
dens_data1 = l1$dens_data
data1 = l1$data

# Probability
pdata2 = l2$pdata
dagg2 = l2$dagg
dens_data2 = l2$dens_data
data2 = l2$data


# Predictions -------------------------------------------------------------

ytop = 22

text1 <- data.frame(x = -18,
                    y = ytop * (1 - 1 / ytop),
                    lab = c("(a)", "(b)"),
                    foraging.behaviour = c("Gulpers", "Mashers"))

text_diff <- data.frame(
  x = c(15, -10),
  y = c(12.5, 12.5),
  name = c("bill > fruit", "bill < fruit"),
  foraging.behaviour = "Mashers"
)



prate <-
  ggplot(pdata1, aes(sizediff, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(aes(sizediff, ymax = dens_scaled, ymin = 0), data = dens_data1,
              color = NA, alpha = 0.2, inherit.aes = F,
              fill = viridis(1, option = viri_op, begin = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(color = NA, alpha = 0.3, fill = viridis(1, option = viri_op)) +
  geom_line(color = viridis(1, option = viri_op), linewidth = 0.5) +
  geom_point(aes(sizediff, fruitstime), data = dagg1, inherit.aes = F,
             alpha = 0.4, size = 1.8) +
  geom_text(aes(x, y, label = lab), data = text1, inherit.aes = F,
            size = 11/.pt) +
  geom_text(aes(x, y, label = name), data = text_diff, inherit.aes = F,
            alpha = 0.7, size = 10/.pt) +
  facet_wrap(vars(foraging.behaviour), axes = "all") +
  nice_theme() +
  theme(axis.title.x = element_blank()) +
  ylab("Consumption rate (fruits / min)") +
  scale_y_continuous(limits = c(0, ytop),
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(-20, 30))
prate

text2 <- data.frame(x = -18,
                    y = 1 * (1 - 1 / ytop),
                    lab = c("(c)", "(d)"),
                    foraging.behaviour = c("Gulpers", "Mashers"))

pprob <-
  ggplot(pdata2, aes(sizediff, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(aes(sizediff, ymax = dens_scaled, ymin = 0), data = dens_data2,
              color = NA, alpha = 0.2, inherit.aes = F,
              fill = viridis(1, option = viri_op, begin = 0.3)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_ribbon(color = NA, alpha = 0.3, fill = viridis(1, option = viri_op)) +
  geom_line(color = viridis(1, option = viri_op), linewidth = 0.5) +
  geom_point(aes(sizediff, cons_bin), data = dagg2, inherit.aes = F,
             alpha = 0.4, size = 1.8) +
  geom_text(aes(x, y, label = lab), data = text2, inherit.aes = F,
            size = 11/.pt) +
  facet_wrap(vars(foraging.behaviour), axes = "all") +
  nice_theme() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Consumption probability") +
  xlab("Bill-fruit size difference (mm)") +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(-20, 30))
pprob

pp <- ggarrange2(prate, pprob, nrow = 2)

ggsave("figures/predictions_both.tiff", plot = pp,
       width = 1961, height = round(1961 * 0.9), units = "px",
       dpi = 300)

# Width: 945 (single column), 1476 (1.5 column) or 1961 (double column) pixels
# (at 300 dpi). Resolution: 300-600 dpi.
# Size: <50 MB (for exceptions on file size, see below).

# Residuals ---------------------------------------------------------------

data1$Model <- "Consumption rate"
data2$Model <- "Consumption probability"

data_res <- rbind(
  data1[, c("res", "sizediff", "foraging.behaviour", "Model")],
  data2[, c("res", "sizediff", "foraging.behaviour", "Model")]
)

data_res$Model <- factor(data_res$Model,
                         levels = c("Consumption rate",
                                    "Consumption probability"))
data_res$model_tit <- "Model"

ggplot(data_res, aes(sizediff, res)) +
  geom_point(alpha = 0.2, size = 1.8) +
  geom_smooth(method = "gam", method.args = list(family = mgcv::betar()),
              formula = y ~ s(x, bs = "cr", k = 20),
              alpha = 0.3, linewidth = 0.5,
              fill = viridis(1, option = "C", begin = 0.2),
              color = viridis(1, option = "C", begin = 0.2)) +
  facet_nested(rows = vars(model_tit, Model), cols = vars(foraging.behaviour),
               axes = "all") +
  ylab("Residuals") +
  xlab("Bill-fruit size difference (mm)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  # facet_wrap(vars(foraging.behaviour), axes = "all") +
  nice_theme() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 10))

ggsave("figures/residuals_both.tiff",
       width = 1961, height = round(1961 * 0.9), units = "px",
       dpi = 300)
ggsave("figures/residuals_both.png",
       width = 1961, height = round(1961 * 0.9), units = "px",
       dpi = 300)
