library(tidyverse)
library(viridis)
theme_set(theme_light())

dtext <- data.frame(
  x = c(1.3, 1.3),
  y = c(1.05, 0.75),
  name = c("Gulpers", "Mashers")
)

dtext2 <- data.frame(
  x = c(0.75, -0.75),
  y = c(0.25, 0.25),
  name = c("bill > fruit", "bill < fruit")
)

ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed") +
  geom_function(fun = function(x) plogis(-1.5 + 10 * x), xlim = c(-1.8, 2),
                color = viridis(1)) +
  geom_function(fun = function(x) plogis(2 + 1.5 * x) * 0.7, xlim = c(-1.8, 2),
                color = viridis(1, begin = 0.5)) +
  geom_text(data = dtext, mapping = aes(x, y, label = name, color = name),
            size = 3.5) +
  scale_color_viridis(discrete = TRUE, end = 0.5) +
  geom_text(data = dtext2, mapping = aes(x, y, label = name), alpha = 0.7,
            size = 3) +
  # scale_x_continuous(breaks = c(0)) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0.005, 0.005)) +
  theme_light() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  ylab("Consumption rate") +
  xlab("Size difference (bill - fruit)")

ggsave("figures/consumption_expected.tiff",
       width = 945, height = round(945 * 0.85), units = "px",
       dpi = 300)

# Width: 945 (single column), 1476 (1.5 column) or 1961 (double column) pixels
# (at 300 dpi). Resolution: 300-600 dpi.
# Size: <50 MB (for exceptions on file size, see below).