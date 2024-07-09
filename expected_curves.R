library(tidyverse)
library(viridis)
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
            size = 10 /.pt) +
  scale_color_viridis(discrete = TRUE, end = 0.5) +
  geom_text(data = dtext2, mapping = aes(x, y, label = name), alpha = 0.7,
            size = 9 /.pt) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0.005, 0.005)) +
  nice_theme() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ylab("Consumption rate") +
  xlab("Bill-fruit size difference")

ggsave("figures/expected_curves.tiff",
       width = 945, height = round(945 * 0.85), units = "px",
       dpi = 300)

# Width: 945 (single column), 1476 (1.5 column) or 1961 (double column) pixels
# (at 300 dpi). Resolution: 300-600 dpi.
# Size: <50 MB (for exceptions on file size, see below).