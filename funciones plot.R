curve(plogis(1.5 + 10 * x) * 0.91, from = -2, to = 2, ylim = c(0, 1),
      ylab = "Probabilidad de consumo",
      xlab = "Tamaño pico - tamaño fruto (cm)")
text(-1.5, 0.8, "Tragador")

curve(plogis(5 + 9 * x) * 0.83, add = TRUE, col = "red")
text(-1.5, 0.6, "Masticador", col = "red")

curve(ifelse(x > 0, 0.6 * (1 - exp(-6 * x)), 0), add = TRUE, col = "forestgreen")
text(-1.5, 0.4, "Uno muy quisquilloso", col = "forestgreen")

abline(v = 0, lty = 2, col = "gray10")
