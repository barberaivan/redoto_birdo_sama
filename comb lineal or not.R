# simulamos datos de peso de aves y de frutos
n <- 150
pajaro <- abs(rnorm(n))
fruto <- abs(rnorm(n))

# si simplemente restamos, estamos multiplicando cada variable por 1
delta_simple <- 1 * pajaro - 1 * fruto

# pero otra opción sería estimar coeficientes distintos de 1 para cada
# variable antes de restar
b_pajaro <- abs(rnorm(1))
b_fruto <- abs(rnorm(1))

delta_coef <- b_pajaro * pajaro - b_fruto * fruto

# si delta_coef es una combinación lineal de delta_simple, no tiene sentido
# hacerlo. En ese caso, la recta debería ajustar perfectamente a los puntos
# (sin variabilidad)
par(mfrow = c(1, 2))
plot(delta_coef ~ delta_simple, pch = 19, col = rgb(0, 0, 0, 0.5),
     main = "Con coeficientes")
abline(lm(delta_coef ~ delta_simple))

# si la transformación es una combinación lineal, la recta ajusta perfecto
# (prueba):
linear_comb <- rnorm(1) + rnorm(1) * delta_simple
plot(linear_comb ~ delta_simple, pch = 19, col = rgb(0, 0, 0, 0.5),
     main = "Combinación lineal")
abline(lm(linear_comb ~ delta_simple))
par(mfrow = c(1, 1))

# o sea, puede ser sensato estimar coeficientes >0 para cada peso antes de
# restar