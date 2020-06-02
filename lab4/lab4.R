size1 <- 20
size2 <- 60
size3 <- 100

sizes <- c(size1, size2, size3)

a_1 = -4
b_1 = 4

a_2 = 6
b_2 = 14

#эмпирические функции распределения

#наша задача - выяснить, сколько z_i попадают в конкретный промежуток

empirical <- function(func1, func2, n, stringname, a, b, x)
{
  selection <- func1 #распределение
  sort(selection) 
  number = (1+3.322*log(n))*(n/4)
  st = (b-a)/number
  array = seq(a, b-st, st) #массив промежутков
  n_i <- vector()
  n_i[1] = 0
  for(i in 2:number){
    arr <- selection[selection < array[i]]
    n_i[i] <- length(arr) / n
  }
  y <- func2
  png(width = 456, height = 456, filename = paste(stringname, "distribution, n =", n,".png"))
  plot(array, n_i, type = "s", main = paste(stringname, "distribution, n =", n), ylim=c(0,1), col = "blue", xlab = "x", ylab = "F(x)")
  lines(x, y, col="red")
  dev.off()
}

#ядерные оценки плотности распределения

KDE <- function(func1, func2, n, stringname, a, b, x, limit)
{
  selection <- func1
  y <- func2
  png(width = 600, height = 600, filename = paste(stringname, "distribution KDE, n =", n,".png"))
  plot(density(selection, adjust=1/2), xlim = c(a, b), col="green", main = paste(stringname, "distribution, n =", n), xlab = "x", ylab = "f(x)", lwd = 2, ylim=c(0,limit))
  lines(x, y, col="black", lwd = 2)
  lines(density(selection, adjust=1), col="blue", lwd = 2)
  lines(density(selection, adjust=2), col = "red", lwd = 2)
  legend("topright", pch = 19, col= c("black", "green", "blue", "red"), legend = c("Prob destiny func", "0.5h", "h", "2h"))
  dev.off()
}

x <- runif(1000, a_1, b_1)
x <- sort(x)
x_1 <- runif(1000, a_2, b_2)
x_1 <- sort(x_1)
x_1 <- round(x_1, 0.1)

#получение ядерных оценок плотности распределения

for (size in sizes)
{
  KDE(rnorm(size, 0, 1), dnorm(x), size, "Normal", a_1, b_1, x, 0.8)
  KDE(rcauchy(size, 0, 1), dcauchy(x), size, "Cauchy", a_1, b_1, x, 0.6)
  KDE(rlaplace(size, 0, 1/sqrt(2)), dlaplace(x, 0, 1/sqrt(2)), size, "Laplace", a_1, b_1, x, 1)
  KDE(rpois(size, 10), dpois(x_1, 10), size, "Poisson", a_2, b_2, x_1, 0.2)
  KDE(runif(size, -sqrt(3), sqrt(3)), dunif(x, -sqrt(3), sqrt(3)), size, "Uniform", a_1, b_1, x, 0.6)
}

#получение эмпирических функций распределения

for (size in sizes)
{
  empirical(rnorm(size, 0, 1), pnorm(x, 0, 1), size, "Normal", a_1, b_1, x)
  empirical(rcauchy(size, 0, 1), pcauchy(x, 0, 1), size, "Cauchy", a_1, b_1, x)
  empirical(rlaplace(size, 0, 1/sqrt(2)), plaplace(x, 0, 1/sqrt(2)), size, "Laplace", a_1, b_1, x)
  empirical(rpois(size, 10), ppois(x_1, 10), size, "Poisson", a_2, b_2, x_1)
  empirical(runif(size, -sqrt(3), sqrt(3)), punif(x, -sqrt(3), sqrt(3)), size, "Uniform", a_1, b_1, x)
}