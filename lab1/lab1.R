#нормальное распределение

normal <- function(number){
  x<-rnorm(number, mean = 0, sd = 1)
  b <- 1 + 3.322*log(number)
  png(width =784, height = 456, filename = paste("Normal distribution, n =", number,".png"))
  hist(x, probability = TRUE, main = paste("Normal distribution, n=", number), breaks = b, xlab = "normalNumbers")
  xx<-seq(min(x)-1, max(x)+1, 0.0001)
  lines(xx, dnorm(xx), col = "red")
  dev.off()
}

# распределение Коши

koshi <- function(number){
  x<-rcauchy(number, location = 0, scale = 1)
  b <- 1 + 3.322*log10(number)
  png(width =784, height = 456, filename = paste("Koshi distribution, n =", number,".png"))
  hist(x, probability = TRUE, main = paste("Koshi distribution, n=", number), breaks = b, xlab = "cauchyNumbers", ylim=c(0,0.4))
  xx<-seq(min(x)-1, max(x)+1, 0.1)
  lines(xx, dcauchy(xx), col = "red")
  dev.off()
}

# распределение Лапласа

laplace <- function(number){
  x<-rlaplace(number, location = 0, scale = 1/sqrt(2))
  b <- 1 + 3.322*log10(number)
  png(width = 784, height = 456, filename = paste("Laplace distribution, n =", number,".png"))
  hist(x, probability = TRUE, main = paste("Laplace distribution, n=", number), breaks = b, xlab = "laplaceNumbers", ylim=c(0,1))
  xx<-seq(min(x)-1, max(x)+1, 0.001)
  lines(xx, dlaplace(xx, 0, 1/sqrt(2)), col = "red")
  dev.off()
}

# распределение Пуассона

poissons <- function(number){
  x<-rpois(n = number, lambda = 10)
  b <- 1 + 3.322*log10(number)
  png(width = 784, height = 456, filename = paste("Poisson distribution, n =", number,".png"))
  hist(x, probability = TRUE, main = paste("Poisson distribution, n=", number), breaks = b, xlab = "poisNumbers")
  xx<-seq(min(x)-1, max(x)+1, 1)
  lines(xx, dpois(xx, lambda = 10), col = "red")
  dev.off()
}

#равномерное распределение

uniform <- function(number){
  x<-runif(number, -sqrt(3), sqrt(3))
  b <- 1 + 3.322*log10(number)
  png(width = 784, height = 456, filename = paste("Uniform distribution, n =", number,".png"))
  hist(x, probability = TRUE, main = paste("Uniform distribution, n=", number), breaks = b, xlab = "unifNumbers")
  xx<-seq(min(x)-1, max(x)+1, 0.1)
  lines(xx, dunif(xx, -sqrt(3), sqrt(3)), col = "red")
  dev.off()
}

#задание размера выборки

n<-c(10, 50, 1000)

for(number in n)
  normal(number)
  koshi(number)
  laplace(number)
  poissons(number)
  uniform(number)
  
  