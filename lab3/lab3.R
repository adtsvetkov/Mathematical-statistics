size1 <- 20
size2 <- 100

#создание боксплота Тьюки

boxplots <- function(x, n1, y, n2, name)
{
  x<-append(x, y)
  arr <- rep(n1, times = n1)
  arr1 <- rep(n2, times = n2)
  arr <- append(arr, arr1)
  png(width = 534, height = 404, filename = paste(name, "distribution.png"))
  boxplot(x ~ arr, horizontal = TRUE, main = paste(name, "distribution"), xlab = "", ylab = "Sample size", col = "coral")
  dev.off()
}

#Нормальное распределение

x<-rnorm(size1, mean = 0, sd = 1)
y<-rnorm(size2, mean =0, sd = 1)
boxplots(x, size1, y, size2, "Normal")

#Распределение Коши
x<-rcauchy(size1, location = 0, scale = 1)
y<-rcauchy(size2, location = 0, scale = 1)
boxplots(x, size1, y, size2, "Cauchy")

#распределение Лапласа
x<-rlaplace(size1, location = 0, scale = 1/sqrt(2))
y<-rlaplace(size2, location = 0, scale = 1/sqrt(2))
boxplots(x, size1, y, size2, "Laplace")

#распределение Пуассона
x<-rpois(size1, lambda = 10)
y<-rpois(size2, lambda = 10)
boxplots(x, size1, y, size2, "Poisson")

#равномерное распределение
x<-runif(size1, -sqrt(3), sqrt(3))
y<-runif(size2, -sqrt(3), sqrt(3))
boxplots(x, size1, y, size2, "Uniform")

#поиск количества выбросов

find_emisson <- function(randX){
  boxplot(randX)
  numInd <- which(randX %in% boxplot.stats(randX)$out)
  return(length(numInd))
}

#подсчет средней доли выбросов

calculate <- function(size)
{
  Normal <- vector()
  Cauchy <- vector()
  Laplace <- vector()
  Pois <- vector()
  Uniform <- vector()
  
  for (i in 1:1000)
  {
    x1 <- rnorm(size, mean = 0, sd = 1)
    Normal <- append(Normal, find_emisson(x1)/size)
    x2 <- rcauchy(size, location = 0, scale = 1)
    Cauchy <- append(Cauchy, find_emisson(x2)/size)
    x3 <- rlaplace(size, location = 0, scale = 1/sqrt(2))
    Laplace <- append(Laplace, find_emisson(x3)/size)
    x4 <- rpois(size, lambda = 10)
    Pois <- append(Pois, find_emisson(x4)/size)
    x5 <- runif(size, -sqrt(3), sqrt(3))
    Uniform <- append(Uniform, find_emisson(x5)/size)
  }
  
  print(paste("Normal distribution: size =", size, "result =", sum(Normal)/1000))
  print(paste("Cauchy distribution: size =", size, "result =", sum(Cauchy)/1000))
  print(paste("Laplace distribution: size =", size, "result =", sum(Laplace)/1000))
  print(paste("Poisson distribution: size =", size, "result =", sum(Pois)/1000))
  print(paste("Uniform distribution: size =", size, "result =", sum(Uniform)/1000))
}

mas <- c(size1, size2)
for (m in mas)
  calculate(m)
  
#теоретические результаты

#расчет границ уса

X1andX2 <- function(q)
{
  return(c(q[1]-(3/2)*(q[2]-q[1]), q[2] + (3/2)*(q[2]-q[1])))
}

theoretical_calculation <- function()
{
  data <- data.frame()
  
  #нормальное распределение
  
  x1 <- qnorm(c(0.25, 0.75), 0, 1)
  x1_whisker <- X1andX2(x1)
  x1_p = pnorm(x1_whisker[1], 0, 1) + 1 - pnorm(x1_whisker[2], 0, 1)
  data <- rbind(data, c(x1[1], x1[2], x1_whisker[1], x1_whisker[2], x1_p), stringsAsFactors = FALSE)
  
  #распределение Коши
  
  x2 <- qcauchy(c(0.25, 0.75), 0, 1)
  x2_whisker <- X1andX2(x2)
  x2_p = pcauchy(x2_whisker[1], 0, 1) + 1 - pcauchy(x2_whisker[2], 0, 1)
  data <- rbind(data, c(x2[1], x2[2], x2_whisker[1], x2_whisker[2], x2_p), stringsAsFactors = FALSE)
  
  #распределение Лапласа
  
  x3 <- qlaplace(c(0.25, 0.75), 1/sqrt(2))
  x3_whisker <- X1andX2(x3)
  x3_p = plaplace(x3_whisker[1], 1/sqrt(2)) + (1 - plaplace(x3_whisker[2], 1/sqrt(2)))
  data <- rbind(data, c(x3[1], x3[2], x3_whisker[1], x3_whisker[2], x3_p), stringsAsFactors = FALSE)
  
  #распределение Пуассона
  
  x4 <- qpois(c(0.25, 0.75), 10)
  x4_whisker <- X1andX2(x4)
  x4_p = ppois(x4_whisker[1], 10) + 1 - ppois(x4_whisker[2], 10)
  data <- rbind(data, c(x4[1], x4[2], x4_whisker[1], x4_whisker[2], x4_p), stringsAsFactors = FALSE)
  
  #равномерное распределение
  
  x5 <- qunif(c(0.25, 0.75), -sqrt(3), sqrt(3))
  x5_whisker <- X1andX2(x5)
  x5_p = punif(x5_whisker[1], -sqrt(3), sqrt(3)) + 1 - punif(x5_whisker[2], -sqrt(3), sqrt(3))
  data <- rbind(data, c(x5[1], x5[2], x5_whisker[1], x5_whisker[2], x5_p), stringsAsFactors = FALSE)
  
  colnames(data)<-c("Q1","Q3", "X1", "X2", "P_в")
  rownames(data)<- c("Normal distribution", "Cauchy distribution", "Laplace distribution", "Poisson distribution", "Uniform distribution")
  
  return(data)
}

data <- theoretical_calculation()
write.csv(data, "C:\\table.csv")
