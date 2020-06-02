gamma = 0.95 #параметр надежности
alpha = 1-gamma #уровень значимости

quantile1 = 1-alpha/2
quantile2 = alpha/2

#метод максимального правдоподобия

max_likelihood <- function(selection, size)
{
  mu <- mean(selection)
  sigma <- sqrt(sum((selection-mu)^2)/size)
  return(c(mu, sigma))
}

#классические интервальные оценки на основе статистик хи квадрат и
#Стьюдента

expected_value_classical <- function(selection, size)
{
  #квантиль распределения Стьюдента с size - 1 степенями свободы и порядка 1-alpha/2
  t <- qt(quantile1, size-1)
  ms <- max_likelihood(selection, size)
  value <-  ms[2]*t/sqrt(size-1)
  left <- ms[1] - value
  right <- ms[1] + value
  return(c(left, right))
}

deviation_classical <- function(seletion, size)
{
  #квантиль распределения хи-квадрат с size-1 степенями свободы и порядка 1-alpha/2
  t1 <- qchisq(quantile1, size-1)
  #квантиль распределения хи-квадрат с size-1 степенями свободы и порядка alpha/2
  t2 <- qchisq(quantile2, size-1)
  ms <- max_likelihood(selection, size)
  value <- ms[2]*sqrt(size)
  left <- value/sqrt(t1)
  right <- value/sqrt(t2)
  return(c(left, right))
}

#асимптотически нормальные интервальные оценки на основе точечных оценок метода максимального
#правдоподобия

expected_value_asymptotically <- function(selection, size)
{
  #квантиль нормального распределения N(0, 1) порядка 1-alpha/2
  u <- qnorm(quantile1, 0, 1)
  ms <- max_likelihood(selection, size)
  value <- ms[2]*u/sqrt(size)
  left <- ms[1]-value
  right <- ms[1]+value
  return (c(left, right))
}

deviation_asymptotically <- function(selection, size)
{
  m4 <- sum((selection - mean(selection))^4)/size
  ms <- max_likelihood(selection, size)
  e <- m4/ms[2]^4 - 3
  u <- qnorm(quantile1, 0, 1)
  value <- 0.5*u*sqrt((e+2)/size)
  left <- ms[2]*(1-0.5*value)
  right <- ms[2]*(1+0.5*value)
  return(c(left, right))
}

table1 <- data.frame()
table2 <- data.frame()

size1 = 20
size2 = 100

sizes <- c(size1, size2)

for (size in sizes)
{
  selection <- rnorm(size, 0, 1)
  
  helpmatrix1 <- data.frame()
  helpmatrix2 <- data.frame()
  
  mu1 <- expected_value_classical(selection, size)
  sigma1 <- deviation_classical(selection, size)
  helpmatrix1 <- rbind(helpmatrix1, c(paste(round(mu1[1], digits = 2), "< mu <", round(mu1[2], digits = 2)), paste(round(sigma1[1], digits = 2), "< sigma <", round(sigma1[2], digits = 2))), stringsAsFactors = FALSE)
  colnames(helpmatrix1)<-c("mu", "sigma")
  
  mu2 <- expected_value_asymptotically(selection, size)
  sigma2 <- deviation_asymptotically(selection, size)
  helpmatrix2 <- rbind(helpmatrix2, c(paste(round(mu2[1], digits = 2), "< mu <", round(mu2[2], digits = 2)), paste(round(sigma2[1], digits = 2), "< sigma <", round(sigma2[2], digits = 2))), stringsAsFactors = FALSE)
  colnames(helpmatrix2)<-c("mu", "sigma")
  
  new_matr <- data.frame(mu = paste("n =", size), sigma = " ", stringsAsFactors = FALSE)
  
  helpmatrix1 <- rbind(new_matr, helpmatrix1)
  helpmatrix2 <- rbind(new_matr, helpmatrix2)
  
  table1 <- rbind(table1, helpmatrix1)
  table2 <- rbind(table2, helpmatrix2)
}

write.csv(table1, "C:\\classical.csv", row.names = FALSE)
write.csv(table2, "C:\\asymptotically.csv", row.names = FALSE)