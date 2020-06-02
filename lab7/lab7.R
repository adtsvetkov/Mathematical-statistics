#метод максимального правдоподобия

max_likelihood <- function(selection, size)
{
  mu <- mean(selection)
  sigma <- sqrt(sum((selection-mu)^2)/size)
  return(c(mu, sigma))
}

#наша гипотеза - что у нас нормальное распределение

#метод максимального правдоподобия

chi_squared <- function(selection, size)
{
  k <- ceiling(1 + 3.322*log10(size)) #берем промежутки по методу Старджесса
  a <- min(selection)
  b <- max(selection)
  h <- (b-a)/k
  x <- seq(a+h, b-h, h)
  x[1] = -Inf
  x[k-1] = Inf
  xx <- vector()
  for (i in 1:(k-2))
    xx <- append(xx, paste(x[i], ";", x[i+1]))
  xx <- append(xx, "---")
  n_i <- vector()
  p_i <- vector()
  np_i <- vector()
  n_i_np_i <- vector()
  answer <- vector()
  for (i in 1:(k-2))
  {
    #частота попадания элементов в i-ый кусочек
    n_i[i] <- length(selection[selection < x[i+1] & selection >= x[i]])
    #вычисляем вероятности с помощью гипотетической функции распределения
    p_i[i] <- pnorm(x[i + 1], 0, 1)- pnorm(x[i], 0, 1)
    #промежуточные вычисления
    np_i[i] <- size * p_i[i]
    #проверим, что число >=5, иначе объединяем промежутки
    print(np_i[i])
    n_i_np_i[i] <- n_i[i] - np_i[i]
    answer[i] <- (n_i_np_i[i]^2) / np_i[i]
  }
  n_i <- append(n_i, sum(n_i))
  p_i <- append(p_i, sum(p_i))
  np_i <- append(np_i, sum(np_i))
  n_i_np_i <- append(n_i_np_i, sum(n_i_np_i))
  answer <- append(answer, sum(answer))
  vec <- list(xx, n_i, p_i, np_i, n_i_np_i, answer)
  matr <- data.frame()
  for (i in vec)
    matr <- rbind(matr, i, stringsAsFactors = FALSE)
  rownames(matr)<-c("Borders", "n_i", "p_i", "np_i", "n_i_np_i", "answer")
  colnames(matr)<- c("1", "2", "3", "4", "5", "6", "sum")
  mxlklhd <- max_likelihood(selection , size)
  new_matr <- data.frame(Borders = paste("k = ", k-2), n_i = paste("size = ", size), p_i = "alpha = 0.05", np_i = "", n_i_np_i = paste("mu =", mxlklhd[1]), answer = paste("sigma =", mxlklhd[2]), stringsAsFactors = FALSE)
  rownames(new_matr)<- c("")
  new_matr <- rbind(new_matr, t(matr))
  print(new_matr)
  return(new_matr)
}

size = 100

matr <- chi_squared(rnorm(size, 0, 1), size)
write.csv(matr, "C:\\table.csv")

