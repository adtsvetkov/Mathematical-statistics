# выборочное среднее (среднее характеристик положения)

sample_mean <- function (sample, number) {
  return(sum(sample)/number)
}

# выборочная медиана (учитываем, что number всегда четный)

selective_median <- function (sample, number)
{
  return((sample[number/2]+sample[number/2+1])/2)
}

#полусумма экстремальных выборочных элементов

half_sum_elem <- function(sample, number)
{
  return((max(sample)+min(sample))/2)
}

#выборочная квартиль порядка p

quartile <- function(sample, number, p)
{
  return(sample[ceiling(number*p)])
}

#полусумма квартилей

half_sum_quartiles <- function(sample, number)
{
  return((quartile(sample, number, 1/4) + quartile(sample, number, 3/4))/2)
}

#усеченное среднее

truncated_mean <- function(sample, number)
{
  r = ceiling(number/4)
  new_sample <- vector()
  new_sample <- append(new_sample, sample[(r+1):(number-r)])
  return((1/(number-2*r))*sum(new_sample))
}

#оценка дисперсии 

position_characteristics <- function(sample, number){
  new_sample <- vector()
  new_sample <- append(new_sample, sample[1:number]^2)
  return (sum(new_sample)/number - (sample_mean(sample, number)^2))
}

#непосредственно подсчет

calculate <- function(number, name)
{
  x_overlined <- vector()
  med_x <- vector()
  z_R <- vector()
  z_Q <- vector()
  z_tr <- vector()
  for (i in 1:1000)
  {
    if(strcmp(name, "Normal") == TRUE)
      sample <- rnorm(number, mean = 0, sd = 1)
    if(strcmp(name, "Cauchy") == TRUE)
      sample <- rcauchy(number, location = 0, scale = 1)
    if(strcmp(name, "Laplace") == TRUE)
      sample <- rlaplace(number, location = 0, scale = 1/sqrt(2))
    if(strcmp(name, "Poisson") == TRUE)
      sample <- rpois(n = number, lambda = 10)
    if(strcmp(name, "Uniform") == TRUE)
      sample <- runif(number, -sqrt(3), sqrt(3))
    sample <- sort(sample)
    x_overlined <- append(x_overlined, sample_mean(sample, number))
    med_x <- append(med_x, selective_median(sample, number))
    z_R <- append(z_R, half_sum_elem(sample, number))
    z_Q <- append(z_Q, half_sum_quartiles(sample, number))
    z_tr <- append(z_tr, truncated_mean(sample, number))
  }
  vec <- list(x_overlined, med_x, z_R, z_Q, z_tr)
  matr <- data.frame()
  for (i in vec)
    matr <- rbind(matr, c(sample_mean(i, 1000), position_characteristics(i, 1000)), stringsAsFactors = FALSE)
  rownames(matr)<-c("x_overlined", "med_x", "z_R", "z_Q", "z_tr")
  colnames(matr)<- c("E(x)", "D(x)")
  new_matr <- data.frame(x_overlined = paste(name, " distribution"), med_x = paste("n=", number), z_R = "", z_Q = "", z_tr = "", stringsAsFactors = FALSE)
  rownames(new_matr)<- c("")
  new_matr <- rbind(new_matr, t(matr))
  return(new_matr)
}

names_arr <- c("Normal", "Cauchy", "Laplace", "Poisson", "Uniform")
numbers <- c(10, 50, 1000)

final_matrix <- data.frame()

for (name in names_arr)
{
  for (num in numbers)
    final_matrix <- rbind(final_matrix, calculate(num, name))
}

write.csv(final_matrix, "C:\\table.csv")
