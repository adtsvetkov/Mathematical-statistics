size1 = 20
size2 = 60
size3 = 100

mu <- c(0, 0)

sizes = c(size1, size2, size3)

rho1 = 0
rho2 = 0.5
rho3 = 0.9

rho = c(rho1, rho2, rho3)

iterations = 1000

#E(z)

sample_mean <- function (sample, number) {
  return(sum(sample)/number)
}

#E(z^2)

sample_mean_squared <- function(sample, number)
{
  return(sum(sample^2)/number)
}

#D(z)

position_characteristics <- function(sample, number){
  new_sample <- vector()
  new_sample <- append(new_sample, sample[1:number]^2)
  return (sum(new_sample)/number - (sample_mean(sample, number)^2))
}

Pearson <- function(matrix)
{
  return(cor(matrix[, 1], matrix[, 2], method = "pearson"))
}

Spearman <- function(matrix)
{
  return(cor(matrix[, 1], matrix[, 2], method = "spearman"))
}

quadrant_cc <- function(matrix, size)
{
  n14 = matrix[matrix[,1] > median(matrix[,1]),] #правая часть от y
  n1 = n14[n14[,2] > median(matrix[,2]),] # первый квадрант
  n23 = matrix[matrix[,1] < median(matrix[,1]),] #левая часть от y
  n3 = n23[n23[,2] < median(matrix[,2]),] #третий квадрант
  count = length(n1[,1]) + length(n3[,1])
  return((2*count - size) / size)
}

#подсчет дисперсий коэффициентов 

get_tables <- function(name, mu, Sigma, Sigma2, rho, size, iterations)
{
  r <- vector() #pearsons
  r_S <- vector() #spearman
  r_Q <- vector() #quadrant
  for (i in 1:iterations)
  {
    if(strcmp(name, "Mixed") == TRUE)
      sample <- mvrnorm(size, mu = mu, Sigma = Sigma)
    if(strcmp(name, "Normal two-dimentional") == TRUE)
      sample <- (0.9 * mvrnorm(size, mu = mu, Sigma = Sigma) + 0.1 * mvrnorm(size, mu = mu, Sigma = Sigma2))
    r <- append(r, Pearson(sample))
    r_S <- append(r_S, Spearman(sample))
    r_Q <- append(r_Q, quadrant_cc(sample, size))
  }
  vec <- list(r, r_S, r_Q)
  matr <- data.frame()
  for (i in vec)
    matr <- rbind(matr, c(sample_mean(i, iterations), sample_mean_squared(i, iterations), position_characteristics(i, iterations)), stringsAsFactors = FALSE)
  rownames(matr)<-c("r", "r_S", "r_Q")
  colnames(matr)<- c("E(z)", "E(z^2)", "D(z)")
  new_matr <- data.frame(r = paste("rho = ", rho), r_S = paste("n = ", size), r_Q = "", stringsAsFactors = FALSE)
  rownames(new_matr)<- c("")
  new_matr <- rbind(new_matr, t(matr))
  write.csv(new_matr, paste("C:\\", name, "distribution, rho = ", rho, "n = ", size, ".csv"))
}

get_ellipses <- function(mu, Sigma, rho, size)
{
  sample <- mvrnorm(size, mu = mu, Sigma = Sigma)
  require(ellipse)
  confidence.ellipse <- ellipse(Sigma,centre=mu,level=0.99,npoints=100)
  png(width = 456, height = 456, filename = paste("Ellipse, size =", size, ", rho =", rho, ".png"))
  plot(confidence.ellipse, type="l", xlim=c(-3, 3), ylim=c(-3, 3), xlab = "x", ylab = "y", main = paste("size =", size, ", rho =", rho), lwd=2)
  par(new=TRUE)
  plot(sample, axes = FALSE, ann = FALSE, col="green", pch=19)
  dev.off()
}

for (size in sizes)
{
  for (p in rho)
  {
    Sigma1 = rbind(c(1, p), c(p, 1))
    get_ellipses(mu, Sigma1, p, size)
  }
}

p1 = -0.9

for (size in sizes)
{
  for (p in rho)
  {
    Sigma1 = rbind(c(1, p), c(p, 1))
    Sigma2 = rbind(c(100, 100*p), c(100*p, 100))
    get_tables("Mixed", mu, Sigma1, Sigma2, p, size, iterations)
  }
  #задание ковариационных матриц корреляций
  #на диагонали дисперсии
  #внедиагональные элементы \rho*\sqrt(disp)*\sqrt(disp)
  #https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  Sigma1 = rbind(c(1, rho3), c(rho3, 1))
  Sigma2 = rbind(c(100, 100*p1), c(100*p1, 100))
  get_tables("Normal two-dimentional", mu, Sigma1, Sigma2, "NONE", size, iterations)
}

