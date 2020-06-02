a = -1.8
b = 2
h = 0.2

size = 20

e_i <-rnorm(size, 0, 1)
x <- seq(a, b, h)

disturb_1 = 10
disturb_20 = -10

#метод наименьших квадратов

least_square_method <- function(x, y)
{
  x_overlined = mean(x)
  y_overlined = mean(y)
  x_overlined_squared = mean(x^2)им
  y_overlined_squared = mean(y^2)
  xy = mean(x * y)
  b1 = (xy - x_overlined * y_overlined) / (x_overlined_squared - x_overlined^2) #формула Крамера
  b0 = y_overlined - x_overlined * b1 #из первого уравнения системы
  return(c(b0, b1))
}

#метод наименьших модулей

least_module_method <- function(x, y, size)
{
  k_q_n = 1.491
  r_Q = mean(sign(x-median(x))*sign(y-median(y)))
  l = size/4
  j = size - l + 1
  q_y = (y[j] - y[l])/k_q_n
  q_x = (x[j] - x[l])/k_q_n
  b1 = r_Q*q_y/q_x
  b0 = median(y) - b1*median(x)
  return(c(b0, b1))
}

paint <- function(x, y, b1, b2, name)
{
  png(width = 700, height = 400, filename = paste(name, "disturbance.png"))
  plot(x, y, col = "green", lwd=2)
  abline(a = 2, b = 2, col = "red")
  abline(a = b1[1], b = b1[2])
  abline(a = b2[1], b = b2[2], col = "blue")
  legend("topleft", c("Выборка", "Модель", "КНК","КНМ"), col = c("green","red","black","blue"), pch = 19)
  dev.off()
}

#выборка без возмущений

y1 = 2 + 2*x + e_i

#выборка с возмущением

y2 = y1
y2[1] = y2[1] + disturb_1
y2[20] = y2[20] + disturb_20

yy = cbind(y1, y2)

names = c("Without", "With")

for (i in 1:length(names))
{
  b1 = least_square_method(x, yy[,i])
  b2 = least_module_method(x, yy[,i], size)
  print(b1)
  print(b2)
  paint(x, yy[,i], b1, b2, names[i])
}