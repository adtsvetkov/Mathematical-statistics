library(EEM)

#regions are Africa and North
region <- c("Africa", "North")

countFiles = 15

#ways to data #пути до файлов с данными
strAfrica = "C:\\Hexane_extr_Kivu_Lake\\"
strNorth = "C:\\VD_DOM_Permafrost\\"

#names of files of Africa #имена файлов Африки
arrStrAfrica <- c("1.1_70.",
                  "1.2_21",
                  "1.3_68",
                  "1.4_114",
                  "1.5_11",
                  "1.6_37",
                  "2.3_5 (400)",
                  "2.3_5 (600)",
                  "2.3_5",
                  "2.4_7",
                  "3.1_14",
                  "3.2_69",
                  "3.3_15 (600)",
                  "3.4_20(800)",
                  #"4.4_87")
                  "3.4_20")

#names of files of North #имена файлов Севера
arrStrNorth <- c("1701",
                 "1702",
                 "1704",
                 "1706",
                 "1708_1to10",
                 "1708_1to20",
                 "1711",
                 "1712",
                 "1727",
                 "1728",
                 "1729",
                 "1730",
                 "1732",
                 "1733",
                 "1734")


#type of files #окончание пути до файла с данными, тип файла
strEnd = ".txt"

#amino acids
vectAA <- c()
vectAA[1:5] = 0

#names of Amino Acids
namesAA <- c("C","A","M","B","T")

#matrix of Amino Acids (AA) from different regions
aminoAcids <- matrix(data = vectAA, nrow = 2, ncol = 5)
colnames(aminoAcids) = namesAA
rownames(aminoAcids) = region

#diferent deltas from Africa's file and North's file
delta <- c(5, 2)

#start number of all files
startNum = 250

#coordinate of rectangles
crdRect <- matrix(data = c(320, 250, 310, 270, 270, 
                           350, 260, 320, 280, 280,
                           420, 380, 380, 300, 320,
                           480, 480, 420, 320, 350),
                  nrow = 5, ncol = 4)
colnames(crdRect) <- c("Left", "Right", "Down", "Up")
rownames(crdRect) <- namesAA
#View(coordsRect)

#initialization of matrixes

corr_names <- c("K x m1", "K x m2")

myK <- matrix(data = 0, nrow = countFiles, ncol = 2)
colnames(myK) <- region

my_m1 <- matrix(data = 0, nrow = countFiles, ncol = 2)
colnames(my_m1) <- region

my_m2 <- matrix(data = 0, nrow = countFiles, ncol = 2)
colnames(my_m2) <- region

spearman_matr <- matrix(data = 0, nrow = 2, ncol = 2)
colnames(spearman_matr) <- region
rownames(spearman_matr) <- corr_names


#function draw rectangles - the area of peaks of intensity
#функция рисует прямоугольники - области пиков интенсивности

drawRect <- function()
{
  rect(crdRect["C","Left"], crdRect["C","Down"], crdRect["C","Right"], crdRect["C","Up"], col = NA, border = "red") 
  rect(crdRect["A","Left"], crdRect["A","Down"], crdRect["A","Right"], crdRect["A","Up"], col = NA, border = "green") 
  rect(crdRect["M","Left"], crdRect["M","Down"], crdRect["M","Right"], crdRect["M","Up"], col = NA, border = "yellow") 
  rect(crdRect["B","Left"], crdRect["B","Down"], crdRect["B","Right"], crdRect["B","Up"], col = NA, border = "black") 
  rect(crdRect["T","Left"], crdRect["T","Down"], crdRect["T","Right"], crdRect["T","Up"], col = NA, border = "blue")
  legend("topleft", c("C", "A", "M","B", "T"), col = c("red","green","yellow","black", "blue"), pch = 19)
}

AfrAAA <- matrix(data = 0, nrow = countFiles, ncol = 5)
colnames(AfrAAA) <- namesAA

NorAAA <- matrix(data = 0, nrow = countFiles, ncol = 5)
colnames(NorAAA) <- namesAA

#find AA from file and add to table AfrAA or NorAA
findAA <- function(dataWithFluo, name = region[1], itr = 1, AfrAA, NorAA)
{
  a <- matrix(0)
  #считываем матрицу
  dataMatr <- dataWithFluo[[1]]
  #собираем аминокислоты из "нужных" мест матрицы.
  #"нужные" места соответствуют прямоугольникам.
  if (name == region[1]) 
  {
    for (i in 1:5)
    {
      a <- dataMatr[(c(crdRect[namesAA[i],"Down"] - startNum + 1)):(c(crdRect[namesAA[i],"Up"] - startNum + 1)),
                    (c(crdRect[namesAA[i],"Left"] - startNum)/delta[1] + 1):(c(crdRect[namesAA[i],"Right"] - startNum)/delta[1] + 1)]
      vectAA[i] = sum(a)
    }
    AfrAA[itr, ] = vectAA
    return(AfrAA)
  } else if (name == region[2]) {
    for (i in 1:5)
    {
      a <- dataMatr[(c(crdRect[namesAA[i],"Down"] - startNum + 1)):(c(crdRect[namesAA[i],"Up"] - startNum + 1)),
                    (c(crdRect[namesAA[i],"Left"] - startNum)/delta[2] + 1):(c(crdRect[namesAA[i],"Right"] - startNum)/delta[2] + 1)]
      vectAA[i] = sum(a)
    }
    NorAA[itr, ] = vectAA
    return(NorAA)
  } else {
    return(-1)
  }
}

#matrix with wavelengths #матрица длин волн
#Africa :: cutEX = 400:600, cutEM = 600:700
#North :: cutEX = 400:800, cutEM = 600:700

mWavelengts <- matrix(data = c(400,600,600,700,400,800,600,700), nrow = 2, ncol = 4, byrow = TRUE)
colnames(mWavelengts) = c("cutExLeft",
                          "cutExRight",
                          "cutEmLeft",
                          "cutEmRight" );
rownames(mWavelengts) = c("Africa", "North")
#View(mWavelengts)

findK <- function(vec)
{
  return((vec[1] + vec[2])/(vec[4] + vec[5]))
}

find_m1 <- function(vec)
{
  return(vec[3]/(vec[1]+vec[2]))
}

find_m2 <- function(vec)
{
  return(vec[3]/(vec[4]+vec[5]))
}

Spearman <- function(vec1, vec2)
{
  return(cor(vec1, vec2, method = "spearman"))
}

#выборочно строим информацию по файлу
example_foo <- function(AfrAA, NorAA, k)
{
  for (k in 1:countFiles)
  {
    #Africa
    strName = arrStrAfrica[k]
    wholeWay = paste0(strAfrica, strName, strEnd)
    
    #get data
    data <- readEEM(wholeWay)
    
    
    png(width =800, height = 500, filename = paste("Filename =", strName,".png"))
    drawEEM(data, n=1)
    dev.off()
    
    #обрезаем график
    dataCut <- cutEEM(data, 
                      cutEX = mWavelengts["Africa","cutExLeft"]:mWavelengts["Africa","cutExRight"],
                      cutEM = mWavelengts["Africa","cutEmLeft"]:mWavelengts["Africa","cutEmRight"])
    
    png(width =800, height = 500, filename = paste("dataCut =", strName,".png"))
    drawEEM(dataCut, n=1)
    dev.off()
    
    #удаляем лучи рэлеевского рассеяния
    dataWithFluo <- delScattering(dataCut, rep = 0)
    
    png(width =800, height = 500, filename = paste("dataWithFluo =", strName,".png"))
    drawEEM(dataWithFluo, n = 1)
    dev.off()
    
    #обозначаем границы аминокислот на рисунке
    
    png(width =800, height = 500, filename = paste("dataWithFluo =", strName,".png"))
    drawEEM(dataWithFluo, n = 1)
    drawRect()
    dev.off()
    
    AfrAA <- findAA(dataWithFluo, name = region[1], itr = k, AfrAA, NorAA)
    
    #North
    strName = arrStrNorth[k]
    wholeWay = paste0(strNorth, strName, strEnd)
    
    #get data
    data <- readEEM(wholeWay)
    
    png(width = 800, height = 500, filename = paste("Filename =", strName,".png"))
    drawEEM(data, n=1)
    dev.off()
    
    #обрезаем график
    dataCut <- cutEEM(data, 
                      cutEX = mWavelengts["North","cutExLeft"]:mWavelengts["North","cutExRight"],
                      cutEM = mWavelengts["North","cutEmLeft"]:mWavelengts["North","cutEmRight"])
    
    png(width =800, height = 500, filename = paste("dataCut =", strName,".png"))
    drawEEM(dataCut, n=1)
    dev.off()
    
    #удаляем лучи рэлеевского рассеяния
    dataWithFluo <- delScattering(dataCut, rep = 0) 
    
    png(width =800, height = 500, filename = paste("dataWithFluo =", strName,".png"))
    drawEEM(dataWithFluo, n = 1)
    dev.off()
    
    #обозначаем границы аминокислот на рисунке
    
    png(width =800, height = 500, filename = paste("dataWithFluo =", strName,".png"))
    drawEEM(dataWithFluo, n = 1)
    drawRect()
    dev.off()
    
    NorAA <- findAA(dataWithFluo, name = region[2], itr = k, AfrAA, NorAA)
    
    #ищем К, m1 и m2
    
    myK[k,1] = findK(AfrAA[k,])
    myK[k,2] = findK(NorAA[k,])
    
    my_m1[k,1] = find_m1(AfrAA[k,])
    my_m1[k,2] = find_m1(NorAA[k,])
    
    my_m2[k,1] = find_m2(AfrAA[k,])
    my_m2[k,2] = find_m2(NorAA[k,])
  }
  colnames(AfrAA) <- namesAA
  colnames(NorAA) <- namesAA
  
  View(AfrAA)
  View(NorAA)
  
  View(aminoAcids)
  
  #рисуем гистограмму для суммы каждой аминокислоты
  
  png(width =600, height = 500, filename ="Aminoacids.png")
  
  barplot(aminoAcids, beside = TRUE, col = c("orange","blue"), ylim = c(0, max(aminoAcids)*1.1), main = "Aminoacids summary", ylab = "intensity", xlab = "agent")
  legend("topright", c("Africa","North"), fill = c("orange","blue"))
  
  dev.off()
  
  View(myK)
  
  #рисуем гистограмму для параметра K
  
  png(width =600, height = 500, filename ="K-parameter 1.png")
  
  barplot(t(myK), beside = TRUE, col = c("brown","skyblue"), ylim = c(0, max(myK)*1.2), main = "K-parameter", ylab = "value", xlab = "file", names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  legend("topright", c("Africa","North"), fill = c("brown","skyblue"))
  
  dev.off()
  
  View(my_m1)
  
  #рисуем гистограмму для параметра m1
  
  png(width =600, height = 500, filename ="m1-parameter.png")
  
  barplot(t(my_m1), beside = TRUE, col = c("brown","skyblue"), ylim = c(0, max(my_m1)*1.2), main = "m1-parameter", ylab = "value", xlab = "file", names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  legend("topleft", c("Africa","North"), fill = c("brown","skyblue"))
  
  dev.off()
  
  View(my_m2)
  
  #рисуем гистограмму для параметра m2
  
  png(width =600, height = 500, filename ="m2-parameter.png")
  
  barplot(t(my_m2), beside = TRUE, col = c("brown","skyblue"), ylim = c(0, max(my_m2)*1.2), main = "m2-parameter", ylab = "value", xlab = "file", names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  legend("topright", c("Africa","North"), fill = c("brown","skyblue"))
  
  dev.off()

  #рисуем сводную гистограмму параметров для Африки
  
  Africa_par <- cbind(myK[, 1], my_m1[, 1], my_m2[, 1])
  
  png(width = 800, height = 500, filename ="Africa parameters.png")
  
  barplot(t(Africa_par), beside = TRUE, col = c("blue", "green", "yellow"), ylim = c(0, max(Africa_par)*1.2), main = "Africa parameters", ylab = "value", xlab = "file", names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  legend("topleft", c("K","m1", "m2"), fill = c("blue","green", "yellow"))
  
  dev.off() 
  
  #риусем сводную гистограмму параметров для севера
  
  North_par <- cbind(myK[, 2], my_m1[, 2], my_m2[, 2])
  
  png(width = 800, height = 500, filename ="North parameters.png")
  
  barplot(t(North_par), beside = TRUE, col = c("blue", "green", "yellow"), ylim = c(0, max(North_par)*1.1), main = "North parameters", ylab = "value", xlab = "file", names.arg = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  legend("topright", c("K","m1", "m2"), fill = c("blue","green", "yellow"))
  
  dev.off()
  
  #находим коэффициент корреляции Спирмена
  
  spearman_matr[1, 1] <- Spearman(myK[, 1], my_m1[, 1])
  spearman_matr[1, 2] <- Spearman(myK[, 2], my_m1[, 2])
  spearman_matr[2, 1] <- Spearman(myK[, 1], my_m2[, 1])
  spearman_matr[2, 2] <- Spearman(myK[, 2], my_m2[, 2])
  
  View(spearman_matr)
  
}

example_foo(AfrAAA, NorAAA, 1)
