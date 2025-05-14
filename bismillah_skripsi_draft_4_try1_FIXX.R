####### Library #######
library(FactoMineR)
library(stats)
library(datasets)
library(magrittr)
library(dplyr)
library(tidyverse)
library(writexl)
library(ca)
library("factoextra")
library(ggplot2)

####### Input Data ####### 
data <- read.csv(file.choose(),header=TRUE) 
head(data)
####### Tabel Kontingensi ####### 
x1=data[,c(1,2)] #Variabel Bahan bangunan atap rumah terluas
x2=data[,c(1,3)] #Variabel Bahan bangunan dinding rumah terluas
x3=data[,c(1,4)] #Variabel Bahan bangunan lantai rumah terluas
x4=data[,c(1,5)] #Variabel Sumber air utama yang digunakan untuk minum
x5=data[,c(1,6)] #Variabel Sumber air utama yang digunakan rumah tangga untuk mandi/cuci/dll
x6=data[,c(1,7)] #Variabel Kepemilikan fasilitas tempat buang air besar dan siapa saja yang menggunakan
x7=data[,c(1,8)] #Variabel Jenis kloset yang digunakan
x8=data[,c(1,9)] #Variabel Tempat pembuangan akhir tinja
x9=data[,c(1,10)] #Variabel Intensitas tangki septik dikosongkan/dilakukan penyedotan dalam 5 tahun terakhir
x10=data[,c(1,11)] #Variabel Kecukupan luas tempat tinggal per kapita

####### Uji Chi-square ####### 
contingency1 = caconv(x1, from="rpm", to="freq")
cont1 = as.table(contingency1)
names(dimnames(cont1))<-c("K","X1")
cont1
chisq.test(cont1)

contingency2 = caconv(x2, from="rpm", to="freq")
cont2 = as.table(contingency2)
names(dimnames(cont2))<-c("K","X2")
cont2
chisq.test(cont2)

contingency3 = caconv(x3, from="rpm", to="freq")
cont3 = as.table(contingency3)
names(dimnames(cont3))<-c("K","X3")
cont3
chisq.test(cont3)

contingency4 = caconv(x4, from="rpm", to="freq")
cont4 = as.table(contingency4)
names(dimnames(cont4))<-c("K","X4")
cont4
chisq.test(cont4)

contingency5 = caconv(x5, from="rpm", to="freq")
cont5 = as.table(contingency5)
names(dimnames(cont5))<-c("K","X5")
cont5
chisq.test(cont5)

contingency6 = caconv(x6, from="rpm", to="freq")
cont6 = as.table(contingency6)
names(dimnames(cont6))<-c("K","X6")
cont6
chisq.test(cont6)

contingency7 = caconv(x7, from="rpm", to="freq")
cont7 = as.table(contingency7)
names(dimnames(cont7))<-c("K","X7")
cont7
chisq.test(cont7)

contingency8 = caconv(x8, from="rpm", to="freq")
cont8 = as.table(contingency8)
names(dimnames(cont8))<-c("K","X8")
cont8
chisq.test(cont8)

contingency9 = caconv(x9, from="rpm", to="freq")
cont9 = as.table(contingency9)
names(dimnames(cont9))<-c("K","X9")
cont9
chisq.test(cont9)

contingency10 = caconv(x10, from="rpm", to="freq")
cont10 = as.table(contingency10)
names(dimnames(cont10))<-c("K","X10")
cont10
chisq.test(cont10)
# H0 ditolak, seluruh variabel digunakan dalam analisis

# save table
#s# write.table(cont10, file="1. cont10.csv",row.names=T,sep=";")

####### Analisis Korespondensi Multipel ####### 
# Menggunakan Package mjca
mca1 <- mjca(data, nd=NA, lambda = "Burt", reti=TRUE)
summary(mca1) 
plot(mca1)
mca1$sv^2
#--# 2 dim 21.0%
#--# 26 dim 70.1%
#--# 60 dim 100%

## Diperoleh Hasil Sebagai berikut:
## Total Dimensi = 60
## Pada dimensi 2 diperoleh hasil cumulative percentage of variance adalah sebesar 21.0%
## Nilai cumulative percentage of variance ini dinilai masih kurang sehingga perlu ditingkatkan lagi
## Salah satu cara yang dapat dilakukan untuk meningkatkan cumulative percentage of variance adalah dengan melakukan rekatogorisasi dengan menggunakan Elliptical Confidence Regions
## dim26 cum%70.1

####### Confidence Regions ####### 
regions <- function(coord, xcoord, ycoord, col = col){
  t <- seq(0, 2*pi, length = 1000)
  pcoord1 <- coord[1] + xcoord*cos(t)
  pcoord2 <- coord[2] + ycoord*sin(t)
  lines(pcoord1, pcoord2, col = col)
}

ca.regions.exe <- function (N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), M =
                              min(I,J) - 1, region = 2, scaleplot = 1.2) {
  #############################################################
  # #
  # Defining features of the contingency table for CA #
  # #
  #############################################################
  I <- nrow(N) # Number of rows of table
  J <- ncol(N) # Number of columns of table
  Inames <- dimnames(N)[[1]] # Row category names
  Jnames <- dimnames(N)[[2]] # Column category names
  n <- sum(N) # Total number classified in the table
  p <- N *(1/n) # Matrix of joint relative proportions
  Imass <- as.matrix(apply(p, 1, sum))
  Jmass <- as.matrix(apply(p, 2, sum))
  ItJ <- Imass %*% t(Jmass)
  y <- p - ItJ
  dI <- diag(Imass[1:I])
  dJ <- diag(Jmass[1:J])
  Ih <- Imass^-0.5
  Jh <- Jmass^-0.5
  dIh <- diag(Ih[1:I])
  dJh <- diag(Jh[1:J])
  x <- dIh%*%y%*%dJh
  sva <- svd(x)
  a <- dIh%*%sva$u
  b <- dJh%*%sva$v
  dmu <- diag(sva$d) # Diagonal matrix of singular values
  f <- a %*% dmu # Row coordinates for Classical CA
  g <- b %*% dmu # Column coordinates for Classical CA
  dimnames(f)[[1]] <- Inames
  dimnames(g)[[1]] <- Jnames
  Principal.Inertia <- diag(t(f[, 1:min(I-1, J-1)])%*%dI%*%
                              f[, 1:min(I-1,J-1)])
  Total.Inertia <- sum(Principal.Inertia)
  Percentage.Inertia <- (Principal.Inertia/Total.Inertia) * 100
  Total.Perc.Inertia.M <- sum(Principal.Inertia[1:M])
  chisq.val <- qchisq(1-alpha, df = (I - 1) * (J - 1) )
  #############################################################
  # #
  # Construction of correspondence plot #
  # #
  #############################################################
  par(pty = "s")
  plot(0, 0, pch = " ", xlim = scaleplot * range(f[, 1:M], g[, 1:M]),
       ylim = scaleplot * range(f[, 1:M], g[, 1:M]),
       xlab = paste("Principal Axis ", a1, "(", round(Percentage.Inertia[a1],
                                                      digits = 2), "%)"), ylab = paste("Principal Axis ", a2, "(",
                                                                                       round(Percentage.Inertia[a2], digits = 2), "%)"))
  text(f[,1], f[,2], labels = Inames, adj = 0, col = cols[1])
  points(f[, a1], f[, a2], pch = "*", col = cols[1])
  text(g[,1], g[,2], labels = Jnames, adj = 1, col = cols[2])
  points(g[, a1], g[, a2], pch = "#", col = cols[2])
  abline(h = 0, v = 0)
  title(main = paste(100 * (1 - alpha), "% Confidence Regions"))
  #############################################################
  # #
  # Calculating the row and column radii length for a #
  # confidence circle #
  # #
  #############################################################
  radii <- sqrt(qchisq(1 - alpha, 2)/(n * Imass))
  radij <- sqrt(qchisq(1 - alpha, 2)/(n * Jmass))
  #############################################################
  # #
  # Calculating the semi-axis lengths for the confidence #
  # ellipses #
  # #
  #############################################################
  hlax1.row <- vector(mode = "numeric", length = I)
  hlax2.row <- vector(mode = "numeric", length = I)
  hlax1.col <- vector(mode = "numeric", length = J)
  hlax2.col <- vector(mode = "numeric", length = J)
  if (M > 2){
    # Semi-axis lengths for the row coordinates in an optimal plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n*Total.Inertia))*
                                        (1/Imass[i] - sum(a[i, 3:M]^2)))
      hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n*Total.Inertia))*
                                        (1/Imass[i] - sum(a[i, 3:M]^2)))
    }
    # Semi-axis lengths for the column coordinates in an optimal plot
    for (j in 1:J){
      hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j] - sum(b[j, 3:M]^2)))
      hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j] - sum(b[j, 3:M]^2)))
    }
  } else {
    # Semi-axis lengths for the row coordinates in a two-dimensional plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Imass[i]))
      hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Imass[i]))
    }
    # Semi-axis lengths for the column coordinates in a two-dimensional plot
    for (j in 1:J){
      hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j]))
      hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j]))
    }
  }
  #############################################################
  # #
  # Eccentricity #
  # #
  #############################################################
  eccentricity <- sqrt(1-(dmu[a2,a2]/dmu[a1,a1])^2)
  #############################################################
  # #
  # Approximate P-values #
  # #
  #############################################################
  pvalrow <- vector(mode = "numeric", length = I)
  pvalrowcircle <- vector(mode = "numeric", length = I)
  pvalcol <- vector(mode = "numeric", length = J)
  pvalcolcircle <- vector(mode = "numeric", length = J)
  for (i in 1:I){
    # Approximate row P-values from Lebart et al.’s (1984) confidence
    # circles
    pvalrowcircle[i] <- 1 - pchisq(n * Imass[i] * (f[i, 1]^2 + f[i, 2]^2),
                                   df = (I-1)*(J-1))
    # Approximate P-values based on Beh’s (2010) confidence ellipses
    if (M > 2){
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * ((1/Imass[i] -
                                                       sum(a[i, 3:M]^2))^(-1)) * ((f[i, 1]/dmu[1, 1])^2
                                                                                  + (f[i, 2]/dmu[2, 2])^2), df = (I - 1) * (J - 1))
    } else {
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * Imass[i]*((f[i, 1]/
                                                                dmu[1, 1])^2 + (f[i, 2]/dmu[2, 2])^2),
                               df = (I - 1) * (J - 1))
    }
  }
  for (j in 1:J){
    # Approximate row P-values based on Lebart et al.’s (1984)
    # confidence circles
    pvalcolcircle[j] <- 1 - pchisq(n * Imass[i] * (g[j, 1]^2 + g[j, 2]^2),
                                   df = (I - 1) * (J - 1))
    # Approximate P-values based on Beh’s (2010) confidence ellipses
    if (M > 2){
      pvalcol[j] <- 1 - pchisq(n * Total.Inertia * ((1/Jmass[j] -
                                                       sum(b[j, 3:M]^2))^(-1)) * ((g[j, 1]/dmu[1, 1])^2
                                                                                  + (g[j, 2]/dmu[2, 2])^2), df = (I - 1)*(J - 1))
    } else {
      pvalcol[j] <- 1 - pchisq(n * Total.Inertia * Jmass[j]*((g[j,1]/
                                                                dmu[1, 1])^2 + (g[j, 2]/dmu[2, 2])^2),
                               df = (I - 1) * (J - 1))
    }
  }
  summ.name <- c("HL Axis 1", "HL Axis 2", "P-value-ellipse",
                 "P-value-circle")
  if (region == 1){
    row.summ <- cbind(radii, radii, pvalrow, pvalrowcircle)
    col.summ <- cbind(radij, radij, pvalcol, pvalcolcircle)
  } else if (region == 2){
    row.summ <- cbind(hlax1.row, hlax2.row, pvalrow, pvalrowcircle)
    col.summ <- cbind(hlax1.col, hlax2.col, pvalcol, pvalcolcircle)
  }
  dimnames(row.summ) <- list(paste(Inames), paste(summ.name))
  dimnames(col.summ) <- list(paste(Jnames), paste(summ.name))
  #############################################################
  # #
  # Superimposing the confidence regions #
  # #
  #############################################################
  if (region == 1){
    # Superimposing the confidence circles
    symbols(f[,a1], f[,a2], circles = radii, add = T, fg = cols[1])
    symbols(g[,a1], g[,a2], circles = radij, add = T, fg = cols[2])
  } else if (region == 2){
    # Superimposing the confidence ellipses
    for (i in 1:I){
      regions(f[i,], xcoord = hlax1.row[i], ycoord = hlax2.row[i],
              col = cols[1])
    }
    for (j in 1:J){
      regions(g[j,], xcoord = hlax1.col[j], ycoord = hlax2.col[j],
              col = cols[2])
    }
  }
  #############################################################
  # #
  # Summary of output #
  # #
  #############################################################
  if (region == 1){
    list(Row.Summary = round(row.summ, digits = 3), Column.Summary =
           round(col.summ, digits = 3),
         Inertia = Principal.Inertia)
  } else if (region == 2){
    list(Eccentricity = round(eccentricity, digits = 3), Row.Summary
         = round(row.summ, digits = 3), Column.Summary =
           round(col.summ, digits = 3), Percentage.Inertia=Percentage.Inertia, Inertia = Principal.Inertia, 
         P=round(p, digits = 4),
         Imass_rtild=round(Imass, digits=4),dI_Drtild =dI,
         cm=round(Jmass, digits=4),  dJ_Dctild=dJ,
         S= round(x, digits=4),sva =sva,
         f=round(f,digits=4),g=round(g,digits=4))
  }
}

emerson.poly <- function (mj, pj) {
  # Menghitung jumlah kolom
  nc <- length(mj)
  
  # Membuat matriks diagonal dari pj
  Dj <- diag(pj)
  
  # Inisialisasi matriks B dengan ukuran (nc+1) x nc, elemen pertama diatur ke 0 pada baris pertama
  B <- matrix(1, (nc + 1), nc)
  B[1, ] <- 0
  
  # Inisialisasi variabel Sh, Th, dan Vh sebagai NULL
  Sh <- Th <- Vh <- NULL
  
  # Loop untuk mengisi nilai pada matriks B mulai dari indeks ke-3 hingga nc+1
  for (i in 3:(nc + 1)) {
    for (j in 1:nc) {
      # Menghitung Th sebagai hasil kali mj, Dj, dan elemen matriks B kuadrat
      Th[i] <- mj %*% Dj %*% B[i - 1, ]^2
      
      # Menghitung Vh sebagai hasil kali mj, Dj, dan produk elemen matriks B[i-1,] dan B[i-2,]
      Vh[i] <- mj %*% Dj %*% (B[i - 1, ] * B[i - 2, ])
      
      # Menghitung Sh sebagai kebalikan dari akar kuadrat dari operasi tertentu
      Sh[i] <- sqrt(mj^2 %*% Dj %*% B[i - 1, ]^2 - Th[i]^2 - Vh[i]^2)^(-1)
      
      # Menghitung elemen B[i, j] berdasarkan rumus yang diberikan
      B[i, j] <- Sh[i] * ((mj[j] - Th[i]) * B[i - 1, j] - Vh[i] * B[i - 2, j])
    }
  }
  # Transpose matriks B untuk mengatur orientasi baris-kolom
  B <- t(B)
  
  # Buat sub-matriks B1 tanpa kolom pertama dan kedua, dan BT tanpa kolom pertama
  B1 <- B[, -c(1, 2)]
  BT <- B[, -c(1)]
  
  # Mengembalikan list dari matriks B1 dan BT
  list(B = B1, BT = BT)
}

soca.regions.exe <- function (N, a1 = 1, a2 = 2, alpha = 0.05, cols = c(2, 4), M =
                                min(I,J) - 1, scaleplot = 1.2) {
  #############################################################
  # #
  # Defining features of the contingency table for CA #
  # #
  #############################################################
  I <- nrow(N) # Number of rows of table
  J <- ncol(N) # Number of columns of table
  Inames <- dimnames(N)[[1]] # Row category names
  Jnames <- dimnames(N)[[2]] # Column category names
  n <- sum(N) # Total number classified in the table
  
  # Calculate row and column mass matrices
  p <- N *(1/n) # Matrix of joint relative proportions
  Imass <- as.matrix(apply(p, 1, sum))
  Jmass <- as.matrix(apply(p, 2, sum))
  
  # Perform SOCA
  ItJ <- Imass %*% t(Jmass)
  y <- p - ItJ
  dI <- diag(Imass[1:I])
  dJ <- diag(Jmass[1:J])
  Ih <- Imass^-0.5
  Jh <- Jmass^-0.5
  dIh <- diag(Ih[1:I])
  dJh <- diag(Jh[1:J])
  x <- dIh%*%y%*%dJh
  sva <- svd(x)
  
  # Define mj as all column indices
  mj <- 1:J
  Jmas <- as.vector(apply(p, 2, sum))
  b.BMD <- emerson.poly(mj,Jmas)$B
  b.HD <- sqrt(dJ) %*% b.BMD
  z <- t(sva$u) %*% x %*% b.HD
  
  a <- dIh%*%sva$u
  f <- a %*% z # Row coordinates for SOCA 
  g <- dJh %*% b.HD # Column coordinates for SOCA
  dimnames(f)[[1]] <- Inames
  dimnames(g)[[1]] <- Jnames
  
  Principal.Inertia <- diag(t(f[, 1:min(I-1, J-1)])%*%dI%*%
                              f[, 1:min(I-1,J-1)])
  Total.Inertia <- sum(Principal.Inertia)
  Percentage.Inertia <- (Principal.Inertia/Total.Inertia) * 100
  Total.Perc.Inertia.M <- sum(Principal.Inertia[1:M])
  
  # Calculate half-length axes and areas
  dmu <- sqrt(Principal.Inertia)
  chisq.val <- qchisq(1-alpha, df = (I - 1) * (J - 1))
  #############################################################
  # #
  # Construction of correspondence plot #
  # #
  #############################################################
  par(pty = "s")
  plot(0, 0, pch = " ", xlim = scaleplot * range(f[, 1:M], g[, 1:M]),
       ylim = scaleplot * range(f[, 1:M], g[, 1:M]),
       xlab = paste("Principal Axis ", a1, "(", round(Percentage.Inertia[a1],
                                                      digits = 2), "%)"), ylab = paste("Principal Axis ", a2, "(",
                                                                                       round(Percentage.Inertia[a2], digits = 2), "%)"))
  text(f[,1], f[,2], labels = Inames, adj = 0, col = cols[1])
  points(f[, a1], f[, a2], pch = "*", col = cols[1])
  text(g[,1], g[,2], labels = Jnames, adj = 1, col = cols[2])
  points(g[, a1], g[, a2], pch = "#", col = cols[2])
  abline(h = 0, v = 0)
  title(main = paste(100 * (1 - alpha), "% Confidence Regions"))
  #############################################################
  # #
  # Calculating the row and column radii length for a #
  # confidence circle #
  # #
  #############################################################
  radij <- sqrt(qchisq(1 - alpha, 2)/(n * Jmass))
  #############################################################
  # #
  # Calculating the semi-axis lengths for the confidence #
  # ellipses #
  # #
  #############################################################
  hlax1.row <- vector(mode = "numeric", length = I)
  hlax2.row <- vector(mode = "numeric", length = I)
  if (M > 2){
    # Semi-axis lengths for the row coordinates in an optimal plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1] * sqrt((chisq.val/(n*Total.Inertia))*
                                      (1/Imass[i] - sum(a[i, 3:M]^2)))
      hlax2.row[i] <- dmu[2] * sqrt((chisq.val/(n*Total.Inertia))*
                                      (1/Imass[i] - sum(a[i, 3:M]^2)))
    }
  } else {
    # Semi-axis lengths for the row coordinates in a two-dimensional plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1] * sqrt((chisq.val/(n * Total.Inertia))*
                                      (1/Imass[i]))
      hlax2.row[i] <- dmu[2] * sqrt((chisq.val/(n * Total.Inertia))*
                                      (1/Imass[i]))
    }
  }
  #############################################################
  # #
  # Eccentricity #
  # #
  #############################################################
  eccentricity <- sqrt(1-(dmu[a2]/dmu[a1])^2)
  #############################################################
  # #
  # Approximate P-values #
  # #
  #############################################################
  pvalrow <- vector(mode = "numeric", length = I)
  pvalcolcircle <- vector(mode = "numeric", length = J)
  for (i in 1:I){
    # Approximate P-values based on Beh’s (2010) confidence ellipses
    if (M > 2){
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * ((1/Imass[i] -
                                                       sum(a[i, 3:M]^2))^(-1)) * ((f[i, 1]/dmu[1])^2
                                                                                  + (f[i, 2]/dmu[2])^2), df = (I - 1) * (J - 1))
    } else {
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * Imass[i]*((f[i, 1]/
                                                                dmu[1])^2 + (f[i, 2]/dmu[2])^2),
                               df = (I - 1) * (J - 1))
    }
  }
  for (j in 1:J){
    # Approximate row P-values based on Lebart et al.’s (1984)
    # confidence circles
    pvalcolcircle[j] <- 1 - pchisq(n * Imass[i] * (g[j, 1]^2 + g[j, 2]^2),
                                   df = (I - 1) * (J - 1))
  }
  summ.name <- c("HL Axis 1", "HL Axis 2", "P-value")
  row.summ <- cbind(hlax1.row, hlax2.row, pvalrow)
  col.summ <- cbind(radij, radij, pvalcolcircle)
  
  dimnames(row.summ) <- list(paste(Inames), paste(summ.name))
  dimnames(col.summ) <- list(paste(Jnames), paste(summ.name))
  #############################################################
  # #
  # Superimposing the confidence regions #
  # #
  #############################################################
  symbols(g[,a1], g[,a2], circles = radij, add = T, fg = cols[2])
  # Superimposing the confidence ellipses
  for (i in 1:I){
    regions(f[i,], xcoord = hlax1.row[i], ycoord = hlax2.row[i],col = cols[1])
  }
  
  #############################################################
  # #
  # Summary of output #
  # #
  #############################################################
  list(Row.Summary = round(row.summ, digits = 3), Column.Summary =round(col.summ, digits = 3),Percentage.Inertia=Percentage.Inertia, Inertia = Principal.Inertia,
       I=I,f=f,a=a,g=g)
}

##### Analisis Correspondence CR ######
corsp.analysis<-function(matrix)
{
  corsp <- ca(matrix)
  KU <- cacoord(corsp,type = c("principal"), dim = NA)
  KUC <- KU$columns
  KUR <- KU$rows
  distance <- dist(KUC,method="euclidean")
  jarak = round(as.matrix(distance), digits = 4)
  hasil=list(Koordinat.Utama = KUC, Jarak = jarak)
  print(hasil)
}

##### VARIABEL X1: Bahan bangunan atap rumah terluas #####
cr1 = ca.regions.exe((data.matrix(cont1)), region = 2); cr1 #X1.Seng, X1.Bambu
# Maka, perlu dilakukan penggabungan kategori
# melihat jarak terdekat
cors1 = corsp.analysis(cont1)
#==== Gabung Data 1 ====
x1re1 = x1
x1re1$X1[x1re1$X1 == "Bambu"] <- "Bambu;Asbes"
x1re1$X1[x1re1$X1 == "Asbes"] <- "Bambu;Asbes"
#== Tabel Kontingensi ==
contingency1re1 = caconv(x1re1, from="rpm", to="freq")
cont1re1 = as.table(contingency1re1);cont1re1
names(dimnames(cont1re1))<-c("K","X1")
chisq.test(cont1re1)
#== Confidence Regions ==
cr1re1 = ca.regions.exe((data.matrix(cont1re1)), region = 2); cr1re1
# Seluruh kategori karakteristik memiliki kontribusi yang signifikan terhadap kebergantungannya


##### VARIABEL X2: Bahan bangunan dinding rumah terluas #####
cr2 = ca.regions.exe((data.matrix(cont2)), region = 2); cr2 #X2.Tembok, X2.Plesteran anyaman bambu/kawat, K.Jakarta Barat
# Maka, perlu dilakukan penggabungan kategori
# melihat jarak terdekat
cors2 = corsp.analysis(cont2)
#==== Gabung Data 1 ====
x2re1 = x2
x2re1$X2[x2re1$X2 == "Tembok"] <- "Tembok;Plesteran anyaman bambu/kawat"
x2re1$X2[x2re1$X2 == "Plesteran anyaman bambu/kawat"] <- "Tembok;Plesteran anyaman bambu/kawat"
#== Tabel Kontingensi ==
contingency2re1 = caconv(x2re1, from="rpm", to="freq")
cont2re1 = as.table(contingency2re1);cont2re1
names(dimnames(cont2re1))<-c("K","X2")
chisq.test(cont2re1)
#== Confidence Regions ==
cr2re1 = ca.regions.exe((data.matrix(cont2re1)), region = 2); cr2re1
# terdapat kategori karakteristik yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya

# melihat jarak terdekat
cors2.re1 = corsp.analysis(cont2re1)
# nilai jarak terkecil Tembok;Plesteran anyaman bambu/kawat dan Kayu/papan
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 2 ====
x2re2 = x2re1
x2re2$X2[x2re2$X2 == "Tembok;Plesteran anyaman bambu/kawat"] <- "Tembok;Plesteran anyaman bambu/kawat;Kayu/papan"
x2re2$X2[x2re2$X2 == "Kayu/papan"] <- "Tembok;Plesteran anyaman bambu/kawat;Kayu/papan"
#== Tabel Kontingensi ==
contingency2re2 = caconv(x2re2, from="rpm", to="freq")
cont2re2 = as.table(contingency2re2);cont2re2
names(dimnames(cont2re2))<-c("K","X2")
chisq.test(cont2re2)
#== Confidence Regions ==
cr2re2 = ca.regions.exe((data.matrix(cont2re2)), region = 2); cr2re2
# terdapat kategori karakteristik yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya

# melihat jarak terdekat
cors2.re2 = corsp.analysis(cont2re2)
# nilai jarak terkecil Tembok;Plesteran anyaman bambu/kawat;Kayu/papan dan Batang kayu
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 3 ====
x2re3 = x2re2
x2re3$X2[x2re3$X2 == "Tembok;Plesteran anyaman bambu/kawat;Kayu/papan"] <- "Tembok;Plesteran anyaman bambu/kawat;Kayu/papan;Batang kayu"
x2re3$X2[x2re3$X2 == "Batang kayu"] <- "Tembok;Plesteran anyaman bambu/kawat;Kayu/papan;Batang kayu"
#== Tabel Kontingensi ==
contingency2re3 = caconv(x2re3, from="rpm", to="freq")
cont2re3 = as.table(contingency2re3);cont2re3
names(dimnames(cont2re3))<-c("K","X2")
chisq.test(cont2re3)
#== Confidence Regions ==
cr2re3 = ca.regions.exe((data.matrix(cont2re3)), region = 2); cr2re3
# terdapat kategori karakteristik yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya


##### VARIABEL X3: Bahan bangunan lantai rumah terluas #####
cr3 = ca.regions.exe((data.matrix(cont3)), region = 2); cr3 #X3.Keramik, X3.Parket/vinil/karpet, X3.Semen/bata merah, K.Jakarta Barat
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 1 ====
x3re1 = x3
x3re1$X3[x3re1$X3 == "Parket/vinil/karpet"] <- "Parket/vinil/karpet;Marmer/granit"
x3re1$X3[x3re1$X3 == "Marmer/granit"] <- "Parket/vinil/karpet;Marmer/granit"
#== Tabel Kontingensi ==
contingency3re1 = caconv(x3re1, from="rpm", to="freq")
cont3re1 = as.table(contingency3re1);cont3re1
names(dimnames(cont3re1))<-c("K","X3")
chisq.test(cont3re1)
#== Confidence Regions ==
cr3re1 = ca.regions.exe((data.matrix(cont3re1)), region = 2); cr3re1
# terdapat kategori yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya

# melihat jarak terdekat
cors3.re1 = corsp.analysis(cont3re1)
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 2 ====
x3re2 = x3re1
x3re2$X3[x3re2$X3 == "Keramik"] <- "Keramik;Semen/bata merah"
x3re2$X3[x3re2$X3 == "Semen/bata merah"] <- "Keramik;Semen/bata merah"
#== Tabel Kontingensi ==
contingency3re2 = caconv(x3re2, from="rpm", to="freq")
cont3re2 = as.table(contingency3re2);cont3re2
names(dimnames(cont3re2))<-c("K","X3")
chisq.test(cont3re2)
#== Confidence Regions ==
cr3re2 = ca.regions.exe((data.matrix(cont3re2)), region = 2); cr3re2
# terdapat kategori yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya

# melihat jarak terdekat
cors3.re2 = corsp.analysis(cont3re2)
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 3 ====
x3re3 = x3re2
x3re3$X3[x3re3$X3 == "Keramik;Semen/bata merah"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit"
x3re3$X3[x3re3$X3 == "Parket/vinil/karpet;Marmer/granit"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit"
#== Tabel Kontingensi ==
contingency3re3 = caconv(x3re3, from="rpm", to="freq")
cont3re3 = as.table(contingency3re3);cont3re3
names(dimnames(cont3re3))<-c("K","X3")
chisq.test(cont3re3)
#== Confidence Regions ==
cr3re3 = ca.regions.exe((data.matrix(cont3re3)), region = 2); cr3re3
# terdapat kategori yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya

# melihat jarak terdekat
cors3.re3 = corsp.analysis(cont3re3)
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 4 ====
x3re4 = x3re3
x3re4$X3[x3re4$X3 == "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit;Ubin/tegel/teraso"
x3re4$X3[x3re4$X3 == "Ubin/tegel/teraso"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit;Ubin/tegel/teraso"
#== Tabel Kontingensi ==
contingency3re4 = caconv(x3re4, from="rpm", to="freq")
cont3re4 = as.table(contingency3re4);cont3re4
names(dimnames(cont3re4))<-c("K","X3")
chisq.test(cont3re4)
#== Confidence Regions ==
cr3re4 = ca.regions.exe((data.matrix(cont3re4)), region = 2); cr3re4
# Seluruh kategori karakteristik memiliki kontribusi yang signifikan terhadap kebergantungannya

#==== Gabung Data 5 ====
x3re5 = x3re4
x3re5$X3[x3re5$X3 == "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit;Ubin/tegel/teraso"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit;Ubin/tegel/teraso;Kayu/papan"
x3re5$X3[x3re5$X3 == "Kayu/papan"] <- "Keramik;Semen/bata merah;Parket/vinil/karpet;Marmer/granit;Ubin/tegel/teraso;Kayu/papan"
#== Tabel Kontingensi ==
contingency3re5 = caconv(x3re5, from="rpm", to="freq")
cont3re5 = as.table(contingency3re5);cont3re5
names(dimnames(cont3re5))<-c("K","X3")
chisq.test(cont3re5)
#== Confidence Regions ==
cr3re5 = ca.regions.exe((data.matrix(cont3re5)), region = 2); cr3re5
# Seluruh kategori karakteristik memiliki kontribusi yang signifikan terhadap kebergantungannya


##### VARIABEL X4: Sumber air utama yang digunakan untuk minum #####
cr4 = ca.regions.exe((data.matrix(cont4)), region = 2); cr4 #X4.Lainnya, X4.Mata air tak terlindung, X4.Mata air terlindung, X4.Sumur terlindung  
# Maka, perlu dilakukan penggabungan kategori
# melihat jarak terdekat
cors4 = corsp.analysis(cont4)
# nilai jarak terkecil Lainnya dan Mata air tak terlindung
# nilai jarak terkecil Mata air terlindung dan Sumur terlindung
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 1 ====
x4re1 = x4
x4re1$X4[x4re1$X4 == "Lainnya"] <- "Lainnya;Mata air tak terlindung"
x4re1$X4[x4re1$X4 == "Mata air tak terlindung"] <- "Lainnya;Mata air tak terlindung"
#== Tabel Kontingensi ==
contingency4re1 = caconv(x4re1, from="rpm", to="freq")
cont4re1 = as.table(contingency4re1);cont4re1
names(dimnames(cont4re1))<-c("K","X4")
chisq.test(cont4re1)
#== Confidence Regions ==
cr4re1 = ca.regions.exe((data.matrix(cont4re1)), region = 2); cr4re1

#==== Gabung Data 2 ====
x4re2 = x4re1
x4re2$X4[x4re2$X4 == "Mata air terlindung"] <- "Mata air terlindung;Sumur terlindung"
x4re2$X4[x4re2$X4 == "Sumur terlindung"] <- "Mata air terlindung;Sumur terlindung"
#== Tabel Kontingensi ==
contingency4re2 = caconv(x4re2, from="rpm", to="freq")
cont4re2 = as.table(contingency4re2);cont4re2
names(dimnames(cont4re2))<-c("K","X4")
chisq.test(cont4re2)
#== Confidence Regions ==
cr4re2 = ca.regions.exe((data.matrix(cont4re2)), region = 2); cr4re2
# Seluruh kategori memiliki kontribusi yang signifikan terhadap kebergantungannya


##### VARIABEL X5: Sumber air utama yang digunakan rumah tangga untuk mandi/cuci/dll #####
cr5 = ca.regions.exe((data.matrix(cont5)), region = 2); cr5 #X5.Air hujan, X5.Air kemasan bermerk  
# Maka, perlu dilakukan penggabungan kategori
# melihat jarak terdekat
cors5 = corsp.analysis(cont5)
# nilai jarak terkecil Air hujan dan Mata air terlindung sebesar 0.0
# nilai jarak terkecil Air kemasan bermerk dan Sumur bor/pompa
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 1 ====
x5re1 = x5
x5re1$X5[x5re1$X5 == "Air hujan"] <- "Air hujan;Mata air terlindung"
x5re1$X5[x5re1$X5 == "Mata air terlindung"] <- "Air hujan;Mata air terlindung"

#== Tabel Kontingensi ==
contingency5re1 = caconv(x5re1, from="rpm", to="freq")
cont5re1 = as.table(contingency5re1);cont5re1
names(dimnames(cont5re1))<-c("K","X5")
chisq.test(cont5re1)
#== Confidence Regions ==
cr5re1 = ca.regions.exe((data.matrix(cont5re1)), region = 2); cr5re1
# Seluruh kategori memiliki kontribusi yang signifikan terhadap kebergantungannya

#==== Gabung Data 2 ====
x5re2=x5re1
x5re2$X5[x5re2$X5 == "Air isi ulang"] <- "Air isi ulang;Air kemasan bermerk"
x5re2$X5[x5re2$X5 == "Air kemasan bermerk"] <- "Air isi ulang;Air kemasan bermerk"

#== Tabel Kontingensi ==
contingency5re2 = caconv(x5re2, from="rpm", to="freq")
cont5re2 = as.table(contingency5re2);cont5re2
names(dimnames(cont5re2))<-c("K","X5")
chisq.test(cont5re2)
#== Confidence Regions ==
cr5re2 = ca.regions.exe((data.matrix(cont5re2)), region = 2); cr5re2
# Seluruh kategori memiliki kontribusi yang signifikan terhadap kebergantungannya


##### VARIABEL X6: Kepemilikan fasilitas tempat buang air besar dan siapa saja yang menggunakan #####
cr6 = ca.regions.exe((data.matrix(cont6)), region = 2); cr6 #K.Jakarta Barat, K.Jakarta Utara 


##### VARIABEL X7: Jenis kloset yang digunakan #####
cr7 = ca.regions.exe((data.matrix(cont7)), region = 2); cr7 #X7.Leher angsa K.Jakarta Barat, K.Jakarta Utara


##### VARIABEL X8: Tempat pembuangan akhir tinja #####
cr8 = ca.regions.exe((data.matrix(cont8)), region = 2); cr8 #X8.Tangki septik
# Maka, perlu dilakukan penggabungan kategori
#==== Gabung Data 1 ====
x8re1 = x8
x8re1$X8[x8re1$X8 == "Tangki septik"] <- "IPAL;Tangki septik"
x8re1$X8[x8re1$X8 == "IPAL"] <- "IPAL;Tangki septik"
#== Tabel Kontingensi ==
contingency8re1 = caconv(x8re1, from="rpm", to="freq")
cont8re1 = as.table(contingency8re1);cont8re1
names(dimnames(cont8re1))<-c("K","X8")
chisq.test(cont8re1)
#== Confidence Regions ==
cr8re1 = ca.regions.exe((data.matrix(cont8re1)), region = 2); cr8re1
# terdapat kategori karakteristik yang tidak memiliki kontribusi yang signifikan terhadap kebergantungannya
# iterasi selesai kategori karakteristik tidak ada yang bisa digabung lagi


##### VARIABEL X9: Intensitas tangki septik dikosongkan/dilakukan penyedotan dalam 5 tahun terakhir #####
cr9 = soca.regions.exe(data.matrix(cont9)); cr9 #K.Jakarta Barat, K.Jakarta Pusat 


####### Multipel #######
df = as.data.frame(c(x1re1,(x2re3[2]),(x3re5[2]),(x4re2[2]),(x5re2[2]),(x6[2]),(x7[2]),(x8re1[2]),(x9[2]),(x10[2]))) 
head(df)
mca2<-mjca(df, nd=NA, lambda = "Burt", reti=TRUE);mca2
mca2$sv^2
summary(mca2) 
#--# 2 dim 23.7%
#--# 21 dim 71.4%
#--# 46 dim 100%

## Diperoleh Hasil Sebagai berikut:
## Total Dimensi = 46
## 60 -> 46
## Pada dimensi 2 diperoleh hasil cumulative percentage of variance adalah sebesar 23.7% asalnya 21%
## dim21    cum%71.4
#--# ada 59 dari asalnya 73 shg bener sudah di rekategorisasi


##Koordinat Standar 
KS = mca2$colcoord ; KS # Koordinat Standar: Menggambarkan Asosiasi antar kategori
#atau
KS = cacoord(mca2, type=c("standard"), dim=NA)
KSC = KS$columns 
head(KSC)
KS = data.frame(as.matrix(KSC))
dim(KS) #59 kategori - 11 variabel = 48
#write_xlsx(KS, "3.Multipel_1 Koordinat Standar.xlsx")

##Koordinat Utama 
KU<-mca2$colpcoord ; KU # Koordinat Utama
#atau
KU = cacoord(mca2, type=c("principal"), dim=NA)
KUC = KU$columns ; KUC
KU = data.frame(as.matrix(KUC))
dim(KU)
#write_xlsx(KU, "3.Multipel_2 Koordinat Utama.xlsx")

##Matriks Jarak Euclidean
# Untuk Semua Dimensi
distance=dist(KU,method="euclidean") ; distance
distance=as.matrix(distance)
Jarak = data.frame(distance)
head(Jarak)
#--# ada 53 (karakteristik) + 6 (kabkota)=59 dari asalnya 67 (karakteristik) + 6 (kabkota)=73 shg bener sudah di rekategorisasi
#--# dim 46 dari asalnya 60 di 100%
#lihat berdasarkan baris dan per kolom x1 hingga x9 (nilai terkecil)
#write_xlsx(Jarak, "3.Multipel_3 Matriks Jarak bismillah_skripsi_draft_4_try1 FIX.xlsx")

####### Multipel manual #######
##Matriks Indikator 
indikator<-read.csv(file.choose(),header=T) ###read file : data clean matriks indikator 
head(indikator)
Z <- as.matrix(indikator)
head(Z)  
dim(Z)

##Matriks Burt 
B<-t(Z)%*%Z 
head(B)
dim(B)
burt = data.frame(as.matrix(B)) 
#write_xlsx(burt, "4.MultipelManual_2 Matriks Burt.xlsx") 

##Matriks Korespondensi Burt 
g<-sum(B) 
g
P<-B/g 
head(P) 
KB=data.frame(as.matrix(P)) 
dim(KB)
#write_xlsx(KB, "4.MultipelManual_3 Matriks Korespondensi Burt.xlsx") 

##Matriks Residual Standar 
###proporsi kolom Burt (penjumlahan terhadap kolom dari matriks korespondensi)
C <- apply(P, 2, sum)  
C 
prop.C=data.frame(as.matrix(C)) 
dim(prop.C)
write_xlsx(prop.C, "4.MultipelManual_4 Matriks Proporsi Kolom Burt.xlsx") 
###Matriks diagonal C 
diagC <- diag(C) 
diagC 
###matriks residual standar 
eP <- C %*% t(C) 
S = (P-eP)/sqrt(eP); S
head(S) 
res=data.frame(as.matrix(S)) ;dim(res)
#write_xlsx(res, "4.MultipelManual_5 Matriks Residual Standar.xlsx") 

##Inersia manual 
J<-dim(indikator)[2] 
Q<-dim(df[,-1])[2] 
dec<-eigen(S) 
delt<-dec$values[1:(J-Q)] 
expl<-100*(delt/sum(delt)) 
lam<-delt^2 
expl2<-100*(lam/sum(lam)) 
rbind(round(lam,3),round(expl2,1)) 

##Koordinat Standar 
KS<-cacoord(mca2,type = c("standard"), dim = NA) 
KSC<-KS$columns 
head(KSC) 
ks=data.frame(as.matrix(KSC)) 
#write_xlsx(ks, "4.MultipelManual_6 Koordinat Standar.xlsx") 

##Koordinat Utama 
KU<-cacoord(mca2,type = c("principal"), dim = NA) 
KUC<-KU$columns 
head(KUC)
ku=data.frame(as.matrix(KUC)) 
#write_xlsx(ku, "4.MultipelManual_7 Koordinat Utama.xlsx") 

##Matriks Euclidean 
distance=dist(KUC,method="euclidean") 
jarak=data.frame(as.matrix(distance)) 
head(jarak) 
#Menyimpan Matriks Euclidean 
#write_xlsx(jarak, "4.MultipelManual_8 Matriks Jarak Euclidean.xlsx") 

####### JCA #######
#Matriks asosiasi MCA
dataJCA <- mjca(df, nd=2, lambda = "JCA",maxit=1000, reti=TRUE,epsilon=0.0001)
dataJCA$Burt.upd
#View(dataJCA$Burt)


summary(dataJCA) ###inersia total
dataJCA$colinertia
singular=dataJCA$sv
#total inersia
eigen=singular^2
total_inertia=sum(eigen)
diagonal=diag(eigen)
inertia_percentage <- (eigen/total_inertia)*100
cum_inertia_percentage <- cumsum(inertia_percentage);cum_inertia_percentage

plot(dataJCA)
#view(dataJCA$Burt.upd)
#view(dataJCA$indmat)
ind=dataJCA$indmat
#write.csv(ind, file="5.JCA_1  Matriks Indikator Awal.csv", row.names=T)

k=dim(ind)[2]
B=dataJCA$Burt
#write.table(B, file="5.JCA_2 Matriks Burt Akhir.csv", row.names=T,sep=";")
B2=dataJCA$Burt.upd
#write.table(B2, file="5.JCA_3 Matriks Rekonstruksi Burt Akhir.csv", row.names=T,sep=";")
b=sum(B)
P=B/b
dim(P)
#write.table(P, file="5.JCA_4 Matriks Korespondensi Awal.csv", row.names=T,sep=";")
c=apply(P,2,sum)
#write.table(c, file="5.JCA_5 Vektor Proporsi Kolom Awal.csv", row.names=T,sep=";")
Dc=diag(c)
eP <- c %*% t(c)
diagC=diag(c)
sn <- solve(sqrt(diagC)) %*% (P - eP) %*% solve(sqrt(diagC))
S <- (P - eP) / sqrt(eP)
dim(S)
#write.table(S, file="5.JCA_6 Matriks Residual Akhir.csv", row.names=T,sep=";")
#view(sn)
#view(S)

eigen_f = (mca2$sv)^2
diagv = diag(eigen_f)
#view(diage)
dec=eigen(S)
L=k-9
eig_value=dec$values[1:L]
#write.table(eig_value, file="5.JCA_7 Vektor Eigen Awal.csv", row.names=T,sep=";")
diage=diag(eig_value);diage
V=dec$vectors[,1:L]
dim(V)
#write.table(V, file="5.JCA_8 Matriks V Awal.csv", row.names=T,sep=";")
SR = V%*%diage%*%(t(V))

##Koordinat standar
KS_r<-cacoord(dataJCA,type = c("standard"))
KSC_r<-KS_r$columns
#write.table(KSC_r, file="5.JCA_9 Matriks Koordinat Standar Awal.csv", row.names=T,sep=";")
#view(KSC_r)

Kor = dataJCA$colpcoord[1:6,]
#view(Kor)
#write.table(Kor, file="5.JCA_910 Matriks Koordinat Utama Awal.csv", row.names=T,sep=";")
dataJCA$colcoord[1:6,]
distanceK <- dist(Kor,method="euclidean")
jarakK = as.matrix(distanceK)
write.table(jarakK, file="5.JCA_911 Matriks Jarak Euclidean pertama.csv", row.names=T,sep=";")

####### df final #######
df_final=df
head(df_final)
head(df_final$K)
unique(df_final$K)
df_final$K[df_final$K == "Jakarta Barat"]<- "JB"
df_final$K[df_final$K == "Jakarta Pusat"]<- "JP"
df_final$K[df_final$K == "Jakarta Selatan"]<- "JS"
df_final$K[df_final$K == "Jakarta Timur"]<- "JT"
df_final$K[df_final$K == "Jakarta Utara"]<- "JU"
df_final$K[df_final$K == "Kepulauan Seribu"]<- "KS"
head(df_final$K)

df_final$X1 = as.factor(df_final$X1)
head(df_final$X1)
unique(df_final$X1)
df_final$X1 = as.numeric(df_final$X1)
head(df_final$X1)
unique(df_final$X1)
df_final$X1[df_final$X1 == "6"]<- "7"
df_final$X1[df_final$X1 == "5"]<- "6"
df_final$X1[df_final$X1 == "4"]<- "5"
df_final$X1[df_final$X1 == "3"]<- "4"
df_final$X1[df_final$X1 == "2"]<- "3"
df_final$X1[df_final$X1 == "1"]<- "1,2"
head(df_final$X1)

df_final$X2 = as.factor(df_final$X2)
head(df_final$X2)
unique(df_final$X2)
df_final$X2 = as.numeric(df_final$X2)
head(df_final$X2)
unique(df_final$X2)
df_final$X2[df_final$X2 == "3"]<- "5"
df_final$X2[df_final$X2 == "4"]<- "3,4,6,7"
head(df_final$X2)

df_final$X3 = as.factor(df_final$X3)
head(df_final$X3)
unique(df_final$X3)
df_final$X3 = as.numeric(df_final$X3)
head(df_final$X3)
unique(df_final$X3)
df_final$X3[df_final$X3 == "2"]<- "6"
df_final$X3[df_final$X3 == "1"]<- "1,2,3,4,5,7"
head(df_final$X3)

df_final$X4 = as.factor(df_final$X4)
head(df_final$X4)
unique(df_final$X4)
df_final$X4 = as.numeric(df_final$X4)
head(df_final$X4)
unique(df_final$X4)
df_final$X4[df_final$X4 == "7"]<- "8"
df_final$X4[df_final$X4 == "6"]<- "7,9"
df_final$X4[df_final$X4 == "4"]<- "4,6"
head(df_final$X4)

df_final$X5 = as.factor(df_final$X5)
head(df_final$X5)
unique(df_final$X5)
df_final$X5 = as.numeric(df_final$X5)
head(df_final$X5)
unique(df_final$X5)
df_final$X5[df_final$X5 == "9"]<- "11"
df_final$X5[df_final$X5 == "8"]<- "10"
df_final$X5[df_final$X5 == "7"]<- "9"
df_final$X5[df_final$X5 == "6"]<- "7"
df_final$X5[df_final$X5 == "5"]<- "6"
df_final$X5[df_final$X5 == "4"]<- "5"
df_final$X5[df_final$X5 == "3"]<- "4"
df_final$X5[df_final$X5 == "2"]<- "2,3"
df_final$X5[df_final$X5 == "1"]<- "1,8"
head(df_final$X5)

df_final$X6 = as.factor(df_final$X6)
head(df_final$X6)
unique(df_final$X6)
df_final$X6 = as.numeric(df_final$X6)
head(df_final$X6)
unique(df_final$X6)
head(df_final$X6)

df_final$X7 = as.factor(df_final$X7)
head(df_final$X7)
unique(df_final$X7)
df_final$X7 = as.numeric(df_final$X7)
head(df_final$X7)
unique(df_final$X7)
head(df_final$X7)

df_final$X8 = as.factor(df_final$X8)
head(df_final$X8)
unique(df_final$X8)
df_final$X8 = as.numeric(df_final$X8)
head(df_final$X8)
unique(df_final$X8)
df_final$X8[df_final$X8 == "6"]<- "7"
df_final$X8[df_final$X8 == "1"]<- "1,6"
head(df_final$X8)

df_final$X9 = as.factor(df_final$X9)
head(df_final$X9)
unique(df_final$X9)
df_final$X9 = as.numeric(df_final$X9)
head(df_final$X9)
unique(df_final$X9)
head(df_final$X9)

df_final$X10 = as.factor(df_final$X10)
head(df_final$X10)
unique(df_final$X10)
df_final$X10 = as.numeric(df_final$X10)
head(df_final$X10)
unique(df_final$X10)
head(df_final$X10)

str(df_final)

####### JCA final #######
JCA_final <- mjca(df_final, nd =2, lambda = "JCA",maxit=160, reti=TRUE,epsilon=0.0001)
summary(JCA_final) #0.000096 #2d 71,6%
help(mjca)
JCA_final$colpcoord

sv_f = JCA_final$sv
##Koordinat Standar
KS_f<-cacoord(JCA_final,type = c("standard"))
KSC_f<-KS_f$columns
head(KSC_f)
##Koordinat Utama
KU_f<-cacoord(JCA_final,type = c("principal"), dim = NA)
KUC_f<-KU_f$columns
rownames(KUC_f)<-rownames(JCA_final$Burt)
#write.table(KUC_f, file="5.JCA_912 Matriks Koordinat Utama Akhir.csv", row.names=T,sep=";")

####====Plot dari JCA=====####
Keterangan=c(1,1,1,1,1,1,
             2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2)
KUC_new=cbind(KUC_f,Keterangan)
ggplot(as.data.frame(KUC_new), aes(x = Dim1, y = Dim2, label = rownames(KUC_new),shape=factor(Keterangan))) +
  geom_point(aes(colour = factor(Keterangan)), size = 4)+
  geom_text(hjust = 0, vjust = 1, size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))

#KUC ZOOM
KUC_p=KUC_f[-c(6,13,18,19,24,26,27,28,30,33,36,39,40,41,44,48,49,50), ] #18
Ket=c(1,1,1,1,1,
      2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2)
KUC_p = cbind(KUC_p,Ket)

ggplot(as.data.frame(KUC_p), aes(x = Dim1, y = Dim2, label = rownames(KUC_p),shape=factor(Ket))) +
  geom_point(aes(colour = factor(Ket)), size = 4)+
  geom_text(hjust = 1, vjust = 1, size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))+
  labs(shape = "Keterangan", colour = "Keterangan")


#MATRIKS JARAK EUCLIDEAN untuk melihat kebergantungan
K_f<-KUC_f[1:59,1:2]
distance_f=dist(K_f,method="euclidean")
df_distance=as.matrix(distance_f)
head(df_distance)
#view(df_distance)
#write.table(df_distance, file="5.JCA_913 df_distance Matriks Jarak Euclidean bismillah_skripsi_draft_4_try1 FIX.csv", row.names=T,sep=",")

##### Nilai Stress Setelah rekategorisasi#####
#Matriks asosiasi MCA
datamca=mjca(df_final, nd=2, lambda = "Burt")
summary(datamca)
dataJCA_ellip <- mjca(df_final, nd=2, lambda = "JCA",maxit=1000,reti=TRUE,epsilon=0.0001)
summary(dataJCA_ellip) ###inersia total
#View(dataJCA$Burt.upd)
#write.table(dataJCA_ellip$Burt.upd, file="6.STRESS_1 JCA Matriks Rekonstruksi Burt Akhir -notes_samain .csv",row.names=T,sep=";")
singular=dataJCA_ellip$sv
# Hitung total inersia
eigen = singular^2
total_inertia <- sum(eigen)
diagonal = diag(eigen)
# Hitung persentase inersia untuk setiap dimensi
inertia_percentage <- (eigen / total_inertia) * 100
cumulative_inertia_percentage <- cumsum(inertia_percentage)
cumulative_inertia_percentage

#View(dataJCA$indmat)
ind=dataJCA$indmat

####### Jarak euclidean ####### 
####### MCA ####### 
#mca<-mjca(data, nd = 2, lambda = c("adjusted", "indicator", "Burt","JCA"))
mca<-mjca(df_final, nd=NA, lambda = "Burt", maxit=0, reti=TRUE)
summary(mca) #inersia total
#Matriks jarak euclidean pada 63 dimensi
KU<-cacoord(mca,type = c("principal"), dim = NA)
KUC<-KU$columns
#write.table(KUC, file="6.STRESS_2 MCA Koordinat Utama -notes_samain .csv",row.names=T,sep=";")

dim(KUC) #59 48
K_mca_rek<-KUC[1:59,1:46]
write.table(K_mca_rek, file="6.STRESS_3 MCA Koordinat Utama.46dim -notes_samain.csv",row.names=T,sep=";")

distance_mca_rek=dist(K_mca_rek,method="euclidean")
jarak_rek_mca<-as.matrix(distance_mca_rek)
write.table(jarak_rek_mca, file="6.STRESS_4 MCA Matriks Jarak Euclidean.46dim -notes_samain.csv",row.names=T,sep=";")

jarak1kolom<-as.matrix(c(jarak_rek_mca))
head(jarak1kolom)
dim(jarak1kolom)
#write.table(jarak1kolom, file="6.STRESS_5 MCA Kolom 1 Jarak Euclidean.46dim.csv",row.names=T,sep=";")
# Mendapatkan elemen-elemen dari segitiga bawah
lower_trianglemca <- jarak[lower.tri(jarak, diag = TRUE)]
#write.table(lower_trianglemca, file="6.STRESS_6 MCA SEGITIGA Jarak Euclidean.csv",row.names=T,sep=";")

####### JCA ####### 
#dataJCA <- mjca(df_final, nd=2, lambda = "JCA",reti=TRUE,epsilon=0.0001)
#Matriks jarak euclidean pada 2 dimensi 2 iterasi
KU2<-cacoord(dataJCA_ellip,type = c("principal"), dim = NA)
KUC2<-KU2$columns
dim(KUC2)
K2<-KUC2[1:59,1:2]
write.table(K2, file="6.STRESS_7 JCA koordinat Utama 2d-notes_samain.csv",row.names=T,sep=";")

distance2=dist(K2,method="euclidean")
jarak2<-as.matrix(distance2)
jarak2
write.table(jarak2, file="6.STRESS_8 JCA Matriks Jarak Euclidean 2d-notes_samain.csv",row.names=T,sep=";")

jarak2kolom1<-as.matrix(c(jarak2))
write.table(jarak2kolom1, file="6.STRESS_9 JCA Kolom 1 Jarak Euclidean.2dim.csv",row.names=T,sep=";")
S_stress2 = sum((jarak-jarak2)^2)/sum(jarak^2);S_stress2 
#---# #0.8

##### STOP ######
