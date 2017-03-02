# Originally from http://rfunctions.blogspot.fi/2015/03/bivariate-maps-bivariatemap-function.html
library(classInt)
library(dismo)
library(maps)
library(raster)
library(rgdal)
library(sp)
library(XML)

#' Create a color matrix
#'
#' @examples
#' # You can specify the number of quantiles, colors and labels of your color
#' # matrix:
#' col.matrix <- colmat(nquantiles = 10, upperleft = "blue", upperright = "yellow",
#'                      bottomleft = "green", bottomright = "red",
#'                      xlab = "My x label", ylab = "My y label")
#'
#' # But let's use this simple code. You can change "nquantiles" to generate color
#' # matrices with different color schemes. For example, change it to 4 to produce
#' # a 4x4 color scheme.
#' col.matrix <- colmat(nquantiles = 4)
#'
colmat <- function(nquantiles = 10,
                   upperleft = rgb(0,150,235, maxColorValue = 255),
                   upperright = rgb(130,0,80, maxColorValue = 255),
                   bottomleft = "grey",
                   bottomright = rgb(255,230,15, maxColorValue = 255),
                   xlab = "x label", ylab = "y label") {

  my.data <- seq(0,1,.01)
  my.class <- classIntervals(my.data, n = nquantiles, style = "quantile")
  my.pal.1 <- findColours(my.class, c(upperleft, bottomleft))
  my.pal.2 <- findColours(my.class, c(upperright, bottomright))
  col.matrix <- matrix(nrow = 101, ncol = 101, NA)

  for (i in 1:101) {
    my.col <- c(paste(my.pal.1[i]), paste(my.pal.2[i]))
    col.matrix[102 - i,] <- findColours(my.class, my.col)
  }

  plot(c(1, 1), pch = 19, col = my.pal.1, cex = 0.5, xli = c(0,1),
       ylim = c(0,1), frame.plot = FALSE, xlab = xlab, ylab = ylab,
       cex.lab = 1.3)

  for (i in 1:101) {
    col.temp <- col.matrix[i - 1,]
    points(my.data, rep((i - 1) / 100, 101), pch = 15, col = col.temp, cex = 1)}
  seqs <- seq(0, 100, (100 / nquantiles))
  seqs[1] <- 1
  col.matrix <- col.matrix[c(seqs), c(seqs)]
}

#' Plot bivariate map
#'
bivariate.map <- function(rasterx, rastery, colormatrix = NULL, nquantiles = 10) {

  if (is.null(colormatrix)) {
    colormatrix <- colmat(nquantiles)
  }

  quanmean <- getValues(rasterx)
  temp <- data.frame(quanmean, quantile = rep(NA, length(quanmean)))
  brks <- with(temp, quantile(temp, na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))))
  r1 <- within(temp, quantile <- cut(quanmean, breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE))
  quantr <- data.frame(r1[,2])

  quanvar <- getValues(rastery)
  temp <- data.frame(quanvar, quantile = rep(NA, length(quanvar)))
  brks <- with(temp, quantile(temp,na.rm = TRUE,
                              probs = c(seq(0, 1, 1 / nquantiles))))
  r2 <- within(temp, quantile <- cut(quanvar, breaks = brks,
                                     labels = 2:length(brks),
                                     include.lowest = TRUE))
  quantr2 <- data.frame(r2[, 2])

  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  col.matrix2 <- colormatrix
  cn <- unique(colormatrix)

  for (i in 1:length(col.matrix2)) {
    ifelse(is.na(col.matrix2[i]),
           col.matrix2[i] <- 1,
           col.matrix2[i] <- which(col.matrix2[i] == cn)[1])
  }

  cols <- numeric(length(quantr[, 1]))
  for (i in 1:length(quantr[, 1])) {
    a <- as.numeric.factor(quantr[i, 1])
    b <- as.numeric.factor(quantr2[i, 1])
    cols[i] <- as.numeric(col.matrix2[b, a])
  }
  r <- rasterx
  r[1:length(r)] <- cols
  return(r)
}
