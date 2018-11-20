## Load package
library(igraph)
# http://igraph.org/r/doc/

library(corrplot)

# https://github.com/cwatson/brainGraph
# http://kateto.net/network-visualization

rm(list=ls())
setwd("E:/UBA/2018-II/DM en Ciencia y Tecnología/Ciclos Sueño")


install.packages("igraph");

download.file("http://moreno.ss.uci.edu/beach.dat", destfile = "windsurfers.dat")
ws <- read.table("windsurfers.dat", skip = 7)
dim(ws)

ws.obs <- as.matrix(ws[1:43, ])
ws.per <- as.matrix(ws[44:86, ])