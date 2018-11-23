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

ws.obs.red <- graph.adjacency(ws.obs, mode="undirected", diag=FALSE, weighted = T)
plot(ws.obs.red)

ws.per.red <- graph.adjacency(ws.per, mode="undirected", diag=FALSE, weighted = T)
plot(ws.per.red)


hist(ws.per[lower.tri(ws.per)], main = "Histograma Interacciones Percibidas")

umbral <-0.5
ws.per.2 <-ws.per
ws.per.2[which(ws.per.2 <= umbral)] <- 0

ws.per.red <- graph.adjacency(ws.per.2, mode="undirected", diag=FALSE, weighted = T)
plot(ws.per.red)

summary(ws.obs.red)
summary(ws.per.red)

#Número nodos y aristas
vcount(ws.obs.red)
ecount(ws.obs.red)

#Número nodos y aristas
vcount(ws.per.red)
ecount(ws.per.red)

#Nodos y aristas
V(ws.obs.red)
E(ws.obs.red)

V(ws.per.red)
E(ws.per.red)

#¿Loops o aristas múltiples?
is.simple(ws.obs.red)
is.simple(ws.per.red)

#¿Completamente conectado?
is.connected(ws.obs.red)
is.connected(ws.per.red)

#diametro red
diameter(ws.obs.red)
get.diameter(ws.obs.red)

diameter(ws.per.red)
get.diameter(ws.per.red)

#Densidad del grafo, relacion entre aristas y nodos
graph.density(ws.obs.red)
graph.density(ws.per.red)

#Coeficiente de clustering global
head(transitivity(ws.obs.red, type = "local"))
head(transitivity(ws.per.red, type = "local"))

#Red de interacciones observadas
transitivity(ws.obs.red, type = "global")

#Red de interacciones percibidas
transitivity(ws.per.red, type = "global")

#Comparar dos coeficientes de clustering
par(mfrow = c(1,2))
hist(transitivity(ws.obs.red, type = "local"), 
     main = "Observada",breaks = seq(0.2, 1, 0.1), xlab = "coefs. de clustering")
hist(transitivity(ws.per.red, type = "local"), 
     main = "Percibida",breaks = seq(0.2, 1, 0.1), xlab = "coefs. de clustering")
par(mfrow = c(1,1))


#Grados de entrada y salida
degree(ws.obs.red)
degree(ws.per.red)

sort(degree(ws.obs.red), decreasing = T)

#Empezamos a analizar que sucede a nivel de grados

qplot(degree(ws.per.red), degree(ws.obs.red))

#Correlación entre el grado que predice el grupo para cada individuo y el grado observado
cor(degree(ws.per.red), degree(ws.obs.red))


#distribuciones de grados y diustribuciones acumuladas
head(degree.distribution(ws.obs.red), 15)
head(degree.distribution(ws.obs.red, cumulative = T))

par(mfrow = c(1,2))
plot(degree.distribution(ws.obs.red),
                         xlab = "Grados", ylab = "Proporción de nodos", type = "h", main ="Observadas")
plot(degree.distribution(ws.per.red),
      xlab = "Grados", ylab = "Proporción de nodos", type = "h", main ="Percibidas")
par(mfrow = c(1,1))


#Prueba mundo pequeño
ws.obs.red.plf <- power.law.fit(degree(ws.obs.red)) 
ws.per.red.plf <- power.law.fit(degree(ws.per.red)) 


#Cálculo de asortividad -> en ambos casos los valores sugieren que no hay asociaciones preferenciales
#entre nodos de un alto grado y por otrolado de bajo grado
assortativity.degree(ws.obs.red)
assortativity.degree(ws.per.red)

#Medidas de centralidad
#Intermedicion
head(sort(betweenness(ws.obs.red), decreasing = T))
head(sort(betweenness(ws.per.red), decreasing = T))

#Cercania
head(sort(closeness(ws.obs.red), decreasing = T))
head(sort(closeness(ws.per.red), decreasing = T))

#Centralidad de autovectores
head(sort(eigen_centrality(ws.obs.red)$vector, decreasing = T))
head(sort(eigen_centrality(ws.per.red)$vector, decreasing = T))

























