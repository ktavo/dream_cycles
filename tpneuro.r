# https://github.com/cwatson/brainGraph
# http://kateto.net/network-visualization

rm(list=ls())
setwd("E:/UBA/2018-II/DM en Ciencia y Tecnología/Ciclos Sueño")

## Load package
library(igraph)
# http://igraph.org/r/doc/

library(corrplot)

########################Generación Redes Promedio########################
N1 <- read.csv("N1promedio.csv",header=FALSE)

#############################N1#############################
N1_suj1 <- read.csv("DataSujetos/N1_suj1.csv",header=FALSE)
N1_suj2 <- read.csv("DataSujetos/N1_suj2.csv",header=FALSE)
N1_suj3 <- read.csv("DataSujetos/N1_suj3.csv",header=FALSE)
N1_suj4 <- read.csv("DataSujetos/N1_suj4.csv",header=FALSE)
N1_suj5 <- read.csv("DataSujetos/N1_suj5.csv",header=FALSE)
N1_suj6 <- read.csv("DataSujetos/N1_suj6.csv",header=FALSE)
N1_suj7 <- read.csv("DataSujetos/N1_suj7.csv",header=FALSE)
N1_suj8 <- read.csv("DataSujetos/N1_suj8.csv",header=FALSE)
N1_suj9 <- read.csv("DataSujetos/N1_suj9.csv",header=FALSE)
N1_suj10 <- read.csv("DataSujetos/N1_suj10.csv",header=FALSE)
N1_suj11 <- read.csv("DataSujetos/N1_suj11.csv",header=FALSE)
N1_suj12 <- read.csv("DataSujetos/N1_suj12.csv",header=FALSE)
N1_suj13 <- read.csv("DataSujetos/N1_suj13.csv",header=FALSE)
N1_suj14 <- read.csv("DataSujetos/N1_suj14.csv",header=FALSE)
N1_suj15 <- read.csv("DataSujetos/N1_suj15.csv",header=FALSE)
N1_suj16 <- read.csv("DataSujetos/N1_suj16.csv",header=FALSE)
N1_suj17 <- read.csv("DataSujetos/N1_suj17.csv",header=FALSE)
N1_suj18 <- read.csv("DataSujetos/N1_suj18.csv",header=FALSE)
#############################END N1#############################



mean_maker <- function(suj1, suj2, suj3, suj4, suj5, suj6, suj7, suj8, suj9,suj10, suj11, suj12, suj13, suj14, suj15, suj16,suj17, suj18) {
  mean_result <- (suj1 + suj2 + suj3 + suj4 + suj5 + suj6 + suj7 + suj8 + suj9 + suj10 + suj11 + suj12 + suj13 + suj14 + suj15 + suj16 + suj17 + suj18)/18
  return(mean_result)
}

N1_suj1 <- as.matrix(N1_suj1)
N1_suj2 <- as.matrix(N1_suj2)


#mean_vector <- c(N1_suj1,N1_suj2)

generated_N1 <- mean_maker(N1_suj1, N1_suj2, N1_suj3, N1_suj4, N1_suj5, N1_suj6, N1_suj7, N1_suj8, N1_suj9, N1_suj10, N1_suj11, N1_suj12, N1_suj13, N1_suj14, N1_suj15, N1_suj16, N1_suj17, N1_suj18)
generated_N1 <-  round(generated_N1,digits=6)

########################FIN Generación Redes Promedio########################


N1 <- read.csv("N1promedio.csv",header=FALSE)
#N2 <- read.csv("N2promedio.csv",header=FALSE)
#N3 <- read.csv("N3promedio.csv",header=FALSE)
#W <- read.csv("Wpromedio.csv",header=FALSE)
aal <- read.csv("aal_extended.csv", header = F)
aalnames <- aal[,2] 

##
N1 <- as.matrix(N1)
N2 <- as.matrix(N2)
N3 <- as.matrix(N3)
W <- as.matrix(W)

colnames(N1) <- aalnames 
# colnames(N2) <- aalnames 
# colnames(N3) <- aalnames 
# colnames(W) <- aalnames 
# 
rownames(N1) <- aalnames 
# rownames(N2) <- aalnames 
# rownames(N3) <- aalnames 
# rownames(W) <- aalnames 

## Cij
corrplot(N1, is.corr=TRUE, title = "N1")#, order="hclust")

## Aij
# Aij = 1, if Cij > ro. ro is selected to keep delta (the link density) constant.
N = dim(N1)[1]
Nmaxlinks = N*(N-1)
##############################Delta calculation##############################
n = 1000
delta = n/Nmaxlinks
print(delta)
##############################Delta calculation##############################


diag(N1)<-0
tmp<-sort(as.vector(N1),decreasing = TRUE)
ro = tmp[n]
N1b = (N1>ro)

corrplot(N1b, is.corr=TRUE, title = "N1 thresholded", outline=FALSE)#, order="hclust")

netN1 <- graph.adjacency(N1b,mode="undirected",diag = FALSE)
V(netN1)$media <- aalnames
plot(netN1)

vcount(netN1)
ecount(netN1)

is.simple(netN1)
is.connected(netN1)

# El diametro del componente conexo
diameter(netN1) 

# Número de aristas divido número de posibles aristas
graph.density(netN1) 

mean(degree(netN1))

plot( degree.distribution(netN1, cumulative = TRUE))

d=degree(netN1)
h=hist(d,breaks = seq(min(d)-0.5,max(d)+0.5,1))
plot(log10(h$mids),log10(h$density))

netN1.plf = power.law.fit(degree(netN1))
netN1.plf$KS.p
netN1.plf$xmin
netN1.plf$alpha

# Centralidad
head( sort( degree(netN1), decreasing = T) )
head( sort( betweenness(netN1), decreasing = T) )
head( sort( closeness(netN1), decreasing = T) )
head( sort( eigen_centrality(netN1)$vector, decreasing = T) )

# Clustering
netN1.cl.eb <- cluster_edge_betweenness(netN1, directed = F, merges = T)
plot(netN1, vertex.color = netN1.cl.eb$membership)

modularity(netN1,netN1.cl.eb$membership)

source("jk.modularity.R")
N = dim(N1)[1]
Nmaxlinks = N*(N-1)
nlist = seq(100,2000,100)
dlist = array(data=NA, dim=length(nlist))
W.mlist = array(data=NA, dim=length(nlist))
N1.mlist = array(data=NA, dim=length(nlist))
N2.mlist = array(data=NA, dim=length(nlist))
N3.mlist = array(data=NA, dim=length(nlist))
k = 0
for (n in nlist) {
  k = k+1
  dlist[k] = n/Nmaxlinks
  #W.mlist[k] = jk.modularity(W,n)
  N1.mlist[k] = jk.modularity(N1,n)
  #N2.mlist[k] = jk.modularity(N2,n)
  #N3.mlist[k] = jk.modularity(N3,n)
}

df <- data.frame(dlist,W.mlist,N1.mlist,N2.mlist,N3.mlist)

library(ggplot2)

ggplot(df, aes(dlist)) +                    # basic graphical object
  geom_line(aes(y=W.mlist), colour="black") +  # first layer
  geom_line(aes(y=N1.mlist), colour="green")  + # second layer
  geom_line(aes(y=N2.mlist), colour="blue")  + # second layer
  geom_line(aes(y=N3.mlist), colour="red")  # second layer



########################Taller Práctico######################################
hist(N1[lower.tri(N1)], main = "Histograma Relaciones N1")

#Conteo de Nodos y aristas
vcount(netN1)
ecount(netN1)

#Nodos y aristas
V(netN1)
E(netN1)

#¿Loops o aristas múltiples?
is.simple(netN1)

#¿Completamente conectado?
is.connected(netN1)


##########################Red Con umbral##########################
#umbral <-0.7377# Umbral N1 por defecto
umbral <-0.5
N1_umbral <-N1
N1_umbral[which(N1_umbral <= umbral)] <- 0
netN1_umbral <- graph.adjacency(N1_umbral, mode="undirected", diag=FALSE, weighted = T)
plot(netN1_umbral)
#Conteo de Nodos y Aristas con umbral
vcount(netN1_umbral)
ecount(netN1_umbral)
##########################Red Con umbral##########################

#Diametro Red
diameter(netN1)
get.diameter(netN1)
diameter(netN1_umbral)
get.diameter(netN1_umbral)


#Densidad del grafo, relacion entre aristas y nodos
graph.density(netN1)
graph.density(netN1_umbral)

#Coeficiente de clustering local 
head(transitivity(netN1, type = "local"))
head(transitivity(netN1_umbral, type = "local"))

#Coeficiente de clustering global
transitivity(netN1, type = "global")
transitivity(netN1_umbral, type = "global")

#########################Comparación coeficintes de clustering#########################
par(mfrow = c(1,2))
hist(transitivity(netN1, type = "local"), 
     main = "N1", xlab = "coefs. de clustering")
hist(transitivity(netN1_umbral, type = "local"), 
     main = "N1 Umbral", xlab = "coefs. de clustering")
par(mfrow = c(1,1))
######################### FIN Comparación coeficintes de clustering#########################

#Grados de entrada y salida
degree(netN1)
degree(netN1_umbral)

sort(degree(netN1), decreasing = T)
sort(degree(netN1_umbral), decreasing = T)

#Empezamos a analizar que sucede a nivel de grados
qplot(degree(netN1), degree(netN1_umbral))

#Correlación entre el grado que predice el grupo para N1 y N1 con umbral
cor(degree(netN1), degree(netN1_umbral))


#distribuciones de grados y diustribuciones acumuladas
head(degree.distribution(netN1), 15)
head(degree.distribution(netN1_umbral), 15)

head(degree.distribution(netN1, cumulative = T))
head(degree.distribution(netN1_umbral, cumulative = T))

#########################Proporción de Nodos#########################
par(mfrow = c(1,2))
plot(degree.distribution(netN1),
     xlab = "Grados", ylab = "Proporción de nodos", type = "h", main ="N1")
plot(degree.distribution(netN1_umbral),
     xlab = "Grados", ylab = "Proporción de nodos", type = "h", main ="N1_Umbral")
par(mfrow = c(1,1))
######################### FIN Proporción de Nodos#########################

#Cálculo de asortividad -> en ambos casos los valores sugieren que no hay asociaciones preferenciales
#entre nodos de un alto grado y por otrolado de bajo grado
assortativity.degree(netN1)
assortativity.degree(netN1_umbral)

#Medidas de centralidad
#Intermedicion
head(sort(betweenness(netN1), decreasing = T))
head(sort(betweenness(netN1_umbral), decreasing = T))

#Cercania
head(sort(closeness(netN1), decreasing = T))
head(sort(closeness(netN1_umbral), decreasing = T))

#Centralidad de autovectores
head(sort(eigen_centrality(netN1)$vector, decreasing = T))
head(sort(eigen_centrality(netN1_umbral)$vector, decreasing = T))


###############################Clustering EB#########################################
net1.cl.eb <- cluster_edge_betweenness(netN1, directed = F, merges = T) 
net1.umbral.cl.eb <- cluster_edge_betweenness(netN1_umbral, directed = F, merges = T) 

par(mfrow = c(1,2))
plot(netN1, vertex.color = net1.cl.eb$membership)
plot(netN1_umbral, vertex.color = net1.umbral.cl.eb$membership)
par(mfrow = c(1,1))

par(mfrow = c(1,2))
plot(netN1, vertex.color = net1.cl.eb$membership, vertex.label=NA, vertex.size=5)
plot(netN1_umbral, vertex.color = net1.umbral.cl.eb$membership, vertex.label=NA, vertex.size=5)
par(mfrow = c(1,1))


############################### END Clustering EB#########################################


###############################Louvain#########################################
net1.cl.lo <- cluster_louvain(netN1, weights = E(netN1)$weight)
net1.umbral.cl.lo <- cluster_louvain(netN1_umbral, weights = E(netN1_umbral)$weight)

par(mfrow = c(1,2))
plot(netN1, vertex.color = netN1$membership)
plot(netN1_umbral, vertex.color = netN1_umbral$membership)
par(mfrow = c(1,1))
############################### END Louvain#########################################

#Membership comparison
table(net1.cl.lo$membership, net1.cl.lo$membership)
plot(ws.obs.red, vertex.color = net1.cl.lo$membership)
#Plot for merbenship with circles and squares
plot(ws.obs.red, vertex.color = ws.obs.red.cl.lo$membership, vertex.shape = ifelse
     (ws.per.red.cl.lo$membership == 1, "circle", "square") )
##############################################################



