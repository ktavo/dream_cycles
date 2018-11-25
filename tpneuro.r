## Load package
library(igraph)
# http://igraph.org/r/doc/

library(corrplot)

# https://github.com/cwatson/brainGraph
# http://kateto.net/network-visualization

rm(list=ls())
setwd("E:/UBA/2018-II/DM en Ciencia y Tecnología/Ciclos Sueño")

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

diameter(netN1) # El diametro del componente conexo

graph.density(netN1) # Número de aristas divido número de posibles aristas

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



##############################################################
hist(N1[lower.tri(N1)], main = "Histograma Interacciones Percibidas")


#Conteo de Nodos y aristas
vcount(netN1)
ecount(netN1)

##########################Red Con umbral##########################
#umbral <-0.7377# Umbral N1 por defecto
umbral <-0.7
N1_umbral <-N1
N1_umbral[which(N1_umbral <= umbral)] <- 0
netN1_umbral <- graph.adjacency(N1_umbral, mode="undirected", diag=FALSE, weighted = T)
plot(netN1_umbral)
#Conteo de Nodos y Aristas con umbral
vcount(netN1_umbral)
ecount(netN1_umbral)
##########################Red Con umbral##########################

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



net1.cl.eb <- cluster_edge_betweenness(netN1, directed = F, merges = T) 
net1.umbral.cl.eb <- cluster_edge_betweenness(netN1_umbral, directed = F, merges = T) 

plot(netN1, vertex.color = net1.cl.eb$membership)
plot(netN1_umbral, vertex.color = net1.umbral.cl.eb$membership)




#Louvain
net1.cl.lo <- cluster_louvain(netN1, weights = E(netN1)$weight)
plot(netN1, vertex.color = netN1$membership)
##############################################################



