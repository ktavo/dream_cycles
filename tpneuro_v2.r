# https://github.com/cwatson/brainGraph
# http://kateto.net/network-visualization

rm(list=ls())
setwd("E:/UBA/2018-II/DM en Ciencia y Tecnología/Ciclos Sueño")

## Load package
library(igraph)
# http://igraph.org/r/doc/

library(corrplot)


W <- read.csv("Wpromedio.csv",header=FALSE)
N1 <- read.csv("N1promedio.csv",header=FALSE)
N2 <- read.csv("N2promedio.csv",header=FALSE)
N3 <- read.csv("N3promedio.csv",header=FALSE)

#Adjusting W, N2 and N3 as N1 for 116x116
W <- W[-1, ]
W[1] <- NULL
N2 <- N2[-1, ]
N2[1] <- NULL
N3 <- N3[-1, ]
N3[1] <- NULL


aal <- read.csv("aal_extended.csv", header = F)
aalnames <- aal[,2] 
aalModule <- aal[,4] 


##
N1 <- as.matrix(N1)
N2 <- as.matrix(N2)
N3 <- as.matrix(N3)
W <- as.matrix(W)

class(W) <- "numeric"
class(N2) <- "numeric"
class(N3) <- "numeric"


colnames(N1) <- aalnames 
colnames(N2) <- aalnames 
colnames(N3) <- aalnames 
colnames(W) <- aalnames 
# 
rownames(N1) <- aalnames 
rownames(N2) <- aalnames 
rownames(N3) <- aalnames 
rownames(W) <- aalnames 


colnames(N1) <- aalModule 
colnames(N2) <- aalModule 
colnames(N3) <- aalModule 
colnames(W) <- aalModule 
# 
rownames(N1) <- aalModule 
rownames(N2) <- aalModule 
rownames(N3) <- aalModule 
rownames(W) <- aalModule 



############################## Delta calculation##############################
N = dim(N1)[1]
Nmaxlinks = N*(N-1)
#n = 7000
#delta = n/Nmaxlinks
#print(delta)

delta = 0.038
n = Nmaxlinks*delta
############################## END Delta calculation##############################

############################# Transformación en grafo no pesado ############################# 
diag(W)<-0
tmp.W<-sort(as.vector(W),decreasing = TRUE)
ro.W = tmp.W[n]
Wb.W = (W>ro.W)

diag(N1)<-0
tmp.N1<-sort(as.vector(N1),decreasing = TRUE)
ro.N1 = tmp.N1[n]
N1b = (N1>ro.N1)

diag(N2)<-0
tmp.N2<-sort(as.vector(N2),decreasing = TRUE)
ro.N2 = tmp.N2[n]
N2b = (N2>ro.N2)

diag(N3)<-0
tmp.N3<-sort(as.vector(N3),decreasing = TRUE)
ro.N3 = tmp.N3[n]
N3b = (N3>ro.N3)

############################# END Transformación en grafo no pesado ############################# 


############################# Gráfica Nodos Aristas ############################# 

netW <- graph.adjacency(Wb.W,mode="undirected",diag = FALSE, weighted = T)
V(netW)$media <- aalModule
V(netW)$color <- aalModule
plot(netW, vertex.label.color = "black", vertex.size=5, main = "W")


netN1 <- graph.adjacency(N1b,mode="undirected",diag = FALSE, weighted = T)
V(netN1)$media <- aalModule
V(netN1)$color <- aalModule
plot(netN1, vertex.label.color = "black", vertex.size=5, main = "N1")


netN2 <- graph.adjacency(N2b,mode="undirected",diag = FALSE, weighted = T)
V(netN2)$media <- aalModule
V(netN2)$color <- aalModule
plot(netN2, vertex.label.color = "black", vertex.size=5, main = "N2")


netN3 <- graph.adjacency(N3b,mode="undirected",diag = FALSE, weighted = T)
V(netN3)$media <- aalModule
V(netN3)$color <- aalModule
plot(netN3, vertex.label.color = "black", vertex.size=5, main = "N3")

############################# END Gráfica Nodos Aristas ############################# 


####################################  Numero de nodos y aristas ####################################
vcount(netW)
ecount(netW)
#vcount(netN1)
ecount(netN1)
#vcount(netN2)
ecount(netN2)
#vcount(netN3)
ecount(netN3)
####################################  END Numero de nodos y aristas ####################################

#################################### Diametro ####################################

# El diametro del componente conexo
diameter(netW) 
diameter(netN1) 
diameter(netN2) 
diameter(netN3) 

#################################### END Diametro ####################################


######################## Taller Práctico ######################################
par(mfrow = c(1,2))
hist(W[lower.tri(W)], main = "Histograma Relaciones W")
hist(N1[lower.tri(N1)], main = "Histograma Relaciones N1")
par(mfrow = c(1,1))

par(mfrow = c(1,2))
hist(N2[lower.tri(N2)], main = "Histograma Relaciones N2")
hist(N3[lower.tri(N3)], main = "Histograma Relaciones N3")
par(mfrow = c(1,1))

#Conteo de Nodos y aristas
vcount(netW)
ecount(netW)

#Nodos y aristas
V(netN1)
E(netN1)

######################### Comparación coeficintes de clustering#########################
par(mfrow = c(1,2))
hist(transitivity(netW, type = "local"), 
     main = "W", xlab = "coefs. de clustering")
hist(transitivity(netN1, type = "local"), 
     main = "N1", xlab = "coefs. de clustering")
par(mfrow = c(1,1))


par(mfrow = c(1,2))
hist(transitivity(netN2, type = "local"), 
     main = "N2", xlab = "coefs. de clustering")
hist(transitivity(netN3, type = "local"), 
     main = "N3", xlab = "coefs. de clustering")
par(mfrow = c(1,1))
######################### FIN Comparación coeficintes de clustering#########################



###############################Clustering EB#########################################

###########Punto 1###########
####Cluster edge betweenness es Girvan-Newman#########
netw.cl.eb <- cluster_edge_betweenness(netW, directed = F, merges = F) 
net1.cl.eb <- cluster_edge_betweenness(netN1, directed = F, merges = F) 
net2.cl.eb <- cluster_edge_betweenness(netN2, directed = F, merges = F) 
net3.cl.eb <- cluster_edge_betweenness(netN3, directed = F, merges = F) 

#V(netN2)$media <- aalnames


#Comparativa entre W y N1 Corregidas
par(mfrow = c(1,2))
plot(netW, vertex.color = netw.cl.eb$membership, vertex.label=aalModule, vertex.size=3, main = "W")
plot(netN1, vertex.color = net1.cl.eb$membership, vertex.label=aalModule, vertex.size=3, main = "N1")
par(mfrow = c(1,1))

#Comparativa entre W y N2 Corregida
par(mfrow = c(1,2))
plot(netN2, vertex.color = net2.cl.eb$membership, vertex.label=aalModule, vertex.size=5, main = "W")
plot(netN3, vertex.color = net3.cl.eb$membership, vertex.label=aalModule, vertex.size=5, main = "N2")
par(mfrow = c(1,1))

#Comparativa entre W y N3 Corregida
par(mfrow = c(1,2))
plot(netN2, vertex.color = net2.cl.eb$membership, vertex.label=aalModule, vertex.size=1, main = "W")
plot(netN3, vertex.color = net3.cl.eb$membership, vertex.label=aalModule, vertex.size=1, main = "N3")
par(mfrow = c(1,1))
##################################Comparativa 4 estados##################################
par(mfrow = c(2,2))
plot(netW, vertex.color = netw.cl.eb$membership, vertex.label=NA, vertex.size=5, main = "W")
plot(netN1, vertex.color = net1.cl.eb$membership, vertex.label=NA, vertex.size=5, main = "N1")
plot(netN2, vertex.color = net2.cl.eb$membership, vertex.label=NA, vertex.size=5, main = "N2")
plot(netN3, vertex.color = net3.cl.eb$membership, vertex.label=NA, vertex.size=5, main = "N3")
par(mfrow = c(1,1))

##################################END Comparativa 4 estados##################################

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



