#########################################
####Taller_4 (análisis multivariado)#####
#### Análisis de Conglomerados ##########
#########################################
library(tidyverse)
require("cluster")
require("factoextra")
require("dendextend")
require("magrittr")
require("readxl")

E.coli_K.pneu_db <- read.csv2('E.coli_&_K.pneu_db', header=TRUE)

names(E.coli_K.pneu_db)

E.coli_K.pneu_db %>% group_by(Bacteria) %>% count()

### Creación de un df que contiene cada una de las medias de resistencia a los dos antibioóticos
### para los dos organismos.
mean_AB_bact <- E.coli_K.pneu_db %>% 
    group_by(Country) %>%
    summarise(Aminoglycosides = mean(Aminoglycosides), 
              Carbapenems = mean(Carbapenems), 
              Fluoroquinolones = mean(Fluoroquinolones), 
              Cephalos_3er_gen = mean(cephalos_3er_gen))

### df para media de resistencia para cada la bacteria Klebsiella Pneumoniae
mean_AB_E.coli <- E.coli_K.pneu_db %>% 
    filter(Bacteria == 'Escherichia coli') %>%
    group_by(Country) %>%
    summarise(Aminoglycosides = mean(Aminoglycosides), 
              Carbapenems = mean(Carbapenems), 
              Fluoroquinolones = mean(Fluoroquinolones), 
              Cephalos_3er_gen = mean(cephalos_3er_gen))

mean_AB_K.pneum <- E.coli_K.pneu_db %>% 
    filter(Bacteria == 'Klebsiella pneumoniae') %>%
    group_by(Country) %>%
    summarise(Aminoglycosides = mean(Aminoglycosides), 
              Carbapenems = mean(Carbapenems), 
              Fluoroquinolones = mean(Fluoroquinolones), 
              Cephalos_3er_gen = mean(cephalos_3er_gen))
    
view(multi_var4_E.coli)

dim(mean_AB_E.coli)
    
multi_var4_E.coli <- as_tibble(mean_AB_E.coli[, -1])

rownames(multi_var4_E.coli) <- mean_AB_E.coli$Country

View(multi_var4_E.coli)

#####################################
# Calculamos la matriz de distancias#
#####################################

dist_db_E.coli <- get_dist(multi_var4_E.coli, stand = TRUE, method = "euclidean") 
round(dist_db_E.coli,2)

class(dist_db_E.coli)
?get_dist

as.matrix(dist_db_E.coli)[1:5,1:5]


fviz_dist(dist_db_E.coli, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


library(gplots) 

# Probamos con algunos dendogramas
hc1_e.coli = hclust(dist_db_E.coli, method = "single") 
plot(hc1_e.coli, hang = -1)
abline(h=1.5,col="red",lty=2)


hc2_e.coli = hclust(dist_db_E.coli, method = "complete") 
plot(hc2_e.coli, hang = -1)
abline(h=1.5,col="red",lty=2)


hc3_e.coli = hclust(dist_db_E.coli, method = "average")
plot(hc3_e.coli, hang = -1)
abline(h=1.5,col="red",lty=2)

### M?todo de Ward

hc4_e.coli = hclust(dist_db_E.coli, method = "ward.D2")
plot(hc4_e.coli, hang = -1)
abline(h=1.5,col="red",lty=2)



# Mejorando los dendogramas
# Dendograma mejorado para average
windows()
fviz_dend(hc3_e.coli, k = 4, # Cortar en 3 grupos
          cex = 0.5, # tama?o de las etiquetas
          k_colors = "npg", # color de los grupos
          color_labels_by_k = TRUE, # color de las etiquetas de los grupos
          rect = TRUE # Agregar un rect?ngulo a los grupos
)

# Dendograma mejorado para ward.D2
fviz_dend(hc4_e.coli, k = 4, # Cortar en 3 grupos
          cex = 0.5, # tama?o de las etiquetas
          k_colors = "npg", # color de los grupos
          color_labels_by_k = TRUE, # color de las etiquetas de los grupos
          rect = TRUE # Agregar un rect?ngulo a los grupos
)




# N?mero ?ptimo de grupos (gr?fico)
# M?todo de Mojena
mojena = function(hc){
    n_hd = length(hc$height)
    alp_g = 0 ; alpha = hc$height[n_hd:1]
    for(i in 1:(n_hd-1)){
        alp_g[i] = mean(alpha[(n_hd-i+1):1])+1.25*sd(alpha[(n_hd-i+1):1])
        # alp_g[i] = mean(alpha[(n_hd-i+1):1])+3.5*sd(alpha[(n_hd-i+1):1])
    }
    nog = sum(alp_g<= alpha[-n_hd]) + 1
    plot(alpha[-n_hd], pch=20, col=(alp_g>alpha[-n_hd])+1, main = paste("Número óptimo de grupos =",nog),
         ylab = expression(alpha[g]), xlab="Nodos")}

windows()
mojena(hc3_e.coli) ## Calcula el n?mero de conglomerados de m?todos jerarquicos



df_e.coli <- scale(multi_var4_E.coli)
head(df_e.coli, 3)


## Estimaci?n del n?mero ?ptimo de clusters
library(factoextra)
windows()

fviz_nbclust(df_e.coli,kmeans,method = "wss")+
    geom_vline(xintercept = 4, linetype=2) ## Advertencia

set.seed(123) ## Semilla
km.res_e.coli <- kmeans(df_e.coli, 4,nstart = 25)

print(km.res_e.coli)

## Ver la clasificaci?n
km.res_e.coli$cluster

## N?mero de individuos por conglomerados
km.res_e.coli$size
dd_e.coli <- cbind(multi_var4_E.coli, 
                   cluster=km.res_e.coli$cluster) ## Agregar la clasificaci?n
head(dd_e.coli)

## Visualizaci?n conglomerados K-medias
windows()
fviz_cluster(km.res_e.coli, 
             data = df_e.coli, 
             palette="jco", 
             star.plot= TRUE,
             pointsize = 2,
             repel = TRUE,
             ggtheme=theme_minimal()
)


## Ayuda para la interpretaci?n
aggregate(multi_var4_E.coli, 
          by=list(cluster=km.res_e.coli$cluster),mean)
