#########################################
####Taller_4 (análisis multivariado)#####
#### Análisis de Conglomerados ##########
#########################################


###EXTRA CON OTRA BACTERIA###

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

view(mean_AB_K.pneum)



multi_var4_K.pneum <- as_tibble(mean_AB_K.pneum[, -1])

rownames(multi_var4_K.pneum) <- mean_AB_K.pneum$Country

View(multi_var4_K.pneum)

# Calculamos la matriz de distancias

dist_db_K.pneum <- get_dist(multi_var4_K.pneum, stand = TRUE, method = "euclidean") 
round(dist_db_K.pneum,2)

class(dist_db_K.pneum)
?get_dist

as.matrix(dist_db_K.pneum)[1:5,1:5]

windows()
fviz_dist(dist_db_K.pneum, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Probamos con algunos dendogramas
hc1_K.pneum = hclust(dist_db_K.pneum, method = "single") 
windows()
plot(hc1_K.pneum, hang = -1)
abline(h=1.5,col="red",lty=2)

windows()
hc2_K.pneum = hclust(dist_db_K.pneum, method = "complete") 
plot(hc2_K.pneum, hang = -1)
abline(h=1.5,col="red",lty=2)

windows()
hc3_K.pneum = hclust(dist_db_K.pneum, method = "average")
plot(hc3_K.pneum, hang = -1)
abline(h=1.5,col="red",lty=2)

### M?todo de Ward
windows()
hc4_K.pneum = hclust(dist_db_K.pneum, method = "ward.D2")
plot(hc4_K.pneum, hang = -1)
abline(h=1.5,col="red",lty=2)



# Mejorando los dendogramas
# Dendograma mejorado para average
windows()
fviz_dend(hc3_K.pneum, k = 3, # Cortar en 3 grupos
          cex = 0.5, # tama?o de las etiquetas
          k_colors = "npg", # color de los grupos
          color_labels_by_k = TRUE, # color de las etiquetas de los grupos
          rect = TRUE # Agregar un rect?ngulo a los grupos
)

# Dendograma mejorado para ward.D2
fviz_dend(hc4_K.pneum, k = 4, # Cortar en 3 grupos
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
mojena(hc3_K.pneum) ## Calcula el n?mero de conglomerados de m?todos jerarquicos



df_K.pneum <- scale(multi_var4_K.pneum)
head(df_K.pneum, 3)


## Estimaci?n del n?mero ?ptimo de clusters
library(factoextra)
windows()
fviz_nbclust(df_K.pneum,kmeans,method = "wss")+
    geom_vline(xintercept = 3, linetype=2) ## Advertencia

set.seed(123) ## Semilla
km.res_K.pneum <- kmeans(df_K.pneum, 3,nstart = 25)

print(km.res_K.pneum)

## Ver la clasificaci?n
km.res_K.pneum$cluster

## N?mero de individuos por conglomerados
km.res_K.pneum$size
dd_K.pneum <- cbind(multi_var4_K.pneum, 
                   cluster=km.res_K.pneum$cluster) ## Agregar la clasificaci?n
head(dd_K.pneum)

## Visualizaci?n conglomerados K-medias
windows()
fviz_cluster(km.res_K.pneum, 
             data = df_K.pneum, 
             palette="Dark2", 
             star.plot= TRUE, 
             repel = TRUE, 
             pointsize = 2,
             pointshape = 21,
             ggtheme=theme_minimal()
)


## Ayuda para la interpretaci?n
aggregate(multi_var4_K.pneum, 
          by=list(cluster=km.res_K.pneum$cluster),mean)
