##Análisis de correlaciones
library(dplyr)
library(tidyverse)
library("FactoMineR")
library("factoextra")

political_moral16ct_db <- read.csv2('political_moral_16cat.csv', header=TRUE)
View(political_moral16ct_db)

political_moral16ct_db %>% group_by(political_, moral_circ) %>% summarise(count=n())


pol_moral16cat <- political_moral16ct_db %>% gather(Key, political_, -moral_circ) %>% 
    group_by(moral_circ, political_) %>%
    summarise(count = n())  %>%
    spread(political_, count, fill = 0) %>%
    as.data.frame()

view(pol_moral16cat)

nrows_16cat <- pol_moral16cat[, 2:4]

rownames(nrows_16cat) <- db_pol_moral16cat$moral_circ

view(nrows_16cat) ##base de datos de trabajo. 



pol_moral16matrix <- as.table(as.matrix(nrows_16cat))
class(pol_moral16matrix)

###MARGINALES###

x <- addmargins(pol_moral16matrix, c(1,2))

round(prop.table(pol_moral16matrix,2)*100, 2)


library("gplots")

balloonplot(t(pol_moral16cat_db), main ="Circulo moral por disposición política", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)


## Prueba de independencia
chisq_pol_moral <- chisq.test(nrows_16cat)
chisq_pol_moral

## An?lisis de Correspondencias
res.ca_16cat_moral <- CA(nrows_16cat, graph = FALSE)
res.ca_16cat_moral

eig.val_16cat_moral <- get_eigenvalue(res.ca_16cat_moral)
round(eig.val_16cat_moral ,2)


res.ca_16cat_moral$row$contrib

## Grafico de sedimentaci?n
fviz_screeplot(res.ca_16cat_moral) +
    geom_hline(yintercept=33.33, linetype=2, color="red")

## Biplot
windows()
fviz_ca_biplot(res.ca_16cat_moral, 
               repel = TRUE, 
               col.col = "contrib", 
               alpha.row ="contrib",
               labelsize = 3,
               pointsize = 4) +
    scale_color_gradient2(low = "white", mid = "blue",
                          high = "red", midpoint = 25) +
    theme_minimal()

## Biplot asim?trico
windows()
fviz_ca_biplot(res.ca_16cat_moral,
               map="rowprincipal", arrow=c(TRUE,TRUE),
               repel=TRUE)
# Si el ?ngulo entre las dos flechas es agudo,
# entonces hay una asociaci?n fuerte entre las 
# filas y columnas correspondientes

##bipplot por contribuciones
fviz_ca_biplot(res.ca_16cat_moral,
               map="rowprincipal", 
               col.row = "contrib", 
               col.col = "contrib",
               arrow=c(TRUE,TRUE),
               labelsize = 3,
               repel=TRUE) +
    scale_color_gradient2(low="white", mid="blue",
                          high="red", midpoint=10)

fviz_pca_biplot(multi1_PCA, 
                # Individuals
                geom.ind = "point",
                fill.ind = E.coli_K.pneu_db$Region, col.ind = "white",
                pointshape = 21, 
                pointsize = 2,
                palette = "jco",
                labelsize = 3,
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", 
                col.var = "contrib",
                gradient.cols = "RdBu") + 
    labs(fill = "Bacteria", 
         color = "Contribución", 
         alpha = "Contribución") +
    theme(legend.title = 
              element_text(face = "bold")) 



multi_corr_surv <- read.csv2('political_moral_otras_vars.csv', header=TRUE)
view(multi_corr_surv)
rownames(multi_corr_surv) <- c(1:nrow(multi_corr_surv))

multi_corr_surv2 <- multi_corr_surv[30:90, ]

## MCA An?lisis de Correspondencias M?ltiples
res.mca_multisurv <- MCA(multi_corr_surv2, graph = FALSE)

## Visualizaci?n
eig.val <- get_eigenvalue(res.mca_multisurv)
fviz_screeplot(res.mca_multisurv, addlabels=TRUE, ylim=c(0,45))

## Biplot
windows()
fviz_mca_biplot(res.mca_multisurv,, 
                ggtheme=theme_minimal())



####tablas de columnas y filas####

nam_x <- c('Variable', 'Descripcion')
y <- c('toda su familia inmediata',  
       'toda su familia extendida',
       'todos sus amigos cercanos',
       'todos sus amigos, incluyendo los distantes',
       'todos sus conocidos',
       'todas las personas con las que se ha cruzado',
       'todas las personas de su país' ,
       'todas las personas en su continente',
       'todas las personas en todos los continentes',
       'todos los mamíferos',
       'todos los anfibios, reptiles, mamíferos, peces y aves',
       'todos los animales en la tierra, incluyendo bacterias y amibas',
       'todos los animales en el universo, incluyendo formas de vida extraterrestres',
       'todas las cosas vivientes en el universo, incluyendo plantas y arboles',
       'todas las cosas naturales en el universo, incluyendo entidades inertes como las rocas',
       'todas las cosas que existen')

x <- c( "familiar_inmediato", "familiar_extendido",
        "amigos_cercanos", "amigos_inclu_lejano", "todos_conocidos",
        "personas_cruzado", "personas_pais", "personas_continente",
        "personas_mundo", "mamiferos", "animal_cordado", 
        "animal_bact",  "animal_extrat", "todo_biotico", "biotico_abiotico", 
        "todo_existente")
z <- c(1:16)

filas <- data.frame(z, x, y)

colnames(filas) <- nam_x

view(filas)



####df columnas####

w <- c('Variable', 'Descripcion')

t <- c('derecha', 'izquierda', 'centro')

r <- c('Todos aquellos que marcaron entre 7 y 10', 
       'Todos aquellos que marcaron entre 0 y 3', 
       'Todos aquellos que marcaron entre 4 y 6')
k <- c(1:3)

columnas <- data.frame(k, t, r)

colnames(columnas) <- w

view(columnas)
