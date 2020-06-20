#limpiando y arreglando datos antibióticos 
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)

ab_bact_db<- read.table('ECDC_surveillance_data_Antimicrobial_resistance.csv', 
                        header=TRUE, sep = ',', stringsAsFactors = FALSE)
view(ab_bact_db)
str(ab_bact_db)
class(ab_bact_db)


options(digits=10)

#Modificación de la base de datos a Tibble, eliminación de columnas no necesarias. 
#creación de nueva columnas conocidas como Bacteria y Antibiotic, eliminación de simbolos
#por Na o espacios. 
ab_bact_db <- ab_bact_db %>% as_tibble(rm.na=TRUE) %>%
    select(-c('HealthTopic', 'Unit', 'TxtValue')) %>%
    separate('Population', c('Bacteria', 'Antibiotic'), sep = '\\|') %>%
    mutate(NumValue = str_replace_all(NumValue, "-", ""))

# Convertir un columna en númerica_ en este caso NumValue. 
ab_bact_db <- ab_bact_db %>% mutate(NumValue = as.numeric(NumValue)) #%>%
    #filter(is.na(NumValue)) %>% 
    #head(n=10)

#######################
###REVISIÓN DE LA BD###
#######################

names(ab_bact_db)
head(ab_bact_db)


### conteo de categorías en cada  columna y creación de columna 'n' para definir el
#número de elementos por columna o variable. 
ab_bact_db %>% count(Indicator, sort = TRUE, name = 'n')
ab_bact_db %>% count(Antibiotic, sort = TRUE, name = 'n')
ab_bact_db %>% count(Bacteria, sort = TRUE, name = 'n')

### Filtro y conteo de elementos para la varible filtrada. 
ab_bact_db %>% filter(Bacteria == 'Enterococcus faecium') %>%
    count(RegionName, sort = TRUE, name = 'n')

### creación de Función para hacer más eficiente los procesos anteriores. 
n_vars <- function(db, var1, element_var1, var2) {
    require(dplyr) #requiere tener la librería dplyr o tidiverse
    #Se filtra por la primera variable (va1) y el elemento introducido y 
    #luego se cuenta por la variable 2 (var2).
    x <- db %>% filter({{var1}} == element_var1) %>% 
        count({{var2}}, sort = TRUE, name = 'n')
    print(x)
}

### Comprobación de la función anterior. 
n_vars(ab_bact_db, Bacteria, 'Klebsiella pneumoniae' , Indicator)

###Función para obtener los NAs y los valores con simbolo determinado al convertir una
###variable en numerica
not_info<- function (db, s){
    values <- suppressWarnings(as.numeric(x))
    values <- is.na(values) | values == s
    values
}


#Aplicación de la función anterior a una columna determinada, se le adiciona el sum() 
#para contar el número de valores que son NAs y/o son el simbolo o patrón buscado. 
not_info(ab_bact_db$NumValue, '-') %>% sum()

#mutate(NumValue = as.numeric(NumValue)) es otra manera de convertir una columna en númerica. 

which(is.na(ab_bact_db$NumValue))

#########################################
###REORGANIZACIÓN DE LA COLUMNAS EN DB###
#########################################

### reorganizando el nombre de los indicadores (INDICATOR)

###Intento de separar indicador por el '-' (no funciona en este caso )
Ind <- ab_bact_db %>% separate(Indicator, c('Letter', 'explain'), sep = '-', 
                               fill = 'right')

###Código para organizar la base de datos por Indicador y cambio de nombre de columnas
### SE DEFINE UNA NUEVA BASE DE DATOS LLAMADA tidy_by_indic (el % se pasa a decimal)
tidy_by_indic <- ab_bact_db %>% spread(Indicator, NumValue) %>% 
    rename("r_isolates" = "R - resistant isolates", 
           "r_percentage" = "R - resistant isolates, percentage  ", 
           "total_isolates" = "Total tested isolates", 
           "i_isolates" = "I - susceptible, increased exposure isolates", 
           "s_isolates" = "S - susceptible isolates", 
           "compl_age" = "Completeness age", 
           "compl_gender" = "Completeness gender", 
           "penicilin_isolate_percentage" = 
               "Penicillin non-wild-type isolates, percentage") %>%
    mutate(r_percentage = r_percentage/100, 
           penicilin_isolate_percentage = penicilin_isolate_percentage/100, 
           compl_age = compl_age/100, 
           compl_gender = compl_gender/100, 
           s_isolates = total_isolates - (r_isolates+i_isolates)) %>%
    select(-c("penicilin_isolate_percentage", 
              "compl_age", 
              "compl_gender"), 
           everything())

view(tidy_by_indic)

##########################
###EXPLORACIÓN DE LA DB###
##########################

tidy_by_indic %>% 
    filter(Bacteria == 'Enterococcus faecium') %>%
    arrange(RegionName) 

### Aplicación función n_vars para identificar bacterias con resistencia >50%. 
n_vars(tidy_by_indic, Bacteria, 'Escherichia coli', r_percentage) %>% 
    arrange(desc(r_percentage)) %>% filter(r_percentage >= 0.5)

tidy_by_indic %>% count(Bacteria)

### Revisión de cuales son los NO NAs en el % de aislados con resistencia a la penicilina no WT.
tidy_by_indic %>% filter(!is.na(penicilin_isolate_percentage)) %>% view()

### Revisión de resistencia de E. Coli > 50% de los aislados por Region
tidy_by_indic %>% 
    filter(Bacteria == 'Escherichia coli', r_percentage >= 0.5) %>% 
    arrange(RegionName) %>% count(RegionName) %>% view()

### Obtener media de porcentaje de resistencia por cada grupo de bacterias. 
tidy_by_indic %>% group_by(Antibiotic) %>% filter(Bacteria == 'Enterococcus faecium') %>%
    summarise(mean = mean(r_percentage, na.rm = TRUE))

### Quitar los NA del porcentaje de resistencia en E.coli.
tidy_by_indic %>% filter(Bacteria == 'Acinetobacter spp.') %>%
    drop_na(r_percentage) %>% view()

### Quitar los valores = a 0 en el total de los aislados. 
tidy_by_indic %>% filter(Bacteria == 'Escherichia coli', 
                         !total_isolates == 0) %>% view()

tidy_by_indic %>% gather(Antibiotic, key,i_isolates:penicilin_isolate_percentage)

tidy_by_indic %>% group_by(Antibiotic) %>% count()


### ejercicio para detectar una cadena de "combined resistance" en las filas. 
y <-tidy_by_indic %>% 
    filter(str_detect(Antibiotic, 'Combined resistance')) %>% 
    mutate(Antibiotic = recode(Antibiotic, 
                               'Combined resistance (fluoroquinolones, aminoglycosides and carbapenems)' =
                                   'cr_Fluor;Amino;Carba', 
                               'Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)' =
                                   'cr_thr;cepha;fluor;amino', 
                               'Combined resistance (at least three of piperac. and tazob., fluoroq., ceftaz., aminogl. and carbapenems)' =
                                   'cr_atlleast_three'))

################################################################################################
###Organizando la BD por Indicador y antibiótico ####BASE DE DATOS ORGANIZADA PARA TRABAJO######
################################################################################################

percent <- function(x) {x/100} # formula para convertir % en frecuencia

# código para limpiar la bd, cambiar nombres de los elementos en las filas para que sean más
# cortos y poder organizarlos mejor; se une el indicador con el AB y luego se pasan las filas a 
# columnas 'spread' de la unión entre AB e indicador; finalmente se limpia la base a través del
# del control de calidad de 'data quality_completeness age.

tidy_by_indic2 <- ab_bact_db %>%
    mutate(Indicator = recode (Indicator, "R - resistant isolates" = "r_isolates", 
           "R - resistant isolates, percentage  " = "r_percentage", 
           "Total tested isolates" = "total_isolates", 
           "I - susceptible, increased exposure isolates" = "i_isolates" , 
           "S - susceptible isolates" = "s_isolates", 
           "Completeness age" = "compl_age", 
           "Completeness gender" = "compl_gender", 
           "Penicillin non-wild-type isolates, percentage" =
               "penicilin_isolate_percentage")) %>% 
    mutate(Antibiotic = recode(Antibiotic, 
                               'Combined resistance (fluoroquinolones, aminoglycosides and carbapenems)' =
                                   'CR_Fluor;Amino;Carba', 
                               'Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)' =
                                   'CR_thr;cepha;fluor;amino', 
                               'Combined resistance (at least three of piperac. and tazob., fluoroq., ceftaz., aminogl. and carbapenems)' =
                                   'CR_atlleast_three')) %>%
    unite(AB_indicator, Antibiotic, Indicator) %>%
    spread(AB_indicator, NumValue) %>% 
    mutate_at(vars(matches('percentage')), percent) %>% 
    drop_na(`Data quality_compl_age`)

str(tidy_by_indic2)


###Creando una base de datos que solo contiene porcetanges de resistencia, elimina los aislados y
###los controles de calidad. 
db_multiv <- tidy_by_indic2 %>% select(!contains('isolates') & !contains('compl'))

view(db_multiv)

db_multiv %>% group_by(Bacteria) %>% count()

db_multiv %>% filter(Bacteria == 'Escherichia coli') %>% view()
db_multiv %>% filter(Bacteria == 'Klebsiella pneumoniae') %>% view()



### se selecciona una sola bacteria y se eliminan columnas con na y filas con na. 
e.coli_db <- db_multiv %>%  
    filter(Bacteria == 'Escherichia coli') %>% 
    select_if(~!all(is.na(.))) %>%
    select(!contains('CR_')) %>%
    drop_na() %>%
    rename('Aminoglycosides'='Aminoglycosides_r_percentage', 
           'Carbapenems'='Carbapenems_r_percentage',
           'Fluoroquinolones'='Fluoroquinolones_r_percentage',
           'Aminopenicillins' = "Aminopenicillins_r_percentage",
           'cephalos_3er_gen' = "Third-generation cephalosporins_r_percentage",
           'Year' = 'Time',
           'Country' = 'RegionName', 
           'Code_C' = 'RegionCode')

e.coli_db %>% nrow()

k.pneum_db <- db_multiv %>%  
    filter(Bacteria == 'Klebsiella pneumoniae') %>% 
    select_if(~!all(is.na(.))) %>%
    select(!contains('CR_')) %>%
    drop_na() %>%
    rename('Aminoglycosides'='Aminoglycosides_r_percentage', 
           'Carbapenems'='Carbapenems_r_percentage',
           'Fluoroquinolones'='Fluoroquinolones_r_percentage',
           'cephalos_3er_gen' = "Third-generation cephalosporins_r_percentage",
           'Year' = 'Time',
           'Country' = 'RegionName', 
           'Code_C' = 'RegionCode')

k.pneum_db %>% nrow()


two_bact_b <- bind_rows(e.coli_db, k.pneum_db) %>% 
    select(!contains('minopenicillins')) %>% 
    view()
names(two_bact_b)

two_bact_b %>% group_by(Country) %>% count() %>% nrow()
nrow(two_bact_b)


##############################################################################################
#####Scraping para obtener regiones y organizacion de las regiones en la base de datos########
##############################################################################################


library(rvest)
library(tidyverse)
library(stringr)
url <- 'https://www.worldatlas.com/articles/the-four-european-regions-as-defined-by-the-united-nations-geoscheme-for-europe.html'
tab <- read_html(url) %>% html_nodes("table")

ws <- tab %>% html_table(fill = TRUE)

ws <- data.frame(ws) %>% 
    rbind(c('', '','Cyprus', '')) %>% 
    rbind(c('', 'Ireland', '', '')) %>%
    mutate(Eastern.Europe = recode(Eastern.Europe, 
                                   'Czech Republic' = 'Czechia')) %>%
    rename('Eastern' = 'Eastern.Europe', 
           'Northern' = 'Northern.Europe', 
           'Southern' = 'Southern.Europe', 
           'Western' = 'Western.Europe')

ws <- ws %>% gather(Region, Country, c(,1:4)) %>% 
    as_tibble() %>% 
    mutate_all(na_if,"") %>% 
    drop_na()

####Base de datos final para trabajo multivariado####
db_AB_regions <- left_join(two_bact_b, ws) %>% 
    write.csv2(file = 'E.coli_&_K.pneu_db', row.names = F) 
