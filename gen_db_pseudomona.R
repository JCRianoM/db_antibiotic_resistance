###############################################################################################
########## CONSTRUCCIÓN Y COMPROBACIÓN DE FUNCIÓN PARA LIMPIEZ DE BASE DE DATOS################
########## DE ANTIBIÓTICOS ####################################################################
###############################################################################################

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
ab_bact_db <- ab_bact_db %>% mutate(NumValue = as.numeric(NumValue)) 

# código para limpiar la bd, cambiar nombres de los elementos en las filas para que sean más
# cortos y poder organizarlos mejor; se une el indicador con el AB y luego se pasan las filas a 
# columnas 'spread' de la unión entre AB e indicador; finalmente se limpia la base a través del
# del control de calidad de 'data quality_completeness age.
percent <- function(x) {x/100} # formula para convertir % en frecuencia
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

view(tidy_by_indic2)


###Creando una base de datos que solo contiene porcetanges de resistencia, elimina los aislados y
###los controles de calidad. 
db_multiv <- tidy_by_indic2 %>% select(!contains('isolates') & !contains('compl'))

### se selecciona una sola bacteria y se eliminan columnas con na y filas con na. 
ps.aur_db <- db_multiv %>%  
    filter(Bacteria == 'Pseudomonas aeruginosa' &
               '') %>% 
    select_if(~!all(is.na(.))) %>%
    select(!contains('CR_')) %>%
    drop_na() %>%
    rename('Aminoglycosides'='Aminoglycosides_r_percentage', 
           'Carbapenems'='Carbapenems_r_percentage', 
           'Ceftazidime'='Ceftazidime_r_percentage', 
           'Fluoroquinolones'='Fluoroquinolones_r_percentage', 
           'Pip_Tazobactam'='PiperacillinTazobactam_r_percentage', 
           'Year' = 'Time',
           'Country' = 'RegionName', 
           'Code_C' = 'RegionCode')

view(ps.aur_db)
nrow(ps.aur_db)

read.csv('pseudomona_db', header=TRUE, sep = ',') %>% view()

### escritura de función para utilizar organizar cualquier bacteria de base de datos 
#orginal: db_multiv; además se deja en la función la escritura de una archivo csv2 
#para omitir separador.

tidy_db_ab <- function(db, name_bacteria) {
    bact_db <- db %>%  
        filter(Bacteria == name_bacteria) %>% 
        select_if(~!all(is.na(.))) %>%
        select(!contains('CR_')) %>%
        drop_na()
    write.csv2(bact_db, file = name_bacteria, row.names = F)
    print(bact_db)
    
}

### comprobación de la función para otras bacterias
tidy_by_indic2 %>% group_by(Bacteria) %>% count()
tidy_db_ab(db_multiv, 'Pseudomonas aeruginosa')

### comprobación del archivo producido en la carpeta a través de read.csv2
read.csv2('Pseudomonas aeruginosa', header=TRUE) %>% view()
