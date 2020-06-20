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

db_AB_regions <- left_join(two_bact_b, ws)

c %>% filter_all(any_vars(is.na(.))) %>% group_by(Country) %>% count()

view(c)

view(ws)

is.na(c) %>% sum()

