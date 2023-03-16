# author: dlueckin
# date: Thu Jan 26 13:04:55 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# import data -------------------------------------------------------------

data <- fread("intermediate/collected_data_with_transposons.tsv") 
data$perc <- data$rel * 100



# barplot comparing rel between groups: -----------------------------------

data <- data %>% 
    filter(str_detect(LETTER, "X"))

ggplot(data, aes(x = label, y = perc, fill = LETTER)) +
    geom_bar(stat = 'identity', position = "stack")
