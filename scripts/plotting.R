# author: dlueckin
# date: Thu Jan 26 13:04:55 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# positions ---------------------------------------------------------------

positions <- c(
    "Cell cycle control, cell division, chromosome partitioning",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Posttranslational modification, protein turnover, chaperones",
    "Signal transduction mechanisms",
    "Defense mechanisms",
    "Extracellular structures",
    "Nuclear structure",
    "Cytoskeleton",
    
    "Mobilome: prophages, transposons",
    "Intracellular trafficking, secretion, and vesicular transport",
    
    "Energy production and conversion",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    
    "RNA processing and modification",
    "Chromatin structure and dynamics",
    "Translation, ribosomal structure and biogenesis",
    "Transcription",
    "Replication, recombination and repair",
    
    "General function prediction only",
    "Function unknown"
)

# import data -------------------------------------------------------------

data <- fread("intermediate/collected_data.tsv") 
data$perc <- data$rel * 100


# wrapped barplot ---------------------------------------------------------

to_plot <- data %>% filter(label != "microbial")
microbial_df <- data %>% filter(label == "microbial")

to_plot$delta_to_microbial <- 0


for(i in 1:nrow(to_plot)){
    to_plot$delta_to_microbial[i] <- to_plot$perc[i] - microbial_df$perc[microbial_df$LETTER == to_plot$LETTER[i]]
}

to_plot$delta_color <- ifelse(to_plot$delta_to_microbial >= 0, 'red', 'blue')


ggplot(to_plot, aes(x = DESCRIPTION, y = delta_to_microbial, fill = delta_color)) +
    geom_bar(stat = 'identity', color = 'black') +
    geom_hline(yintercept = 0) +
    theme_bw() +
    coord_flip() +
    facet_wrap(~label, ncol = length(unique(to_plot$label)),
               labeller = labeller(origin = 
                                       c("ev" = "EV Producer",
                                         "gta" = "GTA Producer",
                                         "viral" = "Transducer"))) +
    theme(legend.position = "None") +
    xlab('') +
    ylab('') +
    scale_fill_manual(values= c("#5B84B1FF", "#FC766AFF")) +
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)) +
    scale_x_discrete(limits = rev(positions)) 

ggsave(filename = "plots/cog_barplot.png", last_plot(), width = 7, height = 4)
ggsave(filename = "plots/cog_barplot.svg", last_plot(), width = 7, height = 4)





