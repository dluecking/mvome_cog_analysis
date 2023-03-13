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

to_plot$fold_change <- 0
to_plot$fold_color <- ""

for(i in 1:nrow(to_plot)){
    # check if 0
    if(to_plot$rel[i] == 0 | microbial_df$rel[microbial_df$LETTER == to_plot$LETTER[i]] == 0){
        to_plot$fold_change[i] <- 0
        to_plot$fold_color[i] <- "blue"
    
        # check which is larger
    }else if(to_plot$rel[i] >= microbial_df$rel[microbial_df$LETTER == to_plot$LETTER[i]] ){
        to_plot$fold_change[i] <- to_plot$rel[i] / microbial_df$rel[microbial_df$LETTER == to_plot$LETTER[i]] - 1
        to_plot$fold_color[i] <- "red"
    }else{
        to_plot$fold_change[i] <- microbial_df$rel[microbial_df$LETTER == to_plot$LETTER[i]] / to_plot$rel[i] * (-1) + 1
        to_plot$fold_color[i] <- "blue"
    }
}

# plot --------------------------------------------------------------------

ggplot(to_plot, aes(x = DESCRIPTION, y = fold_change, fill = fold_color)) +
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


 
# # testing -----------------------------------------------------------------
# 
# ggplot(to_plot[to_plot$label == "ev"], aes(x = to_plot$perc[to_plot$label == "ev"], y = microbial_df$perc)) +
#     geom_point()
# 
# 
# ggplot(to_plot[to_plot$label == "viral"], aes(x = to_plot$perc[to_plot$label == "viral"], y = microbial_df$perc)) +
#     geom_point()
# 
# ggplot(to_plot[to_plot$label == "gta"], aes(x = to_plot$perc[to_plot$label == "gta"], y = microbial_df$perc)) +
#     geom_point()



