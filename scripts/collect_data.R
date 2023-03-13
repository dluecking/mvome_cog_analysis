# author: dlueckin
# date: Thu Jan 26 12:10:39 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)

# prepare template df -----------------------------------------------------

COG_template_df <- fread("helper_files/template_classifier_count.tsv")
# COG_template_df <- fread("~/template_classifier_count.tsv")


# repeat the DF, since we need to check for the origin (later facet wrap) -
repeat_df <- function(d, n) {
    return(do.call("rbind", replicate(n, d, simplify = FALSE)))
}

labels = c("ev", "gta", "viral", "microbial")

COG_template_df <- repeat_df(COG_template_df, length(labels)) 
COG_template_df$label <- rep(labels, each = nrow(COG_template_df) / length(labels))


# import microbial blast out and combine ------------------------------------

blast_df <- data.table()

files <- list.files("intermediate/blast_out/microbial/")
for(file in files){
    file_long = paste0("intermediate/blast_out/microbial/", file)
    if(file.info(file_long)$size  > 0){
        tmp_df <- fread(file_long, select = c("V2"))
        tmp_df$origin <- paste0(file, "_microbial")
        blast_df <- rbind(blast_df, tmp_df)
    }
    
}
rm(tmp_df)


# import ev_gta_and_viral blastout and combine ----------------------------

files <- list.files("intermediate/blast_out/peDNA/")
for(file in files){
    file_long = paste0("intermediate/blast_out/peDNA/", file)
    if(file.info(file_long)$size  > 0){
        tmp_df <- fread(file_long, select = c("V2"))
        tmp_df$origin <- file
        blast_df <- rbind(blast_df, tmp_df)
    }
}
rm(tmp_df)


# add COG db, based on hit ------------------------------------------------

def <- fread("helper_files/cog-20.def.tab")
# cog needs to be treated with sed -i 's/,,,,,/,,,,/' cog-20.cog.csv
cog <- fread("helper_files/cog-20.cog.csv", fill = T, sep = ",") 

blast_df$COG <- cog$V7[match(blast_df$V2, cog$V3)]


# add actual letter per hit -----------------------------------------------

blast_df$letter <- def$V2[match(blast_df$COG, def$V1)]

# this df is set up, so when we fill in the data later on, we dont have to filter the big blast_df one all the time, but just this one small df
df <- blast_df %>% 
    select(origin, letter) %>% 
    group_by(origin, letter) %>% 
    tally()



# fill in data ------------------------------------------------------------

for(i in 1:nrow(COG_template_df)){
    COG_template_df$COUNT[i] <- df %>% 
        filter(str_detect(string = origin, pattern = COG_template_df$label[i])) %>%     # matches origin
        filter(str_detect(string = letter, pattern = COG_template_df$LETTER[i])) %>%    # matches label
        ungroup() %>%                                                                   # ungroup, otherwise the origin appears in the select
        select(n) %>% 
        summarise(total = sum(n)) %>% 
        pull()
}


# normalise by number of hits for this label ------------------------------
# I bet this is possible in pretty

COG_template_df$rel[COG_template_df$label == "ev"] <- 
    COG_template_df$COUNT[COG_template_df$label == "ev"] / sum(COG_template_df$COUNT[COG_template_df$label == "ev"])

COG_template_df$rel[COG_template_df$label == "viral"] <- 
    COG_template_df$COUNT[COG_template_df$label == "viral"] / sum(COG_template_df$COUNT[COG_template_df$label == "viral"])

COG_template_df$rel[COG_template_df$label == "gta"] <- 
    COG_template_df$COUNT[COG_template_df$label == "gta"] / sum(COG_template_df$COUNT[COG_template_df$label == "gta"])

COG_template_df$rel[COG_template_df$label == "microbial"] <- 
    COG_template_df$COUNT[COG_template_df$label == "microbial"] / sum(COG_template_df$COUNT[COG_template_df$label == "microbial"])

fwrite(COG_template_df, "intermediate/collected_data.tsv")
