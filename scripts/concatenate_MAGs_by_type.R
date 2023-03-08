# author: dlueckin
# date: Wed Jan 25 14:46:14 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(dplyr)

# working directory -------------------------------------------------------
this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this_dir)
print(paste0("Setting wd to: \n ", this_dir))


# top20_df ----------------------------------------------------------------

top20_df <- fread("input-data - top20_combined.tsv") %>% 
    filter(Dataset != "Heligoland") %>% 
    filter(label_after_cov != "remove") %>% 
    filter(label_after_cov != "unclear")


# collect and combine -----------------------------------------------------

for(i in 1:nrow(top20_df)){
    current_label = top20_df$label_after_cov[i]
    current_shortname = top20_df$shortname[i]
    current_station = top20_df$Sample[i]
    
    contig_file <- list.files("/run/user/1000/gvfs/sftp:host=linux-desktop-1.mpi-bremen.de,user=dlueckin/home/dlueckin/projects/misc/three_stations/contigs",
                              pattern = current_shortname, full.names = TRUE)
    system(paste0("cat ", contig_file, " >> contigs/", current_station, "_combined_", current_label, "_contigs.fa"))
}
# 
