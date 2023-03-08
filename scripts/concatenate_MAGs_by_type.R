# author: dlueckin
# date: Wed Jan 25 14:46:14 2023

# libraries ---------------------------------------------------------------
library(data.table)
library(stringr)
library(dplyr)

# get input ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
top20_file <- args[1]


# top20_df ----------------------------------------------------------------

top20_df <- fread(top20_file) %>%
    filter(label_after_cov != "remove") %>%
    filter(label_after_cov != "unclear")


# # collect and combine -----------------------------------------------------

for(i in 1:nrow(top20_df)){
    current_label = top20_df$label_after_cov[i]
    current_shortname = top20_df$shortname[i]
    current_station = top20_df$Sample[i]

    contig_file <- list.files("input/MAGs",
                              pattern = current_shortname, full.names = TRUE)
    system(paste0("cat ", contig_file, " >> input/MAGs_sorted/", current_station, "_combined_", current_label, "_contigs.fa"))
}
