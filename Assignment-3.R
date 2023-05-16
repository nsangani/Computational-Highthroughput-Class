# Author: Neel Sangani
# Assignment 3 - Operon Prediction
# Computational Analysis of High Throughput Biomedical Data


# Set working directory 
wd = "C:/Users/sanga/OneDrive/Desktop/Course Work/Spring 2023/Computational High-throughput/Assignments/Assignment_3"
setwd(wd)
# rm(list = ls())

# Open ppt file + skip 2 lines + tab + header func
file_open <- function(file_name) {
  file <- read.csv(file_name, skip = 2, header = T, sep = '\t')
  return(file)
}

# Split the Location column into start/end col using the ".." separator and cbind to df and sort by strand
split_start_end_sort <- function(df) {
  new_df <- data.frame(start = numeric(length(df$Location)), end = numeric(length(df$Location)))
  locations <- strsplit(df$Location, "\\..")
  for (i in 1:length(locations)) {
    new_df$start[i] <- as.numeric(locations[[i]][1])
    new_df$end[i] <- as.numeric(locations[[i]][2])
  }
  merge = cbind(df,new_df)
  return(merge[order(merge$Strand, decreasing = FALSE), ])
}

# Load the file
E_coli = file_open('E_coli_K12_MG1655.ptt')
B_subtilis = file_open('B_subtilis_168.ptt')
Synechocystis = file_open('Synechocystis_PCC6803_uid159873.ptt')
Halobacterium = file_open('Halobacterium_NRC1.ptt')


calculate_distance <- function(df) {
  df$distance <- c(df[2:nrow(df),"start"] - df[1:nrow(df)-1,"end"], NA)
  return(df)
}

operon_pred <- function(df) {
  result <- data.frame(PID = character(), distance = numeric())
  curr_group <- NULL
  prev_distance <- NA
  
  for (i in seq(nrow(df))) {
    if (!is.na(df$distance[i]) && df$distance[i] < 50) {
      if (is.null(curr_group)) {
        curr_group <- toString(df$PID[i])
      } else {
        curr_group <- paste(curr_group, df$PID[i], sep = ", ")
      }
    } else {
      if (!is.null(curr_group)) {
        result <- rbind(result, data.frame(PID = curr_group, distance = prev_distance))
        curr_group <- NULL
      }
      result <- rbind(result, data.frame(PID = df$PID[i], distance = df$distance[i]))
    }
    prev_distance <- df$distance[i]
  }
  if (!is.null(curr_group)) {
    result <- rbind(result, data.frame(PID = curr_group, distance = prev_distance))
  }
  return(result)
}

# Split the df by sign and perform distance calculation, operon pred, compile
split_predit_merge_sign <- function(df){
  df_pos <- df[df$Strand == '+',]
  df_neg <- df[df$Strand == '-',]
  df_pos <- calculate_distance(df_pos) 
  df_neg <- calculate_distance(df_neg)
  df_operon_pos <- operon_pred(df_pos)
  strand = df_pos[1:nrow(df_operon_pos),"Strand"]
  df_operon_pos <- cbind(df_operon_pos[,-ncol(df_operon_pos)], strand)
  
  df_operon_neg <- operon_pred(df_neg)
  strand = df_neg[1:nrow(df_operon_neg),"Strand"]
  df_operon_neg <- cbind(df_operon_neg[,-ncol(df_operon_neg)], strand)
  
  df_operon = rbind(df_operon_pos, df_operon_neg)
  return(df_operon)
}

species_operon_prediction <- function (species_name) {
  species_name = split_start_end_sort(species_name)
  species_name_operon <- split_predit_merge_sign(species_name)
}

E_coli = species_operon_prediction(E_coli)
B_subtilis = species_operon_prediction(B_subtilis)
Synechocystis = species_operon_prediction(Synechocystis)
Halobacterium = species_operon_prediction(Halobacterium)

nrow(E_coli)
nrow(B_subtilis)
nrow(Synechocystis)
nrow(Halobacterium)

# Part 2 - predict operons from crop microbiome
columns = c("PID","source","V3","start","end","period","Strand","zero","ID")
crop <- read.csv('2088090036.gff', skip = 1, header = T, sep = '\t', col.names = columns ) #23908
crop <- subset(crop, V3 == "CDS") #23908

crop_operon = split_predit_merge_sign(crop)

nrow(crop_operon)








