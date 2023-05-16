# Author: Neel Sangani
# Assignment 4 - Motif prediction via PWM 
# Computational Analysis of High Throughput Biomedical Data

# Set working directory 
wd = "C:/Users/sanga/OneDrive/Desktop/Course Work/Spring 2023/Computational High-throughput/Assignments/Assignment_4"
setwd(wd)
#rm(list = ls())

# load library
library(dplyr)

# Load the count matrix
pcm <- read.csv("argR-counts-matrix.txt", sep='\t', header = F, row.names = 1) 
pcm <- pcm[,-1] #drop first column with "|"
pcm <- pcm + 1 #pseudo count adjusted matrix
# head(pcm)

# Convert counts to frequencies
pfm <- pcm %>%
  mutate_if(is.numeric, ~ ./sum(.))

# Convert frequencies to log-odds scores
pwm <- apply(pfm, 2, function(x) log2(x / 0.25))
rm(pcm,pfm)
head(pwm)


E_coli_gene_seq <- read.csv("E_coli_K12_MG1655.400_50", sep = '\\', header= F, col.names = c('geneID','seq','X'), strip.white=T) #
E_coli_gene_seq <- E_coli_gene_seq[,-3] # drop third column
head(E_coli_gene_seq, n=2)

    
motif_scores = function(x) {
  max_Scores = c()
  for (seq_num in x) {
    scores = c()
    window_size = 18
    #seq_num = 'aacggcagaccaacatcaactgcaagctttacgcgaa' # test
    for (i in 1:(nchar(seq_num) - window_size + 1)) {
      window_seq <- substring(seq_num, i, i + window_size - 1)
      # print(window_seq)
      window_seq_str <- strsplit(x=window_seq,split='')
      score = vector()
      for (i in 1:nchar(window_seq)) {
        score[i] <- pwm[window_seq_str[[1]][i],i]
      }
      # print(sum(score))
      scores <- append(scores, sum(score))
    }
    max_Scores <- append(max_Scores, max(scores))
  }
  return(max_Scores)
}

max_scores <- motif_scores(E_coli_gene_seq$seq)
E_coli_gene_seq = cbind(E_coli_gene_seq, max_scores)
E_coli_gene_seq<- E_coli_gene_seq[order(E_coli_gene_seq$max_scores,decreasing=T),]
head(E_coli_gene_seq[, c("geneID","max_scores")], n=30) # Top 30 geneID with highest score

