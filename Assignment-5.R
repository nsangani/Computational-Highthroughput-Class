# Author: Neel Sangani
# Assignment 3 - Protein Interaction network
# Computational Analysis of High Throughput Biomedical Data

###-----------------------------Libraries installation and Loading---------------#########
install.packages("remotes")

install.packages("rstudioapi")

library("dplyr")
library(rstudioapi)
library(igraph)

####-------------------------Relative Dir Setup-----------------------##########

getwd()
wd = "C:/Users/sanga/OneDrive/Desktop/Course Work/Spring 2023/Computational High-throughput/Assignments/Assignment_5"
setwd(wd)
# Get the name of the directory in which the current file is located.
cur_dir = dirname(getSourceEditorContext()$path)

# Change the working directory to the directory with the current file.
setwd(cur_dir)

# Check that the working directory has changed to the directory with the current file.
getwd()

###-----------------------File Load----------------------------------###########


df <- read.table('Human-PPI.txt',header=TRUE,check.names = TRUE)

df = data.frame(df)

###--------------Using Igraph to generate the graph from PPI file----###########

g = graph_from_data_frame(df,directed=FALSE)

###--------------Calculating the degree------------------------------###########
degree = degree(g)



degree_df = data.frame(table(degree))



degree_df$degree = log10(as.numeric(degree_df$degree))
degree_df$Freq = log10(as.numeric(degree_df$Freq))

####----Degree Distribution Plot---------#######

plot(degree_df,main="Degree_Distribution_Plot",
     xlab="Degree", ylab="Frequency", pch=19)


###-------Calculating the coefficient-------------------------------############

coeff = transitivity(g,type='local',isolates = c("zero"))

coeff_avg = transitivity(g,type='average')


Scale = lm(degree_df$Freq ~ degree_df$degree)

#Since if the value of slope is better 2-3 then it is a scale free struture whereas, I got -1.908 , so the structure is not scale free


###----Loading the protein list---------#####

p1 = read.table('protein-list1.txt')

p2 = read.table('protein-list2.txt')

#---Loops to calculate shortest path from the protein list-----####

l_p1 = list()

for(i in 1:nrow(p1)){
  for(k in 1:nrow(p1)){
    
    
    l = try(shortest.paths(g,(p1[i,1]),to=(p1[k,1])),silent = TRUE)
    
    if(l==0 || l=="Error in as.igraph.vs(graph, to) : Invalid vertex names\n")
    {
      next
    } else {
      
      l_p1 = append(l_p1,l)
    }
  }
  
}

###-----Dataframe with the shortest path from protein list 1-----####

l_p1 = unlist(l_p1)

l_p1 = data.frame(l_p1)

l_p1

l_p2 = list()

for(i in 1:nrow(p2)){
  for(k in 1:nrow(p2)){
    
    
    l = try(shortest.paths(g,(p2[i,1]),to =(p2[k,1])),silent = TRUE)
    
    if(l==0 || l=="Error in as.igraph.vs(graph, v) : Invalid vertex names\n" || l=="Error in as.igraph.vs(graph, to) : Invalid vertex names\n")
    {
      next
    } else {
      
      l_p2 = append(l_p2,l)
    }
    
    
  }
  
}

###-----Dataframe with the shortest path from protein list 2-----####

l_p2 = unlist(l_p2)

l_p2 = data.frame(l_p2)

l_p2

####---path length distributions between the two protein sets using a wilcox test--------#####

wilcox = wilcox.test(as.numeric(l_p1$l_p1), as.numeric(l_p2$l_p2), alternative = "two.sided")



