# Author: Neel Sangani
# Assignment 1
# Computational Analysis of High Throughput Biomedical Data


# Set working directory 
setwd("C:/Users/sanga/OneDrive/Desktop/Course Work/Spring 2023/Computational High-throughput/Assignments/Assignment_1")

# Install and Load library
#install.packages('Hmisc')
#install.packages("PerformanceAnalytics")

#library('Hmisc')
library("PerformanceAnalytics")

# load the data

file.1 = read.csv('matrix1.txt', header = TRUE, sep = "\t") # 1046 x 14
file.1 = file.1[apply(file.1[,-1], 1, function(x) any(x != 0 || is.na(x))), ] # drop all rows with zero values
file.1 = file.1[-c(1,14)] # drop column mirna with catergorical data and column X with NA values
file.1 = as.matrix(file.1)
View(file.1)
colnames(file.1)
dim(file.1) # 540 x 12


file.2 = read.csv('matrix2.txt', header = TRUE, sep = "\t") # 1046 x 14
file.2 = file.2[apply(file.2[,-1], 1, function(x) any(x != 0 || is.na(x))), ] # drop all rows with zero values
file.2 = file.2[-c(1,14)] # drop column mirna with catergorical data and column X with NA values
file.2 = as.matrix(file.2)
View(file.2)
colnames(file.2)
dim(file.2) # 421 x 12


# Matrix.1 Pearson correlation
result.1= cor(file.1, method = "pearson", use = "complete.obs")
round(result.1,2)

# Matrix.1 Heatmap
chart.Correlation(result.1)
heatmap(result.1)

# Matrix.2 Pearson correlation
result.2 = cor(file.2, method = "pearson", use = "complete.obs")
round(result.2,2)

# Matrix.2 Heatmap
chart.Correlation(result.2)
heatmap(result.2)

# Matrix.1 x Matrix.2 Pearson correlation
corr = cor(result.1, result.2, method = 'pearson')
round(corr,2)

# Matrix.1 x Matrix.2 Heatmap
heatmap(corr)


