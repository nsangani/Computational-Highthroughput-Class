# Author: Neel Sangani
# Assignment 2 - Half life of Yeast Transcripts
# Computational Analysis of High Throughput Biomedical Data


# Set working directory 
wd = "C:/Users/sanga/OneDrive/Desktop/Course Work/Spring 2023/Computational High-throughput/Assignments/Assignment_2"
setwd(wd)

# Loading data
file.1 = read.delim('DecayTimecourse.txt', header = TRUE) #6185 x 28

# Preprocessing
colnames(file.1) <- as.character(file.1[1,])
file.1 = file.1[-1,]
clean_file = na.omit(file.1) #113 x 28
View(clean_file)


# Slope of each timestamp using linear regression
slope  <-  function(x){
  coef(lm(I(2:10)~x))[2]
}

clean_file$slope1  <-  
  apply(clean_file[,c(2:10)],
        1,
        slope)

clean_file$slope2  <-  
  apply(clean_file[,c(11:19)],
        1,
        slope)

clean_file$slope3  <-  
  apply(clean_file[,c(20:28)],
        1,
        slope)

# avg-slope across three timestamps
clean_file$average.slope  <-  (clean_file$slope1 + clean_file$slope2 + clean_file$slope3)*(-1)/3

# half-life
clean_file$half_life <- (log(2)/clean_file$average.slope)

# sort highest to lowest half-life
clean_file = clean_file[order(-clean_file$half_life),]

View(clean_file)

# Select top 10% of the genes of Ontology
n = floor(nrow(clean_file)*.1) 

head(clean_file[c(1)], n = n)
tail(clean_file[c(1)], n = n)


# Please refer to the attached images for the ontology results.


