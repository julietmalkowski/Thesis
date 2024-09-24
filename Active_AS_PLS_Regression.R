#PLS Regression using SRT time of Active AS Microbial Community to predict COD/TSS Removal over time

library(phyloseq)
library(vegan)
library(picante)
library(ggplot2)
library(tibble)
library(readxl)
library(dplyr)
library(ape)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(plotly)
library(microbiome)
library(mdatools)

#For reproducibility
set.seed(123)

#load excel file
wwtp1 = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/genus_table.xlsx", sheet = 'Sheet1')
wwtp2 = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/genus_table.xlsx", sheet = 'Sheet2')
wwtp3 = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/genus_table.xlsx", sheet = 'Sheet3')
wwtp4 = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/genus_table.xlsx", sheet = 'Sheet4')
wwtp5 = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/genus_table.xlsx", sheet = 'Sheet5')

wwtp = rbind(wwtp1, wwtp2, wwtp3, wwtp4, wwtp5)
#move column names to row 1
colnames(wwtp) = wwtp[1,]
#remove column 1
wwtp = wwtp[,-1]
#change column names to 'week', 'OTU', 'SRT_AS', 'SRT_AD'
colnames(wwtp) = c("week", "OTU", "SRT_AS", "SRT_AD")
#delete last column
wwtp = wwtp[,-ncol(wwtp)]
#make week column rownames and OTUs columns with values being SRT_AS
pivot_table <- wwtp %>%
  pivot_wider(names_from = OTU, values_from = SRT_AS)
#replace na with 0
pivot_table[is.na(pivot_table)] = 0
#remove columns with all 0s
pivot_table = pivot_table[, colSums(pivot_table != 0) > 0]
#remove columns if they contain values in less than 75% of all rows
pivot_table = pivot_table[, colMeans(pivot_table != 0) > 0.75]
#turn week column numeric
pivot_table$week = as.numeric(pivot_table$week)
#turn into matrix
pivot_table = as.matrix(pivot_table)
#make rownames the week column
rownames(pivot_table) = pivot_table[,1]

metadata = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/AS_metadata.xlsx")
#remove column 1
metadata = metadata[,c(-1,-14,-16,-18)]

input_nutrients = metadata[,c(1,2,3,4,5,6,8,10,12,13,14,15,19)]
input_nutrients = input_nutrients %>% group_by(Week) %>% summarise_all(mean)
#add input nutrients to pivot_table based on week
pivot_table = as.data.frame(pivot_table)
input_nutrients = as.data.frame(input_nutrients)
#change Week to week
input_nutrients = input_nutrients %>% rename(week = Week)
pivot_table = merge(pivot_table, input_nutrients, by = "week")
#make week rownames in pivot_table
rownames(pivot_table) = pivot_table[,1]

metadata = metadata[,-c(1,2,3,4,5,6,8,10,12,13,15,14)]
#groupby week and average metadata
metadata = metadata %>% group_by(Week) %>% summarise_all(mean)
#make a matrix
metadata = as.matrix(metadata)
#make week rownames in metadata
rownames(metadata) = metadata[,1]

#keep only the rows in pivot_table that are in metadata based on week column
metadata = as.data.frame(metadata)
pivot_table = as.data.frame(pivot_table)
metadata = metadata[rownames(pivot_table),]

#remove na from metadata
metadata = metadata[complete.cases(metadata),]
#delete column 1
metadata = metadata[,-1]
pivot_table = pivot_table[,-1]
#as matrix
metadata = as.matrix(metadata)
pivot_table = as.matrix(pivot_table)
#Trying PLS only using relative abundances of the microbial community
#find size of otus- with samples as width and zOTUs as length
samples = dim(pivot_table)[1]
zotus = dim(pivot_table)[2]

#subset otus into two sets- training and testing with trained_samples number used for training
train_otus = as.matrix(pivot_table[1:round(samples*.8, 0),])
test_otus = as.matrix(pivot_table[(round(samples*.8, 0)+1):samples,])

## Impact of Influent Microbial community on core AS community
#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(metadata[1:round(samples*.8, 0),])
test_metadata = as.matrix(metadata[(round(samples*.8, 0)+1):samples,])

train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
train_otus = as.matrix(apply(train_otus, 2, as.numeric))
test_otus = as.matrix(apply(test_otus, 2, as.numeric))


## Plotting PLS Model
# ncomp is number of components
pls_model_tss <- pls(train_otus, train_metadata[,3], 10, x.test = test_otus, y.test = test_metadata[,2])
summary(pls_model_tss)
#plot the PLS model
tss_model = plot(pls_model_tss)
plotPredictions(pls_model_tss$res$cal, ncomp = 1, show.stat = TRUE) 
tss_vip = plotVIPScores(pls_model_tss)

