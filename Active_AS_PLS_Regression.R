#PLS Regression on Active AS Community Analysis

#utilizing abundance values of the active AS community to predict TSS Removal and COD Removal
#(since it was most promising from last code iteration)

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

#load phyloseq object
load(file = "/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/wwtp_raw.RData")

otus = otu_table(wwtp)
otus = as.data.frame(otus)
otus = t(otus)
otus = as.matrix(otus)
#move rownames to first column
otus = cbind(rownames(otus), otus)
otus = as.data.frame(otus)
#move rownames to column
colnames(otus)[1] = "OTU"
#remove rows where OTU column does not contain AS-1 and AS-2
otus = otus[grepl("AS-1|AS-2", otus$OTU),]
#remove V1
otus = otus[,c(-1,-2)]
#metadata
metadata = sample_data(wwtp)
metadata = as.data.frame(metadata)

#add columns from metadata to otus as a predictive measure
metadata <- metadata[match(rownames(otus), rownames(metadata)), ]
#move rownames to column 1
metadata = cbind(rownames(metadata), metadata)
#make column 1 date
colnames(metadata)[1] = "date"
#remove last 3 columns
metadata = metadata[,c(-44,-45,-46)]
#keep only characters 6-16
metadata$date = substr(metadata$date, 6, 16)
metadata[,-1] = sapply(metadata[,-1], as.numeric)
#remove date column
metadata = metadata[,-1]
metadata = aggregate(. ~ substr(rownames(metadata), 6, 16), metadata, mean)
#change column 1 to rownames
rownames(metadata) = metadata[,1]
colnames(metadata)[1] = "date"
metadata = metadata[order(metadata$date),]

metadata = metadata[,-1]

#turn rownames into a column
otus = cbind(rownames(otus), otus)

colnames(otus)[1] = "date"
#remove first 3 characters from date
otus$date = substr(otus$date, 6, nchar(otus$date))
#turn date into date format
otus$date = as.Date(otus$date, format = "%m/%d/%Y")
#turn all other columns except date into floats
otus[,-1] = sapply(otus[,-1], as.numeric)
#remove date column
otus = otus[,-1]
#aggregate based on last 10 characters in rownames
otus = aggregate(. ~ substr(rownames(otus), 6, 16), otus, mean)
#move column 1 to rownames
rownames(otus) = otus[,1]
#change column 1 to date
colnames(otus)[1] = "date"
#organize by date
otus = otus[order(otus$date),]
#remove date column
otus = otus[,-1]


#Trying PLS only using relative abundances of the microbial community
#find size of otus- with samples as width and zOTUs as length
samples= dim(otus)[1]
zotus = dim(otus)[2]

#subset otus into two sets- training and testing with trained_samples number used for training
train_otus = as.matrix(otus[1:round(samples*.8, 0),])
test_otus = as.matrix(otus[(round(samples*.8, 0)+1):samples,])

## TSS Load Removed (column 7)
#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(metadata[1:round(samples*.8, 0),])
test_metadata = as.matrix(metadata[(round(samples*.8, 0)+1):samples,])

## Plotting PLS Model
# ncomp is number of components
pls_model_tss <- pls(train_otus, train_metadata[,7], 10, x.test = test_otus, y.test = test_metadata[,7])
summary(pls_model_tss)
#plot the PLS model
cod_model = plot(pls_model_tss)
plotPredictions(pls_model_tss$res$cal, ncomp = 5, show.stat = TRUE) 
cod_vip = plotVIPScores(pls_model_tss)


## COD Load Removed (column 11)
pls_model_cod <- pls(train_otus, train_metadata[,11], 10, x.test = test_otus, y.test = test_metadata[,11])
summary(pls_model_cod)
#plot the PLS model
cod_model = plot(pls_model_cod)
plotPredictions(pls_model_cod$res$cal, ncomp = 5, show.stat = TRUE) 
cod_vip = plotVIPScores(pls_model_cod)
