#PLS Regression on Active AS Community Analysis
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
load(file = "/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/filtered_as_wwtp_p-values.RData")

otus = otu_table(wwtp)
otus = as.data.frame(otus)
otus = t(otus)
otus = as.matrix(otus)
#move rownames to first column
otus = cbind(rownames(otus), otus)
#name column 1 dates
colnames(otus)[1] = "date"
otus = as.data.frame(otus)
#remove first 3 characters from date
otus$date = substr(otus$date, 4, nchar(otus$date))
#print type of character otus$date
class(otus$date)
#turn date into date format
otus$date = as.Date(otus$date, format = "%m/%d/%Y")
#organize by date
otus = otus[order(otus$date),]
#remove date column
otus = otus[,-1]

#metadata
metadata = sample_data(wwtp)
metadata = as.data.frame(metadata)

#add columns from metadata to otus as a predictive measure
metadata <- metadata[match(rownames(otus), rownames(metadata)), ]

## TSS Load Removed
#keep only variables we want to test: columns 7, 9, 11
tss_load = metadata[,c(7)]

#Trying PLS only using relative abundances of the microbial community
#find size of otus- with samples as width and zOTUs as length
samples= dim(otus)[1]
zotus = dim(otus)[2]

#subset otus into two sets- training and testing with trained_samples number used for training
train_otus = as.matrix(otus[1:round(samples*.8, 0),])
test_otus = as.matrix(otus[(round(samples*.8, 0)+1):samples,])

#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(tss_load[rownames(train_otus),])
test_metadata = as.matrix(tss_load[rownames(test_otus),])

#make sure all dataframes are numeric
train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
train_otu = as.matrix(apply(train_otus, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
test_otu = as.matrix(apply(test_otus, 2, as.numeric))

# ncomp is number of components
pls_model_tss <- pls(train_otu, train_metadata, 10, x.test = test_otu, y.test = test_metadata)

## BOD/CBOD Load Removed (col 9)
bod_load = metadata[,c(9)]

#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(bod_load[rownames(train_otus),])
test_metadata = as.matrix(bod_load[rownames(test_otus),])
#make sure all dataframes are numeric
train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
# ncomp is number of components
pls_model_bod <- pls(train_otu, train_metadata, 10, x.test = test_otu, y.test = test_metadata)


## COD Load Removed (11)
cod_load = metadata[,c(11)]
#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(cod_load[rownames(train_otus),])
test_metadata = as.matrix(cod_load[rownames(test_otus),])
#make sure all dataframes are numeric
train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
# ncomp is number of components
pls_model_cod <- pls(train_otu, train_metadata, 10, x.test = test_otu, y.test = test_metadata)



## Ammonia FE (13)
ammonia_load = metadata[,c(13)]
#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(ammonia_load[rownames(train_otus),])
test_metadata = as.matrix(ammonia_load[rownames(test_otus),])
#make sure all dataframes are numeric
train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
# ncomp is number of components
pls_model_ammonia <- pls(train_otu, train_metadata, 10, x.test = test_otu, y.test = test_metadata)


summary(pls_model_ammonia)
summary(pls_model_cod)
summary(pls_model_bod)
summary(pls_model_tss)


par(mfrow = c(2, 2))
plotVIPScores(pls_model_ammonia)
plotVIPScores(pls_model_cod)
plotVIPScores(pls_model_bod)
plotVIPScores(pls_model_tss)

par(mfrow = c(2, 2))
plotPredictions(pls_model_ammonia$res$cal, ncomp = 2, show.stat = TRUE) 
plotPredictions(pls_model_cod$res$cal, ncomp = 1, show.stat = TRUE) 
plotPredictions(pls_model_bod$res$cal, ncomp = 1, show.stat = TRUE)
plotPredictions(pls_model_tss$res$cal, ncomp = 2, show.stat = TRUE) 


