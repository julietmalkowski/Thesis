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
load(file = "/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/filtered_ad_wwtp_p-values.RData")
#metadata
metadata = sample_data(wwtp)
metadata = as.data.frame(metadata)
#make Week column rownames in metadata
rownames(metadata) = metadata$Week

srts = read_excel("/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/filtered_AS_data_over35_weeks_otu.xlsx")
#remove column 1
srts = srts[,-1]
colnames(srts) = c("Week", "SRT_AD","OTU")
srts <- srts %>%
  pivot_wider(names_from = OTU, values_from = SRT_AD)

metadata = as.data.frame(metadata)
srts = as.data.frame(srts)
rownames(srts) = srts$Week

metadata = metadata[rownames(srts),]

#remove metadata rows with na
metadata = metadata[complete.cases(metadata),]
#put a zero instead of NA
srts[is.na(srts)] = 0


#remove column 1
srts = srts[,-1]
srts = srts[rownames(metadata),]

input_nutrients = metadata[,c(1,2,3,4,5,6,13,14)]
#merge srts and input nutrients based on rownames
#srts = cbind(srts, input_nutrients)


#Trying PLS only using relative abundances of the microbial community
#find size of otus- with samples as width and zOTUs as length
samples= dim(srts)[1]
zotus = dim(srts)[2]

#subset otus into two sets- training and testing with trained_samples number used for training
train_otus = as.matrix(srts[1:round(samples*.8, 0),])
test_otus = as.matrix(srts[(round(samples*.8, 0)+1):samples,])

#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(metadata[rownames(train_otus),])
test_metadata = as.matrix(metadata[rownames(test_otus),])

#gas produced
gas_produced_model <- pls(train_otus, train_metadata[,7], 8, x.test = test_otus, y.test = test_metadata[,7])
summary(gas_produced_model)
#plot the PLS model
tss_model = plot(gas_produced_model)
plotPredictions(gas_produced_model$res$cal, ncomp = 2, show.stat = TRUE) 
tss_vip = plotVIPScores(gas_produced_model)
