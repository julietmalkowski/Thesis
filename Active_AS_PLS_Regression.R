#PLS Regression using Influent Microbial Community to Predict Active AS Microbial Community 
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
as = otus[grepl("AS-1|AS-2", otus$OTU),]
#remove V1
as = as[,c(-1,-2)]

#remove rows where OTU column does not contain AS-1 and AS-2
inf = otus[grepl("Inf", otus$OTU),]
#remove V1
inf = inf[,c(-1,-2)]

#Edit AS
colnames(as)[1] = "date"
#remove first 3 characters from date
as$date = substr(as$date, 6, nchar(as$date))
#turn date into date format
as$date = as.Date(as$date, format = "%m/%d/%Y")
#turn all other columns except date into floats
as[,-1] = sapply(as[,-1], as.numeric)
#remove date column
as = as[,-1]
#aggregate based on last 10 characters in rownames
as = aggregate(. ~ substr(rownames(as), 6, 16), as, mean)
#move column 1 to rownames
rownames(as) = as[,1]
#change column 1 to date
colnames(as)[1] = "date"
as$date = as.Date(as$date, format = "%m/%d/%Y")
#organize by date
as = as[order(as$date),]
#remove date column
as = as[,-1]

#move rownames to first column
inf = cbind(rownames(inf), inf)
#call column 1 date
colnames(inf)[1] = "date"
#remove first 3 characters from date
inf$date = substr(inf$date, 5, nchar(inf$date))
#make date rownames
rownames(inf) = inf$date
inf$date = as.Date(inf$date, format = "%m/%d/%Y")
inf = inf[order(inf$date),]
inf = inf[,-1]

#keep only matching rownames from inf and as
#remove last row from as
as = as[-nrow(as),]
#remove row 4
as = as[c(-4,-11),]

#average every 2 rows in inf
as <- as %>%
  group_by(group = (row_number() - 1) %/% 2) %>%
  summarise(across(where(is.numeric), mean))
#remove group column
as = as[,-1]
#keep columns with title "Zotu2", "Zotu3", "Zotu10", "Zotu5", "Zotu289"
as = subset(as,select = c("Zotu2", "Zotu3", "Zotu10", "Zotu5", "Zotu289"))

#Trying PLS only using relative abundances of the microbial community
#find size of otus- with samples as width and zOTUs as length
samples = dim(inf)[1]
zotus = dim(inf)[2]

#turn as and inf dataframes into numeric matrices
as <- data.frame(lapply(as, function(x) as.numeric(as.character(x))))

#subset otus into two sets- training and testing with trained_samples number used for training
train_otus = as.matrix(inf[1:round(samples*.8, 0),])
test_otus = as.matrix(inf[(round(samples*.8, 0)+1):samples,])

## Impact of Influent Microbial community on core AS community
#keep the rownames in metadata that are the same in train_otus
train_metadata = as.matrix(as[1:round(samples*.8, 0),])
test_metadata = as.matrix(as[(round(samples*.8, 0)+1):samples,])

train_metadata = as.matrix(apply(train_metadata, 2, as.numeric))
test_metadata = as.matrix(apply(test_metadata, 2, as.numeric))
train_otus = as.matrix(apply(train_otus, 2, as.numeric))
test_otus = as.matrix(apply(test_otus, 2, as.numeric))


## Plotting PLS Model
# ncomp is number of components
pls_model_tss <- pls(train_otus, train_metadata, 10, x.test = test_otus, y.test = test_metadata)
summary(pls_model_tss)
#plot the PLS model
tss_model = plot(pls_model_tss)
plotPredictions(pls_model_tss$res$cal, ncomp = 1, show.stat = TRUE) 
tss_vip = plotVIPScores(pls_model_tss)


