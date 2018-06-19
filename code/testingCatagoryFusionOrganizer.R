library(readr)
library(readxl)
library(stringr)
library(dplyr)

addCancers <- function(currentCancers, levels) {
  for (lev in levels) {
    print(lev)
    print(currentCancers)
    if (!lev %in% currentCancers) {
      currentCancers <- c(currentCancers,lev)
    }
  }
  return (currentCancers)
}

DB_Data <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/ChimerDB3.0_ChimerSeq.xlsx")
Jack_Data <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/pancanfus.txt")

factoredCancers_DB <- factor(DB_Data$Cancertype_or_disease)
factoredCancers_Jack <- factor(Jack_Data$Cancer)

cancerCats <- c()
cancerCats <- addCancers(cancerCats, levels(factoredCancers_DB))
cancerCats <- addCancers(cancerCats, levels(factoredCancers_Jack))

cancerTable <- tibble(sample(NA,size=length(cancerCats)),replace=TRUE)
colnames(cancerTable) <- cancerCats

filter(DB_Data, Cancertype_or_disease == cancerCats[1])
for (cat in cancerCats) {
  filteredData <- filter(DB_Data, Cancertype_or_disease == cat)
  paste(filteredData$H_gene,"--",filteredData$T_gene)
}

