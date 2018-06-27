library(readr)
library(readxl)
library(stringr)
library(dplyr)

getFactor <- function(geneACol,geneBCol) {
  return (factor(paste(geneACol,geneBCol,sep="--")))
}

DB_Data <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/ChimerDB3.0_ChimerSeq.xlsx")
Jack_Data <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/pancanfus.txt")
Cos_Data <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/CosmicFusionExport.tsv")
Fuse_Data <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/fusion_table.csv")

#Cleaning up DB_Data to get rid of weird dates
counter <- 1
for (gene in DB_Data$H_gene) {
  if (str_detect(gene,"\\.") && str_locate(gene,"\\.") < 4) {
    dotLoc <- str_locate(gene,"\\.")
    gene <- paste(toupper(substr(gene,dotLoc+1,str_length(gene))),substr(gene,0,dotLoc-1),sep="")
    DB_Data[counter,"H_gene"] <- gene
  }
  counter <- counter+1
}

#Creating GeneA and GeneB variables for Cos_Data
Cos_Data$GeneA <- NA
Cos_Data$GeneB <- NA
counter<-1
for (str in Cos_Data$`Translocation Name`) {
  if (!is.na(str)) {
    sub1 <- unlist(str_split(str,"\\{"))
    sub2 <- unlist(str_split(sub1[2],"\\_"))
    Cos_Data[counter,"GeneA"] <- sub1[1]
    Cos_Data[counter,"GeneB"] <- sub2[length(sub2)]
  }
  counter<-counter+1
}

#Cleaning up Fuse_levels with dates
counter<-1
for (gene in Fuse_Data$Head_gene_symbol) {
  if (str_detect(gene,"\\-") && str_locate(gene,"\\-") < 4) {
    newGene <- str_replace(gene,"\\-","")
    Fuse_Data[counter,"Head_gene_symbol"] <- newGene
  }
  counter<-counter+1
}

DB_factor <- getFactor(DB_Data$H_gene,DB_Data$T_gene)
Jack_factor <- getFactor(Jack_Data$Gene_A,Jack_Data$Gene_B)
Fuse_factor <- getFactor(Fuse_Data$Head_gene_symbol,Fuse_Data$Tail_gene_symbol)
Cos_factor <- getFactor(Cos_Data$GeneA,Cos_Data$GeneB)



fusionList <- levels(factor(c(levels(DB_factor),levels(Jack_factor),levels(Fuse_factor),levels(Cos_factor))))


finalFusionTable <- tibble(Fusions = fusionList)
finalFusionTable$DB <- 0
finalFusionTable$Jackson <- 0
finalFusionTable$Fusion <- 0
finalFusionTable$Cosmic <- 0
finalFusionTable$Count <- 0

counter <- 1
for (fuse in finalFusionTable$Fusions) {
  finalFusionTable[counter,"DB"] <- sum(DB_factor == fuse)
  finalFusionTable[counter,"Jackson"] <- sum(Jack_factor == fuse)
  finalFusionTable[counter,"Fusion"] <- sum(Fuse_factor == fuse)
  finalFusionTable[counter,"Cosmic"] <- sum(Cos_factor == fuse)
  finalFusionTable[counter,"Count"] <- sum(finalFusionTable[counter,"DB"],finalFusionTable[counter,"Jackson"],finalFusionTable[counter,"Fusion"],finalFusionTable[counter,"Cosmic"])
  print(paste(counter,"/",43466))
  counter <- counter + 1
}

View(finalFusionTable)
write.table(finalFusionTable, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList.txt", sep="\t",row.names = FALSE)

