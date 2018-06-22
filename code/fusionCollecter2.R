library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

addData <- function(fusionTable,newData,dfName,geneACol,geneBCol) {
  newData_factor <- factor(paste(geneACol,geneBCol,sep="--"))
  fusionTable<-add_column(fusionTable,newDataStudy = 0, .before = "Count")
  names(fusionTable)[names(fusionTable) == 'newDataStudy'] <- dfName
  newData_level<-levels(newData_factor)
  
  print(paste("----------------------------------",dfName,"----------------------------------"))
  counter <- 1
  for (fuse in newData_level) {
    newCount<-sum(newData_factor==fuse)
    fuseCount<- sum(fusionTable$Fusions==fuse)
    if (fuseCount < 1) {
      if (ncol(fusionTable)>3) {
        fusionTable[nrow(fusionTable)+1,] <- c(fuse,sample(0,ncol(fusionTable)-3,replace=TRUE),newCount,newCount)
      }else {
        fusionTable[nrow(fusionTable)+1,] <- c(fuse,newCount,newCount)
      }
    }else {
      fusionTable[fusionTable$Fusions == fuse,dfName] <- newCount
      fusionTable[fusionTable$Fusions == fuse,"Count"] <- as.numeric(fusionTable[fusionTable$Fusions == fuse,'Count']) + newCount
    }
    if (counter %% 1000 == 0) {
      print (paste(counter,"/",length(newData_level)))
    }
    counter<-counter+1
  }
  print("Done")
  return(fusionTable)
}


DB_Data <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/ChimerDB3.0_ChimerSeq.xlsx")
Jack_Data <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/pancanfus.txt")
Cos_Data <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/CosmicFusionExport.tsv")
Fuse_Data <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/fusion_table.csv")

#Cleaning up DB_Data to standarize dates
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

Fusions <- c("Blank")
FusionTable <- tibble(Fusions,Count=0)

#Data from CHIMERDB30
FusionTable <- addData(FusionTable,DB_Data,"DB",DB_Data$H_gene,DB_Data$T_gene)

FusionTable <- FusionTable[-1,]

#Data from Jackson Lab
FusionTable <- addData(FusionTable,Jack_Data,"Jackson",Jack_Data$Gene_A,Jack_Data$Gene_B)

#Data from Fusion Cancer
FusionTable <- addData(FusionTable,Fuse_Data,"Fusion",Fuse_Data$Head_gene_symbol,Fuse_Data$Tail_gene_symbol)

#Data from Cosmic
FusionTable <- addData(FusionTable,Cos_Data,"Cosmic",Cos_Data$GeneA,Cos_Data$GeneB)

FusionTable<- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt")
ChiTA_Data <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/Maher-Nature-fusions.xls")
TIC_Data <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/allseqs_TICdb.txt", col_names = FALSE)
CG_Data <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/ConjoinG.tsv")

#Data from ChiTaRS
FusionTable <- addData(FusionTable,ChiTA_Data,"ChiTaRS",ChiTA_Data$gene1,ChiTA_Data$gene2)

#Data from TICdb
FusionTable <- addData(FusionTable,TIC_Data,"TICdb",TIC_Data$X1,TIC_Data$X2)

#Data from ConjoinG
FusionTable <- addData(FusionTable,CG_Data,"ConjoinG",CG_Data$`Parent gene 1`,CG_Data$`Parent gene 2`)

View(FusionTable)
write.table(FusionTable, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt", sep="\t",row.names = FALSE)

