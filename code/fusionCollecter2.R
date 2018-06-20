library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

addData <- function(fusionTable,newData,dfName,geneAColName,geneBColName) {
  newData_factor <- factor(paste(newData[,geneAColName],newData[,geneBColName],sep="--"))
  
  fusionTable<-add_column(fusionTable,newDataStudy = 0, .before = "Count")
  names(fusionTable)[names(fusionTable) == 'newDataStudy'] <- dfName
  
  
  for (fuse in levels(newData_factor)) {
    count<-sum(newData_factor==fuse)
    if (!(fuse %in% fusionTable$Fusions)) {
      if (ncol(fusionTable)>3) {
        fusionTable[nrow(fusionTable)+1,] <- list(fuse,sample(0,ncol(fusionTable)-3),count,count)
      }else {
        fusionTable[nrow(fusionTable)+1,] <- list(fuse,count,count)
      }
    }else {
      fusionTable[fusionTable$Fusions == fuse,dfName] <- count
      fusionTable[fusionTable$Fusions == fuse,"Count"] <- fusionTable[fusionTable$Fusions == fuse,"Count"] + count
    }
  }
  return(fusionTable)
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

FusionTable <- tibble(Fusions,Count=0)
FusionTable <- addData(FusionTable,DB_Data,"DB","H_gene","T_gene")
FusionTable <- addData(FusionTable,Jack_Data,"Jackson","Gene_A","Gene_B")
FusionTable <- addData(FusionTable,Fuse_Data,"Fusion","Head_gene_symbol","Tail_gene_symbol")
FusionTable <- addData(FusionTable,Cos_Data,"Cosmic","GeneA","GeneB")

View(FusionTable)
write.table(FusionTable, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt", sep="\t",row.names = FALSE)

