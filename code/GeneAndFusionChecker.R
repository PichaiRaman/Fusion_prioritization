library(readr)
library(stringr)

input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/final-list_candidate-fusion-genes.txt")
geneList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/CancerGeneList.txt")
fuseList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList.txt")

output <- input
output$Cancerous_Gene <- F
output$Cancerous_Gene_Symbol <- NA
output$Cancerous_Fusion <- F

geneRep <- 0
repAdjustedGene <- geneList[geneList$Count >= geneRep,]
fuseRep <- 0
repAdjustedFuse <- fuseList[fuseList$Count >= geneRep,]


for (rowIndex in 1:nrow(input)) {
  gene1 <- input[rowIndex,"Gene_1_symbol(5end_fusion_partner)"]
  gene2 <- input[rowIndex,"Gene_2_symbol(3end_fusion_partner)"]
  fuse <- paste(gene1,gene2,sep="--")
  
  gene1In <- gene1 %in% repAdjustedGene$Gene
  gene2In <- gene2 %in% repAdjustedGene$Gene
  fuseIn <- fuse %in% repAdjustedFuse$Fusions
  
  output[rowIndex,"Cancerous_Gene"] <- gene1In || gene2In
  output[rowIndex,"Cancerous_Fusion"] <- fuseIn
  
  if (gene1In && gene2In) {
    output[rowIndex,"Cancerous_Gene_Symbol"] <- paste(gene1,",",gene2)
  }else if (gene1In) {
    output[rowIndex,"Cancerous_Gene_Symbol"] <- gene1
  }else if (gene2In) {
    output[rowIndex,"Cancerous_Gene_Symbol"] <- gene2
  }
}

write.table(output, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/ProcessedInput_geneAndFusion.txt", sep="\t", row.names = FALSE)
View(output)
