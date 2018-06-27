library(readr)
library(stringr)

listChecker <- function(input,gene1ColName,gene2ColName,outputFilePath,geneRep=1,fuseRep=1) {
  output <- input
  output$Cancerous_Gene <- F
  output$Cancerous_Gene_Symbol <- NA
  output$Cancerous_Fusion <- F
  output$Cancerous_Fusion_Symbol <- NA
  
  repAdjustedGene <- geneList[geneList$Count >= geneRep,]
  repAdjustedFuse <- fuseList[fuseList$Count >= geneRep,]
  
  for (rowIndex in 1:nrow(input)) {
    gene1 <- input[rowIndex,gene1ColName]
    gene2 <- input[rowIndex,gene2ColName]
    fuse <- paste(gene1,gene2,sep="--")
    
    gene1In <- gene1 %in% repAdjustedGene$Gene
    gene2In <- gene2 %in% repAdjustedGene$Gene
    fuseIn <- fuse %in% repAdjustedFuse$Fusions
    
    output[rowIndex,"Cancerous_Gene"] <- gene1In || gene2In
    output[rowIndex,"Cancerous_Fusion"] <- fuseIn
    
    if (gene1In && gene2In) {
      output[rowIndex,"Cancerous_Gene_Symbol"] <- paste(gene1,gene2,sep=",")
    }else if (gene1In) {
      output[rowIndex,"Cancerous_Gene_Symbol"] <- gene1
    }else if (gene2In) {
      output[rowIndex,"Cancerous_Gene_Symbol"] <- gene2
    }
    
    if (fuseIn) {
      output[rowIndex,"Cancerous_Fusion_Symbol"] <- paste(gene1,gene2,sep="--")
    }
  }
  
  write.table(output, outputFilePath, sep="\t", row.names = FALSE)
  return (output)
}

input1 <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/final-list_candidate-fusion-genes.txt")
input2 <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv")
geneList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/CancerGeneList.txt")
fuseList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt")

#Adding gene Symbols to input2
input2$LeftGeneSymbol <- NA
input2$RightGeneSymbol <- NA
modifiedLeft <- unlist(str_split(input2$LeftGene,"\\^"))
addCounter <- 1
for (counter in 1:length(modifiedLeft)) {
  if (counter %% 2 != 0) {
    input2[addCounter,"LeftGeneSymbol"] <- modifiedLeft[counter]
    addCounter<- addCounter+1
  }
}
modifiedRight <- unlist(str_split(input2$RightGene,"\\^"))
addCounter <- 1
for (counter in 1:length(modifiedRight)) {
  if (counter %% 2 != 0) {
    input2[addCounter,"RightGeneSymbol"] <- modifiedRight[counter]
    addCounter<- addCounter+1
  }
}



input1Processed <- listChecker(input1,"Gene_1_symbol(5end_fusion_partner)","Gene_2_symbol(3end_fusion_partner)","/home/nick/Desktop/Fusion_prioritization/Data/Processed/fuse_geneAndFusion.txt")
input2Processed <- listChecker(input2,"LeftGeneSymbol","RightGeneSymbol","/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_geneAndFusion.txt")