library(readr)
library(stringr)

input <- read_tsv("../Fusion_prioritization/star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv")
geneList <- geneList <-read_tsv("../Fusion_prioritization/CancerGeneList.txt")

output <- input
output$Cancerous <- F
output$CancerousGene <- NA

geneRep <- 1
repAdjustedList <- geneList[geneList$Count >= geneRep,]

rowIndex = 1
for (fuse in input$'#FusionName') {
  dash_locus <- str_locate(fuse,"--")
  gene1 <- substr(fuse, 0, dash_locus[1,1]-1)
  gene2 <- substr(fuse, dash_locus[1,2]+1,str_length(fuse))
  
  gene1In <- gene1 %in% repAdjustedList$Gene
  gene2In <- gene2 %in% repAdjustedList$Gene
  if (gene1In && gene2In) {
    output[rowIndex,"Cancerous"] <- T
    output[rowIndex,"CancerousGene"] <- paste(gene1,"--",gene2)
  }else if (gene1In) {
    output[rowIndex,"Cancerous"] <- T
    output[rowIndex,"CancerousGene"] <- paste(gene1)
  }else if (gene2In) {
    output[rowIndex,"Cancerous"] <- T
    output[rowIndex,"CancerousGene"] <- paste(gene2)
  }
  rowIndex <- rowIndex+1
}

View(output)