library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(biomaRt)

source("/home/nick/Desktop/Fusion_prioritization/code/getPFAMDomain.R")

geneFuseChecker <- function(input,geneRep=1,fuseRep=1) {
  geneList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/CancerGeneList.txt")
  fuseList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt")
  
  output <- input
  output$Cancerous_Gene <- F
  output$Cancerous_Gene_Symbol <- NA
  output$Cancerous_Fusion <- F
  output$Cancerous_Fusion_Symbol <- NA
  
  repAdjustedGene <- geneList[geneList$Count >= geneRep,]
  repAdjustedFuse <- fuseList[fuseList$Count >= geneRep,]
  
  for (rowIndex in 1:nrow(input)) {
    gene1 <- input[rowIndex,"H_Gene"]
    gene2 <- input[rowIndex,"T_Gene"]
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
  
  return (output)
}

domainChecker <- function(starFusionCombo) {
  starFusionCombo$Cancer_Domains <- NA
  domainList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/Domain_List.tsv")
  rowIndex <- 1
  for (doms in starFusionCombo$H_Gene_PFAM_All) {
    #print(doms[1])
    doms <- str_split(doms,",")

    for (dom in doms[[1]]) {

      if (is.na(dom)) {
        next
      }
      
      if (dom %in% domainList$PF) {
       # print(dom)
        if (is.na(starFusionCombo[rowIndex,"Sig_Domains"] )) {
          starFusionCombo[rowIndex,"Sig_Domains"] <- dom
        }else {
          starFusionCombo[rowIndex,"Sig_Domains"] <- paste(starFusionCombo[rowIndex,"Sig_Domains"],dom,sep=",")  
        }
      }
    }
    rowIndex<-rowIndex+1
  }
  rowIndex<-1
  for (doms in starFusionCombo$T_Gene_PFAM_All) {
      doms <- str_split(doms,",")
      #print(doms[1])
      for (dom in doms[[1]]) {
        if (is.na(dom)) {
          next
        }
        
        if (dom %in% domainList$PF) {
          #print(dom)
          if (is.na(starFusionCombo[rowIndex,"Sig_Domains"] )) {
            starFusionCombo[rowIndex,"Sig_Domains"] <- dom
          }else {
            starFusionCombo[rowIndex,"Sig_Domains"] <- paste(starFusionCombo[rowIndex,"Sig_Domains"],dom,sep=",")  
          }
        }
      }
    rowIndex<-rowIndex+1
  }
  
  
  for (i in 1:nrow(starFusionCombo)) {
    if (!is.na(starFusionCombo[i, "Sig_Domains"])) {
      list <- str_split(starFusionCombo[i,"Sig_Domains"],",")[[1]]
      starFusionCombo[i, "Sig_Domains"] <- paste(levels(factor(list)),collapse=",")
    }
  }
  
  return(starFusionCombo)
}

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_combo.tsv")
input <- getPFAMDomain(input)
input <- getPFAMDomain(input) #Need to fix this, have to run twice in order for it to work
input <- as.tibble(input)
input <- subset(input, select=-c(H_Gene_PFAM_All.x,T_Gene_PFAM_All.x))
colnames(input)[colnames(input)=="H_Gene_PFAM_All.y"] <- "H_Gene_PFAM_All"
colnames(input)[colnames(input)=="T_Gene_PFAM_All.y"] <- "T_Gene_PFAM_All"

input <- domainChecker(input)

input <- geneFuseChecker(input)

input$Targetable <- FALSE
input$Targetable_gene <- NA
input$Drug_name <- NA
input$Drug_chembl_id <- NA

Target_genes <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/Target_List.tsv")

#Making Target Columns
for (i in 1:nrow(input)) {
  if (input$H_Gene[i] %in% Target_genes$Trgt_Genes) {
    input$Targetable[i] <- TRUE
    
    index <- Target_genes$Trgt_Genes == input$H_Gene[i]
    
    input$Targetable_gene[i] <- input$H_Gene[i]
    input$Drug_name[i] <- Target_genes$Drug_name[index]
    input$Drug_chembl_id[i] <- Target_genes$Drug_chembl_id[index]
  }
  
  if (input$T_Gene[i] %in% Target_genes$Trgt_Genes) {
    input$Targetable[i] <- TRUE
    
    index <- Target_genes$Trgt_Genes == input$T_Gene[i]
    
    if (is.na(input$Targetable_gene[i])) {
      input$Targetable_gene[i] <- input$T_Gene[i]
      input$Drug_name[i] <- Target_genes$Drug_name[index]
      input$Drug_chembl_id[i] <- Target_genes$Drug_chembl_id[index]
    }else {
      input$Targetable_gene[i] <- paste(input$Targetable_gene[i],input$T_Gene[i],sep=",")
      input$Drug_name[i] <- paste(input$Drug_name[i],Target_genes$Drug_name[index],sep=",")
      input$Drug_chembl_id[i] <- paste(input$Drug_chembl_id[i],Target_genes$Drug_chembl_id[index],sep=",")
    }
    
  }
}


