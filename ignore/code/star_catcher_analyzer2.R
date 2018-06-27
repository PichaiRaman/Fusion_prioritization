library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(biomaRt)

source("/home/nick/Desktop/Fusion_prioritization/code/getPFAMDomain.R")

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_combo.tsv")

#Cancer Gene and Fusion Prep
geneList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/CancerGeneList.txt")
fuseList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/FusionList2.txt")

input$Cancerous_Gene <- F
input$Cancerous_Gene_Symbol <- NA
input$Cancerous_Fusion <- F
input$Cancerous_Fusion_Symbol <- NA

geneRep <- 1
fuseRep <- 1
repAdjustedGene <- geneList[geneList$Count >= geneRep,]
repAdjustedFuse <- fuseList[fuseList$Count >= fuseRep,]

#Target Prep
input$Targetable <- FALSE
input$Targetable_gene <- NA
input$Drug_name <- NA
input$Drug_chembl_id <- NA

Target_genes <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/Target_List.tsv")

#PFAM Domain
input <- getPFAMDomain(input)
input <- getPFAMDomain(input) #Need to fix this, have to run twice in order for it to work
input <- as.tibble(input)
input <- subset(input, select=-c(H_Gene_PFAM_All.x,T_Gene_PFAM_All.x))
colnames(input)[colnames(input)=="H_Gene_PFAM_All.y"] <- "H_Gene_PFAM_All"
colnames(input)[colnames(input)=="T_Gene_PFAM_All.y"] <- "T_Gene_PFAM_All"

#Domain Prep
input$Cancer_Domains <- NA
domainList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/Domain_List.tsv")

#Big For Loop
for (i in 1:nrow(input)) {
  #Cancer G+F
  gene1 <- input[i,"H_Gene"]
  gene2 <- input[i,"T_Gene"]
  fuse <- paste(gene1,gene2,sep="--")
  
  gene1In <- gene1 %in% repAdjustedGene$Gene
  gene2In <- gene2 %in% repAdjustedGene$Gene
  fuseIn <- fuse %in% repAdjustedFuse$Fusions
  
  input[i,"Cancerous_Gene"] <- gene1In || gene2In
  input[i,"Cancerous_Fusion"] <- fuseIn
  
  if (gene1In && gene2In) {
    input[i,"Cancerous_Gene_Symbol"] <- paste(gene1,gene2,sep=",")
  }else if (gene1In) {
    input[i,"Cancerous_Gene_Symbol"] <- gene1
  }else if (gene2In) {
    input[i,"Cancerous_Gene_Symbol"] <- gene2
  }
  
  if (fuseIn) {
    output[i,"Cancerous_Fusion_Symbol"] <- paste(gene1,gene2,sep="--")
  }
  
  #Target
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
  
  #Domain Checker
  doms <- str_split(input$H_Gene_PFAM_All[i],",")
  
  for (dom in doms[[1]]) {
    
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      # print(dom)
      if (is.na(input[i,"Cancer_Domains"] )) {
        input[i,"Cancer_Domains"] <- dom
      }else {
        input[i,"Cancer_Domains"] <- paste(input[i,"Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  doms <- str_split(input$T_Gene_PFAM_All[i],",")

  for (dom in doms[[1]]) {
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      #print(dom)
      if (is.na(input[i,"Cancer_Domains"] )) {
        input[i,"Cancer_Domains"] <- dom
      }else {
        input[i,"Cancer_Domains"] <- paste(input[i,"Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  if (!is.na(input[i, "Cancer_Domains"])) {
    list <- str_split(input[i,"Cancer_Domains"],",")[[1]]
    input[i, "Cancer_Domains"] <- paste(levels(factor(list)),collapse=",")
  }
}


