library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(biomaRt)


source("/home/nick/Desktop/Fusion_prioritization/code/getPFAMDomain.R")

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
input <- read.delim("/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_combo.tsv")

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
#input <- getPFAMDomain(input) #Need to fix this, have to run twice in order for it to work
# input <- as.tibble(input)
# input <- subset(input, select=-c(H_Gene_PFAM_All.x,T_Gene_PFAM_All.x))
# colnames(input)[colnames(input)=="H_Gene_PFAM_All.y"] <- "H_Gene_PFAM_All"
# colnames(input)[colnames(input)=="T_Gene_PFAM_All.y"] <- "T_Gene_PFAM_All"

#Domain Prep
input$Cancer_Domains <- NA
domainList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/Processed/Domain_List.tsv")

#Big Apply Function
bigApply <- function(x) {
  #Cancer G+F
  gene1 <- x["H_Gene"]
  gene2 <- x["T_Gene"]
  fuse <- paste(gene1,gene2,sep="--")
  
  gene1In <- gene1 %in% repAdjustedGene$Gene
  gene2In <- gene2 %in% repAdjustedGene$Gene
  fuseIn <- fuse %in% repAdjustedFuse$Fusions
  
  x["Cancerous_Gene"] <- gene1In || gene2In
  x["Cancerous_Fusion"] <- fuseIn
  
  if (gene1In && gene2In) {
    x["Cancerous_Gene_Symbol"] <- paste(gene1,gene2,sep=",")
  }else if (gene1In) {
    x["Cancerous_Gene_Symbol"] <- gene1
  }else if (gene2In) {
    x["Cancerous_Gene_Symbol"] <- gene2
  }
  
  if (fuseIn) {
    x["Cancerous_Fusion_Symbol"] <- paste(gene1,gene2,sep="--")
  }
  
  #Target
  if (x["H_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["H_Gene"]
    
    x["Targetable_gene"] <- x["H_Gene"]
    x["Drug_name"] <- Target_genes$Drug_name[index]
    x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
  }
  
  if (x["T_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["T_Gene"]
    
    if (is.na(x["Targetable_gene"])) {
      x["Targetable_gene"] <- x["T_Gene"]
      x["Drug_name"] <- Target_genes$Drug_name[index]
      x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
    }else {
      x["Targetable_gene"] <- paste(x["Targetable_gene"],x["T_Gene"],sep=",")
      x["Drug_name"] <- paste(x["Drug_name"],Target_genes$Drug_name[index],sep=",")
      x["Drug_chembl_id"] <- paste(x["Drug_chembl_id"],Target_genes$Drug_chembl_id[index],sep=",")
    }
    
  }
  
  #Domain Checker
  doms <- c()
  domList <- str_split(unlist(str_split(x["H_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      # print(dom)
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  doms <- c()
  domList <- str_split(unlist(str_split(x["T_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      #print(dom)
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  if (!is.na(x["Cancer_Domains"])) {
    list <- str_split(x["Cancer_Domains"],",")[[1]]
    x["Cancer_Domains"] <- paste(levels(factor(list)),collapse=",")
  }
  return(x)
}

input <- apply(input, MARGIN = 1, FUN = bigApply)
input <- t(input)
input <- as.tibble(input)
input <- input[,c(5,3,1,6,4,7,8,2,9,as.numeric(10:31))]

write.table(input, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_analyzed.tsv", sep="\t", row.names = FALSE)
