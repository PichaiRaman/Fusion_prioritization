#####################################################
#Purpose: Obtain PFAM ID's associated with output from Fusion Prioritization Pipeline       
#Author: Pichai Raman
#Date: June 25, 2018
#####################################################


#Support function to collapse a tall skinny list
#so that data is short and comma separated (serialized)
collapseData <- function(geneDomain=NULL)
{
	out <- geneDomain %>% group_by(hgnc_symbol) %>% summarise(PFAM=paste(pfam, collapse=","))
	out <- data.frame(out);
	return(out);
}


#Input's 
getPFAMDomain <- function(starFusionCombo=NULL)
{
	#Call library
	require("biomaRt")
	require("tidyverse")

	#Choose Mart & Database
	ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

	#Get data for H Gene
	hGenes <- as.character(starFusionCombo[,"H_Gene"]);
	pfamData = getBM(attributes=c("hgnc_symbol", "pfam" ,"pfam_start", "pfam_end", "transcript_start"), 
             filters="hgnc_symbol",
             values=hGenes, 
             mart=ensembl)
	pfamData <- na.omit(pfamData);
	
	#Add all PFAM Domains
	hGeneAllDomains <- collapseData(pfamData[,c("hgnc_symbol", "pfam")]);
	colnames(hGeneAllDomains)[2] <- "H_Gene_PFAM_All"
	starFusionCombo <- merge(starFusionCombo, hGeneAllDomains, by.x="H_Gene", by.y="hgnc_symbol", all.x=T);

	#Get data for T Gene
	tGenes <- as.character(starFusionCombo[,"T_Gene"]);
	pfamData = getBM(attributes=c("hgnc_symbol", "pfam" ,"pfam_start", "pfam_end", "transcript_start"), 
             filters="hgnc_symbol",
             values=tGenes, 
             mart=ensembl)
	pfamData <- na.omit(pfamData);
	
	#Add all PFAM Domains
	tGeneAllDomains <- collapseData(pfamData[,c("hgnc_symbol", "pfam")]);
	colnames(tGeneAllDomains)[2] <- "T_Gene_PFAM_All"
	starFusionCombo <- merge(starFusionCombo, tGeneAllDomains, by.x="T_Gene", by.y="hgnc_symbol", all.x=T);
	return(starFusionCombo);

}


#starFusionCombo2 <- getPFAMDomain(starFusionCombo)


############NOT USED###########################
#starFusionCombo <- read.delim("../data/processed/star_fusion_combo.tsv");

#Call libraries

#myGene <- "TP53"
#myGene2 <- "PIK3CA"
#myGene3 <- "DDX3X"
#myGenes <- c(myGene, myGene2, myGene3);



	#Add subset of the PFAm Domains
#	geneStartSite <- getBM(attributes=c("hgnc_symbol", "start_position"), 
#             filters="hgnc_symbol",
#             values=hGenes, 
#             mart=ensembl)
#	tmpFusionDat <- merge(starFusionCombo, geneStartSite, by.x="H_Gene", by.y="hgnc_symbol", all.x=T);
#	tmpFusionDat[,"AA"] <- round((tmpFusionDat[,"Left_Breakpoint_Pos"]-tmpFusionDat[,"start_position"])/3);
#	pfamData <- merge(pfamData, starFusionCombo[,""])


