library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

#Input from Star Catcher must be named "star_data.tsv" and that from Fusion Catcher "fusion_catcher_data.tsv"
#star_input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/star_data.tsv")
#catcher_input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/fusion_catcher_data.tsv")
star_input <- star.fusion
catcher_input <- fusion.catcher

#star_names <- c("#FusionName","RightBreakpoint","LeftBreakpoint","JunctionReadCount","SpanningFragCount","FUSION_CDS","FUSION_TRANSL","PROT_FUSION_TYPE")
#catcher_names <- c("Gene_1_symbol(5end_fusion_partner)","Gene_2_symbol(3end_fusion_partner)","Counts_of_common_mapping_reads","Spanning_pairs","Spanning_unique_reads","Fusion_sequence","Predicted_fused_transcripts","Predicted_effect")
colnames(star_input) <- c('Gene1','Gene2','#FusionName','JunctionReadCount','SpanningFragCount','SpliceType','LeftGene','chrA','RightGene','chrB','LargeAnchorSupport','LeftBreakpoint','LeftBreakEntropy','RightBreakpoint','RightBreakEntropy','FFPM','annots','CDS_LEFT_ID','CDS_LEFT_RANGE','CDS_RIGHT_ID','CDS_RIGHT_RANGE','PROT_FUSION_TYPE','FUSION_MODEL','FUSION_CDS','FUSION_TRANSL','PFAM_LEFT','PFAM_RIGHT','Sample','Sample_Type')
colnames(catcher_input) <- c('Gene_1_symbol(5end_fusion_partner)','Gene_2_symbol(3end_fusion_partner)','Fusion_description','Counts_of_common_mapping_reads','Spanning_pairs','Spanning_unique_reads','Longest_anchor_found','Fusion_finding_method','chrA','chrB','Fusion_point_for_gene_1(5end_fusion_partner)','Fusion_point_for_gene_2(3end_fusion_partner)','Exon_1_id(5end_fusion_partner)','Exon_2_id(3end_fusion_partner)','Fusion_sequence','Predicted_effect','Predicted_fused_transcripts','Predicted_fused_proteins','Sample','Sample_Type','left_seq','right_seq','perc.identity','Low_complexity')

#Cleaning star BreakPoints
r_raw <- unlist(str_split(star_input$RightBreakpoint,"chr"))
l_raw <- unlist(str_split(star_input$LeftBreakpoint,"chr"))
rowCounter<-1
for (i in 1:length(r_raw)) {
  if (i %% 2 == 0) {
    star_input[rowCounter,"RightBreakpoint"] = r_raw[i]
    rowCounter<-rowCounter+1
  }
}
rowCounter<-1
for (i in 1:length(l_raw)) {
  if (i %% 2 == 0) {
    star_input[rowCounter,"LeftBreakpoint"] = l_raw[i]
    rowCounter<-rowCounter+1
  }
}

#Cleaning Catcher Fusion_Type
catcher_input[catcher_input$Predicted_effect == "in-frame","Predicted_effect"] <- "INFRAME"
catcher_input[catcher_input$Predicted_effect == "out-of-frame","Predicted_effect"] <- "FRAMESHIFT"
catcher_input$Predicted_effect <- toupper(catcher_input$Predicted_effect)


star_fusions <- as.character(star_input$`#FusionName`)
catcher_fusions <- as.character(paste(catcher_input$`Gene_1_symbol(5end_fusion_partner)`,catcher_input$`Gene_2_symbol(3end_fusion_partner)`,sep="--"))
fuseList <- unique(c(star_fusions,catcher_fusions))

#output <- tibble(Fusion_Name=levels(factor(c(levels(star_fusions),levels(catcher_fusions)))))

#sigColNames_star <- c("LeftBreakpoint","RightBreakpoint","JunctionReadCount","SpanningFragCount","FUSION_CDS","FUSION_TRANSL","PROT_FUSION_TYPE")

output <- tibble(Fusion_Name=NA)
output$H_Gene <- NA
output$T_Gene <-NA
output$Left_Breakpoint <-NA
output$Left_Breakpoint_Chr <-NA
output$Left_Breakpoint_Pos <-NA
output$Left_Breakpoint_Str <-NA
output$Right_Breakpoint<- NA
output$Right_Breakpoint_Chr<- NA
output$Right_Breakpoint_Pos<- NA
output$Right_Breakpoint_Str<- NA
output$Star_Junction_Readcount <- NA
output$Star_Spanning_Fragcount <- NA
output$FC_Common_Mapping_Reads <- NA
output$FC_Spanning_Pairs <- NA
output$FC_Spanning_Unique_Reads <- NA
output$Fusion_Transcript_Sequence <- NA
output$Fusion_Protein_Sequence <- NA
output$Fusion_Type <- NA

#The Big For Loop: Loops through each fusion and checks for repeats
counter<-1
for (fuse in fuseList) {
  
  if (fuse %in% star_fusions && fuse %in% catcher_fusions) {
    
    index_star <- which(star_fusions %in% fuse)
    index_catcher <- which(catcher_fusions %in% fuse)
    if (length(index_star)>1 || length(index_catcher)>1) {
      #Checking for duplicates in Star
      sequences_star <- star_input[index_star[1:length(index_star)],"FUSION_TRANSL"]
      index_star <- index_star[!duplicated(sequences_star)]
      
      #Checking for duplicates in Catcher
      sequences_catcher <- catcher_input[index_catcher[1:length(index_catcher)],"Fusion_sequence"]
      index_catcher <- index_catcher[!duplicated(sequences_catcher)]
      for (i in index_star) {
        newRowNum <-nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- star_input[i,"LeftBreakpoint"]
        output[newRowNum,"Right_Breakpoint"] <- star_input[i,"RightBreakpoint"]
        output[newRowNum,"Star_Junction_Readcount"] <- star_input[i,"JunctionReadCount"]
        output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[i,"SpanningFragCount"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[i,"FUSION_CDS"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[i,"FUSION_TRANSL"]
        output[newRowNum,"Fusion_Type"] <- star_input[i,"PROT_FUSION_TYPE"]
        output[newRowNum,"Catcher"] <- "Star"
      }
      
      for (i in index_catcher) {
        newRowNum <-nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_1(5end_fusion_partner)"]
        output[newRowNum,"Right_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_2(3end_fusion_partner)"]
        output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[i,"Counts_of_common_mapping_reads"]
        output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[i,"Spanning_pairs"]
        output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[i,"Spanning_unique_reads"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[i,"Fusion_sequence"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[i,"Predicted_fused_transcripts"]
        output[newRowNum,"Fusion_Type"] <- catcher_input[i,"Predicted_effect"]
        output[newRowNum,"Catcher"] <- "Fusion_Catcher"
      }
    }else {
      newRowNum <-nrow(output)+1
      locus <- index[1]
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- star_input[locus,"LeftBreakpoint"]
      output[newRowNum,"Right_Breakpoint"] <- star_input[locus,"RightBreakpoint"]
      output[newRowNum,"Star_Junction_Readcount"] <- star_input[locus,"JunctionReadCount"]
      output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[locus,"SpanningFragCount"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[locus,"FUSION_CDS"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[locus,"FUSION_TRANSL"]
      output[newRowNum,"Fusion_Type"] <- star_input[locus,"PROT_FUSION_TYPE"]
      output[newRowNum,"Catcher"] <- "Star, Fusion_Catcher"
    }
    
  }else if (fuse %in% star_fusions) {
    
    index <- which(star_fusions %in% fuse)
    
    if (length(index)>1) {
      sequences <- star_input[index[1:length(index)],"FUSION_TRANSL"]
      index <- index[!duplicated(sequences)]
      for (i in index) {
        newRowNum <- nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- star_input[i,"LeftBreakpoint"]
        output[newRowNum,"Right_Breakpoint"] <- star_input[i,"RightBreakpoint"]
        output[newRowNum,"Star_Junction_Readcount"] <- star_input[i,"JunctionReadCount"]
        output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[i,"SpanningFragCount"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[i,"FUSION_CDS"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[i,"FUSION_TRANSL"]
        output[newRowNum,"Fusion_Type"] <- star_input[i,"PROT_FUSION_TYPE"]
        output[newRowNum,"Catcher"] <- "Star"
      }
    }else {
      newRowNum <-nrow(output)+1
      locus <- index[1]
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- star_input[locus,"LeftBreakpoint"]
      output[newRowNum,"Right_Breakpoint"] <- star_input[locus,"RightBreakpoint"]
      output[newRowNum,"Star_Junction_Readcount"] <- star_input[locus,"JunctionReadCount"]
      output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[locus,"SpanningFragCount"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[locus,"FUSION_CDS"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[locus,"FUSION_TRANSL"]
      output[newRowNum,"Fusion_Type"] <- star_input[locus,"PROT_FUSION_TYPE"]
      output[newRowNum,"Catcher"] <- "Star"
    }
  }else {
    
    index <- which(catcher_fusions %in% fuse)
    if (length(index)>1) {
      sequences <- catcher_input[index[1:length(index)],"Fusion_sequence"]
      index <- index[!duplicated(sequences)]
      for (i in index) {
        newRowNum <- nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_1(5end_fusion_partner)"]
        output[newRowNum,"Right_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_2(3end_fusion_partner)"]
        output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[i,"Counts_of_common_mapping_reads"]
        output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[i,"Spanning_pairs"]
        output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[i,"Spanning_unique_reads"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[i,"Fusion_sequence"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[i,"Predicted_fused_transcripts"]
        output[newRowNum,"Fusion_Type"] <- catcher_input[i,"Predicted_effect"]
        output[newRowNum,"Catcher"] <- "Fusion_Catcher"
      }
    }else {
      newRowNum <-nrow(output)+1
      locus <- index[1]
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- catcher_input[locus,"Fusion_point_for_gene_1(5end_fusion_partner)"]
      output[newRowNum,"Right_Breakpoint"] <- catcher_input[locus,"Fusion_point_for_gene_2(3end_fusion_partner)"]
      output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[locus,"Counts_of_common_mapping_reads"]
      output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[locus,"Spanning_pairs"]
      output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[locus,"Spanning_unique_reads"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[locus,"Fusion_sequence"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[locus,"Predicted_fused_transcripts"]
      output[newRowNum,"Fusion_Type"] <- catcher_input[locus,"Predicted_effect"]
      output[newRowNum,"Catcher"] <- "Fusion_Catcher"
    }
  }
  if (counter %% 1000 == 0) {
    print(paste(counter,"/",length(fuseList)))
  }
  counter<-counter+1
}

#Creating Left and Right Gene Vars
fusions_split <- unlist(str_split(output$Fusion_Name,"--"))
rowCounter <- 1
for (i in 1:length(fusions_split)) {
  if (i%%2==0) {
    output[rowCounter,"H_Gene"] <- fusions_split[i]
  }else {
    output[rowCounter,"T_Gene"] <- fusions_split[i]
    rowCounter<-rowCounter+1
  }
}
output <- output[-1,]

#Cleaning up "."s and changing them to NA
dotClean <- function(x) {
  loci <- str_locate(x,"\\.")[,1] == 1
  loci[is.na(loci)] <- FALSE
  x[!is.na(str_locate(x,"\\."))[,1] & loci] <- NA
  return (x)
}
output <- apply(output,FUN=dotClean,MARGIN=2)

output<-as.tibble(output)

#Adding Breakpoint Chromosone Position and Strand Columns
LBP_List <-str_split(output$Left_Breakpoint,":")
rowCounter<-1
for (item in LBP_List) {
  output[rowCounter,"Left_Breakpoint_Chr"] <- item[1]
  output[rowCounter,"Left_Breakpoint_Pos"] <- item[2]
  output[rowCounter,"Left_Breakpoint_Str"] <- item[3]
  rowCounter<-rowCounter+1
}

RBP_List <-str_split(output$Right_Breakpoint,":")
rowCounter<-1
for (item in RBP_List) {
  output[rowCounter,"Right_Breakpoint_Chr"] <- item[1]
  output[rowCounter,"Right_Breakpoint_Pos"] <- item[2]
  output[rowCounter,"Right_Breakpoint_Str"] <- item[3]
  rowCounter<-rowCounter+1
}
View(output)
output <- subset(output, select=-c(Right_Breakpoint,Left_Breakpoint))

write.table(output, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_combo.tsv", sep="\t", row.names = FALSE)