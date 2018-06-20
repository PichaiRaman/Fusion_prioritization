# Fusion_prioritization


Pipeline to prioritize gene fusions coming from the STAR-Fusion workflow. Data can be found here


https://drive.google.com/drive/folders/1iZ1on1exCrN-8k-5C1rC8QnCoG7wmcS2

CancerGeneCompilation.R - Composed CancerGeneList.txt which consists of a list of cancerous genes and a count of the studies in which they were found. The studies can be found in the comments of CancerGeneCompilation.R

basicGeneChecker.R - Takes input from star fusion data and then returns that data with two extra rows (Cancerous and CancerousGene). Based on CancerGeneList.txt

fusionCollector.R - Composed FusionList.txt, a compilation of cancerous fusions and the a count of the studies in which they were found. The studies can be found in the comments of fusionCollector.R

GeneAndFusionChecker.R - Uses FusionList.txt and CancerGeneList.txt to process an input of left and right genes of a fusion to produce an output file with added columns of whether one or two of the partner genes are canerous and whether the fusion is cancerous.

