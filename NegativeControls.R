#Negative Control for Pathway Enrichement Analysis

#-------------------------------------------------------------------------------------------------------------------
#                Random Genes - controls for no pathway passing the FDR (makes sure 
#                        there would NOT be any enriched pathway)
#-------------------------------------------------------------------------------------------------------------------


library(rtracklayer)
library(dplyr)

gtf <- import('~/Desktop/Mus_musculus_wsbeij.WSB_EiJ_v1.113.gtf')
head(gtf)

genes <- gtf[gtf$type=='gene']

colnames(mcols(gtf))

head(genes)


gene_id <- unique(mcols(genes)$projection_parent_gene)
gene_id
set.seed(42) #to enable my random results reproducible 
random_genes <- sample(gene_id, 2978 )
gene_ids_no_version <- sub("\\.\\d+$", "", random_genes)
length(gene_ids_no_version)
gene_ids_no_version
write.table(gene_ids_no_version, "random_mouse_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


allgenes <- read.table('AllGenes.txt')
nrow(allgenes)

#-------------------------------------------------------------------------------------------------------------------
#                7 Not upregulated miRNAs - controls for MAPK not being enriched 
#-------------------------------------------------------------------------------------------------------------------

negativecontrol <- readRDS('~/Desktop/NegativeControlToLoad.RDS')
negativecontrol
length(negativecontrol)

set.seed(44)
random_miRNAs <- sample(negativecontrol, 7 )
random_miRNAs

#[1] "mmu-miR-5119"      "mmu-miR-125b-1-3p" "mmu-miR-218-5p"    "mmu-miR-18a-5p"    "mmu-miR-671-5p"   
#[6] "mmu-miR-5126"      "mmu-miR-425-5p"   

#got these into the DIANA with the same filters I applied for the initial part of my experiment
#gave me a list of 2721 possible gene targets which is similar to the number of the genes i got for my first analysis which were 2978 genes
library(readxl)
GeneTarget <- read_excel('~/Desktop/SevenMiRNAsNegativeControl.xlsx')
GeneTarget
GeneTarget$gene_ensembl_id
unlist(GeneTarget$gene_ensembl_id)

GeneTarget$gene_ensembl_id <- gsub('"', '', GeneTarget$gene_ensembl_id)
NegativeControl7miRNAs <- GeneTarget$gene_ensembl_id

NegativeControl7miRNAs
class(NegativeControl7miRNAs)

write.table(NegativeControl7miRNAs, file = "~/Desktop/NegativeControl7miRNAs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
