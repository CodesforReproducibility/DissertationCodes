#Processing mouse GTF annotation

library(rtracklayer)
library(GenomicRanges)

#These files came straight from the Gencode Website (this has all of the genome annotation)
gtfFile <- 'gencode.vM36.primary_assembly.annotation.gtf.gz'
gtfFile_tRNAs <- 'gencode.vM36.tRNAs.gtf.gz'




#Im converting this file into a Bed because I need it as a bed file to use liftover on it
library(rtracklayer)

# Import miRBase GFF3 file (replace with the correct filename)
mirbase_gff <- import("mmu.gff3")
mirbase_gff

mirbase_gff <- mirbase_gff[mcols(mirbase_gff)$type == "miRNA"] #of course I have to remove primary transcripts duh 



mirbase_bed <- data.frame(
  Chr = as.character(seqnames(mirbase_gff)),  # Chromosome
  Start = start(mirbase_gff) - 1,  # Convert to 0-based for BED format
  End = end(mirbase_gff),  # End stays the same
  Name = mirbase_gff$Name,  # miRNA Name
  Score = rep(".", length(mirbase_gff)),  # Use "." instead of "0"
  Strand = as.character(strand(mirbase_gff))  # Include strand info
)

# Save BED file
write.table(mirbase_bed, "miRBase_mm10.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, na=".")

file_miRNAs2 <- 'miRBase_mm39.bed' 

#now we import all the data
gr <- import(gtfFile) #this is in a GRanges object
trnas <- import(gtfFile_tRNAs)
matMirs2 <- import(file_miRNAs2)

gr #normal GR object 
trnas #again another normal GR object 
matMirs2 #this one too but the I have fixed the name check the code!


#Selecting the columns that are useful for my file (this applied to the gr and then to the trnas)
mcols(gr) <- mcols(gr)[,c(2,6,7)]
gr

mcols(trnas) <- mcols(trnas)[,c(2,7,8)]
trnas

#extracting all the miRNA annotation from the GENCODE annotation file
gr
table(gr$gene_type) #but the thing is we can just take out miRNA because it might be redudant 

gr[gr$gene_name == 'Mir122'] 
#it tells me that theres a gene, a transcript and an exon (redudant)
#we should really only use one (i guess in this case lets go with gene cause mature miRNA comes from gene like coding like for sure - more else but i think we can assume this)


mirs <- gr[gr$gene_type=='miRNA' & gr$type == 'gene'] #omg this is a fun new way to write it slay add a little &
mirs


#maybe is there still some redundacy?
sum(width(reduce(mirs)))
sum(width(mirs))

overlapWidth <- sum(width(mirs)) - sum(width(reduce(mirs)))
#there are 64 positions where there is some overlap which means that they are redundant we wanna remove it 
#how to do it? new function lets find where these elements actually overlap

overPos <- findOverlaps(mirs, mirs)
overPos
length(overPos) #so here we get the majority of these reads mapping to the same position and overlapping with itself which im guessing is like self overlaps so not really helpful either way its like this is inflated
#basically we have self-overlaps - lets select where the positions arent equal

overPosNotEq <- overPos[queryHits(overPos) != subjectHits(overPos)] #with this we have selected the positions that are not equal
mirsOver <- mirs[unique(queryHits(overPosNotEq))] #unique positions that are not equal and we are trying to determine them
mirsOver #this gives us three miRNAs that overlap with each other - we can either remove them (removed true overlaps) or reduce them (remove redudancy but keep them) 

#and we just confirm that these are the same we found before which was that there was an overlap of 64 nt
sum(width(mirsOver))-sum(width(reduce(mirsOver))) 
#also side note I have not removed the redundant miRNA btw - because like its still inside of my mirs but i know that mirsOver is only the overlaping 


#Something to note - we only have the coordinates for miRNA hairpin and not for the individual mature miRNAs 
hist(width(mirs), breaks = 20, xlim = c(0,200))
#see the width lies around 70-100 nt (this is not miRNA)

#looking at the longest ones - they have a GM gene name which suggest these are all predicted miRNAs - so they are not like the ones we want!
mirs[width(mirs)>170]

#what we actually want is the ones that have a name starting with Mir cause like those are true miRNAs i think
mirs[grep('Mir', mirs$gene_name)]


#getting mature miRNA coordinates from miRBase 
matMirs2 

barplot(table(width(matMirs2)), main = 'miRNA sizes in miRBase', xlab='nucleotides') #honestly size distribution looks fine all at the 21-22

#check for overlaps
sum(width(matMirs2)) - sum(width(reduce(matMirs2)))


#So now we will want to removed the miRNAs that we have collected - but which have overlaps with these miRBase mature miRNA
mirs <- subsetByOverlaps(mirs, matMirs2, invert = TRUE)

#this is combining out two miRNA GRanges metadata but we need to make it compatible
mirs #(ive lost 13 sequences)
matMirs2


#so whatever is different it needs to be arranged
matMirs2$gene_name <- matMirs2$name
matMirs2$type <- 'mature_miRNA' #look ive added new columns
matMirs2$gene_type <- 'miRNA'
matMirs2

mcols(matMirs2) <- mcols(matMirs2)[,c('type', 'gene_type', 'gene_name')]
mirs <- c(matMirs2, mirs)
mirs

#other ncRNA of interest
#we are gonna be focusing on ribosomal rRNA because it can be frequently found on the results of sequencing
#in the samples where RNA is more degraded we will have more maps reading to rRNA
gr

rrnas <- gr[gr$gene_type == 'rRNA']
rrnas # we also see the same 3 fold redundancy

rrnas <- rrnas[rrnas$type == 'gene'] #so we pick the gene only 
rrnas

#lets check for redudancy 
sum(width(rrnas)) - sum(width(reduce(rrnas)))
#zero redundacy so ready to use


#there is a possibility that rRNAs overlap with miRNAs and in the case that they do exist we should want to remove it
overMR <- findOverlaps(mirs, rrnas)
overMR
#okay so the next step is to remove them:
mirs[queryHits(overMR)] #so we first extract them 

mirs <- mirs[-queryHits(overMR)] #and then we remove them
#just to check we have removed them:
findOverlaps(mirs, rrnas)


#there are other kinds of ncRNA like tRNAs
overMT <- findOverlaps(mirs, trnas)
overMT

mirs[queryHits(overMT)]
#i get 9 sequences
trnas[subjectHits(overMT)]
#in this i can discard either but too keep as many miRNAs possible lets keep the mirs
trnas <- trnas[-subjectHits(overMT)]
findOverlaps(mirs, trnas) #perf its removed!!

#we also should check the overlaps between tRNAs and rRNAs
overRT <- findOverlaps(rrnas, trnas)
overRT
trnas <- trnas[-subjectHits(overRT)]

#protein coding genes lets consider them as a general category
prot <- gr[gr$gene_type == 'protein_coding']
prot #theres clearly a lot of redundancy and also we dont intros cause they cannot produce a mRNA

prot <- prot[prot$type == 'exon']
prot

#what we will be doing is removing redundacy with the reduce cause like its like this is protein coding so like not sooo important maybe
sum(width(prot))
protNR <- reduce(prot)
sum(width(protNR)) #reduction mostly due tot he all the redundant exons present

#to recover some of the metadata:
origPos <- findOverlaps(protNR, prot, select = 'first')
protNR$gene_name <- prot[origPos]$gene_name
protNR

#okay so in theory they are all like reduced by like in themselves as in no reduction in the data but they might overlap with the rest of our data:
findOverlaps(protNR, mirs) 
findOverlaps(protNR, rrnas)
findOverlaps(protNR, trnas)

#so there are quite a few overlaps so we need to remove them:
protNR <- subsetByOverlaps(protNR, mirs, invert = TRUE)
protNR <- subsetByOverlaps(protNR, rrnas, invert = TRUE)
protNR <- subsetByOverlaps(protNR, trnas, invert = TRUE)

#just to be sure it worked:
findOverlaps(protNR, mirs)


#what do we do with everything else - we might want to keep the rest as a single category so lets call them other
otherNR <- reduce(gr)
otherNR <- subsetByOverlaps(otherNR, mirs, invert = TRUE)
otherNR <- subsetByOverlaps(otherNR, rrnas, invert = TRUE)
otherNR <- subsetByOverlaps(otherNR, trnas, invert = TRUE)
otherNR <- subsetByOverlaps(otherNR, protNR, invert = TRUE)

otherNR #and because we reduce it we have no metadata - but it doesnt really matter cause our data should be mapping here like if it does its not good cause its not miRNA
otherNR$gene_name <- paste0('other', seq_along(otherNR))
otherNR


mirs
#now we just have to make our saf files and then we can just concatenate our results:
mirDF <- data.frame(
  'GeneID' = mirs$gene_name,
  'Chr' = seqnames(mirs),
  'Start' = start(mirs),
  'End' = end(mirs),
  'Strand' = strand(mirs)
)

mirDF$type <- 'miRNA' #add a column
head(mirDF)

rrnas
rrnasDF <- data.frame(
  'GeneID' = rrnas$gene_name,
  'Chr' = seqnames(rrnas),
  'Start' = start(rrnas),
  'End' = end(rrnas),
  'Strand' = strand(rrnas)
)

rrnasDF$type <- 'rRNA'
rrnasDF

trnas
trnasDF <- data.frame(
  'GeneID' = trnas$gene_name,
  'Chr' = seqnames(trnas),
  'Start' = start(trnas),
  'End' = end(trnas),
  'Strand' = strand(trnas)
)

trnasDF$type <- 'tRNA'



protNR
protDF <- data.frame(
  'GeneID' = protNR$gene_name,
  'Chr' = seqnames(protNR),
  'Start' = start(protNR),
  'End' = end(protNR),
  'Strand' = strand(protNR)
)

protDF$type <- 'protein-coding'
protDF



otherNR
otherDF <- data.frame(
  'GeneID' = otherNR$gene_name,
  'Chr' = seqnames(otherNR),
  'Start' = start(otherNR),
  'End' = end(otherNR),
  'Strand' = strand(otherNR)
)

otherDF$type <- 'other'
otherDF

#okay they are READYYYYYY lets combine them yay!
allDF <- rbind(mirDF, rrnasDF, trnasDF, protDF, otherDF)
table(allDF$type)

write.table(allDF, file = gzfile('annotationforall.SAF.gz'), quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)



