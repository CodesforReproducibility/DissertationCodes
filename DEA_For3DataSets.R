#in this script i will be analysing the 3 data sets i have all together 

library(IRanges)
library(Rsubread)
library("edgeR")
library("statmod")
library(limma)

#first i will bring in my saf annotation file (it should be in the desktop)
#get my saf annotation file - in the desktop
saf_data_all_of_genome <- read.table("annotationforall.SAF.gz", header=TRUE, sep="\t", stringsAsFactors=FALSE)
head(saf_data_all_of_genome)

#next set is to get my bam files - i joining all of them so im thinking im gonnna make a new folder so i can list all do them at the same time
bam_files <- list.files(path = "./full_analysis_of_all/", pattern = "*.bam", full.names = TRUE)
length(bam_files) #wanna make sure the 36 files are here

#now the next step is to run feature counts on the files 
counting <- featureCounts(
  files = bam_files,                    # Replace with your actual BAM file
  annot.ext = saf_data_all_of_genome,   # My saf file with the annotated genome
  isGTFAnnotationFile = FALSE,          # its cause its a bam file 
  useMetaFeatures = TRUE,
  strandSpecific = 1,                   # recommended for miRNA in illumina strand specific 1
  nthreads = 4                          # Use multiple cores if available (im keeping at 4 do not change)
)


head((counting$counts)) #just to looking at our counts and they are looking slay as hell period

#i this its really important to save this data set so i dont have to keeping running features counts all the time
#now lets save the matrix
counts_matrix <- counting$counts
head(counts_matrix)

#im saving it as a R document 
saveRDS(counts_matrix, "miRNA_counts_of_ALL_data.rds")


#i can then load it in if i need to GET the data in from FeatureCounts
counts_matrix <- readRDS('~/Desktop/miRNA_counts_of_ALL_data.rds') #its the readRDS function duh
#this is is not needed tho only to get the data in
head(counts_matrix)

dim(counts_matrix) #ive got my 36 samples purio all looking good also i have 92330 rows


# Now we are going to collapse the data
byLane <- split(colnames(counts_matrix), sub("_L1|_L2","",colnames(counts_matrix)))
counts <- matrix(data=0, nrow=nrow(counts_matrix), ncol=length(byLane))
rownames(counts) <- rownames(counts_matrix)
colnames(counts) <- names(byLane)
for (n in names(byLane)) {
  print(n)
  counts[,n] <- rowSums(counts_matrix[,byLane[[n]],drop=FALSE])
}
 
dim(counts) #okay im sure it collapsed correctly cause - 18 files now and still the same 92330 rows
head(counts, 12) #wow this code is acc crazy

#anyways - so i am gonna filter all those miRNAs right away cause i know i have to or else ill get so much counts i acc dont want

mirna_counts <- counts[grepl("^mmu-", rownames(counts)), ] #Okay so just to let u know this babe is getting only the miRNAs from the counts 
dim(mirna_counts) #okay damn i dropped so many of them crazy (i now have 1966 miRNAs)
round(colSums(mirna_counts) / 1e6, 1)

#alright so like now lets filter the lowly expressed miRNAs
# Define a threshold (e.g., keep miRNAs expressed at CPM ≥ 1 in at least 2 samples)
# also im deciding to do cpm greater or equal to 1 cause it will remove any samples that have lowly expressed miRNAs
keep_mirna <- rowSums(cpm(mirna_counts) >= 1) >= 2  
print(table(keep_mirna))
filtered_cpm_mirna <- mirna_counts[keep_mirna, ]

dim(filtered_cpm_mirna) #so now u have 589 possible miRNAs damn so little but more than before okay slay maybe this is a good sign okay eat it up

colnames(filtered_cpm_mirna) #see what they are called 
group <- factor(gsub("_r\\d+\\.for_counting.bam$", "", colnames(filtered_cpm_mirna)))
print(group) #now they are divided into the 6 conditions that i acc have
table(group) #and i can see that i have three replicated in each condition

#now we make a dge list for this whole thing
dge <- DGEList(counts=filtered_cpm_mirna, group=group) #just to be aware by data has not been normalized 
dge

#define the colors by the 6 conditions we have which is the 6 groups
colors <- rainbow(length(levels(group))) 
names(colors) <- levels(group)
colors #and they are defined i think this is crazy

#lets get to the analysis - first: an MDS plot lets see where these samples are at
par(mfrow=c(1,1))
#then we plot it
#pdf('MDSplotnotnormalizsed.pdf')
plotMDS(dge, col = colors[group], main = "MDS Plot not yet normalized")
legend("topright", legend = levels(group), col = colors, pch = 19, title = "Groups")
#dev.off()

#pdf('MDSplotnotnormalizsed_points.pdf')
plotMDS(dge, 
        col = colors[group],
        main = "MDS Plot not yet normalized",
        pch = 16)
legend("topright", legend = levels(group), col = colors, pch = 19, title = "Groups",  xpd = TRUE, inset = c(0.0009, -0.05) )
#dev.off()


#DATA Normalized !!!!!! 
dge <- calcNormFactors(dge)
dge$samples

#pdf('MDSplotnormalized_points.pdf')
plotMDS(dge, 
        col = colors[group],
        main = "MDS Plot of all Datasets Normalized",
        pch = 16)
legend("topright", legend = levels(group), col = colors, pch = 19, title = "Groups",  xpd = TRUE, inset = c(0.0009, -0.05) )
#dev.off()

#check the samples
par(mfrow=c(1,2))

for (i in c(3,4,5)) {
  plotMD(cpm(dge, log=TRUE), column=i)
  grid(col = "blue")
  abline(h=0, col='red', lty=2, lwd=2)
}

#then we are making the design so we can make the BCV plot 
design <- model.matrix(~0+dge$samples$group)
colnames(design) <- levels(dge$samples$group)
design


#pdf('bcvplot.pdf')
par(mar = c(6, 6, 4, 2))
dge <- estimateDisp(dge, design = design, robust = TRUE) #this is to estimate dispersion so to calculate the plot BCV
par(mfrow=c(1,1))
plotBCV(dge, main='Analysis of Full Data Set BCV plot') 
#dev.off()

contrasts <- makeContrasts(
  'FullDose'= 'FD_imp_Flag - FD_res_AGO2',
  'HalfDose_CD45.1' ='HD_imp_Flag - HD_res_AGO2',
  'HalfDose_CD45.2'= 'HD_imp_AGO2 - HD_res_Flag',
  levels = design
)

contrasts
contrast_fulldose <- contrasts[,1]
contrast_fulldose

contrast_halfdoseCD45.1 <- contrasts[,2]
contrast_halfdoseCD45.1

contrast_halfdoseCD45.2 <- contrasts[,3]
contrast_halfdoseCD45.2

#-------------------------------------------------------------------------------------------------------------------
#             For Full Dose
#-------------------------------------------------------------------------------------------------------------------



fit <- glmFit(dge, dispersion = dge$tagwise.dispersion) #explain why I used the tagwise dispersion instead of the remember that it  
lrt <- glmLRT (fit, contrast = contrast_fulldose)
topTags(lrt)

decidetestFD <- decideTests(lrt, adjust.method="BH", p.value=0.05, lfc=0)
table(decidetestFD) #19 miRNAs upregulates

#The smear plot
deGenes <- rownames(lrt)[decidetestFD != 0]
plotSmear(lrt, de.tags = deGenes, main='Plot Smear for CD45.1 Full Dose') #here we have a plot with the genes that are more significant

topTable <- topTags(lrt, n=Inf)$table
topTable

upregulated_miRNAs_for_FullDose <- topTable[topTable$logFC > 0 & topTable$FDR < 0.05, ]
print(upregulated_miRNAs_for_FullDose)



important_miRNAs_for_Full_dose <- c("mmu-miR-3068-3p", "mmu-miR-3535", 
                      "mmu-miR-126a-5p", "mmu-miR-199a-3p", "mmu-miR-205-5p", 
                      "mmu-miR-5099", "mmu-miR-5121")

topTable$negLogFDR <- -log10(topTable$FDR)
# Extract the required values for the important miRNAs
miRNA_positions <- topTable[rownames(topTable) %in% important_miRNAs_for_Full_dose, c("logFC", 'negLogFDR')]


miRNA_positions



upregulated_miRNAs_for_FullDose <- topTable[topTable$logFC > 0 & topTable$FDR < 0.05, ]
upregulated_miRNAs_for_FullDose
nonsig_miRNAs_FullDose <- rownames(decidetestFD)[decidetestFD == 0]
downregulate_miRNAs_for_FullDose <- topTable[topTable$logFC < 0 & topTable$FDR < 0.05, ]
downregulate_miRNAs_for_FullDose



#pdf('~/Desktop/volcanoplotFD.pdf')
plot(topTable$logFC, topTable$negLogFDR, pch=19, col="black", 
     main="Volcano Plot for CD45.1 Donor Cells Full Dose", 
     xlim=c(-10, 12), ylim=c(0, 7), xlab="logFC", ylab="-log10(FDR)")


points(upregulated_miRNAs_for_FullDose$logFC, upregulated_miRNAs_for_FullDose$negLogFDR, pch=19, col="red", cex=1)
points(downregulate_miRNAs_for_FullDose$logFC, downregulate_miRNAs_for_FullDose$negLogFDR, pch=19, col="blue", cex=1)

text(miRNA_positions$logFC, miRNA_positions$negLogFDR, 
     labels = rownames(miRNA_positions), 
     pos = 4,        # Position to the right of the points
     col = 'red',   # Label color
     cex = 0.8,      # Label size
     offset = 0.1)   # Adjust this to move the label closer to the point

legend("topleft", legend=c("Not significant", "Upregulated", "Downregulated"),
       col=c("black", "red", "blue"), pch=19, cex=0.8)
#dev.off()
miRNA_positions

#-------------------------------------------------------------------------------------------------------------------
#               For Half Dose CD45.1 
#-------------------------------------------------------------------------------------------------------------------


#For Half Dose CD45.1
fit <- glmFit(dge, dispersion = dge$tagwise.dispersion) #explain why I used the common dispersion instead of the trend 
lrt <- glmLRT (fit, contrast = contrast_halfdoseCD45.1)
topTags(lrt)

decidetestHDcd45.1 <- decideTests(lrt, adjust.method="BH", p.value=0.05, lfc=0)
table(decidetestHDcd45.1) #11 miRNAs upregulated

#for the negative control
#decidetestHDcd45.1
#fornegativecontrol <- saveRDS(decidetestHDcd45.1, 'ForNegativeControl.RDS')

#The smear plot
deGenes <- rownames(lrt)[decidetestHDcd45.1 != 0]
plotSmear(lrt, de.tags = deGenes, main='Plot Smear for CD45.1 Half Dose') #here we have a plot with the genes that are more significant

topTable <- topTags(lrt, n=Inf)$table
topTable

upregulated_miRNAs_for_HalfDosecd45.1 <- topTable[topTable$logFC > 0 & topTable$FDR < 0.05, ]
print(upregulated_miRNAs_for_HalfDosecd45.1)

important_miRNAs_for_halfDose <- c("mmu-miR-3068-3p", "mmu-miR-3535",
                      "mmu-miR-126a-5p", "mmu-miR-199a-3p", "mmu-miR-205-5p", 
                      "mmu-miR-5099", "mmu-miR-5121")
topTable$negLogFDR <- -log10(topTable$FDR)


miRNA_positions_HD <- topTable[rownames(topTable) %in% important_miRNAs_for_halfDose, c("logFC", 'negLogFDR')]



miRNA_positions_HD






upregulated_miRNAs_for_HD <- topTable[topTable$logFC > 0 & topTable$FDR < 0.05, ]
upregulated_miRNAs_for_HD
nonsig_miRNAs_HD <- rownames(decidetestFD)[decidetestFD == 0]
downregulate_miRNAs_for_HD <- topTable[topTable$logFC < 0 & topTable$FDR < 0.05, ]
downregulate_miRNAs_for_HD

topTable$negLogFDR <- -log10(topTable$FDR)


#pdf('~/Desktop/cd45.1.pdf')
plot(topTable$logFC, topTable$negLogFDR, pch=19, col="black", 
     main="Volcano Plot for CD45.1 Donor Cells Half Dose", 
     xlim=c(-10, 15), ylim=c(0, 6), xlab="logFC", ylab="-log10(FDR)")


points(upregulated_miRNAs_for_HD$logFC, upregulated_miRNAs_for_HD$negLogFDR, pch=19, col="red", cex=1)
points(downregulate_miRNAs_for_HD$logFC, downregulate_miRNAs_for_HD$negLogFDR, pch=19, col="blue", cex=1)

text(miRNA_positions_HD$logFC, miRNA_positions_HD$negLogFDR, 
     labels = rownames(miRNA_positions_HD), 
     pos = 4,        # Position to the right of the points
     col = 'red',   # Label color
     cex = 0.8,      # Label size
     offset = 0.05)   # Adjust this to move the label closer to the point

legend("topleft", legend=c("Not significant", "Upregulated", "Downregulated"),
       col=c("black", "red", "blue"), pch=19, cex=0.8)
#dev.off()
#-------------------------------------------------------------------------------------------------------------------
#                7 Not upregulated miRNAs - controls for MAPK not being enriched 
#-------------------------------------------------------------------------------------------------------------------

#negativecontrol <- rownames(lrt)[decidetestHDcd45.1 == 0]
#negativecontrol
#length(negativecontrol)
#saveRDS(negativecontrol, '~/Desktop/NegativeControlToLoad.RDS')

#-------------------------------------------------------------------------------------------------------------------
#                Done
#-------------------------------------------------------------------------------------------------------------------

#For Half Dose CD45.2
fit <- glmFit(dge, dispersion = dge$tagwise.dispersion) #explain why I used the common dispersion instead of the trend 
lrt <- glmLRT (fit, contrast = contrast_halfdoseCD45.2)
topTags(lrt)

decidetestHDcd45.2 <- decideTests(lrt, adjust.method="BH", p.value=0.05, lfc=0)
table(decidetestHDcd45.2) #still 0 miRNAs

#The smear plot
deGenes <- rownames(lrt)[decidetestHDcd45.2 != 0]
plotSmear(lrt, de.tags = deGenes, main='Plot Smear for CD45.2 Half Dose') #here we have a plot with the genes that are more significant

topTable <- topTags(lrt, n=Inf)$table

upregulated_miRNAs_for_HalfDosecd45.2 <- topTable[topTable$logFC > 0 & topTable$FDR < 0.05, ]
print(upregulated_miRNAs_for_HalfDosecd45.2)


