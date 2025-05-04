library(dplyr)
list.files(path='~/Desktop/pathway_gene_hits/', pattern = '*.txt')
length(list.files(path='~/Desktop/pathway_gene_hits/', pattern = '*.txt')) #i have 38 files


#this includes all possible gene targets from all the miRNA hits I got 
merged_data <- readRDS('~/Desktop/AlltheDATAmerged.rds')
merged_data

#starting with the pathway analysis for GF-RTK-(GRB2 no hits)-SOS-Ras
#first we get all the DAVID gene hits (like what matched)
david_results <- read.delim("MAPK_AllGene_Hits.txt", header = TRUE, sep = "\t") 
david_results$Gene.Name


#im opening a data frame to add it everything
miRNA_hits_df <- data.frame(
  gene_symbol = character(),
  gene_ensembl_id = character(),
  mirna_name = character(),
  match_source = character(),
  stringsAsFactors = FALSE
)



PathwayGenes <- list.files(path='~/Desktop/pathway_gene_hits/', pattern = '*.txt', full.names = TRUE)


for (gene in PathwayGenes) {
  gene_data <- read.delim(gene, header = TRUE, sep = "\t")
  intersect_gene_data <- intersect(gene_data$Gene.Name, david_results$Gene.Name)
  Abbreviation_names <- sub(".*\\((.*)\\)", "\\1", intersect_gene_data)
  for (abbreviation in Abbreviation_names) {
    results_on_hits <- merged_data[merged_data$gene_symbol == abbreviation, ]
    results_on_hits <- results_on_hits %>%
      mutate(match_source = gene)  # Directly assign the current file as the match source
    if (ncol(miRNA_hits_df) != ncol(results_on_hits)) {
      # Align column names in miRNA_hits_df with results_on_hits (this should match your structure)
      colnames(miRNA_hits_df) <- colnames(results_on_hits)
    }
    
    miRNA_hits_df <- rbind(miRNA_hits_df, results_on_hits)
  }
    
  }

head(miRNA_hits_df)

length(unique(miRNA_hits_df$gene_symbol))

length(miRNA_hits_df$gene_symbol)

saveRDS(miRNA_hits_df, 'MAPKmiRNAhits.RDS')


#_________________________________________________________________________________________________________________________________________________
#                             Looking at my data
#__________________________________________________________________________________________________________________________________________________



miRNA_hits_df <- readRDS('MAPKmiRNAhits.RDS')

# Extract the gene name from the match_source column
miRNA_hits_df$match_source <- sub(".*/(.*)_All_Genes\\.txt", "\\1", miRNA_hits_df$match_source)

# View the updated data
head(miRNA_hits_df)

unique(miRNA_hits_df$match_source)
length(unique(miRNA_hits_df$match_source)) #38 


(sum(miRNA_hits_df$match_source=='RTK')) #same as other SLAY

unique(miRNA_hits_df$mirna_name)

miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-5099',] #i got a hit (what)
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-199a-3p',]
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-205-5p',]
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-3535',]
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-5121',]
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-126a-5p',]
miRNA_hits_df[miRNA_hits_df$mirna_name=='mmu-miR-3068-3p',]

nrow(miRNA_hits_df[miRNA_hits_df$match_source=='GF',])


#_________________________________________________________________________________________________________________________________________________
#                             Making a plot for the gene targets 
#__________________________________________________________________________________________________________________________________________________



df <- data.frame(
  mirna_name = character(),
  counts = numeric(),
  stringsAsFactors = FALSE
)


mirna_name <- c(miRNA_hits_df$mirna_name)

for (mirna in unique(miRNA_hits_df$mirna_name)) {
  counts <- nrow(miRNA_hits_df[miRNA_hits_df$mirna_name == mirna, ])
  print(paste(mirna, "→", counts))
  df <- rbind(df, data.frame(mirna_name = mirna, counts = counts))
}

df

library(ggplot2)

pdf('total gene targets in MAPK pathway.pdf')
ggplot(df, aes(x = mirna_name, y = counts, fill = mirna_name)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = counts), vjust = -0.3) + 
  theme_minimal() +
  labs(x = "miRNA", y = "Number of Gene Targets per miRNA", title = "Gene Targets of miRNA in MAPK Signaling Pathway") +
  scale_fill_manual(values = c("skyblue", "orange", 'violetred2', "sienna2",'slateblue4', 'royalblue2', 'peachpuff4'  )) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", )
dev.off()
