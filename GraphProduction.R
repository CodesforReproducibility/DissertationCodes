#______________________________________________________________________________________________________________
#                            ENRICHMENT PATHWAY GRAPH
#______________________________________________________________________________________________________________
siganlingpathways <- read.delim('~/Desktop/TrueData_DavidOutput.txt', header=TRUE, sep='\t')
head(siganlingpathways)
colnames(siganlingpathways)

# Assuming siganlingpathways$Term contains values like "mmu04010:MAPK signaling pathway"
siganlingpathways$ID <- sub(":.*", "", siganlingpathways$Term)  # Extract the ID before the colon
siganlingpathways$term <- sub(".*:", "", siganlingpathways$Term)  # Extract the term after the colon


# Assuming siganling pathways is already loaded and available

terms <- data.frame(
  ID = siganlingpathways$ID,
  term = siganlingpathways$term,
  FDR = siganlingpathways$FDR,
  genes = siganlingpathways$Genes,
  count = siganlingpathways$Count,
  pvalue = siganlingpathways$PValue
)

# Sorting by FDR
labels = terms[order(terms$FDR, decreasing = TRUE), "term"]
terms$term = factor(terms$term, levels = labels)

terms$FDR_label <- ifelse(terms$FDR <= sort(terms$FDR, decreasing = TRUE)[4], 
                          as.character(terms$FDR), 
                          "")


# Drawing the plot
library(ggplot2)

ggplot(data = terms) +
  geom_bar(aes(y = term, x = count, fill = FDR), stat = 'identity') +  # Swap axes
  scale_fill_gradientn(
    colors = c("red", "orange", "blue"),  # Define the colors for the FDR ranges
    values = scales::rescale(c(0, 0.05, 0.25, 1)),  # Scale the values to map to the colors
    limits = c(0, 1)  # Define the range of values for scaling
  ) +
  ylab("Term") + 
  xlab("Gene count") + 
  ggtitle('KEGG Enrichment Pathway for Gene Targets of Upregulated miRNAs')+
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,200)) +  # Adjust x-axis scaling
  scale_y_discrete(expand = c(0, 0)) +  # Adjust y-axis scaling
  theme_bw()

#______________________________________________________________________________________________________________
#                           NEGATIVE CONTROLS ENRICHMENT PATHWAY GRAPH
#______________________________________________________________________________________________________________

#now for negative controls
negativecontrol1 <- read.delim('~/Desktop/negativecontrol1st.txt')
colnames(negativecontrol1)
negativecontrol1$ID <- sub(":.*", "", negativecontrol1$Term)  # Extract the ID before the colon
negativecontrol1$term <- sub(".*:", "", negativecontrol1$Term)  # Extract the term after the colon


negativeterms <- data.frame(
  ID = negativecontrol1$ID,
  term = negativecontrol1$term,
  FDR = negativecontrol1$FDR,
  genes = negativecontrol1$Genes,
  count = negativecontrol1$Count,
  pvalue = negativecontrol1$PValue
)

labels = negativeterms[order(negativeterms$FDR, decreasing = TRUE), "term"]
negativeterms$term = factor(negativeterms$term, levels = labels)

negativeterms$FDR_label <- ifelse(negativeterms$FDR <= sort(negativeterms$FDR, decreasing = TRUE)[4], 
                          as.character(negativeterms$FDR), 
                          "")
par(mfrow = c(1,1))
#pdf('~/Desktop/RandomGenes.pdf')
ggplot(data = negativeterms) +
  geom_bar(aes(y = term, x = count, fill = FDR), stat = 'identity') +  # Swap axes
  scale_fill_gradientn(
    colors = c("red", "orange", "blue"),  # Define the colors for the FDR ranges
    values = scales::rescale(c(0, 0.05, 0.25, 1)),  # Scale the values to map to the colors
    limits = c(0, 1)  # Define the range of values for scaling
  ) +
  ylab("Term") + 
  xlab("Gene count") + 
  ggtitle('KEGG Enrichment Pathway for Random Genes (Negative Control)')+
  theme(axis.text.x = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,140)) +  # Adjust x-axis scaling
  scale_y_discrete(expand = c(0, 0)) +  # Adjust y-axis scaling
  theme_bw()
#dev.off()

#other negative control (upregulated miRNAs)
negativecontrol2 <- read.delim('~/Desktop/forplotnegativecontrol7mirnas.txt')
colnames(negativecontrol2)
negativecontrol2$ID <- sub(":.*", "", negativecontrol2$Term)  # Extract the ID before the colon
negativecontrol2$term <- sub(".*:", "", negativecontrol2$Term)  # Extract the term after the colon


negativetermsfor2 <- data.frame(
  ID = negativecontrol2$ID,
  term = negativecontrol2$term,
  FDR = negativecontrol2$FDR,
  genes = negativecontrol2$Genes,
  count = negativecontrol2$Count,
  pvalue = negativecontrol2$PValue
)

labels = negativetermsfor2[order(negativetermsfor2$FDR, decreasing = TRUE), "term"]
negativetermsfor2$term = factor(negativetermsfor2$term, levels = labels)

negativetermsfor2$FDR_label <- ifelse(negativetermsfor2$FDR <= sort(negativetermsfor2$FDR, decreasing = TRUE)[4], 
                                  as.character(negativetermsfor2$FDR), 
                                  "")


#pdf('~/Desktop/miRNAnoturegulated.pdf')
ggplot(data = negativetermsfor2) +
  geom_bar(aes(y = term, x = count, fill = FDR), stat = 'identity') +  # Swap axes
  scale_fill_gradientn(
    colors = c("red", "orange", "blue"),  # Define the colors for the FDR ranges
    values = scales::rescale(c(0, 0.05, 0.25, 1)),  # Scale the values to map to the colors
    limits = c(0, 1)  # Define the range of values for scaling
  ) +
  ylab("Term") + 
  xlab("Gene count") + 
  ggtitle('KEGG Enrichment Pathway for Target Genes of miRNAs Not Differentially Expressed') +
  theme_minimal() +  # Use theme_minimal only
  theme(axis.text.x = element_text(color = "black", size = 10),  # Adjust x-axis labels size and rotate them
        axis.text.y = element_text(color = "black", size = 6.5, angle = 0, hjust = 1),  # Adjust y-axis label size
        axis.title.x = element_text(size = 12),  # Adjust x-axis title size
        axis.title.y = element_text(size = 12),  # Adjust y-axis title size
        plot.title = element_text(size = 14)) +  # Adjust plot title size
  scale_x_continuous(expand = c(0, 0), limits = c(0, 55)) +  # Adjust x-axis scaling
  scale_y_discrete(expand = c(0, 0))   # Adjust y-axis scaling
#dev.off()

#______________________________________________________________________________________________________________
#                          NETWORK miRNA-MAPK gene targets 
#______________________________________________________________________________________________________________

install.packages("igraph")
library(igraph)

gene_data <- read_excel('~/Desktop/miRNA_hits_MAPK_signaling_Cascade.xlsx')

# Create an edge list from the dataset
edges <- gene_data[, c("gene_symbol", "mirna_name")]

# Create an igraph object
g <- graph_from_data_frame(edges, directed = TRUE)  # If the relationship is directed (e.g., miRNA regulating genes)

# Plot the graph
plot(g, 
     vertex.size = 10,  # Size of the nodes
     vertex.color = "lightblue",  # Color of nodes
     vertex.label.color = "black",  # Color of the labels
     vertex.label.cex = 0.8,  # Size of the labels
     edge.arrow.size = 0.5,  # Arrow size for directed edges
     edge.width = 1,  # Edge width
     main = "Gene-miRNA Network")

library(stringr)

g <- reverse_edges(g)

layout_type <- layout_with_kk(g, niter = 500)
V(g)$name <- str_replace_all(V(g)$name, "^mmu-", "")


V(g)$color <- ifelse(grepl("^miR", V(g)$name), "pink2", "lightblue")  # miRNA nodes are green, genes are blue
V(g)$shape <- ifelse(grepl("^miR", V(g)$name), "square", "circle")  # miRNA nodes are square, genes are circle



# Install and load the ggraph package
install.packages("ggraph")
library(ggraph)

library(ggraph)
library(igraph)
library(ggplot2)

# Color the nodes based on whether they are miRNA or Gene
V(g)$node_type <- ifelse(grepl("^miR", V(g)$name), "miRNA", "Gene")

# Plotting the network
ggraph(g, layout = 'fr') + 
  # Draw edges with some transparency
  geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
  
  # Add nodes (points) and color them based on type (miRNA or Gene)
  geom_node_point(aes(color = node_type), size = 7) +
  
  # Add text labels to the nodes, using repel to avoid overlap
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  
  # Custom color scale for miRNAs and genes
  scale_color_manual(values = c("miRNA" = "pink2", "Gene" = "lightblue")) +
  
  # Add a title and subtitle to the plot
  ggtitle("Gene-miRNA Network", subtitle = "MAPK Signaling Cascade") +
  
  # Custom theme to remove grid lines and axes
  theme_void()


#______________________________________________________________________________________________________________
#                           SIGNALLING CASCADE MAPK GRAPH
#______________________________________________________________________________________________________________

library(tidyverse) 
library(readxl)
df <- read_excel('~/Desktop/miRNA_hits_MAPK_signaling_Cascade.xlsx')



ENSMUSG00000040118

counts_df <- tibble(gene_ID = character(), count = integer())

df
for (gene in unique(df$gene_ensembl_id)) {
  sums <- df %>% filter(gene_ensembl_id == gene) %>% nrow()
  counts <- print(paste("Gene:", gene, "Count:", sums))
  counts_df <- counts_df %>% add_row(gene_ID = gene, count = sums)
}


df[df$match_source=='c-fos',]

counts_df

entrez_map <- bitr(
  counts_df$gene_ID, # ENSEMBL gene IDs from your data
  fromType = "ENSEMBL",            # Specify the type of IDs you're converting from
  toType = "ENTREZID",             # Specify the type of IDs you're converting to
  OrgDb = org.Mm.eg.db             # Specify the organism database (for mouse, "org.Mm.eg.db")
)

head(entrez_map)


# Scale the counts to a range of 0 to 1
counts_df <- counts_df %>%
  mutate(scaled_count = (count - min(count)) / (max(count) - min(count)))

# Join with ENTREZ IDs
counts_df_with_entrez <- left_join(counts_df, entrez_map, by = c("gene_ID" = "ENSEMBL"))

counts_df_with_entrez


# Create named vector with scaled or log-transformed counts
gene_vector <- setNames(counts_df_with_entrez$scaled_count, counts_df_with_entrez$ENTREZID)  # or use log_count


library(pathview)

pathview(
  gene.data = gene_vector,      # Use the named vector of scaled counts
  gene.idtype = "entrez",       # ENTREZ ID type
  pathway.id = "mmu04010",      # Example pathway (replace with your actual pathway ID)
  species = "mmu",              # Mouse species
  out.suffix = "gene_counts",   # Output suffix for filenames
  low = "papayawhip",           # Color for low values
  mid = "paleturquoise",        # Color for mid values (optional)
  high = "plum",                # Color for high values
  na.col = 'whitesmoke'         # Color for missing values (optional)
)


