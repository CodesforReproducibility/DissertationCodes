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


