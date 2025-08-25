library(biomaRt)
# List all available attributes
listAttributes(ensembl)

# Set up the connection to Ensembl
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for gene synonyms and NCBI Gene IDs
# Fetch the gene information, including Ensembl Gene IDs and Transcript IDs
gene_info <- getBM(
  attributes = c("hgnc_symbol", "external_gene_name", "entrezgene_id", "external_synonym", 
                 "ensembl_gene_id", "ensembl_transcript_id"),
  mart = ensembl
)

# View the first few rows of the data
head(gene_info)

# Get the number of rows (genes) in the gene_info dataset
num_genes <- nrow(gene_info)
num_genes
# Save the gene_info data as a CSV file
write.table(gene_info, file = "human_gene.csv", sep = ",", row.names = FALSE, col.names = TRUE)
