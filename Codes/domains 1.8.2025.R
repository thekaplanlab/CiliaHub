# Required packages
packages <- c("BiocManager", "biomaRt", "tidyverse", "readxl", "writexl", "ggrepel", "forcats", "purrr")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

library(biomaRt)
library(tidyverse)
library(readxl)
library(writexl)
library(ggrepel)
library(forcats)
library(purrr)

# 1. Load ciliary genes (assuming the Excel file contains ciliary genes)
cilia_genes <- read_excel("/Users/sebihacevik/Downloads/all gene 1.8.2025.xlsx")[[1]] %>%
  na.omit() %>%
  unique()
message("Loaded ", length(cilia_genes), " ciliary genes.")

# 2. Connect to Ensembl
connect_ensembl <- function() {
  mirrors <- c("useast", "uswest", "asia")
  for (m in mirrors) {
    cat("Trying Ensembl mirror:", m, "...\n")
    try({
      ensembl <- useEnsembl("genes", "hsapiens_gene_ensembl", mirror = m)
      cat("Connected to", m, "mirror.\n")
      return(ensembl)
    }, silent = TRUE)
  }
  cat("All live mirrors failed — trying archive version.\n")
  ensembl <- useEnsembl("genes", "hsapiens_gene_ensembl", version = 110)
  return(ensembl)
}
ensembl <- connect_ensembl()

# 3. Fetch all human genes from Ensembl
all_human_genes <- getBM(
  attributes = c("hgnc_symbol"),
  mart = ensembl
)$hgnc_symbol %>%
  na.omit() %>%
  unique()
message("Fetched ", length(all_human_genes), " human genes.")

# 4. Function to fetch InterPro domains
get_domains <- function(gene_list, filter_type = "hgnc_symbol") {
  df <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "interpro", "interpro_description"),
    filters = filter_type,
    values = gene_list,
    mart = ensembl
  ) %>%
    drop_na(interpro, interpro_description) %>%
    filter(interpro_description != "")
  return(df)
}

# 5. Fetch domains for ciliary genes and all human genes
cilia_domains <- get_domains(cilia_genes)
all_domains <- get_domains(all_human_genes)
message("Fetched ", nrow(cilia_domains), " ciliary domain entries and ", nrow(all_domains), " total domain entries.")

# 6. Summarize domain counts
cilia_domain_counts <- cilia_domains %>%
  count(interpro, interpro_description, name = "CiliaGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)
all_domain_counts <- all_domains %>%
  count(interpro, interpro_description, name = "TotalGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

# 7. Merge domain counts
domain_comparison <- full_join(cilia_domain_counts, all_domain_counts, by = c("DomainID", "DomainDescription")) %>%
  replace_na(list(CiliaGeneCount = 0, TotalGeneCount = 0))

# 8. Perform Fisher’s exact test for enrichment
domain_comparison <- domain_comparison %>%
  mutate(
    NonCiliaGeneCount = TotalGeneCount - CiliaGeneCount,
    NonCiliaNonDomain = length(all_human_genes) - TotalGeneCount,
    PValue = map2_dbl(CiliaGeneCount, TotalGeneCount, ~ {
      matrix <- matrix(c(.x, .y - .x, length(cilia_genes) - .x, length(all_human_genes) - .y - (length(cilia_genes) - .x)), nrow = 2)
      fisher.test(matrix)$p.value
    }),
    AdjPValue = p.adjust(PValue, method = "BH")
  )

# 9. Filter for significant ciliary-specific domains
significant_domains <- domain_comparison %>%
  filter(CiliaGeneCount > 0, AdjPValue < 0.05) %>%
  arrange(AdjPValue)

# 10. Save results
write_xlsx(significant_domains, "CiliarySpecific_Domains.xlsx")

# 11. Plot top significant domains
top_domains <- significant_domains %>%
  arrange(desc(CiliaGeneCount)) %>%
  slice_head(n = 20)
top_domains
ggplot(top_domains, aes(x = CiliaGeneCount, y = fct_reorder(DomainDescription, CiliaGeneCount))) +
  geom_col(fill = "skyblue") +
  labs(
    title = "Top 20 Domains Enriched in Ciliary Genes",
    x = "Number of Ciliary Genes",
    y = "InterPro Domain"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))

ggsave("Top20_CiliarySpecific_Domains.pdf", width = 10, height = 6)