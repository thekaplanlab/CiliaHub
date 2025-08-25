############################################################
# Mouse Ciliary Gene Phenotype Analysis
# Author: [Your Name]
# Date: [YYYY-MM-DD]
# Description:
#   Fetches phenotypes for mouse genes from IMPC,
#   detects ciliopathy-related phenotypes,
#   generates outputs (Excel, plots) for publication and GitHub.
############################################################

# ================================
# 0. Load Libraries
# ================================
library(readxl)
library(dplyr)
library(httr)
library(jsonlite)
library(stringr)
library(tidyr)
library(writexl)
library(ggplot2)

# ================================
# 1. Configuration
# ================================
input_gene_file <- "data/known_ciliary_genes_mouse.xlsx"  # change for your repo
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

# Ciliopathy-related keywords
ciliopathy_terms <- sort(unique(tolower(c(
  "cyst", "hydrocephalus", "polydactyly", "situs", "brain ventricle", "flagellum",
  "axoneme", "cilium", "renal", "kidney", "fibrosis", "dysgenesis", "hydrops",
  "left-right axis", "left-right patterning", "retinal degeneration",
  "congenital heart defect", "biliary", "hedgehog signaling", "infertility",
  "hydrops fetalis", "shortened long bones", "cerebral anomalies",
  "coronary and vascular anomalies", "facial anomalies", "ophthalmic anomalies",
  "nasal anomalies", "cognitive anomalies", "skeletal anomalies",
  "respiratory anomalies", "neural anomalies", "hormonal anomalies",
  "digestive anomalies", "reproductive anomalies", "aural anomalies",
  "liver anomalies", "organ anomalies", "leber congenital amaurosis",
  "nystagmus", "coloboma", "hearing loss", "respiratory infections",
  "bronchiectasis", "sinusitis", "otitis media", "situs inversus",
  "heterotaxy", "congenital heart", "nephronophthisis", "polycystic kidney",
  "hepatic fibrosis", "bile duct", "short rib", "limb shortening",
  "vertebral anomalies", "obesity", "diabetes insipidus", "joubert", "ataxia",
  "seizures", "developmental delay", "ciliopathy", "primary ciliary dyskinesia",
  "bardet-biedl", "meckel", "senior-loken", "oro-facial-digital",
  "retinitis pigmentosa", "flagellar", "anosmia", "cone-rod dystrophy",
  "cystic kidney", "abnormal auditory brainstem response"
))))

# ================================
# 2. Helper Functions
# ================================

# Capitalize first letter of each gene name
capitalize_gene <- function(gene) {
  ifelse(is.na(gene) | gene == "", gene,
         paste0(toupper(substring(gene, 1, 1)), tolower(substring(gene, 2))))
}

# Fetch all phenotypes for a given mouse gene from IMPC API
get_all_phenotypes <- function(gene_symbol) {
  page <- 0
  rows_per_page <- 100
  all_phenos <- c()
  
  repeat {
    url <- paste0(
      "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?",
      "q=marker_symbol:", gene_symbol,
      "&start=", page * rows_per_page,
      "&rows=", rows_per_page,
      "&wt=json"
    )
    
    resp <- GET(url)
    if (status_code(resp) != 200) break
    
    data <- fromJSON(content(resp, "text"), flatten = TRUE)
    phenos <- data$response$docs$mp_term_name
    
    if (length(phenos) == 0) break
    all_phenos <- c(all_phenos, phenos)
    
    if (length(phenos) < rows_per_page) break
    page <- page + 1
  }
  
  unique(all_phenos)
}

# Check if any ciliopathy term is present in phenotypes
has_ciliopathy <- function(phenos) {
  phenos_lower <- tolower(paste(phenos, collapse = ", "))
  any(sapply(ciliopathy_terms, function(term) grepl(term, phenos_lower, fixed = TRUE)))
}

# ================================
# 3. Load Data
# ================================
genes_df <- read_excel(input_gene_file)
genes <- unique(capitalize_gene(genes_df$MouseGene)) # adjust column name if needed

# ================================
# 4. Fetch Phenotypes
# ================================
phenotype_results <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  phenotype_results[[gene]] <- get_all_phenotypes(gene)
}

# ================================
# 5. Create Summary Table
# ================================
phenotype_df <- data.frame(
  Gene = names(phenotype_results),
  Phenotypes = sapply(phenotype_results, paste, collapse = ", "),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Ciliopathy_Flag = ifelse(sapply(phenotype_results, has_ciliopathy), "YES", "NO")
  )

write_xlsx(phenotype_df, file.path(output_dir, "phenotype_results_with_ciliopathy_flag.xlsx"))

# ================================
# 6. Create Gene ?? Phenotype Matrix
# ================================
phenotype_matrix <- matrix(FALSE, nrow = length(ciliopathy_terms), ncol = length(genes),
                           dimnames = list(ciliopathy_terms, genes))

for (gene in genes) {
  phenos <- tolower(paste(phenotype_results[[gene]], collapse = ", "))
  for (term in ciliopathy_terms) {
    if (grepl(term, phenos, fixed = TRUE)) {
      phenotype_matrix[term, gene] <- TRUE
    }
  }
}

phenotype_matrix_df <- cbind(Phenotype = rownames(phenotype_matrix), as.data.frame(phenotype_matrix))
write_xlsx(phenotype_matrix_df, file.path(output_dir, "phenotype_matrix.xlsx"))

# ================================
# 7. Plot: Ciliopathy Phenotypes per Gene (Top 150)
# ================================
top_genes <- phenotype_matrix_df %>%
  select(-Phenotype) %>%
  summarise(across(everything(), ~sum(. == TRUE, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Count") %>%
  arrange(desc(Count)) %>%
  slice_head(n = 150) %>%
  pull(Gene)

phenotype_long <- phenotype_matrix_df %>%
  pivot_longer(cols = -Phenotype, names_to = "Gene", values_to = "Present") %>%
  filter(Present == TRUE & Gene %in% top_genes)

p <- ggplot(phenotype_long, aes(x = Gene, y = Phenotype)) +
  geom_point(color = "steelblue", size = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Ciliopathy Phenotypes per Gene in Mouse",
       x = "Gene", y = "Phenotype")

ggsave(file.path(output_dir, "ciliopathy_phenotypes_per_gene.png"),
       plot = p, width = 20, height = 8, dpi = 300)
