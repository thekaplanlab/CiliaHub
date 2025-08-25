# ==========================================================
# PubMed Scraper for Cilia/Flagella Gene Associations
# ==========================================================
# This script retrieves PubMed articles related to genes of interest 
# with search terms including cilia, cilium, transition zone, basal body, 
# ciliogenesis, flagella, and flagellum. 
# ==========================================================

# -----------------------------
# 1. Required Packages
# -----------------------------
library(rvest)
library(dplyr)
library(writexl)
library(readxl)
library(rentrez)

# -----------------------------
# 2. Function to Fetch PubMed Data
# -----------------------------
get_pubmed_data <- function(gene, max_results = 50) {
  query <- paste0(
    "https://pubmed.ncbi.nlm.nih.gov/?term=",
    gene, 
    "+AND+(cilia+OR+cilium+OR+transition+zone+OR+basal+body+OR+ciliogenesis+OR+flagella+OR+flagellum)",
    "&sort=date&size=", max_results
  )
  
  tryCatch({
    page <- read_html(query)
    
    paper_titles <- page %>%
      html_nodes(".docsum-title") %>%
      html_text() %>%
      trimws()
    
    pmids <- page %>%
      html_nodes(".docsum-pmid") %>%
      html_text() %>%
      trimws()
    
    if (length(paper_titles) == 0 || length(pmids) == 0) {
      cat("No results found for gene:", gene, "\n")
      return(data.frame(Gene = gene, PMID = NA, Title = "No results", stringsAsFactors = FALSE))
    }
    
    if (length(paper_titles) != length(pmids)) {
      cat("Mismatch in titles and PMIDs for gene:", gene, "\n")
    }
    
    data.frame(
      Gene = gene,
      PMID = pmids[1:length(paper_titles)],
      Title = paper_titles,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("Error occurred for gene:", gene, "\n")
    return(data.frame(Gene = gene, PMID = NA, Title = "Error", stringsAsFactors = FALSE))
  })
}

# -----------------------------
# 3. Load Gene List
# -----------------------------
file_path <- "data/human_protein_coding_genes.xlsx"  # <- update with your file path
genes_table <- read_excel(file_path)
genes <- na.omit(genes_table$Locus)

# -----------------------------
# 4. Fetch PubMed Data
# -----------------------------
results_list <- lapply(genes, function(g) get_pubmed_data(g, max_results = 200))
results_df <- bind_rows(results_list)

# -----------------------------
# 5. Remove Duplicates
# -----------------------------
unique_results <- distinct(results_df, Gene, PMID, .keep_all = TRUE)

# Save initial results
write_xlsx(unique_results, "results/pubmed_cilia_flagella_results_raw.xlsx")

# -----------------------------
# 6. Clean Results
# -----------------------------
cleaned_results <- unique_results %>%
  filter(
    !is.na(PMID),
    tolower(trimws(Title)) != "no results",
    trimws(PMID) != ""
  )

# Save cleaned results
write_xlsx(cleaned_results, "results/pubmed_cili_flagella_cleaned.xlsx")
write.csv(cleaned_results, "results/pubmed_cilia_flagella_cleaned.csv", row.names = FALSE)

# ==========================================================
# End of Script
# ==========================================================
