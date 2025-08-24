library(rvest)
library(dplyr)
library(writexl)
library(readxl)
library(rentrez)

# Modified function to include "flagella" in the query
get_pubmed_data <- function(gene, max_results = 50) {
  query <- paste0(
    "https://pubmed.ncbi.nlm.nih.gov/?term=",
    gene, "+AND+(cilia+OR+cilium+OR+transition+zone+OR+basal+body+OR+ciliogenesis+OR+flagella+OR+flagellum)",
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

# ??? Step 1: Read gene file (new path and sheet)
file_path <- "C:\\Users\\lenovo\\Downloads\\unique_results_24482 gene.xlsx"
genes_table <- read_excel(file_path)
genes <- na.omit(genes_table$Locus)

# ??? Step 2: Fetch PubMed data
results_list <- lapply(genes, function(g) get_pubmed_data(g, max_results = 200))
results_df <- bind_rows(results_list)

# ??? Step 3: Remove duplicates
unique_results <- distinct(results_df, Gene, PMID, .keep_all = TRUE)

# ??? Step 4: Save initial results
write_xlsx(unique_results, "C:\\Users\\lenovo\\Documents\\pubmed_flagella_results1.xlsx")

# ??? Step 5: Clean - remove empty and "no results"
cleaned_results <- unique_results %>%
  filter(!is.na(PMID), tolower(trimws(Title)) != "no results", trimws(PMID) != "")

# ??? Step 6: Save cleaned results
write_xlsx(cleaned_results, "C:\\Users\\lenovo\\Documents\\pubmed_flagella_cleaned1.xlsx")
write.csv(cleaned_results, "C:\\Users\\lenovo\\Documents\\pubmed_flagella_cleaned.csv", row.names = FALSE)











