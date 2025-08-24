# --- Setup: Install and load required packages ---
packages <- c("BiocManager", "biomaRt", "tidyverse", "readxl", "writexl", "ggrepel", "forcats", "purrr", "patchwork")

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")

suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggrepel)
  library(forcats)
  library(purrr)
  library(patchwork)
})

# --- Function: Connect to Ensembl with mirror fallback ---
connect_ensembl <- function() {
  mirrors <- c("useast", "uswest", "asia")
  for (m in mirrors) {
    cat("Trying Ensembl mirror:", m, "...\n")
    ensembl <- tryCatch({
      useEnsembl("genes", "hsapiens_gene_ensembl", mirror = m)
    }, error = function(e) {
      message("Failed on mirror ", m, ": ", e$message)
      NULL
    })
    if (!is.null(ensembl)) {
      cat("Connected to", m, "mirror.\n")
      return(ensembl)
    }
  }
  cat("All live mirrors failed — trying archive version 110...\n")
  ensembl <- tryCatch({
    useEnsembl("genes", "hsapiens_gene_ensembl", version = 110)
  }, error = function(e) {
    stop("Failed to connect to Ensembl archive: ", e$message)
  })
  return(ensembl)
}

ensembl <- connect_ensembl()

# --- Load data ---
load_cilia_genes <- function(file_path) {
  if (!file.exists(file_path)) stop("Excel file not found at ", file_path)
  genes <- read_excel(file_path)[[1]] %>%
    na.omit() %>%
    unique()
  message("Loaded ", length(genes), " ciliary genes from file.")
  return(genes)
}

cilia_genes_file <- "/Users/sebihacevik/Downloads/all gene 1.8.2025.xlsx"
cilia_genes <- load_cilia_genes(cilia_genes_file)

# Fetch all human genes
all_human_genes <- getBM(attributes = c("hgnc_symbol"), mart = ensembl)$hgnc_symbol %>%
  na.omit() %>%
  unique()

# Get domains
get_domains <- function(gene_list, filter_type = "hgnc_symbol", mart = ensembl) {
  if (length(gene_list) == 0) {
    warning("Empty gene list supplied to get_domains()")
    return(tibble())
  }
  df <- tryCatch({
    getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "interpro", "interpro_description"),
      filters = filter_type,
      values = gene_list,
      mart = mart
    ) %>%
      drop_na(interpro, interpro_description) %>%
      filter(interpro_description != "")
  }, error = function(e) {
    message("Error fetching domains: ", e$message)
    tibble()
  })
  return(df)
}

cilia_domains <- get_domains(cilia_genes)
all_domains <- get_domains(all_human_genes)

# --- Prepare domain comparison data ---
cilia_domain_counts <- cilia_domains %>%
  count(interpro, interpro_description, name = "CiliaGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

all_domain_counts <- all_domains %>%
  count(interpro, interpro_description, name = "TotalGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

domain_comparison <- full_join(cilia_domain_counts, all_domain_counts, by = c("DomainID", "DomainDescription")) %>%
  replace_na(list(CiliaGeneCount = 0, TotalGeneCount = 0))

# Statistical analysis
total_cilia_genes <- length(cilia_genes)
total_human_genes <- length(all_human_genes)

domain_comparison <- domain_comparison %>%
  mutate(
    NonCiliaGeneCount = TotalGeneCount - CiliaGeneCount,
    NonCiliaNonDomain = total_human_genes - TotalGeneCount,
    PValue_Enrichment = map2_dbl(CiliaGeneCount, TotalGeneCount, ~ {
      matrix <- matrix(c(.x, .y - .x, total_cilia_genes - .x, total_human_genes - .y - (total_cilia_genes - .x)), nrow = 2)
      fisher.test(matrix, alternative = "greater")$p.value
    }),
    AdjPValue_Enrichment = p.adjust(PValue_Enrichment, method = "BH"),
    PValue_Depletion = map2_dbl(CiliaGeneCount, TotalGeneCount, ~ {
      matrix <- matrix(c(.x, .y - .x, total_cilia_genes - .x, total_human_genes - .y - (total_cilia_genes - .x)), nrow = 2)
      fisher.test(matrix, alternative = "less")$p.value
    }),
    AdjPValue_Depletion = p.adjust(PValue_Depletion, method = "BH"),
    FoldEnrichment = (CiliaGeneCount / total_cilia_genes) / (TotalGeneCount / total_human_genes),
    Significance = case_when(
      AdjPValue_Enrichment < 0.05 & FoldEnrichment > 1 ~ "Significantly Enriched",
      AdjPValue_Depletion < 0.05 & FoldEnrichment < 1 ~ "Significantly Depleted",
      TRUE ~ "Not Significant"
    ),
    log2FE = log2(pmax(FoldEnrichment, 1e-6)),
    negLog10P = -log10(pmin(AdjPValue_Enrichment, AdjPValue_Depletion, na.rm = TRUE))
  )

# --- Figure 2A: All genes vs ciliary genes (Donut plot) ---
genes_df <- tibble(
  Category = c("Ciliary Genes", "Non-Ciliary Genes"),
  Count = c(total_cilia_genes, total_human_genes - total_cilia_genes)
) %>%
  mutate(
    Fraction = Count / sum(Count),
    ymax = cumsum(Fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2,
    label = paste0(Count)
  )

fig2a <- ggplot(genes_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(aes(x = 3.5, y = label_pos, label = label), color = "black", size = 5) +
  coord_polar(theta = "y") +
  xlim(2, 4) +
  scale_fill_manual(values = c("#86A3C3", "#B6CEC7")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 12)
  )

# --- Figure 2B: All domains vs ciliary domains (Donut plot) ---
domains_df <- tibble(
  Category = c("Domains in Ciliary Genes", "Domains NOT in Ciliary Genes"),
  Count = c(n_distinct(cilia_domains$interpro), 
            n_distinct(all_domains$interpro) - n_distinct(cilia_domains$interpro))
) %>%
  mutate(
    Fraction = Count / sum(Count),
    ymax = cumsum(Fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2,
    label = paste0(Count)
  )

fig2b <- ggplot(domains_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(aes(x = 3.5, y = label_pos, label = label), color = "black", size = 5) +
  coord_polar(theta = "y") +
  xlim(2, 4) +
  scale_fill_manual(values = c("#D8E0BB", "#B6CEC7")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 12)
  )

# --- Figure 2C: Distribution of Domains ---
stacked_df <- bind_rows(
  all_domain_counts %>%
    mutate(Group = "All Human", Count = TotalGeneCount) %>%
    select(DomainID, DomainDescription, Group, Count),
  cilia_domain_counts %>%
    mutate(Group = "Ciliary", Count = CiliaGeneCount) %>%
    select(DomainID, DomainDescription, Group, Count)
) %>%
  replace_na(list(Count = 0))

total_counts_df <- stacked_df %>%
  group_by(DomainID, DomainDescription) %>%
  summarise(TotalCount = sum(Count), .groups = "drop")

top_domains <- total_counts_df %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 20) %>%
  pull(DomainID)

stacked_df_top <- stacked_df %>%
  filter(DomainID %in% top_domains) %>%
  left_join(total_counts_df, by = c("DomainID", "DomainDescription")) %>%
  mutate(DomainDescription = factor(
    DomainDescription,
    levels = total_counts_df %>%
      arrange(desc(TotalCount)) %>%
      pull(DomainDescription)
  ))

fig2c <- ggplot(stacked_df_top, aes(x = Count, y = DomainDescription, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Ciliary" = "#86A3C3", "All Human" = "#B6CEC7")) +
  labs(x = "Number of Genes with Domain", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

# --- Figure 2D: Top20 Enriched Domains ---
top_enriched <- domain_comparison %>%
  filter(Significance == "Significantly Enriched") %>%
  arrange(desc(FoldEnrichment)) %>%
  slice_head(n = 20)

fig2d <- ggplot(top_enriched, aes(x = FoldEnrichment, y = fct_reorder(DomainDescription, FoldEnrichment))) +
  geom_col(fill = "#86A3C3") +
  labs(x = "Fold Enrichment", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 9)
  )

# --- Figure 2E: Top20 Deprived Domains ---
top_deprived <- domain_comparison %>%
  filter(CiliaGeneCount == 0) %>%
  arrange(desc(TotalGeneCount)) %>%
  slice_head(n = 20)

fig2e <- ggplot(top_deprived, aes(x = TotalGeneCount, y = fct_reorder(DomainDescription, TotalGeneCount))) +
  geom_col(fill = "#B6CEC7") +
  labs(x = "Number of Non-Ciliary Genes", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 9)
  )

# --- Figure 2F: Volcano Plot ---
fig2f <- ggplot(domain_comparison, aes(x = log2FE, y = negLog10P)) +
  geom_point(aes(color = Significance, alpha = Significance, shape = Significance),
             size = 3.5, position = position_jitter(width = 0.05, height = 0)) +
  scale_color_manual(values = c(
    "Significantly Enriched" = "#0072B2",
    "Significantly Depleted" = "#E69F00",
    "Not Significant" = "grey70"
  )) +
  scale_alpha_manual(values = c(
    "Significantly Enriched" = 1,
    "Significantly Depleted" = 1,
    "Not Significant" = 0.3
  )) +
  scale_shape_manual(values = c(
    "Significantly Enriched" = 16,
    "Significantly Depleted" = 17,
    "Not Significant" = 1
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = expression(Log[2]~Fold~Enrichment), y = expression(-Log[10]~P~Value)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  )

# --- Save all figures ---
ggsave("Figure2A.png", fig2a, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("Figure2B.png", fig2b, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("Figure2C.png", fig2c, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Figure2D.png", fig2d, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Figure2E.png", fig2e, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Figure2F.png", fig2f, width = 8, height = 6, dpi = 300, bg = "white")

# --- Save session info ---
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



# Count cilia-specific domains (significantly enriched)
cilia_specific_count <- domain_comparison %>%
  filter(AdjPValue_Enrichment < 0.05, FoldEnrichment > 1) %>%
  nrow()

# Count cilia-deprived domains (significantly depleted)
cilia_deprived_count <- domain_comparison %>%
  filter(AdjPValue_Depletion < 0.05, FoldEnrichment < 1) %>%
  nrow()

cat("Number of cilia-specific domains (enriched):", cilia_specific_count, "\n")
cat("Number of cilia-deprived domains:", cilia_deprived_count, "\n")

# Get complete lists with gene names
# 1. Cilia-specific domains (enriched)
enriched_domains_with_genes <- domain_comparison %>%
  filter(AdjPValue_Enrichment < 0.05, FoldEnrichment > 1) %>%
  left_join(
    cilia_domains %>%
      select(interpro, interpro_description, hgnc_symbol) %>%
      group_by(interpro, interpro_description) %>%
      summarise(Genes = paste(unique(hgnc_symbol), collapse = ", "), 
                by = c("DomainID" = "interpro", "DomainDescription" = "interpro_description")
      ) %>%
      arrange(desc(FoldEnrichment)) %>%
      select(DomainID, DomainDescription, FoldEnrichment, AdjPValue_Enrichment, Genes)
    
    # 2. Cilia-deprived domains (absent or depleted)
    deprived_domains_with_genes <- domain_comparison %>%
      filter(AdjPValue_Depletion < 0.05 | CiliaGeneCount == 0) %>%
      left_join(
        all_domains %>%
          select(interpro, interpro_description, hgnc_symbol) %>%
          group_by(interpro, interpro_description) %>%
          summarise(Genes = paste(unique(hgnc_symbol), collapse = ", "), 
                    by = c("DomainID" = "interpro", "DomainDescription" = "interpro_description")
          ) %>%
          arrange(desc(TotalGeneCount)) %>%
          select(DomainID, DomainDescription, TotalGeneCount, CiliaGeneCount, AdjPValue_Depletion, Genes)
        
        # Print top 10 of each
        cat("\nTop 10 cilia-specific domains:\n")
        print(enriched_domains_with_genes %>% head(10))
        
        cat("\nTop 10 cilia-deprived domains:\n")
        print(deprived_domains_with_genes %>% head(10))
        
        # Save full lists
        write_xlsx(list(
          "Enriched_Domains" = enriched_domains_with_genes,
          "Deprived_Domains" = deprived_domains_with_genes
        ), "Cilia_Domains_With_Genes.xlsx")
        
        # Visualization of counts
        count_data <- tibble(
          Category = c("Cilia-Specific", "Cilia-Deprived"),
          Count = c(cilia_specific_count, cilia_deprived_count)
        )
        
        ggplot(count_data, aes(x = Category, y = Count, fill = Category)) +
          geom_col() +
          geom_text(aes(label = Count), vjust = -0.5) +
          scale_fill_manual(values = c("#86A3C3", "#B6CEC7")) +
          labs(title = "Domain Enrichment/Depletion in Ciliary Genes",
               x = "", y = "Number of Domains") +
          theme_minimal()