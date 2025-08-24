# --- Setup: Install and load required packages ---
packages <- c("BiocManager", "biomaRt", "tidyverse", "readxl", "writexl", "ggrepel", "forcats", "purrr")

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

# --- Function: Load ciliary genes from Excel ---
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

# --- Fetch all human genes (HGNC symbols) ---
message("Fetching all human genes from Ensembl...")
all_human_genes <- getBM(attributes = c("hgnc_symbol"), mart = ensembl)$hgnc_symbol %>%
  na.omit() %>%
  unique()
message("Fetched ", length(all_human_genes), " human genes.")

# Basic counts (replace these with your actual variables)
total_human_genes <- length(all_human_genes)
total_human_genes
total_cilia_genes <- length(cilia_genes)
total_cilia_genes

unique_domains_all <- all_domains %>% distinct(interpro) %>% nrow()
unique_domains_all
unique_domains_cilia <- cilia_domains %>% distinct(interpro) %>% nrow()
unique_domains_cilia

percent_domains_cilia <- round(unique_domains_cilia / unique_domains_all * 100, 1)
percent_domains_cilia

# Domains absent from ciliary genes
domains_not_in_cilia <- domain_comparison %>% filter(CiliaGeneCount == 0)
domains_not_in_cilia
num_domains_not_in_cilia <- nrow(domains_not_in_cilia)
num_domains_not_in_cilia

# Domains significantly enriched in ciliary genes
domains_enriched <- domain_comparison %>%
  filter(AdjPValue_Enrichment < 0.05, FoldEnrichment > 1)
num_domains_enriched <- nrow(domains_enriched)
num_domains_enriched
# Barplot: Number of domains present vs absent in ciliary genes
library(ggplot2)
bar_df <- tibble(
  Category = c("Domains in Ciliary Genes", "Domains NOT in Ciliary Genes"),
  Count = c(unique_domains_cilia, num_domains_not_in_cilia)
)

p_bar <- ggplot(bar_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(title = "Representation of Protein Domains in Ciliary Genes",
       y = "Number of Domains", x = NULL) +
  theme_minimal(base_size = 14)

print(p_bar)

# Barplot: Number of enriched vs not enriched domains
enrich_df <- tibble(
  Category = c("Significantly Enriched Domains", "Other Domains"),
  Count = c(num_domains_enriched, unique_domains_cilia - num_domains_enriched)
)

p_enrich <- ggplot(enrich_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("darkgreen", "gray70")) +
  labs(title = "Significantly Enriched Domains in Ciliary Genes",
       y = "Number of Domains", x = NULL) +
  theme_minimal(base_size = 14)

print(p_enrich)


# --- Function: Retrieve InterPro domains for given gene list ---
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

message("Fetching InterPro domains for ciliary genes...")
cilia_domains <- get_domains(cilia_genes)

message("Fetching InterPro domains for all human genes...")
all_domains <- get_domains(all_human_genes)
all_domains
message("Fetched ", nrow(cilia_domains), " ciliary domain entries and ", nrow(all_domains), " total domain entries.")

# --- Summarize domain counts ---
cilia_domain_counts <- cilia_domains %>%
  count(interpro, interpro_description, name = "CiliaGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

all_domain_counts <- all_domains %>%
  count(interpro, interpro_description, name = "TotalGeneCount") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

# Merge counts; fill NAs with 0
domain_comparison <- full_join(cilia_domain_counts, all_domain_counts, by = c("DomainID", "DomainDescription")) %>%
  replace_na(list(CiliaGeneCount = 0, TotalGeneCount = 0))

# --- Statistical analysis: Fisher's exact test for enrichment and depletion ---
total_cilia_genes <- length(cilia_genes)
total_human_genes <- length(all_human_genes)
total_cilia_genes
domain_comparison <- domain_comparison %>%
  mutate(
    NonCiliaGeneCount = TotalGeneCount - CiliaGeneCount,
    NonCiliaNonDomain = total_human_genes - TotalGeneCount,
    PValue_Enrichment = map2_dbl(CiliaGeneCount, TotalGeneCount, ~ {
      # Contingency table:
      #       Domain    Not Domain
      # CiliaGeneCount       x           Cilia genes without domain
      matrix <- matrix(c(.x, .y - .x, total_cilia_genes - .x, total_human_genes - .y - (total_cilia_genes - .x)), nrow = 2)
      fisher.test(matrix, alternative = "greater")$p.value
    }),
    AdjPValue_Enrichment = p.adjust(PValue_Enrichment, method = "BH"),
    PValue_Depletion = map2_dbl(CiliaGeneCount, TotalGeneCount, ~ {
      matrix <- matrix(c(.x, .y - .x, total_cilia_genes - .x, total_human_genes - .y - (total_cilia_genes - .x)), nrow = 2)
      fisher.test(matrix, alternative = "less")$p.value
    }),
    AdjPValue_Depletion = p.adjust(PValue_Depletion, method = "BH"),
    FoldEnrichment = (CiliaGeneCount / total_cilia_genes) / (TotalGeneCount / total_human_genes)
  )

# --- Identify top 20 cilia-specific domains (significant enrichment) ---
cilia_specific_domains <- domain_comparison %>%
  filter(CiliaGeneCount > 0, AdjPValue_Enrichment < 0.05) %>%
  arrange(AdjPValue_Enrichment) %>%
  slice_head(n = 20)

# Get gene names for cilia-specific domains
cilia_domain_genes <- cilia_domains %>%
  filter(interpro %in% cilia_specific_domains$DomainID) %>%
  select(interpro, interpro_description, hgnc_symbol) %>%
  distinct() %>%
  group_by(interpro, interpro_description) %>%
  summarise(CiliaGeneNames = paste(hgnc_symbol, collapse = ", "), .groups = "drop") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

# Merge gene names with cilia-specific domains
cilia_specific_domains <- left_join(cilia_specific_domains, cilia_domain_genes, by = c("DomainID", "DomainDescription")) %>%
  select(DomainID, DomainDescription, CiliaGeneCount, TotalGeneCount, FoldEnrichment, AdjPValue_Enrichment, CiliaGeneNames)

cilia_specific_domains

# --- Identify top 20 deprived domains (absent in ciliary genes) ---
deprived_domains <- domain_comparison %>%
  filter(CiliaGeneCount == 0) %>%
  arrange(desc(TotalGeneCount)) %>%
  slice_head(n = 20) %>%
  select(DomainID, DomainDescription, TotalGeneCount)
deprived_domains 
# --- Get gene names for all domains ---
all_domain_genes <- all_domains %>%
  select(interpro, interpro_description, hgnc_symbol) %>%
  distinct() %>%
  group_by(interpro, interpro_description) %>%
  summarise(AllGeneNames = paste(hgnc_symbol, collapse = ", "), .groups = "drop") %>%
  rename(DomainID = interpro, DomainDescription = interpro_description)

# Merge all domain gene names and add cilia gene names (fill NAs)
domain_comparison_with_genes <- left_join(domain_comparison, all_domain_genes, by = c("DomainID", "DomainDescription")) %>%
  left_join(cilia_domain_genes, by = c("DomainID", "DomainDescription")) %>%
  mutate(CiliaGeneNames = replace_na(CiliaGeneNames, ""))

# --- Output files (edit paths as needed) ---
output_folder <- "./"  # Change to your desired output folder, e.g. "/Users/sebihacevik/Documents/"

write_xlsx(cilia_specific_domains, file.path(output_folder, "Top20_CiliaSpecific_Domains_WithGenes.xlsx"))
write_xlsx(deprived_domains, file.path(output_folder, "Top20_Deprived_Domains.xlsx"))
write_xlsx(domain_comparison_with_genes, file.path(output_folder, "All_Domains_WithGenes.xlsx"))

# --- Plots ---

library(tidyverse)

# Total gene counts
total_human_genes <- length(all_human_genes)
total_cilia_genes <- length(cilia_genes)
total_human_genes
total_cilia_genes
# Domain counts (unique InterPro domains)
unique_domains_all <- all_domains %>% 
  distinct(interpro) %>%
  nrow()

unique_domains_cilia <- cilia_domains %>%
  distinct(interpro) %>%
  nrow()

# Prepare data for percentage stacked bar plot
combined_df <- tibble(
  Type = rep(c("Genes", "Domains"), each = 2),
  Group = c("All Human", "Ciliary", "All Human", "Ciliary"),
  Count = c(total_human_genes, total_cilia_genes, unique_domains_all, unique_domains_cilia)
)

# Calculate percentages per Type
combined_df <- combined_df %>%
  group_by(Type) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# Plot
p_percent <- ggplot(combined_df, aes(x = Type, y = Percent, fill = Group)) +
  geom_col(position = "fill", color = "black", size = 0.3) +  # position="fill" normalizes to 100%
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Show % on y axis
  scale_fill_manual(values = c("All Human" = "grey70", "Ciliary" = "dodgerblue")) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_fill(vjust = 0.5),
            color = "black",
            size = 5) +
  labs(
    title = "Percentage of Genes and Domains: Human vs Ciliary Genes",
    x = NULL,
    y = "Percentage (%)",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 13),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_percent)



library(dplyr)
library(forcats)

# After filtering to top domains, join the TotalCount info back:
# Create stacked_df with domain counts for both groups
# 1️⃣ Build stacked_df from your existing counts
stacked_df <- bind_rows(
  all_domain_counts %>%
    mutate(Group = "All Human", Count = TotalGeneCount) %>%
    select(DomainID, DomainDescription, Group, Count),
  cilia_domain_counts %>%
    mutate(Group = "Ciliary", Count = CiliaGeneCount) %>%
    select(DomainID, DomainDescription, Group, Count)
)

# 2️⃣ Replace NAs with 0
stacked_df <- stacked_df %>%
  replace_na(list(Count = 0))

# 3️⃣ Calculate total counts per domain
total_counts_df <- stacked_df %>%
  group_by(DomainID, DomainDescription) %>%
  summarise(TotalCount = sum(Count), .groups = "drop")

# 4️⃣ Select top domains
top_domains <- total_counts_df %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 20) %>%
  pull(DomainID)

# 5️⃣ Filter stacked_df for top domains and order factor levels
stacked_df_top <- stacked_df %>%
  filter(DomainID %in% top_domains) %>%
  left_join(total_counts_df, by = c("DomainID", "DomainDescription")) %>%
  mutate(DomainDescription = factor(
    DomainDescription,
    levels = total_counts_df %>%
      arrange(desc(TotalCount)) %>%
      pull(DomainDescription)
  ))


# Now you can plot:
p_stack <- ggplot(stacked_df_top, aes(x = Count, y = DomainDescription, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Ciliary" = "dodgerblue", "All Human" = "grey70")) +
  labs(
    title = "Distribution of Domains in Ciliary vs All Human Genes",
    x = "Number of Genes with Domain",
    y = "InterPro Domain",
    fill = "Gene Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))


print(p_stack)


ggsave(file.path(output_folder, "Domain_comparisons.pdf"), p_stack, width = 9, height = 6, device = cairo_pdf)




# Plot 1: Fold enrichment for top 20 cilia-specific domains
p1 <- ggplot(cilia_specific_domains, aes(x = FoldEnrichment, y = fct_reorder(DomainDescription, FoldEnrichment))) +
  geom_col(fill = "skyblue") +
  labs(
    title = "Top 20 Cilia-Specific Domains by Fold Enrichment",
    x = "Fold Enrichment",
    y = "InterPro Domain"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))
p1
ggsave(file.path(output_folder, "Top20_CiliaSpecific_Domains_FoldEnrichment.pdf"), p1, width = 10, height = 6, device = cairo_pdf)

# Plot 2: Gene counts for top 20 deprived domains
p2 <- ggplot(deprived_domains, aes(x = TotalGeneCount, y = fct_reorder(DomainDescription, TotalGeneCount))) +
  geom_col(fill = "#B6CEC7") +
  labs(
    title = "Top 20 Domains Deprived in Ciliary Genes",
    x = "Number of Non-Ciliary Genes",
    y = "InterPro Domain"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))
p2 
ggsave(file.path(output_folder, "Top20_Deprived_Domains.pdf"), p2, width = 10, height = 6, device = cairo_pdf)

# Plot 3: Volcano plot showing enrichment vs. depletion
# Create a new categorical variable for significance and direction
domain_comparison <- domain_comparison %>%
  mutate(
    Significance = case_when(
      AdjPValue_Enrichment < 0.05 & FoldEnrichment > 1 ~ "Significantly Enriched",
      AdjPValue_Depletion < 0.05 & FoldEnrichment < 1 ~ "Significantly Depleted",
      TRUE ~ "Not Significant"
    )
  )



summary(domain_comparison$FoldEnrichment)
summary(domain_comparison$AdjPValue_Enrichment)
summary(domain_comparison$AdjPValue_Depletion)
table(domain_comparison$Significance, useNA = "ifany")

# Volcano plot with clearer significance groups and colors
library(dplyr)
library(ggplot2)

# Avoid -Inf by replacing 0 or negative FE with a tiny number
domain_comparison <- domain_comparison %>%
  mutate(
    log2FE = log2(pmax(FoldEnrichment, 1e-6)),
    negLog10P = -log10(pmin(AdjPValue_Enrichment, AdjPValue_Depletion, na.rm = TRUE))
  )

library(ggplot2)

# Thresholds
fc_thresh <- 1
pval_thresh <- 0.05


colnames(domain_comparison)
table(domain_comparison$Significance)

# Extract domain lists by Significance
enriched_domains <- domain_comparison$DomainDescription[domain_comparison$Significance == "Significantly Enriched"]
depleted_domains <- domain_comparison$DomainDescription[domain_comparison$Significance == "Significantly Depleted"]
not_significant_domains <- domain_comparison$DomainDescription[domain_comparison$Significance == "Not Significant"]

# Print to check
print("Significantly Enriched Domains:")
print(enriched_domains)

print("Significantly Depleted Domains:")
print(depleted_domains)

print("Not Significant Domains:")
print(not_significant_domains)

writeLines(enriched_domains, "enriched_domains.txt")
writeLines(depleted_domains, "depleted_domains.txt")
writeLines(not_significant_domains, "not_significant_domains.txt")


p_pub <- ggplot(domain_comparison, aes(x = log2FE, y = negLog10P)) +
  geom_point(aes(color = Significance, alpha = Significance, shape = Significance),
             size = 3.5, position = position_jitter(width = 0.05, height = 0)) +
  scale_color_manual(values = c(
    "Significantly Enriched" = "#0072B2",  # Blue
    "Significantly Depleted" = "#E69F00",  # Orange
    "Not Significant" = "grey70"
  )) +
  scale_alpha_manual(values = c(
    "Significantly Enriched" = 1,
    "Significantly Depleted" = 1,
    "Not Significant" = 0.3
  )) +
  scale_shape_manual(values = c(
    "Significantly Enriched" = 16,  # filled circle
    "Significantly Depleted" = 17,  # filled triangle
    "Not Significant" = 1           # hollow circle
  )) +
  coord_cartesian(
    xlim = c(min(domain_comparison$log2FE, na.rm = TRUE) - 1,
             max(domain_comparison$log2FE, na.rm = TRUE) + 1)
  ) +
  labs(
    title = "Domain Enrichment vs. Depletion in Ciliary Genes",
    x = expression(Log[2]~Fold~Enrichment),
    y = expression(-Log[10]~Adjusted~P~Value),
    color = "Significance",
    alpha = "Significance",
    shape = "Significance"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dotted", color = "black") +
  geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

library(ggrepel)

top_labels <- domain_comparison %>%
  filter(Significance != "Not Significant") %>%
  arrange(desc(abs(log2FE))) %>%
  slice_head(n = 10)  # label top 15 by effect size

p_pub +
  geom_text_repel(
    data = top_labels,
    aes(label = DomainDescription),
    size = 3,
    max.overlaps = 10,
    box.padding = 0.3,
    point.padding = 0.5,
    segment.color = 'grey50'
  )



print(p_pub)
# Save high-quality outputs
ggsave("Volcano_Domain_Enrichment_Publication.pdf", p_pub, width = 8, height = 6)
ggsave("Volcano_Domain_Enrichment_Publication.png", p_pub, width = 8, height = 6, dpi = 600)

# --- Save session info for reproducibility ---
writeLines(capture.output(sessionInfo()), file.path(output_folder, "sessionInfo.txt"))

message("Analysis complete. Output files saved in: ", normalizePath(output_folder))

library(dplyr)
library(ggplot2)
library(forcats)

# Combine top enriched and depleted domains
top_domains <- domain_comparison %>%
  filter(Significance != "Not Significant") %>%
  arrange(Significance, desc(abs(log2FE))) %>%
  slice_head(n = 20) %>%
  pull(DomainID)

# Prepare data for plotting counts side by side
counts_df <- bind_rows(
  cilia_domain_counts %>% filter(DomainID %in% top_domains) %>% mutate(Group = "Ciliary", Count = CiliaGeneCount),
  all_domain_counts %>% filter(DomainID %in% top_domains) %>% mutate(Group = "All Human", Count = TotalGeneCount)
) %>%
  select(DomainID, DomainDescription, Group, Count) %>%
  distinct()

counts_df$DomainDescription <- factor(counts_df$DomainDescription,
                                      levels = counts_df %>%
                                        filter(Group == "All Human") %>%
                                        arrange(Count) %>%
                                        pull(DomainDescription))

ggplot(counts_df, aes(x = Count, y = fct_rev(DomainDescription), fill = Group)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Ciliary" = "dodgerblue", "All Human" = "grey70")) +
  labs(
    title = "Gene Counts in Top Enriched and Depleted Domains",
    x = "Number of Genes with Domain",
    y = "Domain Description",
    fill = "Gene Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))



library(tidyr)

# Create a gene-domain presence dataset (domain by gene)
# For simplicity, summarize by counts of genes with domain in each group
gene_domain_presence <- bind_rows(
  cilia_domains %>% mutate(Group = "Ciliary"),
  all_domains %>% mutate(Group = "All Human")
)

# Filter to top domains (enriched + depleted)
top_domain_ids <- domain_comparison %>%
  filter(Significance != "Not Significant") %>%
  arrange(Significance, desc(abs(log2FE))) %>%
  slice_head(n = 20) %>%
  pull(DomainID)

gene_domain_presence <- gene_domain_presence %>%
  filter(interpro %in% top_domain_ids) %>%
  group_by(Group, interpro, interpro_description) %>%
  summarise(GeneCount = n_distinct(hgnc_symbol), .groups = "drop")

gene_domain_presence$interpro_description <- factor(gene_domain_presence$interpro_description,
                                                    levels = gene_domain_presence %>%
                                                      filter(Group == "All Human") %>%
                                                      arrange(GeneCount) %>%
                                                      pull(interpro_description))

ggplot(gene_domain_presence, aes(x = Group, y = fct_rev(interpro_description), size = GeneCount, color = Group)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Ciliary" = "dodgerblue", "All Human" = "grey70")) +
  labs(
    title = "Dotplot of Gene Counts for Top Enriched and Depleted Domains",
    x = "Gene Group",
    y = "Domain Description",
    size = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))


# Example: Load or create domain_function_df with columns: DomainID, FunctionCategory
# domain_function_df <- read_csv("domain_function_annotations.csv")

# Dummy example:
domain_function_df <- tibble(
  DomainID = c("IPR000123", "IPR000456", "IPR000789"),  # replace with real DomainIDs
  FunctionCategory = c("Signal Transduction", "Structural", "Enzymatic Activity")
)

# Merge
domain_comparison <- domain_comparison %>%
  left_join(domain_function_df, by = "DomainID") %>%
  mutate(FunctionCategory = replace_na(FunctionCategory, "Unknown"))

# Example: Color volcano plot by FunctionCategory
ggplot(domain_comparison, aes(x = log2FE, y = negLog10P, color = FunctionCategory)) +
  geom_point(alpha = 0.7) +
  labs(title = "Domain Enrichment Volcano Plot Colored by Function Category") +
  theme_minimal()


#Summary table of top enriched domains with key genes
library(knitr)
library(kableExtra)

# Show nicely in RMarkdown or console
cilia_specific_domains %>%
  select(DomainDescription, CiliaGeneCount, TotalGeneCount, FoldEnrichment, AdjPValue_Enrichment, CiliaGeneNames) %>%
  arrange(AdjPValue_Enrichment) %>%
  slice_head(n = 10) %>%
  kable(caption = "Top 10 Enriched Ciliary Domains with Genes") %>%
  kable_styling(full_width = FALSE)

#Report domain overlaps with known ciliary protein complexes

# Example: known_ciliary_complexes_df with DomainID and ComplexName columns
known_ciliary_complexes_df <- tibble(
  DomainID = c("IPR001234", "IPR005678"),  # replace with your known complex domain IDs
  ComplexName = c("IFT Complex", "BBSome")
)

# Annotate domain_comparison
domain_comparison <- domain_comparison %>%
  left_join(known_ciliary_complexes_df, by = "DomainID")

# Show which enriched domains belong to known complexes
enriched_with_complex <- domain_comparison %>%
  filter(Significance == "Significantly Enriched", !is.na(ComplexName)) %>%
  select(DomainDescription, FoldEnrichment, AdjPValue_Enrichment, ComplexName)

print(enriched_with_complex)

# IFT-A complex (6 proteins)
IFT_A_genes <- c("IFT144", "IFT140", "IFT121", "IFT120", "IFT43")

# IFT-B complex (16 proteins)
IFT_B_genes <- c(
  "IFT172", "IFT88", "IFT81", "IFT80", "IFT74", "IFT70", "IFT56", "IFT54", 
  "IFT57", "IFT52", "IFT46", "IFT38", "IFT27", "IFT25", "IFT22", "IFT20"
)

# BBSome complex (8 proteins)
BBSome_genes <- c("BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS8", "BBS9", "BBIP1")

get_complex_domains <- function(genes, mart, complex_name) {
  getBM(
    attributes = c("hgnc_symbol", "interpro", "interpro_description"),
    filters = "hgnc_symbol",
    values = genes,
    mart = mart
  ) %>%
    distinct(interpro, interpro_description) %>%
    filter(!is.na(interpro)) %>%
    mutate(Complex = complex_name)
}

IFT_A_domains <- get_complex_domains(IFT_A_genes, ensembl, "IFT-A")
IFT_B_domains <- get_complex_domains(IFT_B_genes, ensembl, "IFT-B")
BBSome_domains <- get_complex_domains(BBSome_genes, ensembl, "BBSome")

complex_domains <- bind_rows(IFT_A_domains, IFT_B_domains, BBSome_domains) %>%
  distinct(interpro, Complex)

domain_comparison <- domain_comparison %>%
  left_join(complex_domains, by = c("DomainID" = "interpro"))

enriched_complex_domains <- domain_comparison %>%
  filter(Significance == "Significantly Enriched", !is.na(Complex)) %>%
  select(DomainID, DomainDescription, FoldEnrichment, AdjPValue_Enrichment, Complex)

print(enriched_complex_domains)

library(ggrepel)

p_pub +
  geom_point(data = filter(domain_comparison, !is.na(Complex)), 
             aes(shape = Complex), color = "red", size = 3, alpha = 0.8) +
  geom_text_repel(data = filter(domain_comparison, !is.na(Complex) & Significance == "Significantly Enriched"),
                  aes(label = DomainDescription),
                  size = 3,
                  box.padding = 0.3,
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  scale_shape_manual(values = c(`IFT-A` = 15, `IFT-B` = 17, BBSome = 18))


library(gggenes)
ggplot(cilia_domains, aes(xmin = start, xmax = end, y = hgnc_symbol, 
                          fill = interpro_description)) +
  geom_gene_arrow() +
  facet_wrap(~interpro_description, scales = "free_y", ncol = 1)
