# ================================
# 1. Load Libraries
# ================================
library(STRINGdb)
library(readxl)
library(dplyr)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

# ================================
# 2. Functions
# ================================

# Function: Load gene lists
load_gene_lists <- function(known_path, new_path) {
  known_genes <- read_excel(known_path)[[1]] %>% na.omit() %>% unique()
  new_genes   <- read_excel(new_path)[[1]] %>% na.omit() %>% unique()
  list(known = known_genes, new = new_genes)
}

# Function: Map to STRING IDs
map_genes_to_string <- function(gene_list, species = 9606) {
  string_db <- STRINGdb$new(version = "11.5", species = species, score_threshold = 400)
  mapped <- string_db$map(data.frame(gene = gene_list), "gene", removeUnmappedRows = TRUE)
  mapped
}

# Function: Download and load STRING PPI data
load_string_ppi <- function() {
  url <- "https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz"
  file <- basename(url)
  if (!file.exists(file)) {
    message("Downloading STRING PPI data...")
    download.file(url, file)
  }
  fread(file)
}

# Function: Build PPI network for known & new genes
build_network <- function(ppi_data, known_mapped, new_mapped) {
  all_ids <- c(known_mapped$STRING_id, new_mapped$STRING_id)
  focused <- ppi_data[ppi_data$protein1 %in% all_ids & ppi_data$protein2 %in% all_ids, ]
  
  # Combine all mapped with gene symbols
  all_mapped <- rbind(known_mapped, new_mapped) %>% distinct(STRING_id, .keep_all = TRUE)
  
  # Replace STRING IDs with gene symbols
  focused$protein1 <- all_mapped$gene[match(focused$protein1, all_mapped$STRING_id)]
  focused$protein2 <- all_mapped$gene[match(focused$protein2, all_mapped$STRING_id)]
  
  # Create igraph object
  g <- graph_from_data_frame(focused[, .(protein1, protein2)], directed = FALSE)
  
  # Tag groups (Known/New)
  V(g)$group <- ifelse(V(g)$name %in% known_mapped$gene, "Known", "New")
  V(g)$color <- ifelse(V(g)$group == "Known", "#4682B4", "#8B0000")
  
  list(graph = g, edges = focused)
}

# Function: Plot network (fancy layout with labels)
plot_network <- function(graph, output_png, title = "Ciliary Gene PPI Network") {
  g_tbl <- as_tbl_graph(graph)
  ggraph(g_tbl, layout = "fr") +
    geom_edge_link(color = "gray60", alpha = 0.6) + # slightly darker edges
    geom_node_point(aes(color = factor(group)), size = 5) +
    geom_node_text(
      aes(label = name, color = factor(group)),
      repel = TRUE,
      size = 5,
      fontface = "bold"
    ) +
    scale_color_manual(values = c("Known" = "#4682B4", "New" = "#8B0000")) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 22
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    ) -> p
  
  ggsave(output_png, p, width = 14, height = 12, dpi = 300)
  cat("Saved network to", output_png, "\n")
}

# Function: Extract and plot subnetwork for specific genes
extract_subnetwork <- function(graph, target_genes, output_png) {
  target_nodes <- V(graph)[name %in% target_genes]
  neighbors_nodes <- unique(unlist(neighborhood(graph, order = 1, nodes = target_nodes)))
  subg <- induced_subgraph(graph, vids = neighbors_nodes)
  
  plot_network(subg, output_png, title = paste("Subnetwork:", paste(target_genes, collapse = ", ")))
  subg
}

# ================================
# 3. Main Pipeline
# ================================

# Paths to gene lists (update paths)
known_path <- "C:/Users/lenovo/Desktop/merge/known/known gens c elegans_and domain1 and kegg and pathway and go annotations1 and mice and human phenotype1 and tissue and other localization info and seprated tissue.xlsx"
new_path   <- "C:/Users/lenovo/Desktop/new genes with gene annotations 23.7.2025 - Kopya (2).xlsx"


# Load gene lists
genes <- load_gene_lists(known_path, new_path)

# Map to STRING
known_mapped <- map_genes_to_string(genes$known)
new_mapped   <- map_genes_to_string(genes$new)

# Load STRING PPI
ppi_data <- load_string_ppi()

# Build full PPI network
network <- build_network(ppi_data, known_mapped, new_mapped)

# Save edge list with gene names
write.csv(network$edges, "known_new_interactions_gene_names.csv", row.names = FALSE)
cat("Saved interaction list to known_new_interactions_gene_names.csv\n")

# Plot full network
plot_network(network$graph, "ppi_network_known_new_genes.png")

# ================================
# 4. Subnetwork Example (for NEK5 and others)
# ================================
subnet <- extract_subnetwork(network$graph, c("NEK5", "TMEM145", "ZC2HC1A", "ADAMTS20"),
                             "subnetwork_selected_genes.png")

# ================================
# 5. Automatic Family-based Subnetwork Analysis
# ================================

# Define keyword patterns for families
family_patterns <- list(
  Kinases  = c("NEK", "MAPK", "CDK", "PLK", "PK", "CAMK"),
  Dyneins  = c("^DNAH", "^DNAL", "^DNIC"),
  IFT      = c("^IFT"),
  BBSome   = c("^BBS"),
  CEP      = c("^CEP"),
  MKS_NPHP = c("MKS", "NPHP")
)

# Function to detect family members based on gene name patterns
find_family_genes <- function(graph, patterns) {
  all_genes <- V(graph)$name
  matches <- unlist(lapply(patterns, function(p) grep(p, all_genes, value = TRUE, ignore.case = TRUE)))
  unique(matches)
}

# Iterate through each family and extract subnetworks
for (family in names(family_patterns)) {
  family_genes <- find_family_genes(network$graph, family_patterns[[family]])
  
  if (length(family_genes) > 0) {
    message("Found ", length(family_genes), " genes for family: ", family)
    subnet <- extract_subnetwork(network$graph, family_genes,
                                 paste0("subnetwork_", family, ".png"))
    
    # Export edge list for the subnetwork
    sub_edges <- as_data_frame(subnet, what = "edges")
    write.csv(sub_edges, paste0("subnetwork_", family, "_edges.csv"), row.names = FALSE)
  } else {
    message("No matches found for family: ", family)
  }
}




library(ggplot2)
library(ggrepel)

# Example data
set.seed(123)
df <- data.frame(
  Gene = paste0("Gene", 1:20),
  PC1 = rnorm(20),
  PC2 = rnorm(20),
  Category = rep(c("Known", "New"), each = 10)
)

# Plot
plot_network <- function(graph, output_png, title = "Ciliary Gene PPI Network", highlight_genes = NULL) {
  g_tbl <- as_tbl_graph(graph) %>%
    mutate(degree = centrality_degree())
  
  ggraph(g_tbl, layout = "fr") +
    geom_edge_link(color = "gray60", alpha = 0.6) +
    geom_node_point(aes(color = factor(group), size = degree)) +
    geom_node_text(
      aes(label = ifelse(is.null(highlight_genes) | name %in% highlight_genes, name, ""),
          color = factor(group)),
      repel = TRUE,
      size = 5.5,
      fontface = "bold"
    ) +
    scale_color_manual(values = c("Known" = "#4682B4", "New" = "#8B0000")) +
    scale_size(range = c(3, 8), guide = "none") +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    ) -> p
  
  ggsave(output_png, p, width = 14, height = 12, dpi = 300)
  cat("Saved network to", output_png, "\n")
}




kinase_subg <- extract_family_subnetworks(g, edges, kinase_patterns, "Kinases", output_prefix)
dynein_subg <- extract_family_subnetworks(g, edges, dynein_patterns, "Dyneins", output_prefix)
ift_subg    <- extract_family_subnetworks(g, edges, ift_patterns, "IFT", output_prefix)
bbsome_subg <- extract_family_subnetworks(g, edges, bbsome_patterns, "BBSome", output_prefix)
module_subg <- extract_family_subnetworks(g, edges, module_patterns, "CEP_MKS_NPHP", output_prefix)

