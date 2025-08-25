# =========================================
# Ciliopathy Phenotype Analysis using Open Targets
# =========================================

# -----------------------------
# 0. Load necessary packages
# -----------------------------
library(readxl)
library(dplyr)
library(stringr)
library(httr)
library(jsonlite)
library(pbapply)
library(tidyr)
library(ggplot2)
library(writexl)
library(igraph)
library(ggraph)
library(RColorBrewer)

# -----------------------------
# 1. Load gene list
# -----------------------------
gene_list <- read_excel("input your file")

# Keep only the first Ensembl ID per gene
gene_list <- gene_list %>%
  mutate(Gene_ID = sapply(str_split(ensembl_gene_id, ","), function(x) str_trim(x[1]))) %>%
  select(Gene, Gene_ID) %>%
  filter(!is.na(Gene_ID) & Gene_ID != "")

# -----------------------------
# 2. Query Open Targets API
# -----------------------------
get_open_targets_traits <- function(ensembl_id, threshold = 0.5) {
  query <- paste0('{
    target(ensemblId: "', ensembl_id, '") {
      associatedDiseases {
        rows {
          disease { name }
          score
        }
      }
    }
  }')
  
  url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  tryCatch({
    res <- POST(url, body = list(query = query), encode = "json", timeout(10))
    if (status_code(res) != 200) return(NA)
    parsed <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
    diseases <- parsed$data$target$associatedDiseases$rows
    if (is.null(diseases)) return(NA)
    filtered <- diseases %>% filter(score >= threshold)
    if (nrow(filtered) == 0) return(NA)
    paste(unique(filtered$disease$name), collapse = "; ")
  }, error = function(e) NA)
}

# Apply with progress bar
open_targets_traits <- pbsapply(gene_list$Gene_ID, get_open_targets_traits, threshold = 0.5)

# Add results to gene list and save
gene_list$OpenTargets_Diseases <- open_targets_traits
write_xlsx(gene_list, "human_open_targets_results_threshold_0.5.xlsx")

# -----------------------------
# 3. Define ciliopathy-related terms
# -----------------------------
ciliopathy_terms <- tolower(c(
  "abnormal kidney morphology","small kidney","enlarged kidney","abnormal urinary bladder morphology",
  "increased circulating creatinine level","decreased circulating calcium level","ocular phenotypes",
  "microphthalmia","cataract","abnormal retina morphology","abnormal lens morphology",
  "abnormal retina vasculature morphology","abnormal retina blood vessel morphology","eye hemorrhage",
  "abnormal vitreous body morphology","visual impairments","retinal degeneration","structural eye anomalies",
  "abnormal brain morphology","abnormal cranium morphology","abnormal midbrain development",
  "abnormal neural tube closure","abnormal embryo turning","abnormal head shape","hydrocephaly",
  "persistence of hyaloid vascular system","male infertility","female infertility","enlarged testis",
  "small testis","abnormal testis morphology","abnormal uterus morphology","hydrometra",
  "reproductive system abnormalities","skeletal and craniofacial phenotypes","abnormal sternum morphology",
  "abnormal spine curvature","kyphosis","short tibia","abnormal facial morphology","skeletal dysplasias",
  "thoracic and craniofacial anomalies","cardiovascular phenotypes","abnormal heart morphology",
  "enlarged heart","increased heart weight","shortened qrs complex duration","prolonged qrs complex duration",
  "thick ventricular wall","dilated aorta","hepatic and pancreatic phenotypes","abnormal liver morphology",
  "increased liver weight","increased circulating alanine transaminase level",
  "increased circulating aspartate transaminase level","increased circulating bilirubin level",
  "abnormal pancreas morphology","polydactyly","obesity","developmental delay","intellectual disability",
  "hypogonadism","abnormal auditory brainstem response","cyst","hydrocephalus","polydactyly",
  "situs","brain ventricle","flagellum","axoneme","cilium","renal","kidney","fibrosis","dysgenesis",
  "hydrops","left-right axis","left-right patterning","retinal degeneration","congenital heart defect",
  "biliary","hedgehog signaling","infertility","hydrops fetalis","shortened long bones","cerebral anomalies",
  "coronary and vascular anomalies","facial anomalies","ophthalmic anomalies","nasal anomalies",
  "cognitive anomalies","skeletal anomalies","respiratory anomalies","neural anomalies","hormonal anomalies",
  "renal anomalies","digestive anomalies","reproductive anomalies","aural anomalies","liver anomalies",
  "organ anomalies","leber congenital amaurosis","nystagmus","coloboma","hearing loss","respiratory infections",
  "bronchiectasis","sinusitis","otitis media","situs inversus","heterotaxy","congenital heart",
  "nephronophthisis","polycystic kidney","hepatic fibrosis","bile duct","short rib","limb shortening",
  "vertebral anomalies","obesity","diabetes insipidus","joubert","ataxia","seizures","developmental delay",
  "ciliopathy","primary ciliary dyskinesia","bardet-biedl","meckel","senior-loken","oro-facial-digital",
  "retinitis pigmentosa","flagellar","anosmia","cone-rod dystrophy","cystic kidney"
))
ciliopathy_terms <- sort(unique(ciliopathy_terms))

# -----------------------------
# 4. Create phenotype matrix
# -----------------------------
phenotype_results <- setNames(as.list(gene_list$OpenTargets_Diseases), gene_list$Gene)
genes <- names(phenotype_results)

phenotype_matrix <- matrix(FALSE, nrow = length(ciliopathy_terms), ncol = length(genes),
                           dimnames = list(ciliopathy_terms, genes))

for (gene in genes) {
  phenos <- tolower(paste(phenotype_results[[gene]], collapse = ", "))
  for (term in ciliopathy_terms) {
    if (str_detect(phenos, fixed(term))) {
      phenotype_matrix[term, gene] <- TRUE
    }
  }
}

phenotype_df <- as.data.frame(phenotype_matrix) %>%
  tibble::rownames_to_column("Phenotype")

write_xlsx(phenotype_df, "phenotype_matrix_genes_by_ciliopathy_human55.xlsx")

# -----------------------------
# 5. Dot plot: Phenotypes per gene
# -----------------------------
long_df <- phenotype_df %>%
  pivot_longer(-Phenotype, names_to = "Gene", values_to = "Present") %>%
  filter(Present == TRUE)

dotplot <- ggplot(long_df, aes(x = Phenotype, y = Gene)) +
  geom_point(color = "steelblue", size = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7, face = "bold", color = "black"),
    axis.text.y = element_text(size = 4, face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Ciliopathy Phenotypes per Gene in Human", x = "Phenotype", y = "Gene")

ggsave("Ciliopathy_dotplot_human_swapped1.pdf", plot = dotplot, width = 5, height = 10)
ggsave("Ciliopathy_dotplot_human_swapped.png", plot = dotplot, width = 5, height = 10, dpi = 300)

# -----------------------------
# 6. Bipartite Network: Genes vs Phenotypes
# -----------------------------
# Prepare data
phenotype_matrix <- read_xlsx("input phenotype_matrix file")
phenotype_matrix <- as.data.frame(phenotype_matrix)
rownames(phenotype_matrix) <- phenotype_matrix$Phenotype
phenotype_matrix$Phenotype <- NULL

long_df <- phenotype_matrix %>%
  tibble::rownames_to_column("Phenotype") %>%
  pivot_longer(-Phenotype, names_to = "Gene", values_to = "Present") %>%
  filter(Present == TRUE)

# Nodes
nodes_phenotype <- data.frame(name = unique(long_df$Phenotype), type = TRUE)
nodes_gene <- data.frame(name = unique(long_df$Gene), type = FALSE)
nodes <- rbind(nodes_phenotype, nodes_gene)

# Colors
phenotype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(nodes_phenotype))
nodes$color <- ifelse(nodes$type, phenotype_colors, "steelblue")

# Edges
edges <- long_df %>% select(from = Phenotype, to = Gene)
phenotype_colors_edges <- setNames(phenotype_colors, nodes_phenotype$name)
edges$color <- phenotype_colors_edges[edges$from]

# Graph object
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# Plot bipartite network
set.seed(123)
random_genes <- sample(nodes_gene$name, 10)

p <- ggraph(g, layout = "bipartite") +
  geom_edge_link(aes(color = color), alpha = 0.5) +
  geom_node_point(aes(color = I(color)), size = 2) +
  geom_node_text(aes(label = ifelse(type | name %in% random_genes, name, NA)), 
                 color = "black", fontface = "bold", size = 3, repel = TRUE) +
  theme_void() +
  labs(title = "Bipartite Network: Genes and Phenotypes") +
  guides(edge_color = "none")

ggsave("Bipartite_Network_Genes_Phenotypes.png", plot = p, width = 20, height = 10, dpi = 300)

