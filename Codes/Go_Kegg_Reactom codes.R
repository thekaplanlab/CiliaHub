# =========================================
# Enrichment Analizi: Reactome, KEGG ve GO
# =========================================

# -----------------------------
# 1. Gerekli paketler
# -----------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(writexl)

# -----------------------------
# 2. Gen listesini y??kle
# -----------------------------
gene_file <- "C:/Users/lenovo/Desktop/all genes with gene annotations 23.7.2025.xlsx"
gene_list <- read_excel(gene_file)
gene_symbols <- unique(na.omit(gene_list$Gene))  # Gene s??tunundan benzersiz liste

# -----------------------------
# 3. Gene Symbol -> Entrez ID d??n??????m??
# -----------------------------
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first") %>%
  na.omit() %>%
  unique()

# -----------------------------
# 4. Zenginle??tirme analizleri
# -----------------------------
reactome_result <- enrichPathway(gene = entrez_ids, organism = "human", readable = TRUE)
kegg_result    <- enrichKEGG(gene = entrez_ids, organism = "hsa")
go_bp          <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
go_cc          <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)
go_mf          <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)

# DataFrame'e ??evir
reactome_df <- as.data.frame(reactome_result)
kegg_df     <- as.data.frame(kegg_result)
go_bp_df    <- as.data.frame(go_bp)
go_cc_df    <- as.data.frame(go_cc)
go_mf_df    <- as.data.frame(go_mf)

# -----------------------------
# 5. Dotplot fonksiyonu
# -----------------------------
plot_dot <- function(df, title, top_n = 20) {
  top_df <- df %>% arrange(p.adjust) %>% head(top_n)
  
  if(!"FoldEnrichment" %in% colnames(top_df)) {
    top_df$FoldEnrichment <- (top_df$Count / sum(top_df$Count)) / 
      sapply(top_df$BgRatio, function(x) {
        vals <- strsplit(x, "/")[[1]]; as.numeric(vals[1])/as.numeric(vals[2])
      })
  }
  
  ggplot(top_df, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = title,
         x = "Fold Enrichment", y = "Pathway / GO Term",
         color = "-log10(adj p-value)", size = "Gene Count") +
    theme_minimal(base_size = 13) +
    theme(axis.line = element_line(size = 1.2, color = "black"),
          axis.ticks = element_line(size = 1.1, color = "black"),
          panel.grid.minor = element_blank())
}

# -----------------------------
# 6. Dotplotlar?? olu??tur
# -----------------------------
p_reactome <- plot_dot(reactome_df, "All Genes Reactome Pathway Enrichment")
p_kegg     <- plot_dot(kegg_df, "All Genes KEGG Pathway Enrichment")
p_go_bp    <- plot_dot(go_bp_df, "All Genes GO Biological Process")
p_go_cc    <- plot_dot(go_cc_df, "All Genes GO Cellular Component")
p_go_mf    <- plot_dot(go_mf_df, "All Genes GO Molecular Function")

# G??rselle??tir
print(p_reactome)
print(p_kegg)
print(p_go_bp)
print(p_go_cc)
print(p_go_mf)

# -----------------------------
# 7. Pathway-gene e??lemeleri
# -----------------------------
gene_to_pathway <- reactome_df %>%
  dplyr::select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(GeneSymbol = geneID)

gene_pathway_count <- gene_to_pathway %>%
  group_by(GeneSymbol) %>%
  summarise(PathwayCount = n()) %>%
  arrange(desc(PathwayCount))

# -----------------------------
# 8. Sonu??lar?? kaydet
# -----------------------------
write_xlsx(list(
  Reactome = reactome_df,
  KEGG = kegg_df,
  GO_BP = go_bp_df,
  GO_CC = go_cc_df,
  GO_MF = go_mf_df,
  Gene_to_Pathway = gene_to_pathway,
  Pathway_Counts = gene_pathway_count
), "C:/Users/lenovo/Documents/All Genes_Enrichment_Results_All.xlsx")
