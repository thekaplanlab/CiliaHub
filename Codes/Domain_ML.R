library(biomaRt)
library(tidyverse)
library(readxl)
library(caret)
library(pROC)

# Load data and connect to Ensembl
cilia_genes <- read_excel("/Users/sebihacevik/Downloads/all gene 1.8.2025.xlsx")[[1]] %>%
  na.omit() %>% unique()

ensembl <- useEnsembl("genes", "hsapiens_gene_ensembl")
all_human_genes <- getBM(attributes = "hgnc_symbol", mart = ensembl)$hgnc_symbol %>%
  na.omit() %>% unique()

# Fetch domains
get_domains <- function(genes) {
  getBM(attributes = c("hgnc_symbol", "interpro"), 
        filters = "hgnc_symbol", values = genes, mart = ensembl) %>%
    drop_na() %>% filter(interpro != "")
}

# Create binary matrix
gene_domains <- bind_rows(
  get_domains(cilia_genes) %>% mutate(IsCiliary = 1),
  get_domains(setdiff(all_human_genes, cilia_genes)) %>% mutate(IsCiliary = 0)
) %>%
  distinct(hgnc_symbol, interpro, .keep_all = TRUE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = interpro, values_from = present, values_fill = 0)

# Prepare features and split data
X <- gene_domains %>% select(-hgnc_symbol, -IsCiliary)
y <- factor(ifelse(gene_domains$IsCiliary == 1, "Ciliary", "NonCiliary"))

set.seed(123)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)

# Random Forest training with time estimation
cat("Starting Random Forest training...\n")
start_time <- Sys.time()

rf_model <- train(X[train_idx,], y[train_idx], method = "rf", 
                  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE),
                  tuneLength = 3)

end_time <- Sys.time()
training_time <- difftime(end_time, start_time, units = "mins")
cat("Training completed in:", round(training_time, 2), "minutes\n")

# Quick evaluation
pred <- predict(rf_model, X[-train_idx,])
probs <- predict(rf_model, X[-train_idx,], type = "prob")

print(confusionMatrix(pred, y[-train_idx]))
print(paste("AUC:", round(auc(roc(y[-train_idx], probs$Ciliary)), 3)))

# Top domains (fixed tibble warning)
top_domains <- varImp(rf_model)$importance %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Domain") %>%
  arrange(desc(Overall)) %>%
  slice_head(20)

top_domains %>%
  ggplot(aes(reorder(Domain, Overall), Overall)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Top Predictive Domains", x = "Domain", y = "Importance")