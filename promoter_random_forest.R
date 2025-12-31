## =========================
## Promoter classification: k=4 k-mers + Random Forest (R)
## =========================

## 0) Packages (install once if missing)
# install.packages(c("data.table","caret","randomForest","pROC","ggplot2","reshape2"))

library(data.table)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)
library(reshape2)

set.seed(123)

## 1) Load data
file_path <- "C:/Users/enock_p22oyv9/OneDrive/Desktop/dataset.csv"
dt <- fread(file_path)

## 2) Validate + clean
stopifnot(all(c("id", "sequence", "label") %in% names(dt)))

dt[, sequence := toupper(sequence)]
dt <- dt[grepl("^[ACGT]+$", sequence)]   # keep only valid DNA bases

## Labels (clear class names)
y <- factor(dt$label, levels = c(0, 1), labels = c("Non-Promoter", "Promoter"))

## 3) k-mer encoding (k=4)
make_kmer_matrix <- function(seqs, k = 4, normalize = TRUE) {
  alphabet <- c("A","C","G","T")
  kmers <- alphabet
  for (i in 2:k) kmers <- as.vector(outer(kmers, alphabet, paste0))
  kmers <- sort(kmers)
  kmer_index <- setNames(seq_along(kmers), kmers)
  
  n <- length(seqs); p <- length(kmers)
  X <- matrix(0, nrow = n, ncol = p)
  colnames(X) <- kmers
  
  for (i in seq_len(n)) {
    s <- seqs[i]; L <- nchar(s)
    parts <- substring(s, 1:(L-k+1), 1:(L-k+1) + k - 1)
    tab <- table(parts)
    X[i, kmer_index[names(tab)]] <- as.integer(tab)
    if (normalize) X[i, ] <- X[i, ] / sum(X[i, ])
  }
  X
}

X <- make_kmer_matrix(dt$sequence, k = 4, normalize = TRUE)

## 4) Train/test split (stratified)
idx <- createDataPartition(y, p = 0.80, list = FALSE)
X_train <- X[idx, , drop = FALSE]
X_test  <- X[-idx, , drop = FALSE]
y_train <- y[idx]
y_test  <- y[-idx]

## 5) Train Random Forest
set.seed(123)
rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  mtry = floor(sqrt(ncol(X_train))),
  importance = TRUE
)

## 6) Predict + evaluate (Confusion Matrix + AUC)
rf_pred_class <- predict(rf_model, X_test, type = "class")
rf_pred_prob  <- predict(rf_model, X_test, type = "prob")[, "Promoter"]

cm <- confusionMatrix(rf_pred_class, y_test, positive = "Promoter")
print(cm)

rf_roc <- roc(response = y_test, predictor = rf_pred_prob, levels = c("Non-Promoter","Promoter"))
cat("AUC:", as.numeric(auc(rf_roc)), "\n")

## 7) Top 20 important k-mers (MeanDecreaseGini)
imp <- importance(rf_model)
imp_vec <- imp[, "MeanDecreaseGini"]
top20 <- sort(imp_vec, decreasing = TRUE)[1:20]
print(top20)

## 8) Plot 1: Confusion Matrix heatmap
cm_table <- table(Actual = y_test, Predicted = rf_pred_class)
cm_df <- melt(cm_table)

ggplot(cm_df, aes(x = Predicted, y = Actual, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), size = 6) +
  scale_fill_gradient(low = "lightblue", high = "deepskyblue3") +
  labs(title = "Confusion Matrix", x = "Predicted", y = "Actual") +
  theme_minimal(base_size = 14)

## 9) Plot 2: Top 20 k-mers bar plot
top20_df <- data.frame(
  kmer = factor(names(top20), levels = names(top20)),
  importance = as.numeric(top20)
)

ggplot(top20_df, aes(x = kmer, y = importance)) +
  geom_col() +
  labs(title = "Top 20 Most Important K-mers", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
