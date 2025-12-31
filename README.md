# promoter-prediction-using-kmers
Machine learning–based promoter prediction using k-mer encoding and Random Forest
# Promoter Prediction Using Machine Learning

This project investigates whether machine-learning models can identify promoter regions directly from raw DNA sequence data.

Using a labeled *E. coli* promoter dataset (1,000 sequences), DNA sequences were transformed into numerical features via **4-mer (k-mer) frequency encoding**, and a **Random Forest classifier** was trained to distinguish promoter from non-promoter regions.

---

## Dataset
- Organism: *Escherichia coli*
- Samples: 1,000 DNA sequences
- Labels:
  - `1` → Promoter
  - `0` → Non-Promoter
- Input: Fixed-length DNA sequences (A, C, G, T)

---

## Methods
- **Sequence encoding:** Overlapping 4-mer frequency vectors (256 features)
- **Model:** Random Forest
- **Evaluation:** Stratified train/test split (80/20), confusion matrix, ROC–AUC
- **Interpretation:** Feature importance (Mean Decrease Gini)

---

## Results
- Accuracy: ~96%
- ROC–AUC: ~0.97
- Sensitivity (Promoter): 100%
- Specificity (Non-Promoter): ~92%

Top-ranking k-mers were **AT-rich motifs** such as `AAAA`, `TATA`, `TAAA`, and `ATAA`, consistent with known bacterial promoter elements (e.g., the −10 Pribnow box).

---

## Repository Contents
- `scripts/` – R code for data processing, modeling, and visualization
- `figures/` – Confusion matrix and k-mer importance plots
- `data/` – Promoter dataset (if permitted)

---

## Tools
- R
- randomForest
- caret
- pROC
- ggplot2

---

## Author
Enock Kumi Ackaah  
Bioinformatics & Statistical Genomics
