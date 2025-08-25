library(ape)
library(phylolm)
library(dplyr)

setwd("/Users/juliechen/Julie_Chen/warnecke lab/mcmc_motif_enrichment/")

# Load tree and data
tree <- read.tree("gtdb/bacteria_tree.nwk")
filename = "motif_enrichment_summary_bacteria_k6.csv"
df <- read.csv(filename)

# Choose covariate column (example: "nap")
covariate <- "nap"

# Clean species names
df$species <- gsub("_", "", df$species)
df$species <- tolower(df$species)
tree$tip.label <- gsub("_", "", tree$tip.label)
tree$tip.label <- tolower(tree$tip.label)

# Match tree and dataframe
subset_df <- df %>% filter(species %in% tree$tip.label)
tree <- drop.tip(tree, setdiff(tree$tip.label, subset_df$species))
rownames(subset_df) <- subset_df$species

# Grab motif columns: everything from tRNA_scan_motif onwards
start_idx <- match("tRNA_scan_motif", colnames(subset_df))
motif_cols <- colnames(subset_df)[start_idx:length(colnames(subset_df))]

# Drop any motif columns that end with "kmer"
motif_cols <- motif_cols[!grepl("kmer$", motif_cols)]

# Function to safely run phylolm and extract slope + p-value
run_phylolm <- function(response, model_type) {
  formula <- as.formula(paste(response, "~", covariate))
  fit <- tryCatch(
    phylolm(formula, data = subset_df, phy = tree, model = model_type),
    error = function(e) return(NULL)
  )
  if (is.null(fit)) return(NULL)
  
  coef_summary <- summary(fit)$coefficients
  coef <- coef_summary[2, "Estimate"]
  pval <- coef_summary[2, "p.value"]
  r2 <- summary(fit)$r.squared
  adj_r2 <- summary(fit)$adj.r.squared
  aic <- AIC(fit)
  loglik <- logLik(fit)
  
  return(data.frame(
    motif = response,
    model = model_type,
    coef = coef,
    pval = pval,
    R2 = r2,
    Adj_R2 = adj_r2,
    AIC = aic,
    logLik = as.numeric(loglik),
    row.names = NULL
  ))
}

# Loop over motifs and models
results_list <- list()
for (motif in motif_cols) {
  for (model_type in c("BM", "OUrandomRoot")) {
    res <- run_phylolm(motif, model_type)
    if (!is.null(res)) results_list <- append(results_list, list(res))
  }
}

# Combine into one dataframe
results_df <- bind_rows(results_list)
write.csv(results_df, file = sub("\\.csv$", "_phylolm_results.csv", filename), quote=FALSE, row.names = FALSE)

library(dplyr)
library(ComplexHeatmap)

# Start with results_df from previous step
# Split results into separate dataframes
results_BM <- results_df %>% filter(model == "BM")
results_OU <- results_df %>% filter(model == "OUrandomRoot")

# Drop unwanted motifs
drop_motifs <- c("ARCHAEAL_FIXED_min", "ARCHAEAL_REGEX_min")
results_BM <- results_BM %>% filter(!motif %in% drop_motifs)
results_OU <- results_OU %>% filter(!motif %in% drop_motifs)

# Rename motifs
rename_map <- c(
  "ARCHAEAL_TATA_avg_top3" = "Most depleted TATA-like (top 3)",
  "tRNA_scan_motif"        = "tRNA promoter learned motif",
  "ARCHAEAL_TATA_min"      = "Most depleted TATA-like",
  "overall_min"            = "Most depleted overall"
)

results_BM$motif <- recode(results_BM$motif, !!!rename_map)
results_OU$motif <- recode(results_OU$motif, !!!rename_map)

# Compute signed log p-value
results_BM <- results_BM %>%
  mutate(signed_logp = -log10(pval) * sign(coef))

results_OU <- results_OU %>%
  mutate(signed_logp = -log10(pval) * sign(coef))

# --- Heatmap preparation ---
# Pivot to wide format (motifs x model)
# --- BM results ---
BM_mat <- results_BM %>%
  select(motif, signed_logp) %>%
  distinct() %>%
  column_to_rownames("motif") %>%
  as.matrix()

BM_pval_mat <- results_BM %>%
  select(motif, pval) %>%
  distinct() %>%
  column_to_rownames("motif") %>%
  as.matrix()

# --- OU results ---
OU_mat <- results_OU %>%
  select(motif, signed_logp) %>%
  distinct() %>%
  column_to_rownames("motif") %>%
  as.matrix()

OU_pval_mat <- results_OU %>%
  select(motif, pval) %>%
  distinct() %>%
  column_to_rownames("motif") %>%
  as.matrix()

# --- Combine BM and OU side by side ---
combined_mat <- cbind(
  BM = BM_mat[rownames(BM_mat), 1], 
  OUrandomRoot = OU_mat[rownames(BM_mat), 1]
)

pval_mat <- cbind(
  BM = BM_pval_mat[rownames(BM_pval_mat), 1],
  OUrandomRoot = OU_pval_mat[rownames(BM_pval_mat), 1]
)

rownames(combined_mat) <- rownames(BM_mat)
rownames(pval_mat) <- rownames(BM_pval_mat)

# get all motif names in your matrix
all_motifs <- rownames(combined_mat)

# put your 4 renamed motifs first, in desired order
priority_motifs <- c(
  "tRNA promoter learned motif",
  #"Most depleted TATA-like",
  #"Most depleted TATA-like (top 3)",
  "Most depleted overall"
)

# get the rest of the motifs, excluding those 4
other_motifs <- setdiff(all_motifs, priority_motifs)

# combine into one vector for row_order
row_order <- c(priority_motifs, sort(other_motifs))


# --- Heatmap ---
Heatmap(combined_mat,
        name = "signed -log10(p)",
        col = colorRamp2(c(-0.4, 0), c("blue", "white")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = row_order,
        rect_gp = gpar(col = "white", lwd = 1),  # white borders, line width
        cell_fun = function(j, i, x, y, w, h, fill) {
          p <- pval_mat[i, j]
          if (!is.na(p)) {
            if (p < 0.01) {
              grid.text("**", x, y, gp = gpar(fontsize = 14))
            } else if (p < 0.05) {
              grid.text("*", x, y, gp = gpar(fontsize = 14))
            } else if (p < 0.1) {
              grid.text(".", x, y, gp = gpar(fontsize = 14))
            }
          }
        })
