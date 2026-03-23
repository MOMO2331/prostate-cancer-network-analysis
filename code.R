# ==============================================================================
# Reproducible Network Analysis
# Title: Psychological Distress, Decision‑Making, and Health Behaviors
#        Among Prostate Cancer Survivors with Urinary Incontinence:
#        A Multilevel Network Analysis
#
# This script performs network analysis using the EBICglasso method and
# mgm for sensitivity. It includes data preprocessing, descriptive tables,
# redundancy checks, centrality, stability, and bootstrapping.
#
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS 26.2
#
# Key package versions:
#   bootnet      1.6
#   qgraph       1.9.8
#   mgm          1.2-15
#   psych        2.5.3
#   corrplot     0.95
#   readxl       1.4.5
#   dplyr        1.1.4
#   ggplot2      4.0.1
#   networktools 1.6.0
#   DiagrammeR   1.0.11
#
# All file paths are placeholders – please adjust them to your local setup.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load required packages and set up directories
# ------------------------------------------------------------------------------

# Function to install and load libraries if not already present
load_libraries <- function() {
  required_packages <- c(
    "bootnet", "dplyr", "qgraph", "readxl", "corrplot", "mgm",
    "psych", "ggplot2", "reshape2", "networktools"
  )
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

load_libraries()

# Set directories (change these to your own paths)
data_dir   <- "./data"          # folder containing raw data
output_dir <- "./results"       # folder for all output files

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load and preprocess data
# ------------------------------------------------------------------------------

# Read the cleaned dataset
df <- readxl::read_excel(file.path(data_dir, "cleaned_data.xlsx"))

# ------------------------------------------------------------------------------
# 2. PHQ‑9 analysis
# ------------------------------------------------------------------------------

# Select PHQ‑9 items
phq_vars <- c(
  "phq_interest_num",        # PHQ1  Anhedonia
  "phq_depressed_num",       # PHQ2  Sad mood
  "phq_sleep_num",           # PHQ3  Sleep problems
  "phq_energy_num",          # PHQ4  Low energy
  "phq_appetite_num",        # PHQ5  Appetite changes
  "phq_selfesteem_num",      # PHQ6  Guilt/worthlessness
  "phq_concentration_num",   # PHQ7  Concentration
  "phq_psychomotor_num",     # PHQ8  Psychomotor changes
  "phq_suicidal_num"         # PHQ9  Suicidal ideation
)

df_phq <- df[, phq_vars]
df_phq_z <- scale(df_phq) %>% as.data.frame()

# Set column names (short labels for plotting)
short_phq <- paste0("PHQ", 1:9)
colnames(df_phq_z) <- short_phq

# ------------------------------------------------------------------------------
# 2.1 Informative and redundancy checks (PHQ‑9)
# ------------------------------------------------------------------------------

# Standard deviations (informativeness)
info_sd_phq <- apply(df_phq_z, 2, sd, na.rm = TRUE)
print("PHQ-9: Standard deviations (informativeness):")
print(info_sd_phq)

# Correlation matrix (redundancy)
cor_phq <- cor(df_phq_z, use = "pairwise.complete.obs")
cat("\nPHQ-9: Correlation matrix (first 6 rows):\n")
print(round(cor_phq, 3)[1:6, 1:6])

# Redundant pairs (|r| > 0.25)
redundant_phq <- which(abs(cor_phq) > 0.25 & abs(cor_phq) < 1, arr.ind = TRUE)
cat("\nPHQ-9: Redundant pairs (|r| > 0.25):\n")
print(redundant_phq)

# Plot correlation matrix
pdf(file.path(output_dir, "PHQ9_correlation_matrix.pdf"), width = 8, height = 8)
corrplot::corrplot(cor_phq, method = "color", type = "upper",
                   tl.cex = 0.8, addCoef.col = "black", diag = FALSE,
                   title = "PHQ‑9 Symptom Correlations")
dev.off()

# ------------------------------------------------------------------------------
# 2.2 Descriptive statistics (Table S1)
# ------------------------------------------------------------------------------

desc_phq <- psych::describe(df_phq)
table_S1 <- data.frame(
  Variable = rownames(desc_phq),
  Mean     = round(desc_phq$mean, 2),
  SD       = round(desc_phq$sd, 2),
  Min      = desc_phq$min,
  Max      = desc_phq$max,
  Skewness = round(desc_phq$skew, 2),
  Kurtosis = round(desc_phq$kurtosis, 2)
)
write.csv(table_S1, file.path(output_dir, "Table_S1_PHQ9_distribution.csv"),
          row.names = FALSE)

# Histograms of each item
df_long_phq <- reshape2::melt(df_phq)
p <- ggplot2::ggplot(df_long_phq, aes(x = value)) +
  geom_histogram(bins = 20) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of PHQ‑9 Items")
ggsave(file.path(output_dir, "Figure_S2_PHQ9_histograms.png"), p,
       width = 10, height = 8)

# ------------------------------------------------------------------------------
# 2.3 EBICglasso network estimation
# ------------------------------------------------------------------------------

network_phq <- bootnet::estimateNetwork(
  df_phq_z,
  default = "EBICglasso",
  corMethod = "spearman"
)

adj_phq <- bootnet::getWmat(network_phq)
write.csv(adj_phq, file.path(output_dir, "PHQ9_EBICglasso_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 2.4 mgm sensitivity analysis (Gaussian model)
# ------------------------------------------------------------------------------

p_phq <- ncol(df_phq_z)
mgm_phq <- mgm::mgm(
  data = as.matrix(df_phq_z),
  type = rep("g", p_phq),
  level = rep(1, p_phq),
  lambdaSel = "CV",
  ruleReg = "OR",
  pbar = TRUE
)

# Predictability (R²) for mgm
pred_phq <- mgm::predict(mgm_phq, data = as.matrix(df_phq_z), errorCon = "R2")
pred_phq_df <- pred_phq$error
colnames(pred_phq_df) <- c("Node", "R2")
write.csv(pred_phq_df, file.path(output_dir, "PHQ9_mgm_predictability.csv"),
          row.names = FALSE)

# mgm adjacency matrix
adj_mgm_phq <- mgm_phq$pairwise$wadj
write.csv(adj_mgm_phq, file.path(output_dir, "PHQ9_mgm_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 2.5 Visualize networks (EBICglasso + mgm)
# ------------------------------------------------------------------------------

# Define node labels and colors
legend_phq <- c(
  "PHQ1: Anhedonia",
  "PHQ2: Sad mood",
  "PHQ3: Sleep problems",
  "PHQ4: Low energy",
  "PHQ5: Appetite changes",
  "PHQ6: Guilt/Worthlessness",
  "PHQ7: Poor concentration",
  "PHQ8: Psychomotor changes",
  "PHQ9: Suicidal ideation"
)
node_colors_phq <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD"
)

# Generate fixed layout using qgraph
layout_phq <- qgraph::qgraph(adj_phq, layout = "spring", DoNotPlot = TRUE)$layout

# Plot EBICglasso network with predictability pies
pdf(file.path(output_dir, "PHQ9_EBICglasso_network.pdf"), width = 14, height = 10)
par(mar = c(5, 4, 5, 2))
qgraph::qgraph(
  adj_phq,
  layout = layout_phq,
  labels = short_phq,
  color = node_colors_phq,
  edge.color = c("red", "blue"),
  negDashed = TRUE,
  vsize = 8,
  label.cex = 1.1,
  # Add predictability (R²) as pie charts
  pie = pred_phq_df$R2,
  legend = FALSE,
  title = "PHQ‑9 Network (EBICglasso) with Predictability (R²)",
  title.cex = 1.5
)
# Add legend in separate figure region
par(mar = c(5, 1, 5, 1))
plot.new()
plot.window(xlim = c(0, 20), ylim = c(0, 25))
text(1, 24, "PHQ‑9 Symptoms", adj = 0, cex = 1.5, font = 2)
text(1, 22, "Node colors indicate symptoms", adj = 0, cex = 0.9)
text(1, 21, "Pie charts show predictability (R²)", adj = 0, cex = 0.9)
for (i in 1:9) {
  y_pos <- 19 - (i-1) * 2
  points(1, y_pos, pch = 19, col = node_colors_phq[i], cex = 2.5)
  text(3, y_pos, legend_phq[i], adj = 0, cex = 0.9)
}
dev.off()

# Plot mgm network
pdf(file.path(output_dir, "PHQ9_mgm_network.pdf"), width = 10, height = 8)
qgraph::qgraph(
  adj_mgm_phq,
  layout = layout_phq,
  labels = short_phq,
  color = node_colors_phq,
  title = "PHQ‑9 Network (mgm Gaussian)",
  title.cex = 1.5
)
dev.off()

# ------------------------------------------------------------------------------
# 2.6 Centrality, stability, and bootstrap (EBICglasso)
# ------------------------------------------------------------------------------

# Bootstrap edge and centrality stability
set.seed(1234)
boot_phq <- bootnet::bootnet(
  network_phq,
  statistics = c("edge", "strength", "closeness", "betweenness"),
  nBoots = 1000,
  nCores = 4
)

# Summary table
summary_phq <- summary(boot_phq)
write.csv(summary_phq, file.path(output_dir, "PHQ9_bootstrap_summary.csv"),
          row.names = FALSE)

# Plot bootstrap results
pdf(file.path(output_dir, "PHQ9_bootstrap_plots.pdf"), width = 12, height = 10)
plot(boot_phq, "strength", order = "sample")
plot(boot_phq, "closeness", order = "sample")
plot(boot_phq, "betweenness", order = "sample")
plot(boot_phq, "edge", order = "sample")
dev.off()

# Case‑dropping bootstrap for stability
set.seed(1234)
boot_case_phq <- bootnet::bootnet(
  network_phq,
  nBoots = 2500,
  type = "case",
  nCores = 8,
  statistics = c("strength", "expectedInfluence", "betweenness", "closeness", "edge")
)

# Stability coefficients
stability_phq <- bootnet::corStability(boot_case_phq)
stability_df_phq <- data.frame(
  Statistic = c("Betweenness", "Closeness", "Edge", "ExpectedInfluence", "Strength"),
  CS_Coefficient = stability_phq
)
write.csv(stability_df_phq,
          file.path(output_dir, "PHQ9_stability_coefficients.csv"),
          row.names = FALSE)

# Centrality measures (z‑scored)
centrality_phq <- qgraph::centrality_auto(network_phq)
centrality_raw_phq <- as.data.frame(centrality_phq$node.centrality)
centrality_z_phq <- as.data.frame(scale(centrality_raw_phq))
write.csv(centrality_raw_phq,
          file.path(output_dir, "PHQ9_centrality_raw.csv"),
          row.names = TRUE)
write.csv(centrality_z_phq,
          file.path(output_dir, "PHQ9_centrality_z.csv"),
          row.names = TRUE)

# Centrality plot
pdf(file.path(output_dir, "PHQ9_centrality_plot.pdf"), width = 10, height = 8)
qgraph::centralityPlot(network_phq, orderBy = "Strength",
                       scale = "z-scores", include = "all")
dev.off()

# ------------------------------------------------------------------------------
# 3. GAD‑7 analysis
# ------------------------------------------------------------------------------

# Select GAD‑7 items
gad_vars <- c(
  "gad_nervous_num",           # GAD1  Nervous
  "gad_unable_control_num",    # GAD2  Uncontrollable worry
  "gad_excess_worry_num",      # GAD3  Excessive worry
  "gad_cannot_relax_num",      # GAD4  Cannot relax
  "gad_restless_num",          # GAD5  Restless
  "gad_irritable_num",         # GAD6  Irritable
  "gad_fear_something_num"     # GAD7  Afraid
)

df_gad <- df[, gad_vars]
df_gad_z <- scale(df_gad) %>% as.data.frame()

short_gad <- paste0("GAD", 1:7)
colnames(df_gad_z) <- short_gad

# ------------------------------------------------------------------------------
# 3.1 Redundancy checks (GAD‑7)
# ------------------------------------------------------------------------------

info_sd_gad <- apply(df_gad_z, 2, sd, na.rm = TRUE)
print("GAD-7: Standard deviations (informativeness):")
print(info_sd_gad)

cor_gad <- cor(df_gad_z, use = "pairwise.complete.obs")
redundant_gad <- which(abs(cor_gad) > 0.25 & abs(cor_gad) < 1, arr.ind = TRUE)
cat("\nGAD-7: Redundant pairs (|r| > 0.25):\n")
print(redundant_gad)

# Correlation plot
pdf(file.path(output_dir, "GAD7_correlation_matrix.pdf"), width = 7, height = 7)
corrplot::corrplot(cor_gad, method = "color", type = "upper",
                   tl.cex = 0.8, addCoef.col = "black", diag = FALSE,
                   title = "GAD‑7 Symptom Correlations")
dev.off()

# ------------------------------------------------------------------------------
# 3.2 EBICglasso network estimation
# ------------------------------------------------------------------------------

network_gad <- bootnet::estimateNetwork(
  df_gad_z,
  default = "EBICglasso",
  corMethod = "spearman"
)

adj_gad <- bootnet::getWmat(network_gad)
write.csv(adj_gad, file.path(output_dir, "GAD7_EBICglasso_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 3.3 mgm sensitivity analysis (Gaussian)
# ------------------------------------------------------------------------------

p_gad <- ncol(df_gad_z)
mgm_gad <- mgm::mgm(
  data = as.matrix(df_gad_z),
  type = rep("g", p_gad),
  level = rep(1, p_gad),
  lambdaSel = "CV",
  ruleReg = "OR",
  pbar = TRUE
)

pred_gad <- mgm::predict(mgm_gad, data = as.matrix(df_gad_z), errorCon = "R2")
pred_gad_df <- pred_gad$error
colnames(pred_gad_df) <- c("Node", "R2")
write.csv(pred_gad_df, file.path(output_dir, "GAD7_mgm_predictability.csv"),
          row.names = FALSE)

adj_mgm_gad <- mgm_gad$pairwise$wadj
write.csv(adj_mgm_gad, file.path(output_dir, "GAD7_mgm_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 3.4 Visualize networks (GAD‑7)
# ------------------------------------------------------------------------------

legend_gad <- c(
  "GAD1: Nervous",
  "GAD2: Uncontrollable worry",
  "GAD3: Excessive worry",
  "GAD4: Cannot relax",
  "GAD5: Restless",
  "GAD6: Irritable",
  "GAD7: Afraid"
)
node_colors_gad <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                     "#80B1D3", "#FDB462", "#B3DE69")

layout_gad <- qgraph::qgraph(adj_gad, layout = "spring", DoNotPlot = TRUE)$layout

# EBICglasso network with predictability
pdf(file.path(output_dir, "GAD7_EBICglasso_network.pdf"), width = 14, height = 10)
par(mar = c(5, 4, 5, 2))
qgraph::qgraph(
  adj_gad,
  layout = layout_gad,
  labels = short_gad,
  color = node_colors_gad,
  edge.color = c("red", "blue"),
  negDashed = TRUE,
  vsize = 8,
  label.cex = 1.1,
  pie = pred_gad_df$R2,
  legend = FALSE,
  title = "GAD‑7 Network (EBICglasso) with Predictability (R²)",
  title.cex = 1.5
)
par(mar = c(5, 1, 5, 1))
plot.new()
plot.window(xlim = c(0, 20), ylim = c(0, 20))
text(1, 19, "GAD‑7 Symptoms", adj = 0, cex = 1.5, font = 2)
text(1, 17.5, "Node colors indicate symptoms", adj = 0, cex = 0.9)
text(1, 16.5, "Pie charts show predictability (R²)", adj = 0, cex = 0.9)
for (i in 1:7) {
  y_pos <- 15 - (i-1) * 2
  points(1, y_pos, pch = 19, col = node_colors_gad[i], cex = 2.5)
  text(3, y_pos, legend_gad[i], adj = 0, cex = 0.9)
}
dev.off()

# mgm network
pdf(file.path(output_dir, "GAD7_mgm_network.pdf"), width = 10, height = 8)
qgraph::qgraph(
  adj_mgm_gad,
  layout = layout_gad,
  labels = short_gad,
  color = node_colors_gad,
  title = "GAD‑7 Network (mgm Gaussian)",
  title.cex = 1.5
)
dev.off()

# ------------------------------------------------------------------------------
# 3.5 Centrality, stability, bootstrap (GAD‑7)
# ------------------------------------------------------------------------------

set.seed(1234)
boot_gad <- bootnet::bootnet(
  network_gad,
  statistics = c("edge", "strength", "closeness", "betweenness"),
  nBoots = 1000,
  nCores = 4
)
summary_gad <- summary(boot_gad)
write.csv(summary_gad, file.path(output_dir, "GAD7_bootstrap_summary.csv"),
          row.names = FALSE)

pdf(file.path(output_dir, "GAD7_bootstrap_plots.pdf"), width = 12, height = 10)
plot(boot_gad, "strength", order = "sample")
plot(boot_gad, "closeness", order = "sample")
plot(boot_gad, "betweenness", order = "sample")
plot(boot_gad, "edge", order = "sample")
dev.off()

# Case‑dropping
set.seed(1234)
boot_case_gad <- bootnet::bootnet(
  network_gad,
  nBoots = 2500,
  type = "case",
  nCores = 8,
  statistics = c("strength", "expectedInfluence", "betweenness", "closeness", "edge")
)

stability_gad <- bootnet::corStability(boot_case_gad)
stability_df_gad <- data.frame(
  Statistic = c("Betweenness", "Closeness", "Edge", "ExpectedInfluence", "Strength"),
  CS_Coefficient = stability_gad
)
write.csv(stability_df_gad,
          file.path(output_dir, "GAD7_stability_coefficients.csv"),
          row.names = FALSE)

centrality_gad <- qgraph::centrality_auto(network_gad)
centrality_raw_gad <- as.data.frame(centrality_gad$node.centrality)
centrality_z_gad <- as.data.frame(scale(centrality_raw_gad))
write.csv(centrality_raw_gad,
          file.path(output_dir, "GAD7_centrality_raw.csv"),
          row.names = TRUE)
write.csv(centrality_z_gad,
          file.path(output_dir, "GAD7_centrality_z.csv"),
          row.names = TRUE)

pdf(file.path(output_dir, "GAD7_centrality_plot.pdf"), width = 10, height = 8)
qgraph::centralityPlot(network_gad, orderBy = "Strength",
                       scale = "z-scores", include = "all")
dev.off()

# ------------------------------------------------------------------------------
# 4. Combined variables analysis (SDM, DCS, Anxiety, Depression, Loneliness,
#    Urinary Incontinence, Exercise Adherence)
# ------------------------------------------------------------------------------

# Construct composite variables
df_combined <- df %>%
  dplyr::mutate(
    exercise_adherence = pf_exercise_freq_num * pf_exercise_time_num
  ) %>%
  dplyr::select(
    sdm_sum, dcs_total, gad7_total, phq9_total,
    uls_total, ICIQ_total, exercise_adherence
  )

df_combined_z <- scale(df_combined) %>% as.data.frame()
short_comb <- c("SDM", "DCS", "ANX", "DEP", "LON", "UI", "ADH")
colnames(df_combined_z) <- short_comb

# ------------------------------------------------------------------------------
# 4.1 Redundancy checks (combined)
# ------------------------------------------------------------------------------

cor_comb <- cor(df_combined_z, use = "pairwise.complete.obs")
redundant_comb <- which(abs(cor_comb) > 0.25 & abs(cor_comb) < 1, arr.ind = TRUE)
cat("\nCombined variables: Redundant pairs (|r| > 0.25):\n")
print(redundant_comb)

pdf(file.path(output_dir, "Combined_correlation_matrix.pdf"), width = 7, height = 7)
corrplot::corrplot(cor_comb, method = "color", type = "upper",
                   tl.cex = 0.8, addCoef.col = "black", diag = FALSE,
                   title = "Correlation Matrix (Combined Variables)")
dev.off()

# ------------------------------------------------------------------------------
# 4.2 EBICglasso network
# ------------------------------------------------------------------------------

network_comb <- bootnet::estimateNetwork(
  df_combined_z,
  default = "EBICglasso",
  corMethod = "spearman"
)

adj_comb <- bootnet::getWmat(network_comb)
write.csv(adj_comb, file.path(output_dir, "Combined_EBICglasso_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 4.3 mgm sensitivity (Gaussian)
# ------------------------------------------------------------------------------

p_comb <- ncol(df_combined_z)
mgm_comb <- mgm::mgm(
  data = as.matrix(df_combined_z),
  type = rep("g", p_comb),
  level = rep(1, p_comb),
  lambdaSel = "CV",
  ruleReg = "OR",
  pbar = TRUE
)

pred_comb <- mgm::predict(mgm_comb, data = as.matrix(df_combined_z), errorCon = "R2")
pred_comb_df <- pred_comb$error
colnames(pred_comb_df) <- c("Node", "R2")
write.csv(pred_comb_df, file.path(output_dir, "Combined_mgm_predictability.csv"),
          row.names = FALSE)

adj_mgm_comb <- mgm_comb$pairwise$wadj
write.csv(adj_mgm_comb, file.path(output_dir, "Combined_mgm_adjacency.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 4.4 Visualize combined network (EBICglasso with predictability)
# ------------------------------------------------------------------------------

legend_comb <- c(
  "Shared Decision Making",
  "Decision Conflict",
  "Anxiety (GAD-7)",
  "Depression (PHQ-9)",
  "Loneliness (ULS)",
  "Urinary Incontinence",
  "Exercise Adherence"
)
node_colors_comb <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                      "#80B1D3", "#FDB462", "#B3DE69")

layout_comb <- qgraph::qgraph(adj_comb, layout = "spring", DoNotPlot = TRUE)$layout

pdf(file.path(output_dir, "Combined_EBICglasso_network.pdf"), width = 14, height = 10)
par(mar = c(5, 4, 5, 2))
qgraph::qgraph(
  adj_comb,
  layout = layout_comb,
  labels = short_comb,
  color = node_colors_comb,
  edge.color = c("red", "blue"),
  negDashed = TRUE,
  vsize = 8,
  label.cex = 1.1,
  pie = pred_comb_df$R2,
  legend = FALSE,
  title = "Combined Variables Network (EBICglasso) with Predictability (R²)",
  title.cex = 1.5
)
par(mar = c(5, 1, 5, 1))
plot.new()
plot.window(xlim = c(0, 20), ylim = c(0, 20))
text(1, 19, "Variables", adj = 0, cex = 1.5, font = 2)
text(1, 17.5, "Node colors indicate domains", adj = 0, cex = 0.9)
text(1, 16.5, "Pie charts show predictability (R²)", adj = 0, cex = 0.9)
for (i in 1:7) {
  y_pos <- 15 - (i-1) * 2
  points(1, y_pos, pch = 19, col = node_colors_comb[i], cex = 2.5)
  text(3, y_pos, legend_comb[i], adj = 0, cex = 0.9)
}
dev.off()

# mgm network
pdf(file.path(output_dir, "Combined_mgm_network.pdf"), width = 10, height = 8)
qgraph::qgraph(
  adj_mgm_comb,
  layout = layout_comb,
  labels = short_comb,
  color = node_colors_comb,
  title = "Combined Variables Network (mgm Gaussian)",
  title.cex = 1.5
)
dev.off()

# ------------------------------------------------------------------------------
# 4.5 Centrality, stability, bootstrap (combined)
# ------------------------------------------------------------------------------

set.seed(1234)
boot_comb <- bootnet::bootnet(
  network_comb,
  statistics = c("edge", "strength", "closeness", "betweenness"),
  nBoots = 1000,
  nCores = 4
)
summary_comb <- summary(boot_comb)
write.csv(summary_comb, file.path(output_dir, "Combined_bootstrap_summary.csv"),
          row.names = FALSE)

pdf(file.path(output_dir, "Combined_bootstrap_plots.pdf"), width = 12, height = 10)
plot(boot_comb, "strength", order = "sample")
plot(boot_comb, "closeness", order = "sample")
plot(boot_comb, "betweenness", order = "sample")
plot(boot_comb, "edge", order = "sample")
dev.off()

set.seed(1234)
boot_case_comb <- bootnet::bootnet(
  network_comb,
  nBoots = 2500,
  type = "case",
  nCores = 8,
  statistics = c("strength", "expectedInfluence", "betweenness", "closeness", "edge")
)

stability_comb <- bootnet::corStability(boot_case_comb)
stability_df_comb <- data.frame(
  Statistic = c("Betweenness", "Closeness", "Edge", "ExpectedInfluence", "Strength"),
  CS_Coefficient = stability_comb
)
write.csv(stability_df_comb,
          file.path(output_dir, "Combined_stability_coefficients.csv"),
          row.names = FALSE)

centrality_comb <- qgraph::centrality_auto(network_comb)
centrality_raw_comb <- as.data.frame(centrality_comb$node.centrality)
centrality_z_comb <- as.data.frame(scale(centrality_raw_comb))
write.csv(centrality_raw_comb,
          file.path(output_dir, "Combined_centrality_raw.csv"),
          row.names = TRUE)
write.csv(centrality_z_comb,
          file.path(output_dir, "Combined_centrality_z.csv"),
          row.names = TRUE)

pdf(file.path(output_dir, "Combined_centrality_plot.pdf"), width = 10, height = 8)
qgraph::centralityPlot(network_comb, orderBy = "Strength",
                       scale = "z-scores", include = "all")
dev.off()

