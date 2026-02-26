#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "data/gasch2000.txt"
outfile <- if (length(args) >= 2) args[2] else "results/gene_top10.png"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

if (!dir.exists("results")) dir.create("results")

message("Reading: ", infile)
df <- read_tsv(infile, comment = "#", col_names = TRUE, show_col_types = FALSE)

# First column = gene identifiers
genes <- as.character(df[[1]])
expr <- df[-1]

# Remove condition columns named 'name' or 'gweight' if present
cond_names <- colnames(expr)
drop_idx <- which(tolower(cond_names) %in% c("name", "gweight"))
if (length(drop_idx) > 0) {
  message("Dropping condition columns: ", paste(cond_names[drop_idx], collapse = ", "))
  expr <- expr[, -drop_idx, drop = FALSE]
}

# Coerce expression columns to numeric if needed
expr <- expr %>% mutate(across(everything(), ~ as.numeric(.)))

# Simple checks
if (!all(sapply(expr, is.numeric))) stop("Some expression columns could not be coerced to numeric")

vals <- unlist(expr)
vals <- vals[!is.na(vals)]

# Heuristic to decide whether data is already log-scale
is_log <- NA
if (any(vals < 0)) {
  is_log <- TRUE
  message("Negative values detected -> assuming data already on a log scale.")
} else if (max(vals, na.rm = TRUE) > 50) {
  is_log <- FALSE
  message("Large positive values detected -> will apply log2(x + 1).")
} else {
  is_log <- TRUE
  message("Values in moderate range -> assuming already log-scale.")
}

if (!is_log) {
  expr <- log2(expr + 1)
  message("Applied log2(x + 1) transform to expression values.")
}

# Compute CV per gene (sd / mean), guard mean == 0
means <- rowMeans(expr, na.rm = TRUE)
sds <- apply(expr, 1, sd, na.rm = TRUE)
cv <- sds / means
cv[which(means == 0)] <- NA

summary_tbl <- tibble(gene = genes, mean = means, sd = sds, cv = cv)

top10 <- summary_tbl %>% arrange(desc(cv)) %>% filter(!is.na(cv)) %>% slice_head(n = 10)
if (nrow(top10) == 0) stop("No genes with valid CV found")
message("Top genes selected (by CV): ", paste(top10$gene, collapse = ", "))

# Subset expression to top genes and reshape
expr_top <- expr %>% mutate(gene = genes) %>% filter(gene %in% top10$gene)

# Order genes on x-axis by CV descending
gene_levels <- top10$gene
expr_top$gene <- factor(expr_top$gene, levels = gene_levels)

long <- expr_top %>% pivot_longer(-gene, names_to = "condition", values_to = "expression")

# Plot heatmap
p <- ggplot(long, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "PuGn", na.value = "grey50") +
  labs(title = "gene top 10", x = "gene", y = "conditions", fill = "expression") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Times New Roman", size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.92, 0.92),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

ggsave(outfile, plot = p, width = 5, height = 5, dpi = 300, bg = "white")
message("Saved heatmap to: ", outfile)
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "data/gasch2000.txt"
outfile <- if (length(args) >= 2) args[2] else "results/gene_top10.png"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

if (!dir.exists("results")) dir.create("results")

message("Reading: ", infile)
df <- read_tsv(infile, comment = "#", col_names = TRUE, show_col_types = FALSE)

# First column = gene identifiers
genes <- as.character(df[[1]])
expr <- df[-1]

# Remove condition columns named 'name' or 'gweight' if present
cond_names <- colnames(expr)
drop_idx <- which(tolower(cond_names) %in% c("name", "gweight"))
if (length(drop_idx) > 0) {
  message("Dropping condition columns: ", paste(cond_names[drop_idx], collapse = ", "))
  expr <- expr[, -drop_idx, drop = FALSE]
}

# Coerce expression columns to numeric if needed
expr <- expr %>% mutate(across(everything(), ~ as.numeric(.)))

# Simple checks
if (!all(sapply(expr, is.numeric))) stop("Some expression columns could not be coerced to numeric")

vals <- unlist(expr)
vals <- vals[!is.na(vals)]

# Heuristic to decide whether data is already log-scale
is_log <- NA
if (any(vals < 0)) {
  is_log <- TRUE
  message("Negative values detected -> assuming data already on a log scale.")
} else if (max(vals, na.rm = TRUE) > 50) {
  is_log <- FALSE
  message("Large positive values detected -> will apply log2(x + 1).")
} else {
  is_log <- TRUE
  message("Values in moderate range -> assuming already log-scale.")
}

if (!is_log) {
  expr <- log2(expr + 1)
  message("Applied log2(x + 1) transform to expression values.")
}

# Compute CV per gene (sd / mean), guard mean == 0
means <- rowMeans(expr, na.rm = TRUE)
sds <- apply(expr, 1, sd, na.rm = TRUE)
cv <- sds / means
cv[which(means == 0)] <- NA

summary_tbl <- tibble(gene = genes, mean = means, sd = sds, cv = cv)

top10 <- summary_tbl %>% arrange(desc(cv)) %>% filter(!is.na(cv)) %>% slice_head(n = 10)
if (nrow(top10) == 0) stop("No genes with valid CV found")
message("Top genes selected (by CV): ", paste(top10$gene, collapse = ", "))

# Subset expression to top genes and reshape
expr_top <- expr %>% mutate(gene = genes) %>% filter(gene %in% top10$gene)

# Order genes on x-axis by CV descending
gene_levels <- top10$gene
expr_top$gene <- factor(expr_top$gene, levels = gene_levels)

long <- expr_top %>% pivot_longer(-gene, names_to = "condition", values_to = "expression")

# Plot heatmap
p <- ggplot(long, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "PuGn", na.value = "grey50") +
  labs(title = "gene top 10", x = "gene", y = "conditions", fill = "expression") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Times New Roman", size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.92, 0.92),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

ggsave(outfile, plot = p, width = 5, height = 5, dpi = 300, bg = "white")
message("Saved heatmap to: ", outfile)
