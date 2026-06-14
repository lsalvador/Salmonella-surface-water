# ============================================================
# Serovar co-occurrence heatmap (log10 odds ratios), measure association with log10 odds ratios
# - Presence/absence from serovar proportions ( > 0 )
# - Fisher's exact test p-values
# - BH/FDR correction (correct those p-values using Benjamini-Hochberg false discovery rate) 
# - Black outline where FDR < 0.05
# - Grey cells for upper triangle / diagonal
# - Code developed by Sarita Bugalia, May 2026
# ============================================================

library(readxl) # reads Excel files
library(dplyr) # data wrangling
library(ggplot2) # plotting 

rm(list = ls()) # removes all objects from memory
cat("\014")  # clears the console
graphics.off()  # closes all open plots

# ---------------------------
# 1) Load data 
# ---------------------------
df <- read_excel(
  "EGP_CC_Master_Data_File_3.2.26.xlsx",
  sheet = "Overview-Sal presence",
  skip  = 1
)

df$Date  <- as.Date(df$Date) # Converts the Date column into R Date format 
df$Month <- factor(format(df$Date, "%Y-%m")) # Creates a new Month column by formatting each date as "year-month" like "2026-03", then stores it as a factor

df <- df %>% rename_with(make.names, everything()) # Renames all columns so they become syntactically valid R names. For example, spaces or special characters may be replaced with dots

# keep only creek sites
df <- df %>% filter(Location %in% c("RC1", "RC2", "RC3"))

# for single site, use this by changing locations RC1, RC2, and RC3
# df <- df %>% filter(Location %in% c("RC3"))

# Identify Salmonella presence column
sal_col <- grep("^Salmonella.presence", names(df), value = TRUE)
if (length(sal_col) != 1) stop("Salmonella presence column not found.")

# Remove inconsistent rows
df <- df %>% filter(!( .data[[sal_col]] == 1 & Complexity == 0 ))

# Replace placeholders with NA
df <- df %>%
  mutate(across(where(is.character), ~ na_if(na_if(na_if(., "x"), "X"), "")))

# ---------------------------
# 2) Choose serovar columns
# ---------------------------

serovar_cols <- names(df)[43:74]

# Checks that every chosen serovar column actually exists in df
stopifnot(all(serovar_cols %in% names(df)))

# Ensure numeric (Converts all serovar columns to numeric) 
df[serovar_cols] <- lapply(df[serovar_cols], function(x) suppressWarnings(as.numeric(as.character(x))))

# Drop rows missing all serovar info
df <- df %>% filter(rowSums(is.na(across(all_of(serovar_cols)))) < length(serovar_cols))

# ---------------------------
# 3) Build presence/absence matrix
# ---------------------------

# Turns serovar measurements into binary presence/absence values 
# x > 0 becomes TRUE if present, FALSE if not 
# as.integer(...) turns that into 1 and 0; M becomes a data frame of 0/1 values

M <- as.data.frame(lapply(df[serovar_cols], function(x) as.integer(x > 0)))
M <- as.matrix(M) # Converts M from a data frame to a matrix 

# Prevalence (n samples where serovar present)
prev <- colSums(M, na.rm = TRUE)

# Keep only serovars that appear 5 times or more
sero_order <- names(prev[prev >= 5])

# Order serovars by prevalence (descending)
sero_order <- sero_order[order(prev[sero_order], decreasing = TRUE)]

# Labels with counts in parentheses
sero_labels <- paste0(sero_order, " (", prev[sero_order], ")")

# Restrict matrix to selected serovars only
M <- M[, sero_order, drop = FALSE]

# ---------------------------
# 4) Pairwise OR + Fisher test
# ---------------------------
# the following creates all possible pairs of serovars

pairs <- expand.grid(
  A = sero_order,
  B = sero_order,
  stringsAsFactors = FALSE
)

# Adds numeric positions for A and B in the ordered serovar list.
# These indices help determine where each pair sits in the matrix.

pairs$i <- match(pairs$A, sero_order)
pairs$j <- match(pairs$B, sero_order)

# We'll compute for LOWER triangle only (i > j); upper/diag will be grey
# Creates a logical column lower: TRUE for lower-triangle cells and FALSE for diagonal and upper triangle
# This avoids computing the same pair twice and keeps the plot clean 

pairs$lower <- pairs$i > pairs$j

# Defines a function to calculate statistics for one serovar pair

compute_pair <- function(a, b, M) {
  xa <- M[, a] # xa: presence/absence for serovar a and M: presence/absence matrix
  xb <- M[, b] # xb: presence/absence for serovar b
  
  n11 <- sum(xa == 1 & xb == 1, na.rm = TRUE) # Counts rows where both serovars are present
  n10 <- sum(xa == 1 & xb == 0, na.rm = TRUE) # Counts rows where a is present and b is absent
  n01 <- sum(xa == 0 & xb == 1, na.rm = TRUE) # Counts rows where a is absent and b is present
  n00 <- sum(xa == 0 & xb == 0, na.rm = TRUE) # Counts rows where both are absent
  
  # If any cell is zero, division can break. So the code adds 0.5 to every cell. 
  # That is the Haldane–Anscombe correction to avoid division by zero
  OR <- ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
  
  print(OR)
  
  # OR > 1: positive association / co-occurrence and OR < 1: negative association / avoidance
  
  # Takes the base-10 logarithm of the odds ratio. This makes the scale more symmetric:
  # 0 means no association, positive values mean positive association, and negative values mean negative association
  # Log-transform OR to create a symmetric visualization scale around zero
  log10_or <- log10(OR)
  
  # Fisher exact test: Runs Fisher’s exact test on the 2×2 contingency table.
  # matrix(c(n11, n10, n01, n00), nrow = 2): builds the table
  # alternative = "two.sided": tests for any association, positive or negative
  # $p.value: extracts the p-value
  
  p <- fisher.test(matrix(c(n11, n10, n01, n00), nrow = 2), alternative = "two.sided")$p.value
  
  # Returns three values: log10_or: log10 odds ratio; p: Fisher p-value; n11: count where both serovars are present
  list(log10_or = log10_or, p = p, n11 = n11)
}

# Creates empty columns in pairs to store results
pairs$log10_or <- NA_real_
pairs$p_value  <- NA_real_
pairs$n_both   <- NA_integer_

# Loop through all possible serovar pairs and calculate statistics only for lower-triangle pairs
for (k in seq_len(nrow(pairs))) {
  if (!pairs$lower[k]) next
  a <- pairs$A[k]
  b <- pairs$B[k]
  out <- compute_pair(a, b, M)
  pairs$log10_or[k] <- out$log10_or
  pairs$p_value[k]  <- out$p
  pairs$n_both[k]   <- out$n11
}

# FDR over all tested pairs (lower triangle)
# Identify valid tested pairs for multiple-comparison correction 
tested <- which(pairs$lower & !is.na(pairs$p_value))
pairs$fdr <- NA_real_

# Adjust p-values using Benjamini-Hochberg FDR correction
pairs$fdr[tested] <- p.adjust(pairs$p_value[tested], method = "BH")

pairs$log10_or_clamped <- pairs$log10_or

# Clamp plotted values to [-2, 2] so extreme odds ratios do not dominate the color scale
pairs$log10_or_clamped <- pmax(pmin(pairs$log10_or_clamped, 2), -2)

# Factor ordering + labels
# Create ordered y-axis labels that include serovar prevalence counts
pairs$A_lab <- factor(pairs$A, levels = sero_order, labels = sero_labels)

# Create ordered x-axis labels using serovar names only
pairs$B_lab <- factor(pairs$B, levels = sero_order, labels = sero_order)

# ---------------------------
# 5) Plot
# ---------------------------

# ---------------------------
# Add white padding for visualization only
# ---------------------------

# Set this to the largest number of serovars among panels A-D
max_n_serovars <- 15

n_current <- length(sero_order)
n_pad <- max_n_serovars - n_current

if (n_pad > 0) {
  pad_names <- paste0("PAD_", seq_len(n_pad))
} else {
  pad_names <- character(0)
}

# Plotting order includes real serovars + invisible padding levels
sero_order_plot <- c(sero_order, pad_names)

# Labels: real serovars get labels, padding gets blank labels
x_labels_plot <- setNames(
  c(sero_order, rep("", length(pad_names))),
  sero_order_plot
)

y_labels_plot <- setNames(
  c(sero_labels, rep("", length(pad_names))),
  sero_order_plot
)

# Use raw serovar names as plotting factors
Pairwise_cooccurrence <- pairs %>%
  mutate(
    A_plot = factor(A, levels = sero_order_plot),
    B_plot = factor(B, levels = sero_order_plot)
  )

# Keep significant lower-triangle associations for black outline overlay
sig_df <- Pairwise_cooccurrence %>%
  filter(lower, !is.na(fdr), fdr < 0.05)

# Build the heatmap
p <- ggplot(Pairwise_cooccurrence, aes(x = B_plot, y = A_plot, fill = log10_or_clamped)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(
    data = sig_df,
    aes(x = B_plot, y = A_plot),
    fill = NA,
    color = "black",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = x_labels_plot
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = y_labels_plot
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    na.value = "grey80",
    name = "Log10 odds"
  ) +
  guides(fill = "none") +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(x = NULL, y = NULL)

ggsave(
  "Pairwise_cooccurrence_All.tiff",
  plot = p,
  width = 5,
  height = 4.4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
