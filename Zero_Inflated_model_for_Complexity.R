# ============================================================
# ZERO-INFLATED count models for Salmonella Complexity (count outcome)
# Fixed: Location + Wet.Dry + environmental predictors (X)
#
# Fits and compares:
#   1) Zero-inflated Poisson  (ZIP)
#   2) Zero-inflated Negtive Binomial (ZINB; nbinom2)
#   3) Zero-inflated COM-Poisson (ZICOMP)
#
#   - Non-linear functions for Flow, Width, Rain Day Before and Day of
#   - Code developed by Sarita Bugalia, May 2026
# ============================================================

# -------------------------------
# 0) Load Libraries
# -------------------------------
library(readxl) # for importing Excel files
library(dplyr) # for data cleaning and manipulation
library(forcats) # for categorical/factor variable handling
library(caret) # for correlation-based predictor filtering 
library(glmmTMB) # for fitting zero-inflated count models
library(ggeffects) # for predicted/marginal effects 
library(ggplot2) # for plots
library(tidyr) # for reshaping summary tables
library(performance) # for model-fit metrics
library(DHARMa) # for simulation-based residual diagnostics 
library(corrplot)

rm(list = ls()) # removes all objects from memory
cat("\014")  # clears the console
graphics.off()  # closes all open plots 

# -------------------------------
# 0b) Helper functions
# -------------------------------

# Defines a helper function to replace R-safe variable names with original Excel names

relabel_with_original <- function(x, name_key) {
  x_new <- x
  for (nm in names(name_key)) {
    x_new <- gsub(
      pattern = paste0("\\b", nm, "\\b"),
      replacement = name_key[[nm]],
      x = x_new
    )
  }
  x_new
}

# Defines a helper function for readable model-output labels
pretty_term_labels <- function(x, name_key) {
  x <- relabel_with_original(x, name_key)
  
  # ---- Flow
  x <- gsub("\\bFlow_c2\\b", "Centered Flow²", x)
  x <- gsub("\\bFlow_c\\b", "Centered Flow", x)
  
  # ---- Width
  x <- gsub("\\bWidth_log\\b", "Log-transformed width", x)
  
  # ---- Rain (two-part: occurrence + amount)
  x <- gsub("\\bRainAny_Before\\b", "Rainfall occurrence, previous day", x)
  x <- gsub("\\bRainAny_DayOf\\b", "Rainfall occurrence, sampling day", x)
  
  # ---- Rain spline (previous day)
  x <- gsub("\\bRainSpline_Before1\\b", "Rainfall amount (previous day), spline 1", x)
  x <- gsub("\\bRainSpline_Before2\\b", "Rainfall amount (previous day), spline 2", x)
  x <- gsub("\\bRainSpline_Before3\\b", "Rainfall amount (previous day), spline 3", x)
  
  # ---- Rain spline (sampling day)
  x <- gsub("\\bRainSpline_DayOf1\\b", "Rainfall amount (sampling day), spline 1", x)
  x <- gsub("\\bRainSpline_DayOf2\\b", "Rainfall amount (sampling day), spline 2", x)
  x <- gsub("\\bRainSpline_DayOf3\\b", "Rainfall amount (sampling day), spline 3", x)
  
  # # ---- For log-only version (safe to keep)
  # x <- gsub("\\bRainLog_Before\\b", "Rainfall amount (previous day, log-transformed)", x)
  # x <- gsub("\\bRainLog_DayOf\\b", "Rainfall amount (sampling day, log-transformed)", x)
  
  # Assigns the result of this operation to an object for later use 
  x <- gsub("`", "", x)
  
  x
}

# Defines a helper function for readable plot axis labels for GAM smooth plots
pretty_var_name <- function(v, name_key) {
  out <- if (v %in% names(name_key)) name_key[[v]] else v
  out <- pretty_term_labels(out, name_key)
  out
}

# -------------------------------
# 1) Load and basic data cleaning
# -------------------------------
df <- read_excel(
  "EGP_CC_Master_Data_File_3.2.26.xlsx",
  sheet = "Overview-Sal presence",
  skip  = 1
)

# Stores original Excel column names before renaming
original_names <- names(df)

# Creates valid R column names from original names
safe_names <- make.names(original_names, unique = TRUE)

# Assigns R-safe names to the data frame
names(df) <- safe_names

# Creates a lookup table from safe names to original names 
name_key <- setNames(original_names, safe_names)

# Keeps only RC1, RC2, and RC3 samples
df <- df %>% filter(Location %in% c("RC1", "RC2", "RC3"))

# Identify Salmonella presence column
sal_col <- grep("^Salmonella.presence", names(df), value = TRUE)
if (length(sal_col) != 1) stop("Salmonella presence column not found.")

# Remove biologically inconsistent rows
df <- df %>% filter(!(.data[[sal_col]] == 1 & Complexity == 0))

# Sets the extreme Width value of 410 inches at RC3 to NA while retaining the row
df$Width..in.[df$Location == "RC3" & df$Width..in. == 410] <- NA

# Replace placeholders with NA
df <- df %>%
  mutate(across(where(is.character),
                ~ na_if(na_if(na_if(., "x"), "X"), "")))

# Check raw counts by Wet/Dry
cat("\n=== Raw counts by Wet/Dry and Salmonella presence ===\n")
print(table(df$Wet.Dry, df[[sal_col]]))

# Defines categorical/date columns that should not be converted to numeric
non_numeric_cols <- c("Location", "Wet.Dry", "Date")

# Keeps only listed non-numeric columns that exist in the dataset
non_numeric_cols <- intersect(non_numeric_cols, names(df))

# Identifies numeric-like columns for numeric conversion
num_cols_all <- setdiff(names(df), non_numeric_cols)

# Applies numeric conversion to numeric-like columns
df[num_cols_all] <- lapply(df[num_cols_all], function(x)
  suppressWarnings(as.numeric(as.character(x)))
)

# -------------------------------
# 2) Define outcome
# -------------------------------
df$Complexity <- suppressWarnings(as.numeric(df$Complexity))
df_clean <- df # %>% filter(!is.na(Complexity))

# Counts the number of rows in the cleaned data
cat("Rows remaining:", nrow(df_clean), "\n")
if (nrow(df_clean) < 30) stop("Too few rows.")

# Sets Wet/Dry as a factor and specifies Dry as the reference category 
df_clean$Wet.Dry <- factor(df_clean$Wet.Dry, levels = c("Dry", "Wet"))

# Sets Location as a factor and specifies RC1 as the reference category
df_clean$Location <- factor(df_clean$Location, levels = c("RC1", "RC2", "RC3"))

# # -------------------------------
# # 2b) Table: Sample description and summary statistics (optional)
# # -------------------------------
# 
# # Prints a labeled message to the console
# cat("\n=== Sample description ===\n")
# 
# # Builds a table of sample counts and Complexity summary statistics 
# sample_description <- data.frame(
#   Metric = c(
#     "Total observations after cleaning",
#     "Number of Dry observations",
#     "Number of Wet observations",
#     "Number of Salmonella absence observations",
#     "Number of Salmonella presence observations",
#     "Mean Complexity",
#     "SD Complexity",
#     "Median Complexity",
#     "Min Complexity",
#     "Max Complexity"
#   ),
#   Value = c(
#     nrow(df_clean),
#     sum(df_clean$Wet.Dry == "Dry", na.rm = TRUE),
#     sum(df_clean$Wet.Dry == "Wet", na.rm = TRUE),
#     sum(df_clean[[sal_col]] == 0, na.rm = TRUE),
#     sum(df_clean[[sal_col]] == 1, na.rm = TRUE),
#     round(mean(df_clean$Complexity, na.rm = TRUE), 3),
#     round(sd(df_clean$Complexity, na.rm = TRUE), 3),
#     round(median(df_clean$Complexity, na.rm = TRUE), 3),
#     round(min(df_clean$Complexity, na.rm = TRUE), 3),
#     round(max(df_clean$Complexity, na.rm = TRUE), 3)
#   )
# )
# 
# # Prints the sample-description table in tab-delimited form 
# write.table(sample_description, sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # Prints a labeled message to the console 
# cat("\n=== Summary statistics for numeric variables in cleaned data ===\n")
# 
# # Constructs summary statistics for numeric variables
# numeric_summary <- df_clean %>%
#   select(where(is.numeric)) %>%
#   summarise(across(
#     everything(),
#     list(
#       N = ~sum(!is.na(.)),
#       Mean = ~mean(., na.rm = TRUE),
#       SD = ~sd(., na.rm = TRUE),
#       Median = ~median(., na.rm = TRUE),
#       Min = ~min(., na.rm = TRUE),
#       Max = ~max(., na.rm = TRUE)
#     ),
#     .names = "{.col}__{.fn}"
#   )) %>%
#   # Reshapes the summary table from wide to long format
#   pivot_longer(
#     cols = everything(),
#     names_to = c("Variable", "Statistic"),
#     names_sep = "__",
#     values_to = "Value"
#   ) %>%
#   # Reshapes the summary table so each statistic becomes a column  
#   pivot_wider(names_from = Statistic, values_from = Value) %>%
#   mutate(
#     # Applies readable labels to variable names 
#     Variable = pretty_term_labels(Variable, name_key),
#     across(where(is.numeric), ~ round(., 3))
#   )
# 
# # Prints the numeric summary table in tab-delimited form 
# write.table(numeric_summary, sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------
# 3) Build predictor matrix X
# -------------------------------

# Defines candidate predictor columns
predictor_cols <- c("Location", "Wet.Dry", names(df_clean)[7:42])

# # Removes requested predictor names that are absent from the dataset
# predictor_cols <- predictor_cols[predictor_cols %in% names(df_clean)]

# Creates the candidate predictor matrix
X <- df_clean[, predictor_cols, drop = FALSE]

# Assigns the result of this operation to an object for later use
X <- X %>%
  # Begins replacing placeholder character values with NA
  mutate(across(where(is.character), as.factor)) %>%
  # Keeps missing factor values as an explicit Missing category
  mutate(across(where(is.factor), ~ fct_explicit_na(., "Missing")))

# Remove zero-variance predictors (keep Location and Wet.Dry)
keep_cols <- sapply(names(X), function(nm) {
  if (nm %in% c("Location", "Wet.Dry")) return(TRUE)
  v <- X[[nm]]
  length(unique(v[!is.na(v)])) > 1
})
X <- X[, keep_cols, drop = FALSE]

# Prints a labeled message to the console
cat("Predictors after zero-variance filter:", ncol(X), "\n")

# -------------------------------
# 4) Remove correlated numeric predictors (|r| > 0.75)
# -------------------------------

# Extracts numeric predictors for correlation screening 
num_predictors <- X[, sapply(X, is.numeric), drop = FALSE]

if (ncol(num_predictors) > 1) {
  # Computes pairwise correlations among numeric predictors 
  corr_matrix <- cor(num_predictors, use = "pairwise.complete.obs")
  
  # Replaces NA correlations with zero
  corr_matrix[is.na(corr_matrix)] <- 0
  
  # Identifies highly correlated predictors to remove
  high_corr_idx <- findCorrelation(corr_matrix, cutoff = 0.75)
  
  if (length(high_corr_idx) > 0) {
    # Stores names of predictors removed for high correlation 
    removed <- colnames(num_predictors)[high_corr_idx]
    
    # Prints a labeled message to the console 
    cat("Removing correlated predictors:",
        paste(pretty_term_labels(removed, name_key), collapse = ", "),
        "\n")
    
    # Removes highly correlated predictors from X
    X <- X[, !(names(X) %in% removed), drop = FALSE]
    
    # Updates the numeric predictor matrix after removing correlated predictors
    num_predictors <- num_predictors[, -high_corr_idx, drop = FALSE]
  }
}

# Final correlation matrix after removing highly correlated predictors
rem_corr <- cor(num_predictors, use = "pairwise.complete.obs")

# Replaces any remaining NA correlations with zero for plotting
rem_corr[is.na(rem_corr)] <- 0

# Restore original Excel names for plotting
# Replaces safe column names in the correlation plot with original Excel labels
colnames(rem_corr) <- name_key[colnames(rem_corr)]

# Replaces safe row names in the correlation plot with original Excel labels
rownames(rem_corr) <- name_key[rownames(rem_corr)]

# plot 
tiff("correlation_heatmap_plot.tiff",
     width = 6.5, height = 6.5, units = "in",
     res = 600)

# Sets base-R plot margins 
par(mar = c(0.5, 0.5, 0.5, 0.5), oma = c(0, 0, 0, 0))

# Creates the final predictor correlation plot 
corrplot(
  rem_corr,
  method = "ellipse",
  type = "lower",
  diag = TRUE,
  tl.pos = "ld",
  tl.col = "black",
  tl.cex = 0.8,
  cl.pos = "b",
  addCoef.col = "black",
  number.cex = 0.5,
  col = colorRampPalette(c("blue", "white", "red"))(200),
  bg = "white",
  mar = c(0, 0, 0, 0)
)

# Closes the graphics device and writes the image file to disk
dev.off()

# -------------------------------
# 5) Optional visualization (GAM smooths)
# -------------------------------

# Creates a plotting dataset for exploratory GAM plots
plot_df <- data.frame(Complexity = df_clean$Complexity, X)

# Lists selected predictors for exploratory plots
vars_to_plot <- c(
  "Flow",
  "Width..in.",
  "Rain..in..Day.Before",
  "Rain..in..Day.of"
)

# Assigns the result to an object for later use 
vars_to_plot <- vars_to_plot[vars_to_plot %in% names(plot_df)]

# Loops over selected variables and creates one plot per variable 
for (v in vars_to_plot) {
  x <- plot_df[[v]]
  # Checks whether a variable has more than one observed value 
  k_use <- min(4, length(unique(x[is.finite(x)])) - 1)
  # Chooses a small GAM smoothing basis dimension based on unique values
  if (!is.finite(k_use) || k_use < 3) k_use <- 3
  
  # Assigns the result of this operation to an object for later use 
  pretty_x <- ifelse(v %in% names(name_key), name_key[[v]], v)
  
  # Begins creating a ggplot object 
  p <- ggplot(plot_df, aes(x = .data[[v]], y = Complexity)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = k_use), se = FALSE) +
    labs(
      x = pretty_x,
      y = "Complexity"
    ) +
    theme_minimal(base_size = 12)
  
  # Creates a safe output file name for a plot
  file_name <- paste0(gsub("[^A-Za-z0-9]+", "_", v), ".tiff")
  
  # Draws the figure
  tiff(file_name, width = 3, height = 2, units = "in", res = 600)
  print(p)
  dev.off()
}

# ============================================================
# 6) ADD NONLINEAR / TWO-PART PREDICTORS
#     - Flow quadratic (centered)
#     - Width log transform or exponential decay 
#     - Rain Day Before + Day Of as two-part predictors
# ============================================================

# ---- detect column names robustly
flow_col <- grep("^Flow$", names(X), value = TRUE)
if (length(flow_col) == 0) flow_col <- grep("Flow", names(X), value = TRUE)[1]

width_col <- grep("^Width", names(X), value = TRUE)
if (length(width_col) == 0) width_col <- grep("Width", names(X), value = TRUE)[1]

rain_before_col <- grep("^Rain.*Day\\.Before$", names(X), value = TRUE)
if (length(rain_before_col) == 0) rain_before_col <- grep("Rain.*Before", names(X), value = TRUE)[1]

rain_dayof_col <- grep("^Rain.*Day\\.of$", names(X), value = TRUE)
if (length(rain_dayof_col) == 0) rain_dayof_col <- grep("Rain.*Day\\.of", names(X), value = TRUE)[1]

# ---- Flow quadratic (centered)
if (!is.na(flow_col) && flow_col %in% names(df_clean)) {
  Flow_c <- as.numeric(scale(df_clean[[flow_col]], center = TRUE, scale = FALSE))
  # Adds centered Flow to X
  X$Flow_c  <- Flow_c
  # Adds squared centered Flow to X for a quadratic relationship
  X$Flow_c2 <- Flow_c^2
  # Removes raw Flow after creating quadratic Flow terms
  if (flow_col %in% names(X)) X[[flow_col]] <- NULL
} else {
  warning("Flow column not found; skipping Flow quadratic.")
}

# ---- Width log transform
if (!is.na(width_col) && width_col %in% names(df_clean)) {
  # the transformation is log(1+Width)
  X$Width_log <- log1p(df_clean[[width_col]])
  # Removes raw Width after creating the transformed Width variable
  if (width_col %in% names(X)) X[[width_col]] <- NULL
} else {
  warning("Width column not found; skipping Width_log.")
}

# # ---- Uncomment this part for alternative exponential decay function for Width
# # ---- Width exponential-decay transform
# if (!is.na(width_col) && width_col %in% names(df_clean)) {
#   Width_z <- as.numeric(scale(df_clean[[width_col]], center = TRUE, scale = TRUE))
#   X$Width_expdecay <- exp(-Width_z)
# 
#   if (width_col %in% names(X)) X[[width_col]] <- NULL
# } else {
#   warning("Width column not found; skipping Width exponential-decay transform.")
# }

# # ---- Rain two-part (Day Before) log transformed
# if (!is.na(rain_before_col) && rain_before_col %in% names(df_clean)) {
#   rb <- df_clean[[rain_before_col]]
#   X$RainAny_Before <- ifelse(is.na(rb), NA_integer_, as.integer(rb > 0))
#   X$RainLog_Before <- ifelse(is.na(rb), NA_real_, ifelse(rb > 0, rb, 0))
#   if (rain_before_col %in% names(X)) X[[rain_before_col]] <- NULL
# } else {
#   warning("Rain Day Before column not found; skipping RainAny_Before/RainLog_Before.")
# }
# 
# # ---- Rain two-part (Day Of) log transformed
# if (!is.na(rain_dayof_col) && rain_dayof_col %in% names(df_clean)) {
#   rd <- df_clean[[rain_dayof_col]]
#   X$RainAny_DayOf <- ifelse(is.na(rd), NA_integer_, as.integer(rd > 0))
#   X$RainLog_DayOf <- ifelse(is.na(rd), NA_real_, ifelse(rd > 0, rd, 0))
#   if (rain_dayof_col %in% names(X)) X[[rain_dayof_col]] <- NULL
# } else {
#   warning("Rain Day Of column not found; skipping RainAny_DayOf/RainLog_DayOf.")
# }

# # ---- Rain two-part (Day Before): Any rain + quadratic for positive amount
# if (!is.na(rain_before_col) && rain_before_col %in% names(df_clean)) {
#   rb <- df_clean[[rain_before_col]]
# 
#   X$RainAny_Before <- ifelse(is.na(rb), NA_integer_, as.integer(rb > 0))
# 
#   rb_pos <- ifelse(!is.na(rb) & rb > 0, rb, NA_real_)
#   rb_pos_idx <- which(!is.na(rb_pos))
# 
#   if (length(rb_pos_idx) > 0) {
#     rb_center <- rb_pos[rb_pos_idx] - mean(rb_pos[rb_pos_idx], na.rm = TRUE)
# 
#     X$RainLog_Before  <- ifelse(is.na(rb), NA_real_, 0)
#     X$RainLog_Before2 <- ifelse(is.na(rb), NA_real_, 0)
# 
#     X$RainLog_Before[rb_pos_idx]  <- rb_center
#     X$RainLog_Before2[rb_pos_idx] <- rb_center^2
#   } else {
#     X$RainLog_Before  <- ifelse(is.na(rb), NA_real_, 0)
#     X$RainLog_Before2 <- ifelse(is.na(rb), NA_real_, 0)
#   }
# 
#   if (rain_before_col %in% names(X)) X[[rain_before_col]] <- NULL
# } else {
#   warning("Rain Day Before column not found; skipping RainAny_Before/RainLog_Before.")
# }
# 
# # ---- Rain two-part (Day Of): Any rain + quadratic for positive amount
# if (!is.na(rain_dayof_col) && rain_dayof_col %in% names(df_clean)) {
#   rd <- df_clean[[rain_dayof_col]]
# 
#   X$RainAny_DayOf <- ifelse(is.na(rd), NA_integer_, as.integer(rd > 0))
# 
#   rd_pos <- ifelse(!is.na(rd) & rd > 0, rd, NA_real_)
#   rd_pos_idx <- which(!is.na(rd_pos))
# 
#   if (length(rd_pos_idx) > 0) {
#     rd_center <- rd_pos[rd_pos_idx] - mean(rd_pos[rd_pos_idx], na.rm = TRUE)
# 
#     X$RainLog_DayOf  <- ifelse(is.na(rd), NA_real_, 0)
#     X$RainLog_DayOf2 <- ifelse(is.na(rd), NA_real_, 0)
# 
#     X$RainLog_DayOf[rd_pos_idx]  <- rd_center
#     X$RainLog_DayOf2[rd_pos_idx] <- rd_center^2
#   } else {
#     X$RainLog_DayOf  <- ifelse(is.na(rd), NA_real_, 0)
#     X$RainLog_DayOf2 <- ifelse(is.na(rd), NA_real_, 0)
#   }
# 
#   if (rain_dayof_col %in% names(X)) X[[rain_dayof_col]] <- NULL
# } else {
#   warning("Rain Day Of column not found; skipping RainAny_DayOf/RainLog_DayOf.")
# }

# ---- Rain two-part (Day Before): Any rain + spline for positive amount
if (!is.na(rain_before_col) && rain_before_col %in% names(df_clean)) {
  # Stores previous-day rainfall values
  rb <- df_clean[[rain_before_col]]

  # Creates previous-day rainfall occurrence or spline value depending on context
  X$RainAny_Before <- ifelse(is.na(rb), NA_integer_, as.integer(rb > 0))

  # Find positive previous-day rainfall
  rb_pos <- ifelse(!is.na(rb) & rb > 0, rb, NA_real_)
  rb_pos_idx <- which(!is.na(rb_pos))
  # Checks whether a variable has more than one observed value
  rb_n_unique <- length(unique(rb_pos[rb_pos_idx]))

  if (rb_n_unique >= 3) {
    # Sets degrees of freedom for the previous-day rainfall spline
    df_rb <- min(3, rb_n_unique - 1)
    # Creates natural spline basis functions for previous-day rainfall
    rb_ns <- splines::ns(rb_pos[rb_pos_idx], df = df_rb)

    X$RainSpline_Before1 <- ifelse(is.na(rb), NA_real_, 0)
    X$RainSpline_Before2 <- ifelse(is.na(rb), NA_real_, 0)
    X$RainSpline_Before3 <- ifelse(is.na(rb), NA_real_, 0)

    X$RainSpline_Before1[rb_pos_idx] <- rb_ns[, 1]
    if (ncol(rb_ns) >= 2) X$RainSpline_Before2[rb_pos_idx] <- rb_ns[, 2]
    if (ncol(rb_ns) >= 3) X$RainSpline_Before3[rb_pos_idx] <- rb_ns[, 3]
  } else {
    X$RainSpline_Before1 <- ifelse(is.na(rb), NA_real_, ifelse(rb > 0, rb, 0))
  }

  if (rain_before_col %in% names(X)) X[[rain_before_col]] <- NULL
} else {
  warning("Rain Day Before column not found; skipping RainAny_Before/RainSpline_Before.")
}

# ---- Rain two-part (Day Of): Any rain + spline for positive amount, same as above
if (!is.na(rain_dayof_col) && rain_dayof_col %in% names(df_clean)) {
  rd <- df_clean[[rain_dayof_col]]

  X$RainAny_DayOf <- ifelse(is.na(rd), NA_integer_, as.integer(rd > 0))

  rd_pos <- ifelse(!is.na(rd) & rd > 0, rd, NA_real_)
  rd_pos_idx <- which(!is.na(rd_pos))
  rd_n_unique <- length(unique(rd_pos[rd_pos_idx]))

  if (rd_n_unique >= 3) {
    df_rd <- min(3, rd_n_unique - 1)
    rd_ns <- splines::ns(rd_pos[rd_pos_idx], df = df_rd)

    X$RainSpline_DayOf1 <- ifelse(is.na(rd), NA_real_, 0)
    X$RainSpline_DayOf2 <- ifelse(is.na(rd), NA_real_, 0)
    X$RainSpline_DayOf3 <- ifelse(is.na(rd), NA_real_, 0)

    X$RainSpline_DayOf1[rd_pos_idx] <- rd_ns[, 1]
    if (ncol(rd_ns) >= 2) X$RainSpline_DayOf2[rd_pos_idx] <- rd_ns[, 2]
    if (ncol(rd_ns) >= 3) X$RainSpline_DayOf3[rd_pos_idx] <- rd_ns[, 3]
  } else {
    X$RainSpline_DayOf1 <- ifelse(is.na(rd), NA_real_, ifelse(rd > 0, rd, 0))
  }

  if (rain_dayof_col %in% names(X)) X[[rain_dayof_col]] <- NULL
} else {
  warning("Rain Day Of column not found; skipping RainAny_DayOf/RainSpline_DayOf.")
}

# # ---- Rain two-part (Day Of): Any rain + spline for positive amount (df = 2)
# if (!is.na(rain_dayof_col) && rain_dayof_col %in% names(df_clean)) {
#   rd <- df_clean[[rain_dayof_col]]
# 
#   X$RainAny_DayOf <- ifelse(is.na(rd), NA_integer_, as.integer(rd > 0))
# 
#   rd_pos <- ifelse(!is.na(rd) & rd > 0, rd, NA_real_)
#   rd_pos_idx <- which(!is.na(rd_pos))
#   rd_n_unique <- length(unique(rd_pos[rd_pos_idx]))
# 
#   if (rd_n_unique >= 3) {
#     df_rd <- min(2, rd_n_unique - 1)
#     rd_ns <- splines::ns(rd_pos[rd_pos_idx], df = df_rd)
# 
#     X$RainSpline_DayOf1 <- ifelse(is.na(rd), NA_real_, 0)
#     X$RainSpline_DayOf2 <- ifelse(is.na(rd), NA_real_, 0)
# 
#     X$RainSpline_DayOf1[rd_pos_idx] <- rd_ns[, 1]
#     if (ncol(rd_ns) >= 2) X$RainSpline_DayOf2[rd_pos_idx] <- rd_ns[, 2]
#   } else {
#     X$RainSpline_DayOf1 <- ifelse(is.na(rd), NA_real_, ifelse(rd > 0, rd, 0))
#   }
# 
#   if (rain_dayof_col %in% names(X)) X[[rain_dayof_col]] <- NULL
# } else {
#   warning("Rain Day Of column not found; skipping RainAny_DayOf/RainSpline_DayOf.")
# }

# Drop any new zero-variance columns created
keep_cols2 <- sapply(names(X), function(nm) {
  if (nm %in% c("Location", "Wet.Dry")) return(TRUE)
  v <- X[[nm]]
  length(unique(v[!is.na(v)])) > 1
})
X <- X[, keep_cols2, drop = FALSE]


# -------------------------------
# 7) Fit ZERO-INFLATED models (ZIP vs ZINB vs ZICOMP)
# -------------------------------

# Creates the final model data frame containing response and predictors
df_model <- data.frame(
  Complexity = df_clean$Complexity,
  X
)

# Constructs the regression formula automatically from X 
formula_glmm <- as.formula(
  paste("Complexity ~", paste(colnames(X), collapse = " + "))
)

# Specifies an intercept-only zero-inflation model
zi_form <- ~ 1

# alternative if predictor variables need to be included 
# zi_form <- ~ Wet.Dry + Location

# Fits the zero-inflated Poisson model
m_zip   <- glmmTMB(formula_glmm, ziformula = zi_form, data = df_model, family = poisson())
# Fits the zero-inflated negative binomial model
m_zinb  <- glmmTMB(formula_glmm, ziformula = zi_form, data = df_model, family = nbinom2())
# Fits the zero-inflated COM-Poisson model
m_zicom <- glmmTMB(formula_glmm, ziformula = zi_form, data = df_model, family = compois())

# DHARMa diagnostics for all three models 
res_zip   <- DHARMa::simulateResiduals(m_zip,   n = 2000)
res_zinb  <- DHARMa::simulateResiduals(m_zinb,  n = 2000)
res_zicom <- DHARMa::simulateResiduals(m_zicom, n = 2000)

# Tests residual dispersion for all three models 
cat("\n--- DHARMa Dispersion (ZI models) ---\n")
print(DHARMa::testDispersion(res_zip))
print(DHARMa::testDispersion(res_zinb))
print(DHARMa::testDispersion(res_zicom))

# Tests zero inflation whether excess zeros remain 
cat("\n--- DHARMa Zero Inflation (ZI models) ---\n")
print(DHARMa::testZeroInflation(res_zip))
print(DHARMa::testZeroInflation(res_zinb))
print(DHARMa::testZeroInflation(res_zicom))

# Tests for outliers in simulated residuals 
cat("\n--- DHARMa Outliers (ZI models) ---\n")
print(DHARMa::testOutliers(res_zip))
print(DHARMa::testOutliers(res_zinb))
print(DHARMa::testOutliers(res_zicom))

# Computes AIC for three models comparison
cat("\n--- AIC comparison (ZI models) ---\n")
aic_tab <- AIC(m_zip, m_zinb, m_zicom)
print(aic_tab)

# Creates a Word-friendly AIC table
aic_table_word <- data.frame(
  Model = rownames(aic_tab),
  df = aic_tab$df,
  AIC = round(aic_tab$AIC, 3),
  row.names = NULL
)

cat("\n=== AIC comparison of candidate models (TAB-DELIMITED) ===\n")
write.table(aic_table_word, sep = "\t", row.names = FALSE, quote = FALSE)

# Selects the lowest-AIC model name 
best_name <- rownames(aic_tab)[which.min(aic_tab$AIC)]
cat("\nBest ZI model by AIC:", best_name, "\n")

# Stores the selected best model object
m_best <- switch(best_name,
                 "m_zip"   = m_zip,
                 "m_zinb"  = m_zinb,
                 "m_zicom" = m_zicom)

# Stores the residual object for the selected best model
res_best <- switch(best_name,
                   "m_zip"   = res_zip,
                   "m_zinb"  = res_zinb,
                   "m_zicom" = res_zicom)

# # DHARMa residual plots with titles for the best model 
# tiff("res_best.tiff",
#      width = 7, height = 5, units = "in",
#      res = 600)
# 
# # Sets base-R plot margins
# par(oma = c(0, 0, 3, 0))  # add outer top margin
# 
# plot(res_best)
# 
# # Adds a custom title to a base-R plot
# mtext("Best model", side = 3, outer = TRUE, line = 1, cex = 1.2, font = 1, adj = 0.29)
# 
# # Closes the graphics device and writes the image file to disk
# dev.off()
# 
# # plot residual for Poisson model
# tiff("res_zip.tiff",
#      width = 7, height = 5, units = "in",
#      res = 600)
# 
# par(oma = c(0, 0, 3, 0))  # add outer top margin
# 
# plot(res_zip)
# 
# mtext("Poisson model", side = 3, outer = TRUE, line = 1, cex = 1.2, font = 1, adj = 0.26)
# dev.off()

# plot residual for Negative binomial model
tiff("res_zinb.tiff",
     width = 7, height = 5, units = "in",
     res = 600)

par(oma = c(0, 0, 3, 0))  # add outer top margin
plot(res_zinb)
mtext("Negative binomial model", side = 3, outer = TRUE, line = 1, cex = 1.2, font = 1, adj = 0.12)

dev.off()

# plot residual for COM-Poisson model
tiff("res_zicom.tiff",
     width = 7, height = 5, units = "in",
     res = 600)

par(oma = c(0, 0, 3, 0))  # add outer top margin

plot(res_zicom)
mtext("COM-Poisson model", side = 3, outer = TRUE, line = 1, cex = 1.2, font = 1, adj = 0.18)

dev.off()

# Prints the summary of the best model
cat("\n=== BEST MODEL SUMMARY ===\n")
print(summary(m_best))

# Computes ZI/Hurdle R2 
cat("\n=== R² (ZI/Hurdle R2) ===\n")
print(tryCatch(performance::r2(m_best), error = function(e) NA))

# Computes McFadden Pseudo-R2
cat("\n=== McFadden Pseudo-R2 ===\n")
print(tryCatch(performance::r2_mcfadden(m_best), error = function(e) NA))

# -------------------------------
# 8) Effect size tables (conditional + zero-inflation parts)
# -------------------------------

# Conditional (count mean) part, prints the summary of the selected model
coefs_cond <- summary(m_best)$coefficients$cond

# Builds the conditional count-effect table
effect_table_cond <- data.frame(
  Predictor  = rownames(coefs_cond),
  Estimate   = coefs_cond[, "Estimate"],
  # Computes rate ratios by exponentiating count-model coefficients
  RateRatio  = exp(coefs_cond[, "Estimate"]),
  # Computes lower Wald confidence interval bound
  CI_lower   = exp(coefs_cond[, "Estimate"] - 1.96 * coefs_cond[, "Std. Error"]),
  # Computes upper Wald confidence interval bound
  CI_upper   = exp(coefs_cond[, "Estimate"] + 1.96 * coefs_cond[, "Std. Error"]),
  p_value    = coefs_cond[, "Pr(>|z|)"],
  row.names  = NULL
) %>%
 # filter(Predictor != "(Intercept)") %>%
  mutate(
    Predictor = pretty_term_labels(Predictor, name_key)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n=== Conditional (count) effect table (TAB-DELIMITED) ===\n")
write.table(effect_table_cond, sep = "\t", row.names = FALSE, quote = FALSE)

# Zero-inflation (structural zero) part, prints the summary of the selected model
coefs_zi <- summary(m_best)$coefficients$zi

# Builds the zero-inflation effect table
if (!is.null(coefs_zi) && nrow(coefs_zi) > 0) {
  effect_table_zi <- data.frame(
    Predictor  = rownames(coefs_zi),
    Estimate   = coefs_zi[, "Estimate"],
    # Computes odds ratios by exponentiating zero-inflation coefficients
    OddsRatio  = exp(coefs_zi[, "Estimate"]),
    # Computes lower Wald confidence interval bound
    CI_lower   = exp(coefs_zi[, "Estimate"] - 1.96 * coefs_zi[, "Std. Error"]),
    # Computes upper Wald confidence interval bound
    CI_upper   = exp(coefs_zi[, "Estimate"] + 1.96 * coefs_zi[, "Std. Error"]),
    p_value    = coefs_zi[, "Pr(>|z|)"],
    # Converts zero-inflation log-odds to probability
    StructuralZeroProb = plogis(coefs_zi[, "Estimate"]),
    row.names  = NULL
  ) %>%
    mutate(
      Predictor = pretty_term_labels(Predictor, name_key)
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  cat("\n=== Zero-inflation (structural zero) effect table (TAB-DELIMITED) ===\n")
  write.table(effect_table_zi, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  cat("\n(No zero-inflation coefficients found; ziformula may be ~0.)\n")
}
  
# -------------------------------
# 9) Pretty list of predictors in final model
# -------------------------------
model_terms_pretty <- pretty_term_labels(colnames(X), name_key)

cat("\n=== Predictors included in final model ===\n")
print(model_terms_pretty)


# BELOW IS OPTIONAL TO CHECK THE PERFORMENE OF POISSON MODEL 
print(tryCatch(performance::r2(m_zip), error = function(e) NA))
print(tryCatch(performance::r2_mcfadden(m_zip), error = function(e) NA))

coefs_cond_zip <- summary(m_zip)$coefficients$cond

effect_table_cond_zip <- data.frame(
  Predictor  = rownames(coefs_cond_zip),
  Estimate   = coefs_cond_zip[, "Estimate"],
  RateRatio  = exp(coefs_cond_zip[, "Estimate"]),
  CI_lower   = exp(coefs_cond_zip[, "Estimate"] - 1.96 * coefs_cond_zip[, "Std. Error"]),
  CI_upper   = exp(coefs_cond_zip[, "Estimate"] + 1.96 * coefs_cond_zip[, "Std. Error"]),
  p_value    = coefs_cond_zip[, "Pr(>|z|)"],
  row.names  = NULL
) %>%
  # filter(Predictor != "(Intercept)") %>%
  mutate(
    Predictor = pretty_term_labels(Predictor, name_key)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))


cat("\n=== Conditional (count) effect table (TAB-DELIMITED) ===\n")
write.table(effect_table_cond_zip, sep = "\t", row.names = FALSE, quote = FALSE)

# BELOW IS OPTION TO CHECK THE PERFORMENE OF NEGATIVE BINOMIAL MODEL 
print(tryCatch(performance::r2(m_zinb), error = function(e) NA))
print(tryCatch(performance::r2_mcfadden(m_zinb), error = function(e) NA))

coefs_cond_zinb <- summary(m_zinb)$coefficients$cond

effect_table_cond_zinb <- data.frame(
  Predictor  = rownames(coefs_cond_zinb),
  Estimate   = coefs_cond_zinb[, "Estimate"],
  RateRatio  = exp(coefs_cond_zinb[, "Estimate"]),
  CI_lower   = exp(coefs_cond_zinb[, "Estimate"] - 1.96 * coefs_cond_zinb[, "Std. Error"]),
  CI_upper   = exp(coefs_cond_zinb[, "Estimate"] + 1.96 * coefs_cond_zinb[, "Std. Error"]),
  p_value    = coefs_cond_zinb[, "Pr(>|z|)"],
  row.names  = NULL
) %>%
  # filter(Predictor != "(Intercept)") %>%
  mutate(
    Predictor = pretty_term_labels(Predictor, name_key)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n=== Conditional (count) effect table (TAB-DELIMITED) ===\n")
write.table(effect_table_cond_zinb, sep = "\t", row.names = FALSE, quote = FALSE)

