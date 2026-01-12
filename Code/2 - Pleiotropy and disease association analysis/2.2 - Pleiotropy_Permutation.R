# ============================================================================
# PERMUTATION ANALYSIS OF PLEIOTROPY IN LONGEVITY-ASSOCIATED GENES
# ============================================================================
# Description: This script compares ARC-based pleiotropy between longevity-
#              associated genes (GenAge database) and the background of all
#              human protein-coding genes using permutation testing.
# ============================================================================


# ==============================================================================
# 1) LIBRARIES
# ==============================================================================
# Load required packages with startup messages suppressed
suppressPackageStartupMessages({
  library(dplyr)              # Data manipulation
  library(ggplot2)            # Data visualization
  library(ggtext)             # Enhanced text rendering in ggplot (for markdown/HTML in titles)
  library(cowplot)            # Plot arrangement utilities (for plot_grid function)
  library(EnsDb.Hsapiens.v86) # Ensembl database for human genome annotations (protein-coding genes)
})


# ==============================================================================
# 2) EXECUTABLE BLOCK
# ==============================================================================

## (Optional) Set seed for reproducibility
# Ensures that random processes (like permutations) yield the same results each time
set.seed(42)

## --- Load project-specific input data ---
# Adjust file paths as needed for your specific environment
# Load GenAge human genes (aging-related genes in humans)
GenAgeHumArr <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")

# Load GenAge model organism genes (aging-related genes in model organisms)
GenAgeModArr <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")

# Load ARD/ARC gene annotation frame and rename columns for clarity
GenArdArcFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>%
  dplyr::rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)

## --- Define universe of protein-coding genes ---
## If you already have a gene universe (AllGen list), use it directly and skip this section.
## Here we assume you want to build the universe from the Ensembl database.
#AllGen <- GenArdArcFrm$Gen %>% unique()
#AllGenFrm <- data.frame(Gen = AllGen, stringsAsFactors = FALSE)


# Retrieve protein-coding genes from Ensembl database
# Using EnsDb.Hsapiens.v86 with filter for protein-coding biotype
protein_genes <- genes(EnsDb.Hsapiens.v86, 
                       filter = ~ gene_biotype == "protein_coding")

# Convert to dataframe for easier manipulation
protein_genes_df <- as.data.frame(protein_genes)

# Select relevant columns (remove unnecessary metadata)
protein_genes_df <- protein_genes_df[, c("gene_id", "gene_name", "seqnames", 
                                         "start", "end", "strand", "gene_biotype")]

# Extract unique gene symbols to define the protein-coding gene universe
AllGen = protein_genes_df$gene_name %>% unique()

# Create a basic dataframe with all protein-coding genes
AllGenFrm = data.frame(Gen = AllGen)

# Add logical flags indicating whether each gene is in GenAge databases
# AgeHum: TRUE if gene is in human aging database
# AgeMod: TRUE if gene is in model organism aging database
AllGenFrm = AllGenFrm %>% mutate(AgeHum = ifelse(Gen %in% GenAgeHumArr, TRUE, FALSE),
                                 AgeMod = ifelse(Gen %in% GenAgeModArr, TRUE, FALSE))


## --- Alternative approach for GenAge group labeling (commented out) ---
#AllGenFrm <- AllGenFrm %>%
#  mutate(
#    AgeHum = Gen %in% GenAgeHumArr,
#    AgeMod = Gen %in% GenAgeModArr
#  )

## --- Calculate pleiotropy per gene (count of distinct ARCs per gene) ---
# Pleiotropy here is operationalized as the number of distinct ARC terms associated with each gene
ArcPltFrm <- GenArdArcFrm %>%
  dplyr::select(Arc, Gen) %>%          # Select only ARC and Gene columns
  distinct() %>%                       # Remove duplicate ARC-gene pairs
  pull(Gen) %>%                        # Extract gene column as vector
  table() %>%                          # Count frequency of each gene (number of distinct ARCs)
  as.data.frame()                      # Convert to dataframe

# Rename columns for clarity
colnames(ArcPltFrm) <- c("Gen", "Plt")  # Plt = Pleiotropy score

## --- Integrate all data and clean NA values ---
# Merge gene universe with pleiotropy scores (full outer join)
AllPltAgeGenFrm <- merge(AllGenFrm, ArcPltFrm, by = "Gen", all = TRUE)

# Replace NA pleiotropy values with 0 (genes with no associated ARCs)
AllPltAgeGenFrm$Plt[is.na(AllPltAgeGenFrm$Plt)] <- 0

# Quick check of pleiotropy distribution
AllPltAgeGenFrm$Plt %>% table()

# Ensure logical flags don't have NA values (set to FALSE if NA)
AllPltAgeGenFrm$AgeHum[is.na(AllPltAgeGenFrm$AgeHum)] <- FALSE
AllPltAgeGenFrm$AgeMod[is.na(AllPltAgeGenFrm$AgeMod)] <- FALSE

# Data verification summary - useful for debugging and understanding data structure
cat("=== DATA VERIFICATION ===\n")
cat("Total genes:", nrow(AllPltAgeGenFrm), "\n")
cat("Genes AgeHum TRUE:", sum(AllPltAgeGenFrm$AgeHum), "\n")
cat("Genes AgeMod TRUE:", sum(AllPltAgeGenFrm$AgeMod), "\n")
cat("Pleiotropy range:", range(AllPltAgeGenFrm$Plt), "\n\n")

## --- Permutation analysis (two-sided test) ---
# Compare pleiotropy of human aging genes vs. background
# Using label permutation method which maintains group sizes
results_agehum <- compare_pleiotropy(
  df = AllPltAgeGenFrm,
  group_var = "AgeHum",          # Grouping variable: human aging genes
  n_iterations = 10000,          # Number of permutations for null distribution
  alternative = "two.sided",     # Test both tails (could also use "greater" or "less")
  method = "label",              # Permute group labels while keeping group size constant
  seed = 1,                      # For reproducibility
  verbose = TRUE                 # Print progress messages
)

# Compare pleiotropy of model organism aging genes vs. background
results_agemod <- compare_pleiotropy(
  df = AllPltAgeGenFrm,
  group_var = "AgeMod",          # Grouping variable: model organism aging genes
  n_iterations = 10000,
  alternative = "two.sided",
  method = "label",
  seed = 1,
  verbose = TRUE
)

## --- Console summary output ---
# Helper function to print permutation test results in readable format
print_perm_summary <- function(res, nice_name) {
  if (is.null(res)) {
    cat(nice_name, ": no results\n")
    return()
  }
  cat("=== ", nice_name, " ===\n", sep = "")
  cat("Group size (n): ", res$group_size, "\n", sep = "")
  cat("Observed mean: ", signif(res$observed_mean, 4), "\n", sep = "")
  cat("Null mean:      ", signif(res$null_mean, 4), "\n", sep = "")
  cat("P-value (", res$alternative, "): ", signif(res$p_value, 4), "\n\n", sep = "")
}

# Print summaries for both analyses
print_perm_summary(results_agehum, "GenAge.Hum")
print_perm_summary(results_agemod, "GenAge.Mod")

## --- Generate visualization plots ---
# Create histogram showing permutation distribution for human aging genes
p1 <- plot_histogram_perm(results_agehum, display_name = "GenAge.Hum", bins = 30)

# Create histogram showing permutation distribution for model organism aging genes
p2 <- plot_histogram_perm(results_agemod, display_name = "GenAge.Mod", bins = 30)

# Arrange both plots side by side for comparison
plot_grid(p1, p2, ncol = 2)

# ==============================================================================
# 3) FUNCTIONS (at the end)
# ==============================================================================

# Function to compare pleiotropy between a group of interest and background
# Expected dataframe: must contain columns 'Plt' (numeric) and a logical column specified by group_var
# Parameters:
#   - df: Input dataframe
#   - group_var: Name of logical column defining group membership
#   - n_iterations: Number of permutations for null distribution
#   - alternative: "greater", "less", or "two.sided" hypothesis test
#   - method: "label" (permute group labels) or "resample" (resample values)
#   - seed: Optional random seed for reproducibility
#   - verbose: Whether to print progress messages
compare_pleiotropy <- function(df,
                               group_var,
                               n_iterations = 10000,
                               alternative = c("greater", "less", "two.sided"),
                               method = c("label", "resample"),
                               seed = NULL,
                               verbose = TRUE) {
  # Validate alternative hypothesis and method parameters
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  
  # Input validation checks
  if (!("Plt" %in% names(df))) stop("Dataframe must contain 'Plt' column.")
  if (!(group_var %in% names(df))) stop("Specified group_var column does not exist.")
  if (!is.logical(df[[group_var]])) stop("group_var must be logical (TRUE/FALSE).")
  
  # Clean Plt column - remove rows with NA or non-finite values
  bad_plt <- is.na(df$Plt) | !is.finite(df$Plt)
  if (any(bad_plt)) {
    if (verbose) cat("Warning: removing", sum(bad_plt), "rows with non-finite Plt values.\n")
    df <- df[!bad_plt, ]
  }
  
  # Clean group column - remove rows with NA values
  if (any(is.na(df[[group_var]]))) {
    if (verbose) cat("Warning:", sum(is.na(df[[group_var]])), "NA values in", group_var, "- rows removed.\n")
    df <- df[!is.na(df[[group_var]]), ]
  }
  
  # Identify indices of genes in the group of interest
  idx_grp <- df[[group_var]] == TRUE
  group_size <- sum(idx_grp)
  
  # Check if group has any members
  if (group_size == 0) {
    if (verbose) cat("No genes with", group_var, "= TRUE\n")
    return(NULL)
  }
  
  # Set random seed for reproducibility if specified
  if (!is.null(seed)) set.seed(seed)
  
  # Calculate observed statistic: mean pleiotropy in the group of interest
  observed_mean <- mean(df$Plt[idx_grp])
  
  # Generate null distribution via permutation
  B <- n_iterations
  null_distribution <- numeric(B)  # Initialize vector for null distribution
  
  if (method == "label") {
    # Method 1: Permute group labels while maintaining original group size
    n <- nrow(df)
    for (b in seq_len(B)) {
      idx_perm <- rep(FALSE, n)  # Start with all FALSE
      # Randomly assign TRUE to group_size genes
      idx_perm[sample.int(n, size = group_size, replace = FALSE)] <- TRUE
      # Calculate mean pleiotropy for this random group
      null_distribution[b] <- mean(df$Plt[idx_perm])
    }
  } else {
    # Method 2: Resample pleiotropy values directly
    for (b in seq_len(B)) {
      # Randomly sample pleiotropy scores (without replacement)
      null_distribution[b] <- mean(sample(df$Plt, size = group_size, replace = FALSE))
    }
  }
  
  # Calculate mean of null distribution
  null_mean <- mean(null_distribution)
  
  # Calculate p-value with +1 correction (to avoid p=0)
  if (alternative == "greater") {
    # Right-tailed test: probability of observing this high or higher value
    r <- sum(null_distribution >= observed_mean)
    p_val <- (r + 1) / (B + 1)
  } else if (alternative == "less") {
    # Left-tailed test: probability of observing this low or lower value
    r <- sum(null_distribution <= observed_mean)
    p_val <- (r + 1) / (B + 1)
  } else {
    # Two-tailed test: probability of observing this extreme or more extreme value
    r <- sum(abs(null_distribution - null_mean) >= abs(observed_mean - null_mean))
    p_val <- (r + 1) / (B + 1)
  }
  
  # Return comprehensive results list
  list(
    group_name        = group_var,
    group_size        = group_size,
    observed_mean     = observed_mean,
    null_mean         = null_mean,
    null_sd           = sd(null_distribution),
    p_value           = p_val,
    alternative       = alternative,
    method            = method,
    n_iterations      = B,
    null_distribution = null_distribution
  )
}

# Function to create histogram visualization of permutation test results
# Parameters:
#   - results: Output from compare_pleiotropy() function
#   - display_name: Custom name for display (defaults to group_name)
#   - bins: Number of histogram bins
plot_histogram_perm <- function(results,
                                display_name = NULL,
                                bins = 30) {
  # Handle case where no results were returned
  if (is.null(results)) return(ggplot() + theme_void())
  
  # Prepare data frame for plotting
  perm_data <- data.frame(statistic = results$null_distribution)
  
  # Helper function for formatting numeric values
  fmt <- function(x) ifelse(is.na(x), "NA",
                            ifelse(abs(x) < 1e-3, sprintf("%.2e", x), round(x, 3)))
  
  # Calculate summary statistics
  PerMean <- mean(results$null_distribution)
  PerMedian <- median(results$null_distribution)
  
  # Use custom display name if provided, otherwise use group name
  if (is.null(display_name)) display_name <- results$group_name
  
  # Create subtitle text (plain text version without HTML formatting)
  sub_txt <- paste0(
    "P-value (", results$alternative, "): ", fmt(results$p_value),
    " | Null mean: ", fmt(PerMean),
    " | ", display_name, ": ", fmt(results$observed_mean),
    " | n = ", results$group_size
  )
  
  # Create the plot
  ggplot(perm_data, aes(x = statistic)) +
    # Histogram of null distribution
    geom_histogram(aes(y = after_stat(count)), bins = bins,
                   fill = "lightblue", color = "black", alpha = 0.7) +
    # Vertical line for observed statistic (red)
    geom_vline(xintercept = results$observed_mean, color = "red",
               linetype = "dashed", linewidth = 1.1) +
    # Vertical line for null median (blue)
    geom_vline(xintercept = PerMedian, color = "blue",
               linetype = "dashed", linewidth = 1.1) +
    # Labels and titles
    labs(
      title = paste("Permutation test of ARC-Pleiotropy:\n", display_name, " vs null distribution"),
      subtitle = sub_txt,
      x = "Mean pleiotropy under null",
      y = "Frequency"
    ) +
    # Theme and styling
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),  # Using plain element_text (no HTML formatting)
      axis.title.x  = element_text(size = 9),
      axis.title.y  = element_text(size = 9)
    )
}
