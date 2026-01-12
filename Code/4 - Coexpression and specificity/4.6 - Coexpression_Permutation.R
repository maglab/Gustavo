# ==============================================================================
# COEXPRESSION ANALYSIS FOR ARC AND ARD GENES
# ==============================================================================
# This script analyzes coexpression patterns for ARC and ARD related genes
# using permutation tests to assess statistical significance

library(dplyr)
library(purrr)
library(ggplot2)
library(ggtext)
library(cowplot)
library(data.table)

# ==============================================================================
# 1. LOAD AND PREPROCESS GENE MAPPING DATA
# ==============================================================================

# Load Ensembl-HGNC mapping
EnsHgnFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Ranges/HgncEnsembl_hg38_v115.rds")
rownames(EnsHgnFrm) <- EnsHgnFrm$ensembl_gene_id

# Filter and deduplicate gene mapping
FltEnsHgnFrm <- EnsHgnFrm %>% 
  filter(external_gene_name != "")

NonDupEns <- FltEnsHgnFrm$ensembl_gene_id[!duplicated(FltEnsHgnFrm$external_gene_name)]
FltFltEnsHgnFrm <- FltEnsHgnFrm[NonDupEns, ]

# ==============================================================================
# 2. LOAD GENE SETS
# ==============================================================================

load_gene_frame <- function(file_path, meaning_col, set_name = NULL) {
  # Load and format gene frame
  frame <- readRDS(file_path) %>%
    dplyr::rename(Phn = !!sym(meaning_col), Gen = Gene) %>%
    dplyr::select(Phn, Gen)
  
  if (!is.null(set_name)) {
    unique_genes <- frame$Gen %>% unique()
    all_set_frame <- data.frame(Phn = set_name, Gen = unique_genes)
    frame <- rbind(frame, all_set_frame)
  }
  
  return(frame)
}

# High pleiotropy frame
GenArdArcFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>%
  dplyr::rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)
High_Pleiotropy_Frame = Get_high_pleiotropy_frame(GenArdArcFrm, c("Phn","Gen"))

# Load gene sets
GenArcFrm <- load_gene_frame(
  "Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds",
  "ARC_Meaning",
  "AllArc"
)

# Bind diseases frame and high pleiotropy frame
GenArcFrm = rbind(GenArcFrm,High_Pleiotropy_Frame)

GenArdFrm <- load_gene_frame(
  "Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds",
  "ARD_Meaning",
  "AllArd"
)

# Load aging-related genes
GenAgeHum <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
GenAgeMod <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")

# Create aging gene frames
create_aging_frame <- function(genes, prefix) {
  data.frame(Phn = prefix, Gen = genes)
}

GenAgeHumFrm <- create_aging_frame(GenAgeHum, "HumAge")
GenAgeModFrm <- create_aging_frame(GenAgeMod, "ModAge")

# Combine all gene sets
GenArcAgeFrm <- rbind(GenArcFrm, GenAgeHumFrm, GenAgeModFrm)
GenArdAgeFrm <- rbind(GenArdFrm, GenAgeHumFrm, GenAgeModFrm)

# Save combined frames
saveRDS(GenArcAgeFrm, 'Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARC_GenAge.rds')
saveRDS(GenArdAgeFrm, 'Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_GenAge.rds')

# ==============================================================================
# 3. LOAD AND FILTER COEXPRESSION MATRIX (CORREGIDO)
# ==============================================================================

load_coexpression_matrix <- function(gene_list, hgnc_mapping, matrix_path, output_path = NULL) {
  # Load subset of coexpression matrix for specified genes
  
  # Convert gene symbols to Ensembl IDs
  ens_ids <- hgnc_mapping %>%
    filter(external_gene_name %in% gene_list) %>%
    pull(ensembl_gene_id)
  
  # Read column names first
  colnames <- fread(matrix_path, nrows = 0, sep = "\t")
  keep_cols <- which(names(colnames) %in% ens_ids)
  
  # Read only needed columns
  dt <- fread(
    matrix_path,
    sep = "\t",
    select = c(1, keep_cols)
  )
  
  # Rename first column
  setnames(dt, 1, "Ens")
  
  # Create a mapping from Ensembl ID to gene symbol
  ens_to_gene <- setNames(hgnc_mapping$external_gene_name, hgnc_mapping$ensembl_gene_id)
  
  # Rename columns (skip first column which is Ens)
  for (i in 2:ncol(dt)) {
    col_name <- names(dt)[i]
    if (col_name %in% names(ens_to_gene)) {
      gene_name <- ens_to_gene[[col_name]]
      if (!is.na(gene_name) && gene_name != "") {
        setnames(dt, i, gene_name)
      }
    }
  }
  
  # Convert first column from Ensembl ID to gene symbol
  dt$Ens <- ifelse(dt$Ens %in% names(ens_to_gene), 
                   ens_to_gene[dt$Ens], 
                   dt$Ens)
  
  setnames(dt, "Ens", "Gen")
  
  # Convert to data.frame for consistency with rest of code
  dt_df <- as.data.frame(dt)
  
  # Save if output path provided
  if (!is.null(output_path)) {
    saveRDS(dt_df, output_path)
  }
  
  return(dt_df)
}

# Get unique genes for ARC analysis
arc_genes <- unique(c(GenArcFrm$Gen, GenAgeHum, GenAgeMod))

# Load coexpression matrix for ARC genes
coex_matrix <- load_coexpression_matrix(
  gene_list = arc_genes,
  hgnc_mapping = FltFltEnsHgnFrm,
  matrix_path = "Data/Retrieved/Network_Sources/Coexpression/GeneFriends/human_genes_correlation_matrix.tsv",
  output_path = "Data/Generated/Specificity_and_coexpression/Coexpression/CoexpressionMatrix_AgeAndDisease_vs_All.rds"
)

# Take absolute values for correlation coefficients
coex_matrix <- as.data.frame(coex_matrix)
gene_cols <- setdiff(colnames(coex_matrix), "Ens")
coex_matrix[gene_cols] <- abs(coex_matrix[gene_cols])


# ==============================================================================
# 3.1 MAP ENSEMBL IDs TO HGNC GENE SYMBOLS
# ==============================================================================

coex_matrix = readRDS("Data/Generated/Specificity_and_coexpression/Coexpression/CoexpressionMatrix_AgeAndDisease_vs_All.rds")

# Load Ensembl-HGNC mapping from RDS file
EnsHgnFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Ranges/HgncEnsembl_hg38_v115.rds")

# Set row names to Ensembl gene IDs for easy indexing
rownames(EnsHgnFrm) <- EnsHgnFrm$ensembl_gene_id

# Filter out entries without HGNC gene symbols
FltEnsHgnFrm <- EnsHgnFrm %>% filter(external_gene_name != "")

# Identify unique Ensembl IDs (remove duplicates where multiple Ensembl IDs map to same HGNC symbol)
NonDupEns <- FltEnsHgnFrm$ensembl_gene_id[!duplicated(FltEnsHgnFrm$external_gene_name)]

# Create deduplicated mapping dataframe
FltFltEnsHgnFrm <- FltEnsHgnFrm[NonDupEns, ]

# Rename first column of coexpression matrix to "Ens" (temporary name)
colnames(coex_matrix)[1] <- "Ens"

# Get current column names for iteration
current_colnames <- colnames(coex_matrix)

# Map column names from Ensembl IDs to HGNC symbols
# Iterate through all columns except the first one (which contains gene IDs)
for (i in 2:length(current_colnames)) {
  ens_id <- current_colnames[i]
  
  # Check if this Ensembl ID exists in our mapping table
  if (ens_id %in% rownames(FltFltEnsHgnFrm)) {
    
    # Get corresponding HGNC symbol
    hgnc_symbol <- FltFltEnsHgnFrm[ens_id, "external_gene_name"]
    
    # Only rename if HGNC symbol is valid (not NA and not empty)
    if (!is.na(hgnc_symbol) && hgnc_symbol != "") {
      colnames(coex_matrix)[i] <- hgnc_symbol
    }
  }
}

# Map values in the first column from Ensembl IDs to HGNC symbols
# This converts the gene identifiers in rows from Ensembl to HGNC format
coex_matrix$Ens <- ifelse(coex_matrix$Ens %in% rownames(FltFltEnsHgnFrm),
                          FltFltEnsHgnFrm[coex_matrix$Ens, "external_gene_name"],
                          coex_matrix$Ens)

# Rename first column from "Ens" to "Gen" (standardized name for gene column)
colnames(coex_matrix)[1] <- "Gen"

# Display first 10 rows and columns to verify mapping was successful
cat("First 10 rows and columns after mapping:\n")
print(coex_matrix[1:10, 1:10])

# ==============================================================================
# 3.2 APPLY ABS() TO NUMERIC COLUMNS
# ==============================================================================

# Get all column names from the matrix
all_cols <- colnames(coex_matrix)

# Identify all columns except "Gen" (which contains gene names)
non_gen_cols <- setdiff(all_cols, "Gen")

# Identify which of these columns contain numeric data
# sapply() tests each column with is.numeric() function
coex_matrix = as.data.frame(coex_matrix)

numeric_cols <- sapply(coex_matrix[, non_gen_cols, drop = FALSE], is.numeric)

# Create vector of column names that are numeric (for applying abs())
gene_cols <- non_gen_cols[numeric_cols]

# Log message showing how many numeric columns will be processed
cat("\nApplying abs() to", length(gene_cols), "numeric columns...\n")

# Apply absolute value function to all numeric columns
# This converts correlation coefficients to absolute values (0 to 1 range)
coex_matrix[, gene_cols] <- abs(coex_matrix[, gene_cols])

# Display verification of the transformation
cat("\nValues after abs() transformation (first 5 rows and columns):\n")
print(coex_matrix[1:5, 1:5])

# ==============================================================================
# 4. CALCULATE MEAN COEXPRESSION PROFILES (REVISED)
# ==============================================================================

calculate_mean_coexpression <- function(gene_frame, coexpression_matrix, set_name) {
  # Calculate mean coexpression for each gene set
  
  all_results <- data.frame()
  phenotypes <- unique(gene_frame$Phn)
  
  for (phenotype in phenotypes) {
    # Get genes for this phenotype
    phenotype_genes <- gene_frame %>%
      filter(Phn == phenotype) %>%
      pull(Gen) %>%
      intersect(coexpression_matrix$Gen)
    
    if (length(phenotype_genes) == 0) {
      warning(sprintf("No genes found for phenotype: %s", phenotype))
      next
    }
    
    # Subset matrix for this phenotype's genes
    # Ensure we only select columns that exist in the matrix
    available_genes <- phenotype_genes[phenotype_genes %in% colnames(coexpression_matrix)]
    
    if (length(available_genes) == 0) {
      warning(sprintf("No columns found in matrix for phenotype: %s", phenotype))
      next
    }
    
    # Select relevant columns
    cols_to_select <- c("Gen", available_genes)
    subset_matrix <- coexpression_matrix[, cols_to_select, drop = FALSE]
    
    # Calculate mean coexpression
    subset_matrix$Mean <- rowMeans(subset_matrix[, available_genes, drop = FALSE], na.rm = TRUE)
    
    # Store results
    result <- subset_matrix %>%
      select(Gen, Mean) %>%
      filter(!is.na(Gen) & !is.na(Mean)) %>%
      mutate(Set = phenotype, Dis = set_name)
    
    all_results <- rbind(all_results, result)
  }
  
  return(all_results)
}

# Calculate for ARC
arc_coexpression <- calculate_mean_coexpression(GenArcAgeFrm, coex_matrix, "Arc")

# Calculate for ARD
ard_coexpression <- calculate_mean_coexpression(GenArdAgeFrm, coex_matrix, "Ard")

# ==============================================================================
# 4.1 CALCULATE ARC MEAN (ACROSS ALL ARC PHENOTYPES) - MISSING PART
# ==============================================================================

# Calculate mean coexpression across all ARC phenotypes (excluding special sets)
Arc_Mean <- arc_coexpression %>%
  filter(!(Set %in% c("HumAge", "ModAge", "AllArc"))) %>%   # Exclude these sets
  group_by(Gen) %>%
  summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%  # Average per gene
  mutate(Set = "Arc_Mean", Dis = "Arc") %>%
  select(Gen, Mean, Set, Dis) %>%
  as.data.frame()

# Calculate mean coexpression across all ARD phenotypes (excluding special sets)
Ard_Mean <- ard_coexpression %>%
  filter(!(Set %in% c("HumAge", "ModAge", "AllArd"))) %>%   # Exclude these sets
  group_by(Gen) %>%
  summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%  # Average per gene
  mutate(Set = "Ard_Mean", Dis = "Ard") %>%
  select(Gen, Mean, Set, Dis) %>%
  as.data.frame()

# ==============================================================================
# 4.2 COMBINE ALL RESULTS
# ==============================================================================

# Combine all coexpression results
combined_results <- rbind(arc_coexpression, Arc_Mean, ard_coexpression, Ard_Mean)

# Save combined results
saveRDS(combined_results, 
        "Data/Retrieved/Genes_and_diseases/Diseases/GeneCoexpression_with_diseases.rds")

# Display summary
cat("Coexpression profiles calculated:\n")
cat("- ARC phenotypes:", length(unique(arc_coexpression$Set)), "sets\n")
cat("- ARD phenotypes:", length(unique(ard_coexpression$Set)), "sets\n")
cat("- Total unique genes in Arc_Mean:", nrow(Arc_Mean), "\n")
cat("- Total unique genes in Ard_Mean:", nrow(Ard_Mean), "\n")
cat("- Combined results:", nrow(combined_results), "rows\n")

# Verify structure
cat("\nStructure of combined results:\n")
print(str(combined_results))

# Show unique sets
cat("\nUnique sets in combined results:\n")
print(unique(combined_results$Set))

# ==============================================================================
# 5. PLEIOTROPY ANALYSIS
# ==============================================================================

# Calculate pleiotropy scores
#pleiotropy_data <- GenArcAgeFrm %>%
#  filter(!(Phn %in% c("HumAge", "ModAge", "AllArc"))) %>%
#  group_by(Gen) %>%
#  summarise(Count = n(), .groups = "drop") %>%
#  mutate(Phn = ifelse(Count >= 4, "High ARC-Pleiotropy", "Low ARC-Pleiotropy"),
#         Count = NULL)

# Update gene frame with pleiotropy information
#GenArcAgeFrm <- rbind(GenArcAgeFrm, pleiotropy_data)

# Standardize naming
GenArcAgeFrm <- GenArcAgeFrm %>%
  mutate(Phn = case_when(
    Phn == "AllArc" ~ "Diseases",
    Phn == "HumAge" ~ "GenAge.Hum",
    Phn == "ModAge" ~ "GenAge.Mod",
    TRUE ~ Phn
  ))

arc_coexpression <- arc_coexpression %>%
  mutate(Set = case_when(
    Set == "AllArc" ~ "Diseases",
    Set == "HumAge" ~ "GenAge.Hum",
    Set == "ModAge" ~ "GenAge.Mod",
    TRUE ~ Set
  ))

# ==============================================================================
# 6. PERMUTATION TEST FUNCTIONS
# ==============================================================================

perform_cox_permutation_test <- function(cox_frame, com_query, ref_val, len_query, 
                                         n_perm = 10000, test_statistic = "mean", 
                                         seed = 42) {
  # Perform permutation test for coexpression values
  
  if (!is.null(seed)) set.seed(seed)
  
  permuted_stats <- map_dbl(1:n_perm, function(i) {
    sampled_genes <- sample(cox_frame$Gen, size = len_query, replace = FALSE)
    
    permuted_cox <- cox_frame %>%
      filter(Gen %in% sampled_genes) %>%
      pull(cox)
    
    if (test_statistic == "mean") {
      return(mean(permuted_cox, na.rm = TRUE))
    } else if (test_statistic == "median") {
      return(median(permuted_cox, na.rm = TRUE))
    } else {
      stop("test_statistic must be 'mean' or 'median'")
    }
  })
  
  # Calculate two-tailed p-value
  permuted_stats <- permuted_stats[is.finite(permuted_stats)]
  center <- mean(permuted_stats, na.rm = TRUE)
  p_value <- (sum(abs(permuted_stats - center) >= abs(ref_val - center)) + 1) / 
    (length(permuted_stats) + 1)
  
  return(list(
    ComQry = com_query,
    RefVal = ref_val,
    LenQry = len_query,
    real_statistic = ref_val,
    permuted_stats = permuted_stats,
    p_value = p_value,
    n_permutations = n_perm,
    test_statistic = test_statistic
  ))
}

plot_cox_permutation_distribution <- function(permutation_results,n_pvalues=4,len_query) {
  # Plot permutation test results
  
  perm_data <- data.frame(statistic = permutation_results$permuted_stats)
  real_stat <- permutation_results$real_statistic
  p_value <- permutation_results$p_value
  com_query <- permutation_results$ComQry
  
  # Calculate null distribution statistics
  null_mean <- mean(permutation_results$permuted_stats)
  
  # Format values for display
  format_value <- function(x) {
    if (is.na(x)) return("NA")
    ifelse(x < 1e-2, sprintf("%.2e", x), sprintf("%.3f", x))
  }
  
  # Create subtitle with formatted values
  #subtitle <- paste0(
  #  "P-value: ", format_value(p_value),
  #  "<br><span style='color:blue;'>Null mean: ", format_value(null_mean), "</span>",
  #  " | <span style='color:red;'>", com_query, ": ", format_value(real_stat), "</span>"
  #)
  
  pretty_num <- function(x) {
    x <- as.numeric(x)
    ifelse(x == 1, "1", format(x, scientific = TRUE))
  }
  
  padj = p_bonferroni(p_value,n_pvalues) #%>% pretty_num()
  
  subtitle <- paste0(
    "P_adj: ", format_value(padj), " |  n=", len_query,
    "<br><span style='color:blue;'>Null mean: ", format_value(null_mean), "</span>",
    " | <span style='color:red;'>", com_query, ": ", format_value(real_stat), "</span>"
  )
  
  # Create plot
  ggplot(perm_data, aes(x = statistic)) +
    geom_histogram(aes(y = after_stat(count)), bins = 30,
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = real_stat, color = "red",
               linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = null_mean, color = "blue",
               linetype = "dashed", linewidth = 1.2) +
    labs(
      title = paste(com_query, "vs\nNull Distribution"),
      subtitle = subtitle,
      x = "Mean Coexpression Value",
      y = "Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_markdown(size = 8),
      axis.title = element_text(size = 9)
    )
}


# ==============================================================================
# 7. RUN PERMUTATION TESTS
# ==============================================================================

# Define communities to test
communities <- c("GenAge.Hum", "GenAge.Mod", "Diseases", "High ARC-Pleiotropy")

# Prepare coexpression data
coexpression_data <- combined_results %>%
  filter(Set == "Arc_Mean") %>%
  select(Gen, Mean) %>%
  rename(cox = Mean)

# Run permutation tests
permutation_results <- list()

for (community in communities) {
  cat("Processing:", community, "\n")
  
  # Get reference value
  ref_val <- arc_coexpression %>%
    filter(Set == community) %>%
    pull(Mean) %>%
    mean(na.rm = TRUE)
  
  # Get gene set size
  len_query <- GenArcAgeFrm %>%
    filter(Phn == community) %>%
    nrow()
  
  # Run permutation test
  permutation_results[[community]] <- perform_cox_permutation_test(
    cox_frame = coexpression_data,
    com_query = community,
    ref_val = ref_val,
    len_query = len_query,
    n_perm = 10000
  )
  
  # Add plot to results
  permutation_results[[community]]$plot <- plot_cox_permutation_distribution(
    permutation_results=permutation_results[[community]],
    len_query = len_query
  )
}

# ==============================================================================
# 8. SAVE AND VISUALIZE RESULTS
# ==============================================================================

# Save results
saveRDS(permutation_results, 
        "Data/Generated/Permutations/Coexpression/cox_permutation_results.rds")

# Create summary table
create_summary_table <- function(results_list) {
  summary_data <- data.frame()
  
  for (community in names(results_list)) {
    result <- results_list[[community]]
    
    summary_row <- data.frame(
      Community = community,
      Reference_Value = round(result$RefVal, 4),
      Gene_Set_Size = result$LenQry,
      Null_Mean = round(mean(result$permuted_stats), 4),
      Null_SD = round(sd(result$permuted_stats), 4),
      Difference = round(result$RefVal - mean(result$permuted_stats), 4),
      P_Value = result$p_value,
      N_Permutations = result$n_permutations
    )
    
    summary_data <- rbind(summary_data, summary_row)
  }
  
  return(summary_data)
}

summary_table <- create_summary_table(permutation_results)
write.csv(summary_table, 
          "Data/Generated/Permutations/Coexpression/cox_permutation_summary.csv", 
          row.names = FALSE)

# ==============================================================================
# 9. CREATE FINAL VISUALIZATION
# ==============================================================================

# Extract plots
plots <- map(communities, ~ permutation_results[[.]]$plot)

# Create plot grid
plot_grid <- plot_grid(
  plotlist = plots,
  labels = "AUTO",
  label_size = 12,
  ncol = 2
)

# Add title
title_plot <- ggdraw() + 
  draw_label("Permutation Tests of Mean Inter-Set Coexpression with ARC-Related Genes",
             fontface = "bold", x = 0.5, hjust = 0.5, size = 14)

# Combine title and plots
final_plot <- plot_grid(title_plot, plot_grid, ncol = 1, rel_heights = c(0.05, 1))





final_plot <- ggdraw() +
  draw_plot(
    plot_grid(title_plot, plot_grid, ncol = 1, rel_heights = c(0.05, 1)),
    x = 0.05,   # margen izquierdo
    y = 0,
    width = 0.90,  # 1 - (izq + der)
    height = 1
  )




# Display final plot
print(final_plot)


# ==============================================================================
# SESSION INFO
# ==============================================================================

cat("\n=== SESSION INFO ===\n")
sessionInfo()
cat("\nAnalysis completed successfully!\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

p_bonferroni <- function(p, n) {
  p_adj <- p * n
  if (p_adj > 1) p_adj <- 1
  return(p_adj)
}


Get_high_pleiotropy_genes = function(GenArdArcFrm){
  GenPlt = GenArdArcFrm %>% dplyr::select(Gen,Arc) %>% unique() %>% pull(Gen) %>% table()
  HhgPltGen = GenPlt[GenPlt>=4] %>% names()
  return(HhgPltGen)
}

Get_high_pleiotropy_frame = function(GenArdArcFrm,ColNames = c("Arc","Gen")){
  GenPlt = GenArdArcFrm %>% dplyr::select(Gen,Arc) %>% unique() %>% pull(Gen) %>% table()
  HhgPltGen = GenPlt[GenPlt>=4] %>% names()
  HghPltFrm = data.frame("High ARC-Pleiotropy",HhgPltGen)
  colnames(HghPltFrm) = ColNames
  return(HghPltFrm)
}

Get_low_pleiotropy_genes = function(GenArdArcFrm){
  GenPlt = GenArdArcFrm %>% dplyr::select(Gen,Arc) %>% unique() %>% pull(Gen) %>% table()
  LowPltGen = GenPlt[GenPlt<=3] %>% names()
  return(LowPltGen)
}

Get_low_pleiotropy_frame = function(GenArdArcFrm,ColNames = c("Arc","Gen")){
  GenPlt = GenArdArcFrm %>% dplyr::select(Gen,Arc) %>% unique() %>% pull(Gen) %>% table()
  LowPltGen = GenPlt[GenPlt<=3] %>% names()
  LowPltFrm = data.frame("Low ARC-Pleiotropy",LowPltGen)
  colnames(LowPltFrm) = ColNames
  return(LowPltFrm)
}