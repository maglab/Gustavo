# Load required libraries
library(dplyr)
library(purrr)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(cowplot)

################################################################################
# DATA LOADING AND PREPROCESSING
################################################################################

# Load main dataset
GenArcAgeFrm <- readRDS('Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARC_GenAge.rds')

# Remove empty gene entries
GenArcAgeFrm <- GenArcAgeFrm %>%
  filter(Gen != "")

################################################################################
# PLEIOTROPY CALCULATION
################################################################################

# Calculate pleiotropy scores (number of phenotypes per gene)
#GenPlt <- GenArcAgeFrm %>% 
#  filter(!(Phn %in% c("HumAge", "ModAge", "Diseases"))) %>% 
#  distinct(Gen, Phn) %>% 
#  count(Gen, name = "Count") %>% 
#  pull(Count, name = Gen)

# Define high and low pleiotropy genes
#HghPlt <- names(GenPlt[GenPlt >= 4])
#LowPlt <- names(GenPlt[GenPlt < 4])


# High pleiotropy frame
GenArdArcFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>%
  dplyr::rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)
High_Pleiotropy_Frame = Get_high_pleiotropy_frame(GenArdArcFrm, c("Phn","Gen"))
Low_Pleiotropy_Frame = Get_low_pleiotropy_frame(GenArdArcFrm, c("Phn","Gen"))
#GenArcAgeFrm = rbind(GenArcAgeFrm,High_Pleiotropy_Frame)


# Add pleiotropy categories to main dataframe
GenArcAgeFrm <- GenArcAgeFrm %>%
  bind_rows(
    High_Pleiotropy_Frame,
    Low_Pleiotropy_Frame
    #data.frame(Phn = "High_ARC_Pleiotropy", Gen = HghPlt),
    #data.frame(Phn = "Low_ARC_Pleiotropy", Gen = LowPlt)
  ) %>%
  distinct()

################################################################################
# GENE SPECIFICITY DATA PROCESSING
################################################################################

# Load tissue specificity (tau) data
OrgTauFrm <- read.csv("Data/Retrieved/Specificity/Tau_gene_V8.csv")

# Remove NA values in tau
OrgTauFrm <- OrgTauFrm %>% 
  filter(!is.na(tau))

# Load gene ID mapping
EnsHgnFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Ranges/HgncEnsembl_hg38_v115.rds")

# Map Ensembl IDs to gene names
EnsHgnFrm <- EnsHgnFrm %>% 
  as_tibble() %>% 
  select(ensembl_gene_id, external_gene_name) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Create clean tau dataframe with gene names
TauFrm <- OrgTauFrm %>%
  left_join(EnsHgnFrm, by = c("gene_id" = "ensembl_gene_id")) %>%
  rename(Gen = external_gene_name) %>%
  select(Gen, tau) %>%
  filter(!is.na(tau), !is.na(Gen), Gen != "")

# Clean up to save memory
rm(OrgTauFrm)

################################################################################
# CALCULATE MEAN TAU FOR EACH PHENOTYPE CATEGORY
################################################################################

# Get unique phenotype categories
ComArr <- unique(GenArcAgeFrm$Phn)

# Get available genes
GenArr <- unique(TauFrm$Gen)

# Calculate mean tau for each category
ArcTauFrm <- map_df(ComArr, function(ComQry) {
  cat("Processing:", ComQry, "\n")
  
  # Get genes for this category
  GenQry <- GenArcAgeFrm %>% 
    filter(Phn == ComQry) %>% 
    pull(Gen) %>% 
    intersect(GenArr)
  
  # Calculate mean tau
  if (length(GenQry) > 0) {
    mean_tau <- TauFrm %>% 
      filter(Gen %in% GenQry) %>% 
      pull(tau) %>% 
      mean(na.rm = TRUE)
  } else {
    mean_tau <- NA
  }
  
  # Return result
  data.frame(Set = ComQry, tau = mean_tau)
})

# Save results
saveRDS(ArcTauFrm, "Data/Generated/Specificity_and_coexpression/Specificity/MeanSpecificity.rds")

################################################################################
# DATA CLEANING FOR ANALYSIS
################################################################################

# Reload data for consistency
GenArcAgeFrm <- readRDS('Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARC_GenAge.rds')
GenArcAgeFrm <- GenArcAgeFrm %>% 
  filter(Gen != "")

# Load mean specificity results
ArcTauFrm <- readRDS("Data/Generated/Specificity_and_coexpression/Specificity/MeanSpecificity.rds")

# Standardize naming conventions
ArcTauFrm <- ArcTauFrm %>%
  mutate(Set = case_when(
    #Set == "immunological/systemic disorders" ~ "immune_disorders",
    Set == "AllArc" ~ "Diseases",
    Set == "High_ARC_Pleiotropy" ~ "High ARC-Pleiotropy",
    Set == "HumAge" ~ "GenAge.Hum",
    Set == "ModAge" ~ "GenAge.Mod",
    TRUE ~ Set
  ))

GenArcAgeFrm <- GenArcAgeFrm %>%
  mutate(Phn = case_when(
    #Phn == "immunological/systemic disorders" ~ "immune_disorders",
    Phn == "AllArc" ~ "Diseases",
    #Phn == "immune_disorders" ~ "High ARC-Pleiotropy",
    Phn == "HumAge" ~ "GenAge.Hum",
    Phn == "ModAge" ~ "GenAge.Mod",
    TRUE ~ Phn
  ))

# Filter TauFrm for valid entries
TauFrm <- TauFrm %>% 
  filter(!is.na(Gen), Gen != "")

################################################################################
# PERMUTATION TEST FUNCTION
################################################################################

perform_tau_permutation_test <- function(TauFrm,           # Dataframe with Gen and tau
                                         ComQry,           # Group/community to analyze
                                         RefVal,           # Reference value (real tau)
                                         LenQry,           # Size of gene set
                                         n_per = 10000,    # Number of permutations
                                         test_statistic = "mean",
                                         seed = 42) {
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Perform permutations
  permuted_stats <- map_dbl(1:n_per, function(i) {
    # Random sample of genes with same size as query set
    sampled_genes <- sample(TauFrm$Gen, size = LenQry, replace = FALSE)
    
    # Get tau values for sampled genes
    permuted_tau <- TauFrm %>%
      filter(Gen %in% sampled_genes) %>%
      pull(tau)
    
    # Calculate statistic for this permutation
    if (test_statistic == "mean") {
      stat_value <- mean(permuted_tau, na.rm = TRUE)
    } else if (test_statistic == "median") {
      stat_value <- median(permuted_tau, na.rm = TRUE)
    }
    
    return(stat_value)
  })
  
  # Remove infinite values
  permuted_stats <- permuted_stats[is.finite(permuted_stats)]
  
  # Calculate two-tailed p-value
  center <- mean(permuted_stats, na.rm = TRUE)
  B <- length(permuted_stats)
  p_value <- (sum(abs(permuted_stats - center) >= abs(RefVal - center)) + 1) / (B + 1)
  
  # Create results list
  results <- list(
    ComQry = ComQry,
    RefVal = RefVal,
    LenQry = LenQry,
    real_statistic = RefVal,
    permuted_stats = permuted_stats,
    p_value = p_value,
    n_permutations = n_per,
    test_statistic = test_statistic
  )
  
  # Add plot to results
  results$plot <- plot_tau_permutation_distribution(results, LenQry=LenQry)
  
  return(results)
}

################################################################################
# PLOTTING FUNCTION
################################################################################

plot_tau_permutation_distribution <- function(permutation_results,n_pvalues=4,LenQry) {
  
  # Prepare data for plotting
  perm_data <- data.frame(statistic = permutation_results$permuted_stats)
  real_stat <- permutation_results$real_statistic
  p_value <- permutation_results$p_value
  ComQry <- permutation_results$ComQry
  
  # Calculate null distribution statistics
  PerMean <- mean(permutation_results$permuted_stats)
  PerMedian <- median(permutation_results$permuted_stats)
  
  # Format values for display
  format_value <- function(x) {
    ifelse(x < 1e-2, sprintf("%.2e", x), round(x, 3))
  }
  
  TrgVal <- format_value(real_stat)
  NulVal <- format_value(PerMean)
  PVal <- format_value(p_value)
  
  padj = p_bonferroni(p_value,n_pvalues) %>% format_value()
  
  # Create subtitle with HTML formatting
  #sub_txt <- paste0(
  #  "P-value: ", PVal,
  #  "<br><span style='color:blue;'>Null mean: ", NulVal, "</span>",
  #  " | <span style='color:red;'>", ComQry, ": ", TrgVal, "</span>"
  #)
  
  sub_txt <- paste0(
    "P_adj: ", padj, " | n=", LenQry,
    "<br><span style='color:blue;'>Null mean: ", NulVal, "</span>",
    " | <span style='color:red;'>", ComQry, ": ", TrgVal, "</span>"
  )
  
  # Create histogram
  p <- ggplot(perm_data, aes(x = statistic)) +
    geom_histogram(aes(y = after_stat(count)), bins = 30,
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = real_stat, color = "red",
               linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = PerMedian, color = "blue",
               linetype = "dashed", linewidth = 1.2) +
    labs(
      title = paste(ComQry, "vs\nNull Distribution"),
      subtitle = sub_txt,
      x = "Mean Tau Value",
      y = "Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_markdown(size = 8),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  return(p)
}

################################################################################
# COMPREHENSIVE ANALYSIS FUNCTION
################################################################################

run_comprehensive_tau_analysis <- function(TauFrm,           # Dataframe with Gen and tau
                                           ArcTauFrm,        # Dataframe with Set and tau (reference values)
                                           GenArcAgeFrm,     # Dataframe with Phn for LenQry
                                           ComArr,           # Array of communities to analyze
                                           n_per = 10000,    # Number of permutations
                                           test_statistic = "mean",
                                           seed = 42) {
  
  results <- list()
  
  for (ComQry in ComArr) {
    cat("Processing:", ComQry, "\n")
    
    # Get reference value
    RefVal <- ArcTauFrm %>% 
      filter(Set == ComQry) %>% 
      pull(tau)
    
    # Get set size
    LenQry <- GenArcAgeFrm %>% 
      filter(Phn == ComQry) %>% 
      nrow()
    
    # Skip if no data
    if (length(RefVal) == 0 || LenQry == 0) {
      warning(paste("Skipping", ComQry, ": No data available"))
      next
    }
    
    # Execute permutation test
    results[[ComQry]] <- perform_tau_permutation_test(
      TauFrm = TauFrm,
      ComQry = ComQry,
      RefVal = RefVal,
      LenQry = LenQry,
      n_per = n_per,
      test_statistic = test_statistic,
      seed = seed
    )
  }
  
  return(results)
}

################################################################################
# RESULTS SUMMARY FUNCTION
################################################################################

create_tau_summary_table <- function(results_list) {
  summary_data <- map_df(names(results_list), function(ComQry) {
    result <- results_list[[ComQry]]
    
    data.frame(
      Community = ComQry,
      Reference_Value = round(result$RefVal, 4),
      Gene_Set_Size = result$LenQry,
      Null_Mean = round(mean(result$permuted_stats), 4),
      Null_SD = round(sd(result$permuted_stats), 4),
      Difference = round(result$RefVal - mean(result$permuted_stats), 4),
      P_Value = result$p_value,
      N_Permutations = result$n_permutations
    )
  })
  
  return(summary_data)
}

################################################################################
# MAIN ANALYSIS
################################################################################

# Define communities to analyze
ComArr <- c("GenAge.Hum", "GenAge.Mod", "Diseases", "High ARC-Pleiotropy")

#GenArcAgeFrm[GenArcAgeFrm$Phn=="immune_disorders","Phn"] = "High ARC-Pleiotropy"

# Run comprehensive analysis
tau_results <- run_comprehensive_tau_analysis(
  TauFrm = TauFrm,
  ArcTauFrm = ArcTauFrm,
  GenArcAgeFrm = GenArcAgeFrm,
  ComArr = ComArr,
  n_per = 10000,
  seed = 42
)

# Create summary table
summary_table_tau <- create_tau_summary_table(tau_results)

# Save results
saveRDS(tau_results, "Data/Generated/Permutations/Specificity/tau_permutation_results.rds")
write.csv(summary_table_tau, "Data/Generated/Permutations/Specificity/tau_permutation_summary.csv", row.names = FALSE)

################################################################################
# VISUALIZATION
################################################################################

# Create combined plot grid
p <- plot_grid(
  tau_results$GenAge.Hum$plot,
  tau_results$GenAge.Mod$plot,
  tau_results$Diseases$plot,
  tau_results$`High ARC-Pleiotropy`$plot,
  labels = "AUTO",
  label_size = 12,
  ncol = 2
)

# Add main title
title <- ggdraw() + 
  draw_label(
    "Tau Permutation Analysis",
    fontface = "bold",
    size = 13.5,
    x = 0.5,
    hjust = 0.5
  )

# Combine title and plots
final_plot <- plot_grid(
  title,
  p,
  ncol = 1,
  rel_heights = c(0.08, 1),
  align = "v"
) +
  theme(
    plot.margin = margin(-5, 10, 0, 0, "pt")
  )

# Display final plot
print(final_plot)



# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

p_bonferroni <- function(p, n) {
  p_adj <- p * n
  if (p_adj > 1) p_adj <- 1
  return(p_adj)
}

pretty_num <- function(x) {
  x <- as.numeric(x)
  ifelse(x == 1, "1", format(x, scientific = TRUE))
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