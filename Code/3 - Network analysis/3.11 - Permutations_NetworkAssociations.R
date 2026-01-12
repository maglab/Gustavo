# ==============================================================================
# LIBRARIES
# ==============================================================================
library(dplyr)      # Data manipulation
library(tidyr)      # Data tidying
library(purrr)      # Functional programming
library(ggplot2)    # Data visualization
library(igraph)     # Network analysis
library(tibble)     # Modern data frames
library(stringr)    # String manipulation
library(ggtext)     # Enhanced text in plots
library(cowplot)    # Plot arrangement
library(openxlsx)   # Excel file operations

# ==============================================================================
# DATA PREPARATION FUNCTIONS
# ==============================================================================

#' Convert proximity to distance metric
#' @param x Proximity value
#' @return Distance value (1/x - 1)
Prx2Dst <- function(x) {
  (1 / x) - 1
}

#' Normalize key columns (Gen and Ntw) for consistent merging
#' @param df Input dataframe
#' @return Dataframe with trimmed character keys
norm_keys <- function(df) {
  df %>%
    mutate(
      Gen = trimws(as.character(Gen)),
      Ntw = trimws(as.character(Ntw))
    )
}

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

# Load gene-ARD-ARC mapping
GenArdArcFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)

# Calculate ARC pleiotropy (number of ARCs per gene)
ArcPltFrm <- GenArdArcFrm %>% 
  dplyr::select(Gen, Arc) %>%
  unique() %>%
  pull(Gen) %>%
  table() %>% 
  sort() %>%
  rev() %>%
  as.data.frame() %>%
  setNames(c("Gen", "NumArcPlt"))

# Load proximity and interactor data
ArcPrx_ArcInt_Frm <- read.csv("Data/Generated/Networks_and_predictions/Topological_data/Distance_and_centrality/Proximity_and_Interactors.csv") %>% 
  dplyr::rename(
    Gen = Gene,
    GenTyp = GeneType,
    AgeDst = Distance_to_GenAgeHum,
    ModDst = Distance_to_GenAgeMod,
    DisDst = Distance_to_ARCs,
    NumArcInt = Number_Of_Neighbouring_ARCs,
    MeanPrx2Arc = MeanProximity_to_ARCs,
    MeanDst2Arc = MeanDistance_to_ARCs,
    Ntw = Network
  ) %>% 
  merge(ArcPltFrm, by = "Gen", all = TRUE) %>%
  mutate(NumArcPlt = ifelse(is.na(NumArcPlt), 0, NumArcPlt))


# Load random walk association data
ArcArdAss <- "Arc"  # Options: "Arc", "Ard"
ArcRnd_Frm <- read.csv(
  paste("Data/Generated/Networks_and_predictions/Topological_data/RWR/", 
        toupper(ArcArdAss), "/Summaries_and_Multiplex/","RWR_ARC_MeanAssociation_MonoMulti.csv", sep = "")
) %>%
  rename(
    Gen = Gene,
    MeanWlk2Arc = Mean_Score_to_Trait,
    Ntw = Network
  ) %>%
  select(Gen, MeanWlk2Arc, Ntw)

# Normalize keys for consistent merging
ArcRnd2  <- norm_keys(ArcRnd_Frm)
ArcPrx2  <- norm_keys(ArcPrx_ArcInt_Frm)

# Merge proximity and random walk data
ArcPrx_ArcInt_ArcWlk_Frm <- ArcRnd2 %>%
  left_join(ArcPrx2, by = c("Gen", "Ntw"), suffix = c("", ".arcint"))

# ==============================================================================
# CORE ANALYSIS FUNCTIONS
# ==============================================================================

#' Generate permutation test distribution plot
#' @param permutation_results List containing permutation test results
#' @return ggplot object showing null distribution and observed value
plot_permutation_distribution <- function(permutation_results) {
  
  # Extract summary statistics
  PerMean <- mean(permutation_results$permuted_stats, na.rm = TRUE)
  PerMedian <- median(permutation_results$permuted_stats, na.rm = TRUE)
  real_stat <- permutation_results$real_statistic
  p_value <- permutation_results$p_value
  
  # Prepare data for histogram
  perm_data <- data.frame(statistic = permutation_results$permuted_stats)
  x <- perm_data$statistic
  
  # Calculate optimal bin count
  target_bins <- 30
  uniq_vals <- length(unique(na.omit(x)))
  bins <- max(1, min(target_bins, uniq_vals))
  
  # Define text based on score type
  score_texts <- list(
    "Wlk" = list(text = "RWR score to ARCs", title = "RWR score to ARCs"),
    "Prx" = list(text = "Proximity to ARCs", title = "Proximity to ARCs"),
    "Int" = list(text = "number of indirectly connected ARCs", title = "ARC.Interactions to ARCs"),
    "Plt" = list(text = "number of directly associated ARCs", title = "ARC.Pleiotropy to ARCs")
  )
  
  ScrInfo <- score_texts[[permutation_results$score_type]]
  
  # Define group labels
  group_labels <- list(
    "AgeGen" = list(label = "GenAge.Hum", x_label = paste("Mean", ScrInfo$text), y_label = "Count"),
    "ModGen" = list(label = "GenAge.Mod", x_label = paste("Mean", ScrInfo$text), y_label = "Count"),
    "DisGen" = list(label = "Diseases", x_label = paste("Mean", ScrInfo$text), y_label = "Count"),
    "DisNgb" = list(label = "Neighbours", x_label = paste("Mean", ScrInfo$text), y_label = "Count"),
    "HghGen" = list(label = "High.Pleiotropy", x_label = paste("Mean", ScrInfo$text), y_label = "Count")
  )
  
  GrpInfo <- group_labels[[permutation_results$target_group]]
  
  # Format values for display
  format_value <- function(val) {
    if (val < 1e-2) sprintf("%.2e", val) else round(val, 2)
  }
  
  TrgVal <- format_value(real_stat)
  NulVal <- format_value(PerMean)
  PVal <- format_value(p_value)
  
  # Create subtitle with formatted text
  sub_txt <- paste0(
    "P-value: ", PVal,
    "<br><span style='color:blue;'>Null mean: ", NulVal, "</span>",
    " | <span style='color:red;'>", GrpInfo$label, ": ", TrgVal, "</span>"
  )
  
  # Generate plot
  p <- ggplot(perm_data, aes(x = statistic)) +
    geom_histogram(aes(y = after_stat(count)), bins = bins,
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = real_stat, color = "red",
               linetype = "dashed", linewidth = 1.2) +
    geom_vline(xintercept = PerMedian, color = "blue",
               linetype = "dashed", linewidth = 1.2) +
    labs(
      title = paste(GrpInfo$label, " vs Null distribution", sep = ""),
      subtitle = sub_txt,
      x = GrpInfo$x_label,
      y = GrpInfo$y_label
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = ggtext::element_markdown(size = 7),
      axis.title.x = element_text(size = 8, margin = margin(t = 6)),
      axis.title.y = element_text(size = 8, margin = margin(r = 6)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    )
  
  return(p)
}

#' Perform single permutation test
#' @param score_df Dataframe with Gen, Scr, Deg columns
#' @param group_assignment_df Dataframe with Gen, GenTyp columns
#' @param all_genes Vector of all unique genes
#' @param target_group Gene group to analyze (e.g., "AgeGen")
#' @param n_permutations Number of permutations (default: 1000)
#' @param balance_by_degree Whether to balance sampling by degree (default: FALSE)
#' @param test_statistic Statistic to use: "mean" or "median" (default: "mean")
#' @param network Network identifier
#' @param score_type Type of score being analyzed
#' @param i Iteration counter for plotting
#' @param seed Random seed (default: 42)
#' @return List containing permutation test results
perform_permutation_test <- function(score_df, group_assignment_df, all_genes,
                                     target_group, n_permutations = 1000,
                                     balance_by_degree = FALSE, test_statistic = "mean",
                                     network, score_type, i,
                                     association_type,
                                     seed = 42) {
  
  # Extract target group genes
  group_genes <- group_assignment_df %>%
    filter(GenTyp == target_group) %>%
    distinct(Gen) %>%
    pull(Gen)
  
  # Calculate real statistic
  real_scores <- score_df %>%
    filter(Gen %in% group_genes) %>%
    pull(Scr)
  
  
  if (test_statistic == "mean") {
    real_stat <- mean(real_scores, na.rm = TRUE) #%>% Prx2Dst()
  } else if (test_statistic == "median") {
    real_stat <- median(real_scores, na.rm = TRUE) #%>% Prx2Dst()
  }
  
  if(association_type=="Dst"){
    real_stat = real_stat %>% Prx2Dst()
  }
  
  # Prepare degree-stratified sampling if requested
  if (balance_by_degree) {
    degree_data <- score_df %>%
      distinct(Gen, .keep_all = TRUE) %>%
      mutate(degree_stratum = ntile(Deg, 5))  # 5 equal-size groups
  }
  
  # Perform permutations
  if (!is.null(seed)) set.seed(seed)
  
  permuted_stats <- map_dbl(1:n_permutations, function(iter) {
    
    if (balance_by_degree) {
      # Stratified sampling by degree
      sampled_data <- degree_data %>%
        group_by(degree_stratum) %>%
        mutate(
          target_size = max(1, floor(length(group_genes) * n() / nrow(degree_data))),
          actual_size = min(n(), target_size)
        ) %>%
        sample_n(size = first(actual_size), replace = FALSE)
      
      sampled_genes <- sampled_data %>% pull(Gen)
    } else {
      # Simple random sampling
      sampled_genes <- sample(all_genes, size = length(group_genes), replace = FALSE)
    }
    
    # Calculate statistic for permuted sample
    permuted_scores <- score_df %>%
      filter(Gen %in% sampled_genes) %>%
      pull(Scr)
    
    if (test_statistic == "mean") {
      stat_value <- mean(permuted_scores, na.rm = TRUE) #%>% Prx2Dst()
    } else if (test_statistic == "median") {
      stat_value <- median(permuted_scores, na.rm = TRUE) #%>% Prx2Dst()
    }
    
    if(association_type=="Dst"){
      stat_value = stat_value %>% Prx2Dst()
    }
    
    return(stat_value)
  })
  
  # Calculate two-sided p-value
  perm_center <- mean(permuted_stats, na.rm = TRUE)
  B <- length(permuted_stats)
  p_value <- (sum(abs(permuted_stats - perm_center) >= abs(real_stat - perm_center)) + 1) / (B + 1)
  
  # Compile results
  results <- list(
    target_group = target_group,
    n_genes = length(group_genes),
    real_statistic = real_stat,
    permuted_stats = permuted_stats,
    p_value = p_value,
    n_permutations = n_permutations,
    balance_by_degree = balance_by_degree,
    test_statistic = test_statistic,
    network = network,
    score_type = score_type,
    i = i
  )
  
  # Generate plot
  results$plot <- plot_permutation_distribution(results)
  
  return(results)
}

#' Run comprehensive permutation analysis across multiple groups and balance types
#' @param score_df Dataframe with Gen, Scr, Deg columns
#' @param group_assignment_df Dataframe with Gen, GenTyp columns
#' @param all_genes Vector of all unique genes
#' @param group_types Vector of gene groups to analyze
#' @param n_permutations Number of permutations (default: 1000)
#' @param test_statistic Statistic to use (default: "mean")
#' @param network Network identifier
#' @param score_type Type of score being analyzed
#' @param i Iteration counter
#' @param seed Random seed (default: 42)
#' @return Nested list containing all permutation results
run_comprehensive_permutation_analysis <- function(score_df, group_assignment_df, 
                                                   all_genes, group_types,
                                                   n_permutations = 1000,
                                                   test_statistic = "mean",
                                                   network, score_type, i, 
                                                   association_type,
                                                   seed = 42) {
  
  results <- list()
  
  # Test both balanced and unbalanced sampling
  for (balance_type in c(FALSE)) {
    balance_label <- ifelse(balance_type, "Balanced_by_Degree", "Unbalanced")
    
    cat("Running", balance_label, "permutations...\n")
    
    group_results <- list()
    group_type = group_types[1]
    for (group_type in group_types) {
      cat("  Processing", group_type, "...\n")
      
      result <- perform_permutation_test(
        score_df = score_df,
        group_assignment_df = group_assignment_df,
        all_genes = all_genes,
        target_group = group_type,
        n_permutations = n_permutations,
        balance_by_degree = balance_type,
        test_statistic = test_statistic,
        network = network,
        score_type = score_type,
        i = i,
        association_type,
        seed = seed
      )
      
      group_results[[group_type]] <- result
    }
    
    results[[balance_label]] <- group_results
  }
  
  return(results)
}

#' Create comprehensive summary table from permutation results
#' @param results_list Nested list of permutation results
#' @return Dataframe with summary statistics
create_comprehensive_summary <- function(results_list) {
  summary_data <- data.frame()
  
  for (balance_type in names(results_list)) {
    for (group_type in names(results_list[[balance_type]])) {
      result <- results_list[[balance_type]][[group_type]]
      
      summary_row <- data.frame(
        Group = group_type,
        Balance_Type = balance_type,
        N_Genes = result$n_genes,
        Real_Statistic = round(result$real_statistic, 4),
        Permuted_Mean = round(mean(result$permuted_stats), 4),
        Permuted_SD = round(sd(result$permuted_stats), 4),
        P_Value = result$p_value,
        N_Permutations = result$n_permutations,
        Test_Statistic = result$test_statistic
      )
      
      summary_data <- bind_rows(summary_data, summary_row)
    }
  }
  
  return(summary_data)
}

# ==============================================================================
# MAIN ANALYSIS EXECUTION
# ==============================================================================

# Analysis parameters
NumPer <- 10000      # Number of permutations
seed <- 42           # Random seed for 

# Score type options:
# Int = ARC_Interactions
# Prx = Proximity  
# Dst = Distance
# Wlk = Random walk
ScrOpc <- "Dst"      # Score type: "Int", "Dst", "Prx", "Wlk", "Plt"

GenSet <- c("DisGen", "AgeGen", "ModGen", "DisNgb", "HghGen")
NtwArr <- c('PPI', 'COX90', 'COX95', 'KEGG')

# Execute analysis for each network
results <- list()
i <- 1

#max(ArcPrx_ArcInt_ArcWlk_Frm$MeanDst2Arc)

#x = ArcPrx_ArcInt_ArcWlk_Frm$MeanPrx2Arc

# Cantidad de NA
#n_NA <- sum(is.na(x))

# Cantidad de Inf (incluye Inf y -Inf)
#n_Inf <- sum(is.infinite(x))



ArcPrx_ArcInt_ArcWlk_Frm$MeanDst2Arc = (1/ArcPrx_ArcInt_ArcWlk_Frm$MeanPrx2Arc) - 1

for (NtwQry in NtwArr) {
  cat("\n=== Processing network:", NtwQry, "===\n")
  
  # Load network data
  network_file <- paste("Data/Generated/Networks_and_predictions/Networks/", 
                        NtwQry, "/Lists/GraphsList.rds", sep = '')
  Grp <- readRDS(network_file)
  
  # Calculate node degrees
  DegFrm <- igraph::degree(Grp$GenGenNtw) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gen") %>%
    setNames(c("Gen", "Deg"))
  
  all_genes <- DegFrm$Gen
  
  # Select appropriate score column based on analysis type
  ScrFrm <- ArcPrx_ArcInt_ArcWlk_Frm %>% 
    filter(Ntw %in% NtwQry)
  
  # Map score option to column name
  score_columns <- list(
    "Int" = "NumArcInt",
    "Dst" = "MeanPrx2Arc",
    "Prx" = "MeanPrx2Arc",
    "Wlk" = "MeanWlk2Arc",
    "Plt" = "NumArcPlt"
  )
  
  selected_col <- score_columns[[ScrOpc]]
  ScrFrm <- ScrFrm %>% 
    select(Gen, !!selected_col) %>% 
    unique() %>%
    setNames(c("Gen", "Scr"))
  
  # Merge scores with degrees
  score_df <- merge(ScrFrm, DegFrm, by = "Gen")
  
  # Prepare group assignments
  group_assignment_df <- ArcPrx_ArcInt_ArcWlk_Frm %>% 
    filter(Ntw %in% NtwQry) %>%
    select(Gen, GenTyp) %>%
    filter(GenTyp %in% GenSet)
  
  # Run permutation analysis
  results[[NtwQry]] <- run_comprehensive_permutation_analysis(
    score_df = score_df,
    group_assignment_df = group_assignment_df,
    all_genes = all_genes,
    group_types = GenSet,
    n_permutations = NumPer,
    test_statistic = "mean",
    network = NtwQry,
    score_type = ScrOpc,
    i = i,
    association_type = ScrOpc,
    seed = seed
  )
  
  i <- i + 1
}


if(ScrOpc == "Int")
ScrOpcTxt = "Interactors"
if(ScrOpc == "Dst")
  ScrOpcTxt = "Distance"
if(ScrOpc == "Prx")
  ScrOpcTxt = "Proximity"
if(ScrOpc == "Wlk")
  ScrOpcTxt = "Walker"

# Save results
results_file <- paste("Data/Generated/Permutations/Network_Associations/PermutationList_",
                      toupper(ArcArdAss),'_', ScrOpcTxt,"_", NumPer, ".rds", sep = "")
#results_file <- paste("Data/Generated/Permutations/Network_Associations/PrmLst_",
#                      ArcArdAss, ScrOpc, "Dst_", NumPer, ".rds", sep = "")
saveRDS(results, results_file)

# ==============================================================================
# RESULTS PROCESSING AND VISUALIZATION
# ==============================================================================

# Load results if needed
results <- readRDS(results_file)

# Process results into dataframe
#GenSetArr <- c('AgeGen', 'ModGen', 'DisGen', 'DisNgb')
GenSetArr <- c('AgeGen', 'ModGen', 'DisGen', 'HghGen')

ResFrm <- data.frame()

for (NtwQry in NtwArr) {
  for (GenSetQry in GenSetArr) {
    QryLst <- results[[NtwQry]]$Unbalanced[[GenSetQry]]
    
    # Map group names to display names
    group_names <- list(
      "AgeGen" = "GenAge.Hum",
      "ModGen" = "GenAge.Mod",
      "DisGen" = "Diseases",
      'HghGen' = "High ARC-Pleiotropy"
      #"DisNgb" = "Neighbours"
    )
    
    GenSetTxt <- group_names[[GenSetQry]]
    
    # Create summary row
    sResFrm <- data.frame(
      Network = NtwQry,
      GeneSet = GenSetTxt,
      ObservedValue = QryLst$real_statistic,
      NullMean = mean(QryLst$permuted_stats),
      Difference = QryLst$real_statistic - mean(QryLst$permuted_stats),
      pvalue = QryLst$p_value
    )
    
    ResFrm <- rbind(ResFrm, sResFrm)
  }
}

ResFrm$padj <- p.adjust(ResFrm$pvalue, method = "bonferroni")

# Save results to Excel
#output_file <- paste("Data/Generated/Permutations/Network_Associations/PrmFrm_",
#                     ArcArdAss, ScrOpc, "Dst_", NumPer, ".xlsx", sep = "")
output_file <- paste("Data/Generated/Permutations/Network_Associations/PermutationTable_",
                     toupper(ArcArdAss),'_', ScrOpcTxt,"_", NumPer, ".xlsx", sep = "")
write.xlsx(ResFrm, output_file)

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create composite plot for a specific network
#' @param results_list Results list for a specific network
#' @param network_name Name of the network for title
#' @param score_type Type of score for subtitle
#' @return Combined plot object
create_network_summary_plot <- function(results_list, network_name, score_type) {
  
  # Map score type to title text
  score_titles <- list(
    "Int" = "Mean number of interactions with ARCs",
    "Prx" = "Mean Proximity to ARCs",
    "Dst" = "Mean Distance to ARCs",
    "Wlk" = "Mean RWR score to ARCs",
    "Plt" = "Mean Pleiotropy to ARCs"
  )
  
  TitTxt <- score_titles[[score_type]]
  
  # Extract individual plots
  plots <- list(
    results_list$Unbalanced$AgeGen$plot,
    results_list$Unbalanced$ModGen$plot,
    results_list$Unbalanced$DisGen$plot,
    results_list$Unbalanced$HghGen$plot
  )
  
  # Combine plots in a grid
  p <- plot_grid(plotlist = plots, 
                 labels = c("a", "b", "c", "d"), 
                 label_size = 11, 
                 ncol = 2)
  
  # Add main title
  title <- ggdraw() + 
    draw_label(
      paste("Permutation results for ", network_name, ": ", TitTxt, sep = ""),
      fontface = "bold", 
      size = 11,
      x = 0.5, 
      hjust = 0.5
    )
  
  # Combine title and plots
  final_plot <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  
  return(final_plot)
}

# Generate summary plot for KEGG network (example)
kegg_plot <- create_network_summary_plot(results$KEGG, "KEGG", ScrOpc)
print(kegg_plot)

# ==============================================================================
# COMPREHENSIVE SUMMARY GENERATION
# ==============================================================================

# Generate summary tables for all networks
summary_tables <- list()
for (network in NtwArr) {
  summary_tables[[network]] <- create_comprehensive_summary(results[[network]]) %>%
    mutate(Network = network)
}

# Combine all summaries
combined_summary <- bind_rows(summary_tables)

# Adjuste pvalues
combined_summary$P_adj = p.adjust(combined_summary$P_Value, method = "bonferroni")

# Save combined summary

summary_file <- paste("Data/Generated/Permutations/Network_Associations/PermutationResultsSummary_",
                      toupper(ArcArdAss),"_", ScrOpcTxt, "_", NumPer, ".csv", sep = "")
#summary_file <- paste("Data/Generated/Permutations/Network_Associations/PrmResSum_",
#                      ArcArdAss, ScrOpc, "Dst_", NumPer, ".csv", sep = "")
write.csv(combined_summary, summary_file, row.names = FALSE)

# Print summary
print(combined_summary)





