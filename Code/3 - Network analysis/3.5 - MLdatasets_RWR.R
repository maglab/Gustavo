# ==============================================================================
# RANDOM WALK WITH RESTART (RWR) ANALYSIS FOR GENE NETWORKS
# Unified analysis for both ARC and ARD levels - VERSION CORREGIDA
# ==============================================================================

# ==============================================================================
# CONFIGURATION 
# ==============================================================================

# Set GRIN environment (conda libraries)
env_lib <- "C:/Users/Usuario/miniconda3/envs/GRIN/Lib/R/library"
if (!dir.exists(env_lib)) dir.create(env_lib, recursive = TRUE)
.libPaths(c(env_lib, .libPaths()))

# Load required packages
suppressPackageStartupMessages({
  library(igraph)
  library(tidyr)
  library(dplyr)
  library(magrittr)
  library(tibble)
  library(readr)
  library(RandomWalkRestartMH)
})

# ==============================================================================
# USER INPUT: CHOOSE ANALYSIS LEVEL
# ==============================================================================

# Set this variable to either "ARD" or "ARC" to choose analysis level
analysis_level <- "ARC"  # Change to "ARD" for ARD analysis

# Validate user input
if (!analysis_level %in% c("ARD", "ARC")) {
  stop("analysis_level must be either 'ARD' or 'ARC'")
}

message("=========================================================")
message("Starting RWR analysis at ", analysis_level, " level")
message("=========================================================")

# ==============================================================================
# INPUT PATHS 
# ==============================================================================

base_in <- "Data/Generated/Networks_and_predictions/Networks"

# Network graph lists (RDS files containing igraph objects)
input_paths <- list(
  PPI   = file.path(base_in, "PPI",   "Lists", "GraphsList.rds"),
  COX90 = file.path(base_in, "COX90", "Lists", "GraphsList.rds"),
  COX95 = file.path(base_in, "COX95", "Lists", "GraphsList.rds"),
  KEGG  = file.path(base_in, "KEGG",  "Lists", "GraphsList.rds")
)

# ARD/ARC gene annotations
ard_arc_rds <- file.path("Data/Retrieved/Genes_and_diseases/Diseases", "GeneFrame_ARD_ARC.rds")

GenAgeHum = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")



# ==============================================================================
# OUTPUT SETUP 
# ==============================================================================

# Create output directory with analysis level
#output_base <- file.path("D:/Respaldo_PHD/Nature_Data/Data/Generated/Reviews/RWR", analysis_level)
output_base <- file.path("Data/Generated/Networks_and_predictions/Topological_data/RWR", analysis_level)
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# RWR PARAMETERS 
# ==============================================================================

restart_probability <- 0.7  # Restart probability for RWR

# Tau parameters for weighting network layers
tau_single <- 1
tau_multiplex <- c(1, 1, 1, 1)  # Equal weights for all 4 networks

# ==============================================================================
# UTILITY FUNCTIONS - VERSION CORREGIDA
# ==============================================================================

#' Coerce object to igraph if possible
#' @param x Object to coerce
#' @return igraph object
coerce_to_igraph <- function(x) {
  if (inherits(x, "igraph")) return(x)
  stop("Loaded object is not an 'igraph'. Check 'GenGenNtw' or your reading method.")
}

#' Build MultiplexObject for a single network layer
#' @param g igraph object
#' @param layer_name Name for the network layer
#' @param delta Inter-layer coupling parameter
#' @return List with mpo (MultiplexObject) and adjnorm (normalized adjacency matrix)
build_mpo_single <- function(g, layer_name = "L1", delta = 0.5) {
  # 1) Basic validation
  if (!igraph::is_igraph(g)) {
    stop("Expected an 'igraph' object in 'g'.")
  }
  if (is.null(igraph::V(g)$name)) {
    stop("Graph lacks V(g)$name attribute with gene IDs.")
  }
  
  # Ensure unique character vertex names
  igraph::V(g)$name <- as.character(igraph::V(g)$name)
  if (anyDuplicated(igraph::V(g)$name)) {
    stop("Duplicate vertex names in V(g)$name.")
  }
  
  # 2) Build single-layer multiplex
  layers <- list()
  layers[[layer_name]] <- g
  mpo <- create.multiplex(LayersList = layers)
  
  # 3) Compute and normalize supra-adjacency matrix
  adj <- compute.adjacency.matrix(mpo, delta = delta)
  adjnorm <- normalize.multiplex.adjacency(adj)
  
  list(mpo = mpo, adjnorm = adjnorm)
}

#' Build MultiplexObject from multiple network layers
#' @param graph_list_named Named list of igraph objects
#' @param delta Inter-layer coupling parameter
#' @return List with mpo (MultiplexObject) and adjnorm (normalized adjacency matrix)
build_mpo_multiplex <- function(graph_list_named, delta = 0.5) {
  # Get common node universe across all layers
  all_nodes <- sort(unique(unlist(lapply(graph_list_named, 
                                         function(g) igraph::V(g)$name))))
  
  # Ensure each layer contains all nodes (add isolated vertices if missing)
  graph_list_padded <- lapply(graph_list_named, function(g) {
    g <- coerce_to_igraph(g)
    missing_nodes <- setdiff(all_nodes, igraph::V(g)$name)
    
    if (length(missing_nodes) > 0) {
      g <- igraph::add_vertices(g, length(missing_nodes), name = missing_nodes)
    }
    
    # Reorder vertices alphabetically
    g <- igraph::induced_subgraph(g, vids = all_nodes)
    g
  })
  
  # Create multiplex object
  mpo <- create.multiplex(LayersList = graph_list_padded)
  
  # Compute and normalize supra-adjacency matrix
  adj <- compute.adjacency.matrix(mpo, delta = delta)
  adjnorm <- normalize.multiplex.adjacency(adj)
  
  list(mpo = mpo, adjnorm = adjnorm)
}

#' Run RWR for a single Trait (ARD or ARC) - VERSION CORREGIDA
#' @param mpo MultiplexObject
#' @param adjnorm Normalized adjacency matrix
#' @param seeds Vector of seed genes
#' @param restart Restart probability
#' @param tau Layer weights
#' @param Trait_name Name of the Trait (ARD or ARC)
#' @param level Analysis level ("ARD" or "ARC")
#' @return Tibble with RWR results for all genes
run_rwr_one_Trait <- function(mpo, adjnorm, seeds, restart = 0.7, 
                                  tau = 1, Trait_name = NA_character_, 
                                  level = "ARD") {
  node_pool <- mpo$Pool_of_Nodes
  
  # Find seeds present in the network
  seeds_in_network <- intersect(seeds, node_pool)
  
  if (length(seeds_in_network) < 1) {
    warning(sprintf("%s '%s': 0 seeds in network; skipping.", level, Trait_name))
    
    # Devolver tibble con la columna correcta según el nivel
    if (level == "ARD") {
      return(tibble(
        Ard = character(0), 
        Gene = character(0), 
        Score = numeric(0), 
        Rank = numeric(0)
      ))
    } else {
      return(tibble(
        Arc = character(0), 
        Gene = character(0), 
        Score = numeric(0), 
        Rank = numeric(0)
      ))
    }
  }
  
  # Run RWR
  rwr_result <- Random.Walk.Restart.Multiplex(
    x = adjnorm, 
    MultiplexObject = mpo,
    Seeds = seeds_in_network, 
    r = restart, 
    tau = tau, 
    weights = 1
  )
  
  # Process results - CORREGIDO: usar la columna correcta directamente
  result_table <- rwr_result$RWRM_Results %>%
    as_tibble() %>%
    transmute(Gene = NodeNames, Score = Score)
  
  # Fill missing genes with score 0 - CORREGIDO: crear columna correcta directamente
  if (level == "ARD") {
    full_result <- tibble(Gene = node_pool) %>%
      left_join(result_table, by = "Gene") %>%
      mutate(
        Score = coalesce(Score, 0),
        Ard = Trait_name,
        Rank = min_rank(desc(Score))
      ) %>%
      relocate(Ard, Gene, Score, Rank)
  } else {
    full_result <- tibble(Gene = node_pool) %>%
      left_join(result_table, by = "Gene") %>%
      mutate(
        Score = coalesce(Score, 0),
        Arc = Trait_name,
        Rank = min_rank(desc(Score))
      ) %>%
      relocate(Arc, Gene, Score, Rank)
  }
  
  full_result
}

#' Run RWR for all communities (ARDs or ARCs)
#' @param mpo MultiplexObject
#' @param adjnorm Normalized adjacency matrix
#' @param Trait_df Tibble with columns 'Gen' (gene) and Trait column
#' @param restart Restart probability
#' @param tau Layer weights
#' @param level Analysis level ("ARD" or "ARC")
#' @return Long table with RWR results for all communities
run_rwr_all_communities <- function(mpo, adjnorm, Trait_df, restart = 0.7, 
                                    tau = 1, level = "ARD") {
  # Determine Trait column name
  comm_col <- ifelse(level == "ARD", "Ard", "Arc")
  
  # Get unique communities
  communities <- sort(unique(Trait_df[[comm_col]]))
  results_list <- vector("list", length(communities))
  names(results_list) <- communities
  
  for (i in seq_along(communities)) {
    comm <- communities[i]
    seeds <- Trait_df %>% 
      filter(.data[[comm_col]] == comm) %>% 
      pull(Gen) %>% 
      unique()
    
    results_list[[i]] <- run_rwr_one_Trait(
      mpo, adjnorm, seeds, 
      restart = restart, 
      tau = tau, 
      Trait_name = comm,
      level = level
    )
  }
  
  bind_rows(results_list)
}

# ==============================================================================
# DATA LOADING 
# ==============================================================================

# 1) Load network graphs
graph_list <- list()

for (network_name in names(input_paths)) {
  graph_data <- readRDS(input_paths[[network_name]])
  # Assuming object has $GenGenNtw which is an 'igraph'
  graph_list[[network_name]] <- graph_data$GenGenNtw
}

# 2) Load ARD/ARC annotations and select appropriate column
gene_Trait_df <- readRDS(ard_arc_rds) %>%
  rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)

# Select appropriate column based on analysis level
if (analysis_level == "ARD") {
  gene_Trait_df <- gene_Trait_df %>%
    select(Gen, Ard) %>%
    distinct()
} else {
  gene_Trait_df <- gene_Trait_df %>%
    select(Gen, Arc) %>%
    distinct()
}

# ==============================================================================
# MONOPLEX ANALYSIS - VERSION CORREGIDA
# ==============================================================================

message("\nStarting monoplex (individual network) analysis...")

# Initialize storage for monoplex results
monoplex_results <- tibble()

for (network_name in names(graph_list)) {
  message(">>> Processing monoplex network: ", network_name)
  
  # Create output directory for this network
  network_output_dir <- file.path(output_base, paste0(network_name))
  dir.create(network_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Build multiplex object for single network
  mpo_single <- build_mpo_single(
    graph_list[[network_name]], 
    layer_name = network_name, 
    delta = 0.5
  )
  
  # Run RWR for all communities
  Trait_associations <- run_rwr_all_communities(
    mpo_single$mpo, 
    mpo_single$adjnorm, 
    gene_Trait_df, 
    restart = restart_probability, 
    tau = tau_single,
    level = analysis_level
  )
  
  # Calculate Z-scores within each Trait - CORREGIDO
  Trait_col <- ifelse(analysis_level == "ARD", "Ard", "Arc")
  
  Trait_associations <- Trait_associations %>%
    group_by(.data[[Trait_col]]) %>%
    mutate(Score_Z = (Score - mean(Score)) / sd(Score)) %>%
    ungroup()
  
  # Save detailed Trait associations
  output_file <- file.path(
    network_output_dir, 
    paste0("RWR_", analysis_level, "_Association",".csv")
    #paste0("RWR_", analysis_level, "Association_", network_name, ".csv")
  )
  write.csv(Trait_associations, output_file, row.names = FALSE)

  
  # Convert to a wide dataset, apply the GenAge_Hum label and save
  wide_score <- to_wide_by_firstcol(Trait_associations, gene_col = "Gene", value_col = "Score")
  wide_score$AgeGen = ifelse(row.names(wide_score) %in% GenAgeHum, "Age", "NotAge")
  write.csv(wide_score, paste("Data/Generated/Networks_and_predictions/Networks/",network_name,
                              "/Ageing_Prediction/Datasets/RWR.Proximity2",analysis_level,".csv",sep=""))


  # Calculate gene-level summary statistics
  gene_summary <- Trait_associations %>%
    group_by(Gene) %>%
    summarise(
      TopTrait = .data[[Trait_col]][which.max(Score)],
      TopScore = max(Score), 
      .groups = "drop"
    )
  
  
  
  
  # Save gene summary
  summary_file <- file.path(
    network_output_dir, 
    #paste0("RWR_", analysis_level, "GeneSummary_", network_name, ".csv")
    paste0("RWR_", analysis_level, "_GeneSummary", ".csv")
  )
  write.csv(gene_summary, summary_file, row.names = FALSE)
  
  # Calculate mean statistics per gene across communities
  gene_mean_stats <- Trait_associations %>%
    group_by(Gene) %>%
    summarise(
      Mean_Score_to_Trait = mean(Score, na.rm = TRUE),
      Mean_Rank_to_Trait = mean(Rank, na.rm = TRUE),
      Mean_Zscore_to_Trait = mean(Score_Z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Score_to_Trait)) %>%
    mutate(Network = network_name)
  
  # Accumulate results
  monoplex_results <- bind_rows(monoplex_results, gene_mean_stats)
}

# Save combined monoplex results
monoplex_output_file <- file.path(output_base, "Summaries_and_Multiplex", paste0("RWR_", analysis_level, "_MeanAssociation_Monoplex.csv"))
write.csv(monoplex_results, monoplex_output_file, row.names = FALSE)

# ==============================================================================
# MULTIPLEX ANALYSIS 
# ==============================================================================

message("\n>>> Building multiplex network with PPI, COX90, COX95, KEGG...")

# Create output directory for multiplex results
multiplex_output_dir <- file.path(output_base, "Summaries_and_Multiplex")
dir.create(multiplex_output_dir, showWarnings = FALSE, recursive = TRUE)

# Build multiplex object from all networks
mpo_multiplex <- build_mpo_multiplex(graph_list)

# Run RWR for all communities on multiplex
multiplex_associations <- run_rwr_all_communities(
  mpo_multiplex$mpo, 
  mpo_multiplex$adjnorm, 
  gene_Trait_df,
  restart = restart_probability, 
  tau = tau_multiplex,
  level = analysis_level
)

# Calculate Z-scores - CORREGIDO
Trait_col <- ifelse(analysis_level == "ARD", "Ard", "Arc")

multiplex_associations <- multiplex_associations %>%
  group_by(.data[[Trait_col]]) %>%
  mutate(Score_Z = (Score - mean(Score)) / sd(Score)) %>%
  ungroup()

# Save detailed multiplex associations
multiplex_assoc_file <- file.path(multiplex_output_dir, paste0("RWR_", analysis_level, "_Association_Multiplex.csv"))
write.csv(multiplex_associations, multiplex_assoc_file, row.names = FALSE)

# Calculate gene-level summary for multiplex
gene_summary_multiplex <- multiplex_associations %>%
  group_by(Gene) %>%
  summarise(
    TopTrait = .data[[Trait_col]][which.max(Score)],
    TopScore = max(Score), 
    .groups = "drop"
  )

multiplex_summary_file <- file.path(multiplex_output_dir, paste0("RWR_", analysis_level, "_GeneSummary_Multiplex.csv"))
write.csv(gene_summary_multiplex, multiplex_summary_file, row.names = FALSE)

# Calculate mean statistics per gene for multiplex
multiplex_mean_stats <- multiplex_associations %>%
  group_by(Gene) %>%
  summarise(
    Mean_Score_to_Trait = mean(Score, na.rm = TRUE),
    Mean_Rank_to_Trait = mean(Rank, na.rm = TRUE),
    Mean_Zscore_to_Trait = mean(Score_Z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Score_to_Trait)) %>%
  mutate(Network = "MULTIPLEX")

multiplex_mean_file <- file.path(multiplex_output_dir, paste0("RWR_", analysis_level, "_MeanAssociation_Multiplex.csv"))
write.csv(multiplex_mean_stats, multiplex_mean_file, row.names = FALSE)


# ==============================================================================
# SAVE MULTIPLEX DATASETS
# ==============================================================================

GenAge = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")

multiplex_mean_stats_arc = read.csv("Data/Generated/Networks_and_predictions/Topological_data/RWR/ARC/Summaries_and_Multiplex/RWR_ARC_Association_Multiplex.csv")
multiplex_mean_stats_ard = read.csv("Data/Generated/Networks_and_predictions/Topological_data/RWR/ARD/Summaries_and_Multiplex/RWR_ARD_Association_Multiplex.csv")

# Convert to matrix-like mode
ArcScoreMat <- multiplex_mean_stats_arc %>%
  dplyr::select(Gene, Arc, Score) %>%
  pivot_wider(
    names_from  = Arc,
    values_from = Score
  ) %>%
  column_to_rownames("Gene")

ArdScoreMat <- multiplex_mean_stats_ard %>%
  dplyr::select(Gene, Ard, Score) %>%
  pivot_wider(
    names_from  = Ard,
    values_from = Score
  ) %>%
  column_to_rownames("Gene")


ArcScoreMat$AgeGen = ifelse(row.names(ArcScoreMat) %in% GenAge, "Age", "NotAge")
ArdScoreMat$AgeGen = ifelse(row.names(ArdScoreMat) %in% GenAge, "Age", "NotAge")

write.csv(ArcScoreMat,"Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Datasets/RWR.Proximity2ARC.csv")
write.csv(ArdScoreMat,"Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Datasets/RWR.Proximity2ARD.csv")


# ==============================================================================
# COMBINE RESULTS 
# ==============================================================================

# Combine monoplex and multiplex results
combined_results <- bind_rows(monoplex_results, multiplex_mean_stats)

combined_output_file <- file.path(output_base, "Summaries_and_Multiplex", paste0("RWR_", analysis_level, "_MeanAssociation_MonoMulti.csv"))
write.csv(combined_results, combined_output_file, row.names = FALSE)

message("\n=========================================================")
message(analysis_level, " analysis complete!")
message("Results saved to: ", output_base)
message("=========================================================")


# ==============================================================================
# FUNCTIONS
# ==============================================================================

to_wide_by_firstcol <- function(df, gene_col = "Gene", value_col = "Score",
                                fill = NA_real_, agg = NULL) {
  group_col <- names(df)[1]  # <-- primera columna (Arc u otro nombre)
  
  # Si hay duplicados por (grupo, gen), define cómo agregarlos:
  # - agg = min / max / mean / etc.
  if (!is.null(agg)) {
    df <- df %>%
      group_by(across(all_of(c(group_col, gene_col)))) %>%
      summarise(.value = agg(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    value_col2 <- ".value"
  } else {
    value_col2 <- value_col
  }
  
  out <- df %>%
    select(all_of(c(group_col, gene_col, value_col2))) %>%
    pivot_wider(
      names_from  = all_of(group_col),
      values_from = all_of(value_col2),
      values_fill = fill
    ) %>%
    column_to_rownames(gene_col)
  
  out
}
