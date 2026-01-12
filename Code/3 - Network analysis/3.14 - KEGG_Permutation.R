# Load required libraries
library(dplyr)
library(igraph)
library(cowplot)
library(ggplot2)
library(ggtext)
library(openxlsx)

remotes::install_version("markdown", version = "1.12", repos = "https://cran.rstudio.com")
library("markdown")

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load KEGG interaction network
GenGenKeg = readRDS("Data/Retrieved/Network_Sources/KEGG/Other/KEGG_Interactions.rds")

# Load gene sets
GenAge.Hum <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
GenAge.Mod <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")

# Load ARC/ARD gene annotations
GenArdArcFrm <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>%
  dplyr::rename(Cod = Code, Ard = ARD_Meaning, Arc = ARC_Meaning, Gen = Gene)

# Define gene groups
High_Pleiotropy = Get_high_pleiotropy_genes(GenArdArcFrm)

Diseases = GenArdArcFrm %>% 
  pull(Gen) %>% 
  unique()

# ==============================================================================
# COMPUTE KEGG METRICS
# ==============================================================================

head(GenGenKeg)
# Calculate KEGG metrics (distance to leaves and roots)
kegg_metrics_df <- compute_kegg_metrics_no_scc(
  GenGenKeg, 
  from_col = "Gen_from", 
  to_col = "Gen_to"
)

# ==============================================================================
# RUN PERMUTATION TESTS
# ==============================================================================

# Define gene groups for permutation testing
groups <- list(
  GenAge.Hum = GenAge.Hum,
  GenAge.Mod = GenAge.Mod,
  High_Pleiotropy = High_Pleiotropy,
  Diseases = Diseases
)

# Run permutation test suite
suite <- run_perm_suite(
  metrics_df = kegg_metrics_df,
  groups_named_list = groups,
  metrics = c("leaf_mean", "root_mean"),
  n_perm = 10000,
  seed = 42,
  out_dir = "perm_leaf_tests",
  save_plots = TRUE,
  keep_plots = TRUE,
  keep_results = TRUE,
  progress_each = TRUE
)

# ==============================================================================
# CREATE FORMATTED PLOTS WITH ADJUSTED P-VALUES
# ==============================================================================

# Create custom plots with adjusted p-values
create_formatted_plots <- function(suite, distance_type = "root") {
  # Get group names
  group_names <- names(suite$plots)
  plots_list <- list()
  
  for(group_name in group_names) {
    # Extract plot
    p <- suite$plots[[group_name]][[distance_type]]$mean
    
    # Get adjusted p-value from summary table
    p_val_adj <- suite$summary_table %>% 
      filter(metric == paste0(distance_type, "_mean"), 
             group == group_name) %>% 
      pull(p_adjust_BH) %>% 
      format_mixto()
    
    # Extract statistics from plot subtitle
    sub_txt <- p$labels$subtitle
    
    # Extract sample size
    n <- sub(".*n = ([0-9\\.eE-]+).*", "\\1", sub_txt)
    
    # Extract null mean
    null_mean <- sub(".*Null mean: ([0-9\\.eE-]+).*", "\\1", sub_txt)
    
    # Extract observed value
    obs_value <- sub(".*d62728;'>[^:]+:\\s*([0-9eE\\.-]+).*", "\\1", sub_txt)
    
    # Create formatted plot with adjusted p-value
    formatted_plot <- p +
      labs(subtitle = NULL) +  # Remove original subtitle
      labs(
        subtitle = paste0(
          "Padj: ", p_val_adj, " | n = ", n, "<br>",
          "<span style='color:blue;'>Null mean: ", null_mean, "</span>",
          " | ",
          "<span style='color:red;'>", group_name, ": ", obs_value, "</span>"
        )
      )
    
    plots_list[[group_name]] <- formatted_plot
  }
  
  return(plots_list)
}

# Create plots for both distance types
leaf_plots <- create_formatted_plots(suite, "leaf")
root_plots <- create_formatted_plots(suite, "root")

# ==============================================================================
# CREATE FINAL GRID PLOT
# ==============================================================================

# Rotation functions for better visualization
plot_leaf_rotated <- function(p, 
                              xlab = "Mean distance to\nKEGG leaves",
                              ylab = "Frequency under null") {
  p +
    coord_flip() +
    scale_y_reverse(expand = expansion(mult = c(0.12, 0))) +
    scale_x_continuous(position = "top") +
    labs(x = xlab, y = ylab) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold")
    )
}

plot_root_rotated <- function(p, 
                              xlab = "Frequency under null",
                              ylab = "Mean distance to\nKEGG roots") {
  p +
    coord_flip() +
    scale_x_reverse() +
    labs(x = ylab, y = xlab) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold", angle = 90),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    )
}

# Create list of all plots (2 columns: root, leaf for each group)
#plots_list <- list(
#  # Column 1: Root distances
#  suite$plots$GenAge.Hum$root$mean,
#  suite$plots$GenAge.Mod$root$mean,
#  suite$plots$High_Pleiotropy$root$mean,
#  suite$plots$Diseases$root$mean,
#  # Column 2: Leaf distances
#  suite$plots$GenAge.Hum$leaf$mean,
#  suite$plots$GenAge.Mod$leaf$mean,
#  suite$plots$High_Pleiotropy$leaf$mean,
#  suite$plots$Diseases$leaf$mean
#)

# Create list of all plots (2 columns: root, leaf for each group)
plots_list <- list(

  suite$plots$GenAge.Hum$root$mean,
  suite$plots$GenAge.Hum$leaf$mean,
  
  suite$plots$GenAge.Mod$root$mean,
  suite$plots$GenAge.Mod$leaf$mean,
  
  suite$plots$High_Pleiotropy$root$mean,
  suite$plots$High_Pleiotropy$leaf$mean,
  
  suite$plots$Diseases$root$mean,
  suite$plots$Diseases$leaf$mean

)

# Create grid of plots
grid_plots <- plot_grid(
  plotlist = plots_list,
  ncol = 2,
  align = "hv",
  labels = letters[1:length(plots_list)],
  label_size = 11
)

# Create y-axis label
y_label <- ggdraw() +
  draw_label("Count", angle = 90, size = 10, hjust = 0.5, vjust = 1)

# Combine y-label with plots
main_with_label <- plot_grid(
  y_label, grid_plots, ncol = 2, rel_widths = c(0.05, 1)
)

# Create title
title <- ggdraw() + 
  draw_label(
    "Permutation tests: KEGG distances to roots and leaves",
    fontface = "bold", size = 12
  )

# Create caption
caption <- ggdraw() +
  draw_label("Mean distance to roots and leaves", 
             size = 10, hjust = 0.5, vjust = 0)

# Assemble final plot
final_plot <- plot_grid(
  title,
  main_with_label,
  caption,
  ncol = 1,
  rel_heights = c(0.08, 2, 0.06)
)

# Display final plot
#print(final_plot)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save plot
saveRDS(final_plot, "Data/Generated/Permutations/KEGG_Hierarchy/KEGG_Distance_PermutationPlot.rds")

# Save results
write.xlsx(
  suite$summary_table,
  file = "Data/Generated/Permutations/KEGG_Hierarchy/KEGG_Distance_PermutationTable.xlsx"
)


# ==============================================================================
# PRINTING 
# ==============================================================================


# --- IMPROVED DESIGN PARAMETERS ---
row_names_width <- 0.08   # Width for row name labels
counts_width    <- 0.025  # Width for "counts" label
left_total_width <- row_names_width + counts_width
row_gap_height  <- 0.08   # Gap between rows
column_gap_width <- 0.22  # Gap between columns (as fraction of total width)

# Global margins (tight for final composition)
outer_margin <- margin(t = 5, r = 15, b = 0, l = 5, unit = "pt")

# Helper function to remove plot margins
tighten_margins <- function(plot_obj) {
  plot_obj + theme(plot.margin = margin(0, 0, 0, 0))
}

# Generate labels for subplots (a, b, c, ...)
labels_vector <- letters[1:length(plots_list)]

# ========= CREATE MAIN GRID WITH COLUMN GAPS =====
create_row_with_gap <- function(left_plot, right_plot,
                                left_label, right_label,
                                gap_proportion = column_gap_width) {
  cowplot::plot_grid(
    left_plot,
    NULL,  # Empty space for gap
    right_plot,
    ncol = 3,
    rel_widths = c(1, gap_proportion, 1),
    align = "hv",
    labels = c(left_label, "", right_label),
    label_size = 10,
    label_fontface = "bold",
    label_x = 0.02,
    label_y = 0.98
  )
}

# Apply to each row (assuming plots_list has even number of plots)
rows_list <- lapply(seq(1, length(plots_list), by = 2), function(i) {
  create_row_with_gap(
    tighten_margins(plots_list[[i]]),
    tighten_margins(plots_list[[i + 1]]),
    left_label = labels_vector[i],
    right_label = labels_vector[i + 1]
  )
})

# Create empty spacer for row separation
empty_spacer <- ggdraw() + theme_void()

# Combine rows with gaps
grid_plots <- plot_grid(
  rows_list[[1]], empty_spacer,
  rows_list[[2]], empty_spacer,
  rows_list[[3]], empty_spacer,
  rows_list[[4]],
  ncol = 1,
  rel_heights = c(1, row_gap_height, 1, row_gap_height, 1, row_gap_height, 1)
)

# ========= CREATE LEFT SIDE LABELS (Row names + counts) =====
left_labels <- plot_grid(
  ggdraw() + draw_label("GenAge.Hum", angle = 90, x = 0.48, y = 0.5, size = 10, fontface = "bold"),
  ggdraw() + draw_label("GenAge.Mod", angle = 90, x = 0.48, y = 0.5, size = 10, fontface = "bold"),
  ggdraw() + draw_label("High ARC-Pleiotropy", angle = 90, x = 0.48, y = 0.5, size = 10, lineheight = 0.95, fontface = "bold"),
  ggdraw() + draw_label("Diseases", angle = 90, x = 0.48, y = 0.5, size = 10, fontface = "bold"),
  ncol = 1,
  rel_heights = c(1, 1, 1, 1)
)

# "counts" label
counts_column <- ggdraw() +
  draw_label("counts", angle = 90, x = 0.02, y = 0.5, hjust = 0, vjust = 0.5, size = 10)

# Combine left side labels
left_strip <- plot_grid(
  left_labels,
  counts_column,
  ncol = 2,
  rel_widths = c(row_names_width, counts_width)
)

# Combine left strip with main plots
middle_section <- plot_grid(
  left_strip,
  grid_plots,
  ncol = 2,
  rel_widths = c(left_total_width, 1)
)

# ======== CREATE BOTTOM CAPTIONS =====
bottom_captions <- plot_grid(
  ggdraw(),
  ggdraw() + draw_label("Number of edges to roots", x = 0.6, y = 0.9, hjust = 0.5, vjust = 1, size = 10),
  ggdraw() + draw_label("Number of edges to leafs", x = 0.6, y = 0.9, hjust = 0.5, vjust = 1, size = 10),
  ncol = 3,
  rel_widths = c(left_total_width, 1, 1)
)

# ======== CREATE COLUMN HEADERS =====
column_headers <- plot_grid(
  ggdraw(),
  ggdraw() + draw_label("Distance to KEGG roots\n(smaller = more upstream/regulatory)", 
                        x = 0.6, y = 0.5, hjust = 0.5, vjust = 0.5, size = 10, fontface = "bold"),
  NULL,  # Space for column gap
  ggdraw() + draw_label("Distance to KEGG leafs\n(smaller = more downstream/terminal)", 
                        x = 0.6, y = 0.5, hjust = 0.5, vjust = 0.5, size = 10, fontface = "bold"),
  ncol = 4,
  rel_widths = c(left_total_width, 1, column_gap_width, 1)
)

# ======== ASSEMBLE FINAL PLOT =====
final_plot <- plot_grid(
  column_headers,
  middle_section,
  bottom_captions,
  ncol = 1,
  rel_heights = c(0.07, 1, 0.06)
)

# Apply global margins
final_plot <- ggdraw(final_plot) + 
  theme(plot.margin = outer_margin)

# Display final plot
print(final_plot)

# ======== SAVE PLOTS IN MULTIPLE FORMATS =====

# Save as SVG
ggsave(
  filename = "Last_Figures/Supp_Kegg_RootLeaf_8.svg",
  plot = final_plot,
  width = 8.5,
  height = 6.5,
  units = "in",
  device = "svg",
  limitsize = FALSE
)


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Format numbers in scientific notation with 2 decimal places
#' @param x Numeric vector to format
#' @return Character vector with formatted numbers
format_mixto <- function(x) {
  formatC(x, format = "e", digits = 2)
}

#' Create permutation histogram plot
#' @param null_stats Null distribution statistics
#' @param obs_stat Observed statistic
#' @param null_mean Mean of null distribution
#' @param p_two_tailed Two-tailed p-value
#' @param n_group Sample size in group
#' @param title Plot title
#' @param xlab X-axis label
#' @param bins Number of histogram bins
#' @return ggplot object
plot_perm_hist <- function(null_stats, obs_stat, null_mean,
                           p_two_tailed, n_group,
                           title, xlab = "Mean under null",
                           bins = 30, n_pvalues=8) {
  
  # Filter finite values
  ns <- null_stats[is.finite(null_stats)]
  df <- data.frame(null = ns)
  
  #pretty_num <- function(x) {
  #  x <- as.numeric(x)
  #  ifelse(x == 1, "1", format(x, scientific = TRUE))
  #}
  
  pretty_num <- function(x) {
    x <- as.numeric(x)
    ifelse(
      x == 1,
      "1",
      sprintf("%.2e", x)
    )
  }
  
  padj = p_bonferroni(p_two_tailed,n_pvalues) %>% pretty_num()
  
  # Create formatted subtitle
  #subtitle <- sprintf(
  #  "P-value (two-sided): %s   |   n = %d<br><span style='color:#1f77b4;'>Null mean: %.3f</span>   |   <span style='color:#d62728;'>Observed mean: %.3f</span>",
  #  formatC(p_two_tailed, format = "e", digits = 2),
  #  n_group, null_mean, obs_stat
  #)
  
  subtitle <- sprintf(
    "P_adj: %s   |   n = %d<br><span style='color:#1f77b4;'>Null mean: %.3f</span>   |   <span style='color:#d62728;'>Observed mean: %.3f</span>",
    formatC(padj, format = "e", digits = 2),
    n_group, null_mean, obs_stat
  )
  
  # Create plot
  ggplot(df, aes(x = null)) +
    geom_histogram(bins = bins, fill = "#cfe8f3", color = "#cfe8f3") +
    geom_vline(xintercept = null_mean, linetype = "dashed", linewidth = 1, color = "#1f77b4") +
    geom_vline(xintercept = obs_stat,  linetype = "dashed", linewidth = 1, color = "#d62728") +
    labs(
      subtitle = subtitle,
      x = xlab,
      y = "Frequency"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_bw(base_size = 12) +
    theme(
      plot.subtitle = element_markdown(size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

#' Run permutation test for one group and metric
#' @param metrics_df Data frame with gene metrics
#' @param group_genes Vector of gene names in group
#' @param metric_col Name of metric column to test
#' @param n_perm Number of permutations
#' @param seed Random seed
#' @param return_null Whether to return null distribution
#' @param gname Group name for progress reporting
#' @return List with summary statistics and null distribution
perm_test_one <- function(metrics_df, group_genes, metric_col,
                          n_perm = 10000, seed = 42,
                          return_null = TRUE, gname) {
  stopifnot(metric_col %in% names(metrics_df))
  if (!is.null(seed)) set.seed(seed)
  
  # Prepare population data
  dt <- metrics_df[is.finite(metrics_df[[metric_col]]),
                   c("gene", metric_col)]
  pop_genes <- dt$gene
  pop_vals  <- setNames(dt[[metric_col]], dt$gene)
  
  # Intersection with group
  g_in <- intersect(group_genes, pop_genes)
  if (length(g_in) == 0L) {
    warning(sprintf("No intersection for group %s with metric %s", gname, metric_col))
    return(NULL)
  }
  k <- length(g_in)
  
  # Observed statistic
  obs <- mean(pop_vals[g_in])
  
  # Permutation test
  null_stats <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    samp <- sample(pop_genes, k, replace = FALSE)
    null_stats[i] <- mean(pop_vals[samp])
    
    # Progress reporting
    if (i %% 1000 == 0) {
      cat(gname, " - Iteration", i, "of", n_perm, "\n")
      flush.console()
    }
  }
  
  # Calculate statistics
  null_mean <- mean(null_stats)
  null_sd   <- sd(null_stats)
  r_ge <- sum(null_stats >= obs)
  r_le <- sum(null_stats <= obs)
  p_ge <- (r_ge + 1) / (n_perm + 1)
  p_le <- (r_le + 1) / (n_perm + 1)
  p_two <- min(1, 2 * min(p_ge, p_le))
  z <- if (null_sd > 0) (obs - null_mean) / null_sd else NA_real_
  
  # Prepare output
  out <- list(
    summary = data.frame(
      metric       = metric_col,
      n_group      = k,
      n_perm       = n_perm,
      obs_mean     = obs,
      null_mean    = null_mean,
      effect_size  = obs - null_mean,
      null_sd      = null_sd,
      z            = z,
      p_greater    = p_ge,
      p_less       = p_le,
      p_two_tailed = p_two,
      stringsAsFactors = FALSE
    )
  )
  if (return_null) out$null_stats <- null_stats
  out
}

#' Run permutation test suite for multiple groups and metrics
#' @param metrics_df Data frame with gene metrics
#' @param groups_named_list Named list of gene groups
#' @param metrics Vector of metric names to test
#' @param n_perm Number of permutations
#' @param seed Random seed
#' @param out_dir Output directory
#' @param save_plots Whether to save plots
#' @param keep_plots Whether to keep plots in output
#' @param keep_results Whether to keep results in output
#' @param progress_each Whether to show progress for each test
#' @return List with results, summary table, and plots
run_perm_suite <- function(metrics_df,
                           groups_named_list,
                           metrics = c("leaf_mean","root_mean"),
                           n_perm = 10000,
                           seed = 42,
                           out_dir = "perm_results",
                           save_plots = TRUE,
                           keep_plots = TRUE,
                           keep_results = TRUE,
                           progress_each = TRUE) {
  
  # Create output directory if needed
  if (!dir.exists(out_dir) && save_plots) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Helper to parse metric names
  parse_metric <- function(m) {
    parts <- strsplit(m, "_", fixed = TRUE)[[1]]
    side <- parts[1]  # "leaf" or "root"
    opt  <- ifelse(length(parts) > 1, parts[2], "mean")
    list(side = side, opt = opt)
  }
  
  # Initialize containers
  results_hier <- list()
  plots_hier   <- list()
  rows <- list()
  
  # Run tests for each group and metric
  for (gname in names(groups_named_list)) {
    for (mcol in metrics) {
      meta <- parse_metric(mcol)
      message(sprintf(">> %s | %s_%s", gname, meta$side, meta$opt))
      
      res <- perm_test_one(metrics_df, groups_named_list[[gname]], mcol,
                           n_perm = n_perm, seed = seed,
                           return_null = TRUE, gname)
      if (is.null(res)) next
      
      # Create plot title
      plot_title <- ifelse(mcol == "leaf_mean", 
                           "Mean distance to leaves",
                           "Mean distance to roots")
      
      # Create plot
      p <- plot_perm_hist(
        null_stats   = res$null_stats,
        obs_stat     = res$summary$obs_mean,
        null_mean    = res$summary$null_mean,
        p_two_tailed = res$summary$p_two_tailed,
        n_group      = res$summary$n_group,
        title        = sprintf("%s vs null distribution", gname),
        xlab         = "Mean under the null",
        bins         = 30
      )
      
      # Store in hierarchical structure
      if (keep_plots) {
        plots_hier[[gname]][[meta$side]][[meta$opt]] <- p
      }
      if (keep_results) {
        results_hier[[gname]][[meta$side]][[meta$opt]] <- res
      }
      
      # Prepare row for summary table
      row <- transform(res$summary,
                       group = gname,
                       metric = mcol)[, c("group","metric","n_group","n_perm",
                                          "obs_mean","null_mean","effect_size",
                                          "null_sd","z","p_greater","p_less","p_two_tailed")]
      rows[[paste(gname, mcol, sep = "_")]] <- row
    }
  }
  
  # Create summary table with adjusted p-values
  res_tbl <- do.call(rbind, rows)
  rownames(res_tbl) <- NULL
  if (!is.null(res_tbl) && nrow(res_tbl) > 0) {
    res_tbl$p_adjust_BH <- p.adjust(res_tbl$p_two_tailed, method = "BH")
  }
  res_tbl <- res_tbl[order(res_tbl$p_adjust_BH, res_tbl$group, res_tbl$metric), ]
  
  list(results = results_hier, 
       summary_table = res_tbl, 
       plots = plots_hier, 
       out_dir = out_dir)
}

#' Build clean directed graph from edge list
#' @param edge_df Data frame with edges
#' @param from_col Name of source column
#' @param to_col Name of target column
#' @return Clean igraph object
build_clean_graph <- function(edge_df, from_col = "Gen_from", to_col = "Gen_to") {
  stopifnot(is.data.frame(edge_df),
            all(c(from_col, to_col) %in% names(edge_df)))
  
  edges <- unique(edge_df[, c(from_col, to_col)])
  g0 <- graph_from_data_frame(edges, directed = TRUE)
  igraph::simplify(
    g0,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = "first"
  )
}

#' Reverse edges in a graph while preserving all vertices
#' @param g Input graph
#' @return Graph with reversed edges
reverse_edges <- function(g) {
  # Get all vertex names
  all_vertices <- V(g)$name
  
  # Get edge list
  el <- as_edgelist(g, names = TRUE)
  
  if (nrow(el) == 0) {
    # If no edges, return empty graph with same vertices
    return(make_empty_graph(n = length(all_vertices), directed = TRUE) %>% 
             set_vertex_attr("name", value = all_vertices))
  }
  
  # Reverse edges
  el_rev <- el[, 2:1, drop = FALSE]
  
  # Create new graph with ALL original vertices
  # Create vertices data frame first
  vertices_df <- data.frame(name = all_vertices)
  
  # Create graph from reversed edges with all vertices
  gt <- graph_from_data_frame(
    as.data.frame(el_rev), 
    directed = TRUE, 
    vertices = vertices_df
  )
  
  return(gt)
}

#' Summarize distance vector (exclude infinite values)
#' @param v Numeric vector of distances
#' @return Named vector with min, max, mean, and count
summarize_dist_vec <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0L) {
    c(min = NA_real_, max = NA_real_, mean = NA_real_, count = 0L)
  } else {
    c(min = min(v), max = max(v), mean = mean(v), count = length(v))
  }
}

#' Compute KEGG metrics for distance to leaves and roots
#' @param edge_df Edge list data frame
#' @param from_col Source column name
#' @param to_col Target column name
#' @param betw_normalized Whether to normalize betweenness
#' @return Data frame with gene-level metrics
compute_kegg_metrics_no_scc <- function(edge_df,
                                        from_col = "Gen_from",
                                        to_col = "Gen_to",
                                        betw_normalized = FALSE) {
  
  # Build clean graph
  g <- build_clean_graph(edge_df, from_col, to_col)
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(seq_len(vcount(g)))
  }
  
  # Calculate basic node metrics
  deg_in  <- degree(g, mode = "in")
  deg_out <- degree(g, mode = "out")
  deg_all <- degree(g, mode = "all")
  btw     <- betweenness(g, directed = TRUE, normalized = betw_normalized)
  
  # --- DISTANCES TO LEAVES ---
  leaves <- V(g)[deg_out == 0]
  if (length(leaves) == 0) {
    stop("No leaves (out-degree = 0) in the graph.")
  }
  
  D_leaf <- distances(g, v = V(g), to = leaves, mode = "out")
  leaf_stats <- t(apply(D_leaf, 1, summarize_dist_vec))
  leaf_stats <- as.data.frame(leaf_stats)
  colnames(leaf_stats) <- paste0("leaf_", colnames(leaf_stats))
  
  # --- DISTANCES TO ROOTS ---
  roots <- V(g)[deg_in == 0]
  if (length(roots) == 0) {
    stop("No roots (in-degree = 0) in the graph.")
  }
  
  # Use reversed graph for root distances
  gt <- reverse_edges(g)
  root_names <- V(g)[roots]$name
  
  # Verify all root vertices exist in reversed graph
  missing_vertices <- root_names[!root_names %in% V(gt)$name]
  if (length(missing_vertices) > 0) {
    warning("Some root vertices not found in reversed graph: ", 
            paste(missing_vertices, collapse = ", "))
  }
  
  D_root <- distances(gt, v = V(gt)$name, to = roots$name, mode = "out")
  root_stats <- t(apply(D_root, 1, summarize_dist_vec))
  colnames(root_stats) <- paste0("root_", colnames(root_stats))
  
  # Compile results
  out <- data.frame(
    gene = V(g)$name,
    # Leaf metrics
    leaf_min = leaf_stats[, "leaf_min"],
    leaf_max = leaf_stats[, "leaf_max"],
    leaf_mean = leaf_stats[, "leaf_mean"],
    n_leafs_access = as.integer(leaf_stats[, "leaf_count"]),
    # Root metrics
    root_min = root_stats[, "root_min"],
    root_max = root_stats[, "root_max"],
    root_mean = root_stats[, "root_mean"],
    n_roots_access = as.integer(root_stats[, "root_count"]),
    # Node centrality metrics
    deg_in = as.numeric(deg_in[V(g)]),
    deg_out = as.numeric(deg_out[V(g)]),
    deg_all = as.numeric(deg_all[V(g)]),
    betweenness = as.numeric(btw[V(g)]),
    stringsAsFactors = FALSE
  )
  
  # Sort by leaf mean and gene name
  out[order(out$leaf_mean, out$gene), ]
}


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