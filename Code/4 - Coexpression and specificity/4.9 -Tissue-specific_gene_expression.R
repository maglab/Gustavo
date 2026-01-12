# =============================================================================
# TISSUE SPECIFICITY AND PLEIOTROPY ANALYSIS
# =============================================================================
# Analysis of tissue-specific gene expression patterns across disease communities
# =============================================================================

# Load required libraries
library(dplyr)
library(igraph)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(stringr)
library(ComplexHeatmap)
library(viridisLite)
library(circlize)
library(grid)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Visualization settings
PLOT_CONFIG <- list(
  plot_type = "box",               # "violin" or "box"
  alpha = 0.4,
  text_angle = 45,
  text_sizes = list(
    test = 2.5,
    annotation = 3.0,
    axis_label = 8
  ),
  colors = list(
    primary = "#F8696D",
    secondary = "#5A8AC6",
    background = "#CAE2EE"
  )
)

# Plot titles and labels
PLOT_LABELS <- list(
  title = "Tissue Specificity vs. Direct Pleiotropy",
  x_axis = "Number of Indirectly\nConnected Communities",
  y_axis = "Tau Specificity Index"
)

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

load_gene_community_data <- function() {
  # Load gene-community associations
  gene_community_df <- readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
    dplyr::rename(
      gene = Gene,
      phenotype = Code,
      phenotype_meaning = ARD_Meaning,
      community = ARC_Meaning
    )
  
  # Load aging-related gene sets
  human_aging_genes <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
  model_aging_genes <- readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")
  
  # Add aging gene annotations
  human_aging_df <- data.frame(
    gene = human_aging_genes,
    phenotype = "GenAge",
    phenotype_meaning = "GenAge",
    community = "GenAge"
  )
  
  model_aging_df <- data.frame(
    gene = model_aging_genes,
    phenotype = "ModAge",
    phenotype_meaning = "ModAge",
    community = "ModAge"
  )
  
  # Combine all gene annotations
  combined_df <- rbind(gene_community_df, model_aging_df, human_aging_df) %>% 
    unique()
  
  return(combined_df)
}

load_tissue_data <- function() {
  # Load tissue specificity data
  tau_df <- read.csv("Data/Retrieved/Specificity/Tau_gene_V8.csv")
  tissue_disease_df <- read.csv("Data/Retrieved/Specificity/Tissues_Diseses_TestAnnotation.csv")
  ensembl_mapping <- readRDS("Data/Retrieved/Genes_and_diseases/Ranges/HgncEnsembl_hg38_v115.rds")
  
  # Map Ensembl IDs to gene symbols
  rownames(ensembl_mapping) <- ensembl_mapping$ensembl_gene_id
  tau_df$gene <- ensembl_mapping[tau_df$gene_id, "external_gene_name"]
  
  # Clean and filter data
  tau_df <- tau_df %>% 
    filter(!is.na(gene), gene != "")
  
  # Extract tissue columns
  tissue_columns <- setdiff(colnames(tau_df), c("gene_id", "tau", "gene"))
  
  return(list(
    tau = tau_df,
    tissue_disease = tissue_disease_df,
    tissue_columns = tissue_columns
  ))
}

# =============================================================================
# PART 1: TISSUE BOXPLOTS (5x3 grid plots)
# =============================================================================

# Function to create tissue boxplots (from your original code)
make_tissue_boxplots <- function(df,
                                 log_x = FALSE,
                                 drop_outliers = FALSE,
                                 x_limits = NULL,
                                 x_limits_mode = c("zoom", "truncate"),
                                 show_titles = TRUE,
                                 title_size = 9,
                                 x_text_size = 7,
                                 x_axis_title_size = 7,
                                 y_text_size = 7,
                                 y_axis_title_size = 7,
                                 ascending_top = TRUE) {
  
  x_limits_mode <- match.arg(x_limits_mode)
  
  df2 <- df %>%
    dplyr::mutate(
      community_display = dplyr::case_when(
        community == "GenAge" ~ "GenAge.Hum",
        community == "ModAge" ~ "GenAge.Mod",
        community == "immunological/systemic disorders" ~ "Immunological Disorders",
        community == "immunological disorders" ~ "Immunological Disorders",
        TRUE ~ as.character(community)
      ),
      tissue = factor(tissue)
    )
  
  tissue_levels <- levels(df2$tissue)
  if (is.null(tissue_levels)) tissue_levels <- unique(df2$tissue)
  
  # Helper function for inliers
  .inliers <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(v)
    q1  <- stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE)
    q3  <- stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE)
    iqr <- q3 - q1
    v[v >= (q1 - 1.5*iqr) & v <= (q3 + 1.5*iqr)]
  }
  
  plots <- list()
  for (ti in tissue_levels) {
    tissue_data <- df2 %>% dplyr::filter(tissue == ti)
    
    # Scale for ordering
    if (log_x) {
      x_ord <- log10(pmax(tissue_data$expression, .Machine$double.eps))
    } else {
      x_ord <- tissue_data$expression
    }
    
    ord_df <- tibble::tibble(community_display = tissue_data$community_display, x_ord = x_ord)
    
    if (drop_outliers) {
      ord_df <- ord_df %>%
        dplyr::group_by(community_display) %>%
        dplyr::reframe(x_ord = .inliers(x_ord))
    }
    
    # Order by median
    community_order <- ord_df %>%
      dplyr::group_by(community_display) %>%
      dplyr::summarise(med = stats::median(x_ord, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(med, community_display) %>%
      dplyr::pull(community_display)
    
    if (ascending_top) community_order <- rev(community_order)
    
    # Put "Immunological Disorders" first
    if ("Immunological Disorders" %in% community_order) {
      community_order <- c(setdiff(community_order, "Immunological Disorders"),
                           "Immunological Disorders")
    }
    
    tissue_data <- tissue_data %>% 
      dplyr::mutate(community_display = factor(community_display, levels = community_order))
    
    # Colors
    colors_all <- setNames(rep("#CAE2EE", length(community_order)), community_order)
    if ("GenAge.Hum" %in% names(colors_all)) colors_all["GenAge.Hum"] <- "#F8696B"
    if ("GenAge.Mod" %in% names(colors_all)) colors_all["GenAge.Mod"] <- "#F8696B"
    if ("Immunological Disorders" %in% names(colors_all)) colors_all["Immunological Disorders"] <- "#5A8AC6"
    
    # Plot
    p <- ggplot(tissue_data, aes(x = expression, y = community_display, fill = community_display)) +
      geom_boxplot(
        width = 0.7,
        outlier.shape = if (drop_outliers) NA else 16,
        outlier.alpha = if (drop_outliers) NA else 0.4
      ) +
      scale_fill_manual(values = colors_all, guide = "none") +
      labs(
        title = if (show_titles) as.character(ti) else NULL,
        x = "",  # "Expression (RPKM)",
        y = ""   # "Gene set"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title   = element_text(face = "bold", hjust = 0.5, size = title_size),
        axis.text.x  = element_text(size = x_text_size),
        axis.title.x = element_text(size = x_axis_title_size),
        axis.text.y  = element_text(size = y_text_size),
        axis.title.y = element_text(size = y_axis_title_size)
      )
    
    # Scales and limits
    if (log_x) {
      p <- p + scale_x_log10()
      if (!is.null(x_limits)) {
        lims <- sort(as.numeric(x_limits))[1:2]
        lims[1] <- max(lims[1], .Machine$double.eps)
        p <- if (x_limits_mode == "zoom") p + coord_cartesian(xlim = lims)
        else p + scale_x_log10(limits = lims)
      }
    } else {
      if (!is.null(x_limits)) {
        p <- if (x_limits_mode == "zoom") p + coord_cartesian(xlim = x_limits)
        else p + scale_x_continuous(limits = x_limits)
      }
    }
    
    plots[[ti]] <- p
  }
  
  return(plots)
}

# Function to create tissue expression data frame
create_tissue_expression_data <- function(gene_community_df, tissue_data) {
  # Create tissue expression data frame for all communities
  communities <- unique(gene_community_df$community)
  tissue_expr_df <- data.frame()
  
  for (community in communities) {
    # Get community-specific genes
    community_genes <- gene_community_df %>% 
      filter(community %in% !!community) %>% 
      pull(gene)
    
    if (length(community_genes) > 0) {
      # Get tissue expression for these genes
      community_expr <- tissue_data$tau %>% 
        filter(gene %in% community_genes) %>% 
        select(all_of(tissue_data$tissue_columns))
      
      # Reshape to long format
      for (tissue in tissue_data$tissue_columns) {
        tissue_values <- community_expr[[tissue]]
        temp_df <- data.frame(
          tissue = tissue,
          expression = tissue_values,
          community = community
        )
        tissue_expr_df <- rbind(tissue_expr_df, temp_df)
      }
    }
  }
  
  return(tissue_expr_df)
}

# =============================================================================
# PART 2: HEATMAP FUNCTION
# =============================================================================

plot_grouped_heatmaps_equal <- function(
    expression_matrix, 
    group_annotations, 
    tau_values,
    gene_column = "gene",
    label_size = 8, 
    label_rot = 30,
    palette = "viridis",
    center_white = NULL,
    tau_width_mm = 3,
    tau_range = NULL,
    tau_na_fill = "global_median",
    tau_border_lwd = 0.8,
    row_title_size = 10,
    row_title_face = "bold",
    expr_legend_title_text = "Expression\n(RPKM)",
    expr_legend_title_size = 8,
    expr_legend_labels_size = 7,
    tau_legend_title_text = "tau",
    tau_legend_title_size = 8,
    tau_legend_labels_size = 7,
    tau_colname_text = "tau",
    tau_colname_size = 8,
    tau_colname_show_first_only = TRUE,
    tau_colname_rotation = 0,
    padding = grid::unit(c(10, 2, 10, 2), "mm"),
    use_raster = NULL,
    quiet = TRUE,
    normalize_upper = FALSE,
    verbose = TRUE) {
  
  if (quiet && exists("ht_opt")) ht_opt$message = FALSE
  
  stopifnot(all(c("gene", "set") %in% names(group_annotations)))
  stopifnot(all(c("gene", "tau") %in% names(tau_values)))
  
  # Normalization function
  .normalize <- function(x) {
    x <- trimws(as.character(x))
    if (normalize_upper) x <- toupper(x)
    x
  }
  
  # Set row names if gene column exists
  if (!is.null(gene_column) && gene_column %in% colnames(expression_matrix)) {
    expression_matrix[[gene_column]] <- .normalize(expression_matrix[[gene_column]])
    if (any(duplicated(expression_matrix[[gene_column]]))) {
      expression_matrix <- expression_matrix %>%
        dplyr::group_by(.data[[gene_column]]) %>%
        dplyr::summarise(dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
                         .groups = "drop") %>%
        as.data.frame()
    }
    rownames(expression_matrix) <- expression_matrix[[gene_column]]
    expression_matrix[[gene_column]] <- NULL
  } else {
    rn <- rownames(expression_matrix)
    if (is.null(rn)) stop("Expression matrix must have row names as genes.")
    rownames(expression_matrix) <- .normalize(rn)
  }
  
  # Normalize group annotations and tau values
  group_annotations$gene <- .normalize(group_annotations$gene)
  tau_values$gene <- .normalize(tau_values$gene)
  tau_values$tau <- as.numeric(tau_values$tau)
  
  # Filter genes present in expression matrix
  group_annotations <- group_annotations %>% 
    dplyr::filter(gene %in% rownames(expression_matrix))
  
  if (nrow(group_annotations) == 0) {
    stop("No genes from group annotations are present in the expression matrix.")
  }
  
  # Create combined matrix
  combined_matrix <- as.matrix(expression_matrix[group_annotations$gene, , drop = FALSE])
  mode(combined_matrix) <- "numeric"
  combined_matrix[is.infinite(combined_matrix)] <- NA
  
  # Impute for column ordering
  impute_row_mean <- function(mat) {
    if (!anyNA(mat)) return(mat)
    rmean <- rowMeans(mat, na.rm = TRUE)
    for (i in seq_len(nrow(mat))) {
      nas <- is.na(mat[i, ])
      if (any(nas)) mat[i, nas] <- rmean[i]
    }
    mat
  }
  
  combined_imputed <- impute_row_mean(combined_matrix)
  
  # Global column order
  col_hc <- hclust(dist(t(combined_imputed)), method = "complete")
  col_order <- colnames(combined_imputed)[order.dendrogram(as.dendrogram(col_hc))]
  combined_matrix <- combined_matrix[, col_order, drop = FALSE]
  combined_imputed <- combined_imputed[, col_order, drop = FALSE]
  
  # Expression colormap
  exp_values <- combined_matrix[is.finite(combined_matrix)]
  exp_min <- min(exp_values)
  exp_max <- max(exp_values)
  
  pal <- switch(palette,
                viridis = viridisLite::viridis(256),
                magma   = viridisLite::magma(256),
                plasma  = viridisLite::plasma(256),
                inferno = viridisLite::inferno(256),
                cividis = viridisLite::cividis(256),
                viridisLite::viridis(256))
  
  if (is.null(center_white)) {
    exp_col_fun <- circlize::colorRamp2(seq(exp_min, 10, length.out = 256), pal)
  } else {
    mid <- if (identical(center_white, "median")) median(exp_values) else as.numeric(center_white)
    exp_col_fun <- circlize::colorRamp2(c(exp_min, mid, exp_max), c(pal[1], "white", pal[length(pal)]))
  }
  
  # Tau values
  tau_values <- tau_values[tau_values$gene %in% rownames(combined_matrix), , drop = FALSE]
  tau_vec <- setNames(tau_values$tau, tau_values$gene)
  
  if (is.null(tau_range)) {
    tau_min <- min(tau_values$tau, na.rm = TRUE)
    tau_max <- max(tau_values$tau, na.rm = TRUE)
  } else {
    tau_min <- tau_range[1]
    tau_max <- tau_range[2]
  }
  
  tau_mid <- (tau_min + tau_max) / 2
  tau_col_fun <- circlize::colorRamp2(c(tau_min, tau_mid, tau_max),
                                      c("#5A8AC6", "white", "#F8696B"))
  tau_breaks <- pretty(c(tau_min, tau_max), n = 4)
  
  # Create heatmaps per set
  sets <- unique(group_annotations$set)
  heatmap_list <- NULL
  tau_annotation_names <- character(length(sets))
  
  for (i in seq_along(sets)) {
    current_set <- sets[i]
    set_indices <- group_annotations$set == current_set
    
    set_matrix <- combined_matrix[set_indices, , drop = FALSE]
    set_matrix_imputed <- combined_imputed[set_indices, , drop = FALSE]
    
    # Row clustering
    row_clustering <- TRUE
    if (nrow(set_matrix_imputed) >= 2) {
      rhc <- try(hclust(stats::dist(set_matrix_imputed), method = "complete"), silent = TRUE)
      if (!inherits(rhc, "try-error")) row_clustering <- as.dendrogram(rhc)
    }
    
    # Tau values for this set
    set_genes <- rownames(set_matrix)
    set_tau_raw <- tau_vec[set_genes]
    
    if (verbose) {
      n_non_na <- sum(!is.na(set_tau_raw))
      message(sprintf("[%-20s] Non-NA tau: %d/%d", current_set, n_non_na, length(set_tau_raw)))
    }
    
    set_tau <- set_tau_raw
    
    # Handle NA tau values
    if (anyNA(set_tau)) {
      if (identical(tau_na_fill, "none")) {
        # Leave as NA
      } else if (identical(tau_na_fill, "group_median")) {
        set_tau[is.na(set_tau)] <- median(set_tau, na.rm = TRUE)
      } else if (identical(tau_na_fill, "global_median")) {
        set_tau[is.na(set_tau)] <- median(tau_vec, na.rm = TRUE)
      } else if (identical(tau_na_fill, "min")) {
        set_tau[is.na(set_tau)] <- tau_min
      } else if (identical(tau_na_fill, "max")) {
        set_tau[is.na(set_tau)] <- tau_max
      } else if (is.numeric(tau_na_fill) && length(tau_na_fill) == 1L) {
        set_tau[is.na(set_tau)] <- tau_na_fill
      }
    }
    
    # Create annotation name
    annotation_name <- paste0("tau_", make.names(current_set))
    tau_annotation_names[i] <- annotation_name
    
    # Right annotation (tau bar)
    right_annotation <- do.call(ComplexHeatmap::rowAnnotation, c(
      setNames(
        list(ComplexHeatmap::anno_simple(
          set_tau,
          which = "row",
          col = tau_col_fun,
          na_col = "grey85",
          width = grid::unit(tau_width_mm, "mm"),
          border = FALSE,
          gp = grid::gpar(col = NA),
          pt_size = grid::unit(0, "mm")
        )),
        annotation_name
      ),
      list(
        annotation_name_gp = grid::gpar(fontsize = 0),
        show_legend = FALSE
      )
    ))
    
    heatmap_name <- paste0("Expr_", make.names(current_set))
    
    # Create heatmap
    set_heatmap <- ComplexHeatmap::Heatmap(
      set_matrix,
      name = heatmap_name,
      col = exp_col_fun,
      cluster_rows = row_clustering,
      cluster_columns = FALSE,
      column_order = col_order,
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_gp = grid::gpar(fontsize = label_size),
      column_names_rot = label_rot,
      rect_gp = grid::gpar(col = NA),
      border = TRUE,
      row_title = as.character(current_set),
      row_title_gp = grid::gpar(fontface = row_title_face, fontsize = row_title_size),
      right_annotation = right_annotation,
      height = grid::unit(1, "null"),
      use_raster = use_raster,
      show_heatmap_legend = (i == 1),
      heatmap_legend_param = list(
        title = expr_legend_title_text,
        title_gp = grid::gpar(fontsize = expr_legend_title_size),
        labels_gp = grid::gpar(fontsize = expr_legend_labels_size)
      )
    )
    
    heatmap_list <- if (is.null(heatmap_list)) set_heatmap else heatmap_list %v% set_heatmap
  }
  
  # Tau legend
  tau_legend <- ComplexHeatmap::Legend(
    title = tau_legend_title_text,
    title_gp = grid::gpar(fontsize = tau_legend_title_size),
    labels_gp = grid::gpar(fontsize = tau_legend_labels_size),
    col_fun = tau_col_fun,
    at = tau_breaks,
    labels = format(tau_breaks, trim = TRUE)
  )
  
  # Draw the heatmap
  ComplexHeatmap::draw(
    heatmap_list,
    heatmap_legend_side = "right",
    merge_legends = TRUE,
    ht_gap = grid::unit(2, "mm"),
    annotation_legend_list = list(tau_legend),
    padding = padding
  )
  
  # Add vertical separators for expression heatmaps
  for (s in sets) {
    ComplexHeatmap::decorate_heatmap_body(paste0("Expr_", make.names(s)), slice = 1, {
      nc <- ncol(combined_matrix)
      if (nc > 1) {
        xs <- seq(1/nc, 1 - 1/nc, length.out = nc - 1)
        for (x in xs) {
          grid::grid.lines(
            x = grid::unit(c(x, x), "npc"),
            y = grid::unit(c(0, 1), "npc"),
            gp = grid::gpar(col = "black", lwd = 0.6)
          )
        }
      }
    })
  }
  
  # Add border to tau bars
  for (ann_name in tau_annotation_names) {
    ComplexHeatmap::decorate_annotation(ann_name, slice = 1, {
      grid::grid.rect(gp = grid::gpar(col = "black", lwd = tau_border_lwd, fill = NA))
    })
  }
  
  # Add tau label below last panel
  if (!is.null(tau_colname_text) && nzchar(tau_colname_text)) {
    last_ann <- utils::tail(tau_annotation_names, 1)
    ComplexHeatmap::decorate_annotation(last_ann, slice = 1, {
      grid::grid.text(
        tau_colname_text,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0, "npc") - grid::unit(2, "mm"),
        rot = tau_colname_rotation,
        gp = grid::gpar(fontsize = tau_colname_size)
      )
    })
  }
  
  invisible(NULL)
}

# =============================================================================
# MAIN EXECUTION - BOTH ANALYSES
# =============================================================================

cat("Loading data...\n")
gene_community_data <- load_gene_community_data()
tissue_data <- load_tissue_data()

# =============================================================================
# PART 1: GENERATE TISSUE BOXPLOTS (5x3 grids)
# =============================================================================

cat("\n=== PART 1: Generating Tissue Expression Boxplots ===\n")

# Create tissue expression data frame
tissue_expr_data <- create_tissue_expression_data(gene_community_data, tissue_data)

# Filter out specific communities if needed
filtered_tissue_data <- tissue_expr_data %>% 
  dplyr::filter(!(community %in% c("CelAge", "GenDR", "cancer")))

# Create boxplots
cat("Creating tissue boxplots...\n")
tissue_boxplots <- make_tissue_boxplots(
  filtered_tissue_data,
  log_x = FALSE,
  drop_outliers = FALSE
)

# Split into two sets for the 5x3 grid display
tissue_names <- names(tissue_boxplots)
n_tissues <- length(tissue_names)

# First grid: tissues 1-15
if (n_tissues >= 15) {
  cat("\nGrid 1: Tissues 1-15 (5 rows x 3 columns)\n")
  grid1_plots <- tissue_boxplots[1:15]
  grid1_labels <- LETTERS[1:15]
  
  plot_grid1 <- cowplot::plot_grid(
    plotlist = grid1_plots,
    ncol = 3,
    nrow = 5,
    labels = grid1_labels,
    label_size = 12
  )
  
  # Print first grid
  print(plot_grid1)
  cat("First grid displayed (15 tissues)\n")
  
  # Second grid: tissues 16-30
  if (n_tissues >= 30) {
    cat("\nGrid 2: Tissues 16-30 (5 rows x 3 columns)\n")
    grid2_plots <- tissue_boxplots[16:30]
    grid2_labels <- LETTERS[1:15]  # Reset labels for second grid
    
    plot_grid2 <- cowplot::plot_grid(
      plotlist = grid2_plots,
      ncol = 3,
      nrow = 5,
      labels = grid2_labels,
      label_size = 12
    )
    
    # Print second grid
    print(plot_grid2)
    cat("Second grid displayed (15 tissues)\n")
  } else if (n_tissues > 15) {
    cat(sprintf("\nGrid 2: Tissues 16-%d\n", n_tissues))
    grid2_plots <- tissue_boxplots[16:n_tissues]
    grid2_labels <- LETTERS[1:length(grid2_plots)]
    
    plot_grid2 <- cowplot::plot_grid(
      plotlist = grid2_plots,
      ncol = 3,
      labels = grid2_labels,
      label_size = 12
    )
    
    print(plot_grid2)
  }
} else {
  # If fewer than 15 tissues, create a single grid
  cat(sprintf("\nCreating single grid with %d tissues\n", n_tissues))
  plot_single_grid <- cowplot::plot_grid(
    plotlist = tissue_boxplots,
    ncol = 3,
    labels = LETTERS[1:n_tissues],
    label_size = 12
  )
  
  print(plot_single_grid)
}




# Carta vertical
ggsave(
  filename = "Last_Figures/Supp_Tissues_Box1.svg",
  plot = plot_grid1,
  device = "svg",
  width = 8.5,
  height = 11,
  units = "in"
)

# Carta vertical
ggsave(
  filename = "Last_Figures/Supp_Tissues_Box2.svg",
  plot = plot_grid2,
  device = "svg",
  width = 8.5,
  height = 11,
  units = "in"
)


# =============================================================================
# PART 2: GENERATE HEATMAP
# =============================================================================

cat("\n=== PART 2: Generating Grouped Heatmap ===\n")

# Prepare data for heatmap
cat("Preparing heatmap data...\n")

# Get unique genes and their pleiotropy status
gene_pleiotropy <- gene_community_data %>%
  select(gene, community) %>%
  filter(!(community %in% c('GenAge', 'ModAge'))) %>%
  unique() %>%
  group_by(gene) %>%
  summarise(frequency = n()) %>%
  mutate(pleiotropy = ifelse(frequency >= 4, "High_ARC.Pleiotropy", "Low_ARC.Pleiotropy"))

# Create final group annotations
heatmap_annotations <- gene_community_data %>%
  select(gene, community) %>%
  filter(community %in% c("GenAge", "ModAge")) %>%
  mutate(set = case_when(
    community == "GenAge" ~ "GenAge.Hum",
    community == "ModAge" ~ "GenAge.Mod"
  )) %>%
  select(gene, set) %>%
  rbind(
    gene_pleiotropy %>%
      select(gene, set = pleiotropy)
  )

# Prepare expression matrix
expression_matrix <- tissue_data$tau %>%
  select(-c(gene_id)) %>%
  distinct(gene, .keep_all = TRUE)

# Prepare tau values
tau_values <- tissue_data$tau %>%
  select(gene, tau)

# Filter to genes with valid tau values
valid_genes <- tau_values %>%
  filter(!is.na(tau)) %>%
  pull(gene)

expression_matrix <- expression_matrix %>%
  filter(gene %in% valid_genes)

heatmap_annotations <- heatmap_annotations %>%
  filter(gene %in% valid_genes)

tau_values <- tau_values %>%
  filter(gene %in% valid_genes)

# Remove gene column from expression matrix for heatmap
expression_matrix_heatmap <- expression_matrix
rownames(expression_matrix_heatmap) <- expression_matrix_heatmap$gene
expression_matrix_heatmap$gene <- NULL
expression_matrix_heatmap$tau <- NULL

# Generate the heatmap
cat("Generating heatmap...\n")
tryCatch({
  plot_grouped_heatmaps_equal(
    expression_matrix_heatmap,
    heatmap_annotations,
    tau_values,
    gene_col = NULL,  # Already set as row names
    label_size = 6,
    label_rot = 45,
    palette = "viridis",
    tau_range = c(0, 1),
    padding = grid::unit(c(2, 2, 2, 2), "mm")
  )
  cat("Heatmap generated successfully.\n")
}, error = function(e) {
  cat("Error generating heatmap:", e$message, "\n")
})

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Total genes analyzed: %d\n", nrow(gene_community_data)))
cat(sprintf("Communities analyzed: %d\n", length(unique(gene_community_data$community))))
cat(sprintf("Tissues analyzed: %d\n", length(tissue_data$tissue_columns)))
cat("Two outputs generated:\n")
cat("1. Tissue expression boxplots (5x3 grid layout)\n")
cat("2. Grouped heatmap with tau specificity annotations\n")