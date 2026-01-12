# Load required libraries
library(cowplot)
library(ggplot2)
library(dplyr)
library(rlang)

# Remove unused libraries: grid, svglite, scales are not needed for core functionality
# If svglite is needed for saving, keep only when actually saving

# ==============================================================================
# CONFIGURATION AND DATA LOADING
# ==============================================================================

# Set analysis parameters
ArcArdAss = "Arc"
ScrOpc = "Int"  # Options: "Int", "Prx", "Wlk", "Plt", "Dst"
NumPer = 10000

# Define network types to analyze
NtwArr = c('PPI', 'COX90', 'COX95', 'KEGG')

if(ScrOpc == "Dst"){
  ScrOpcTxt = "Distance"
}
if(ScrOpc == "Int"){
  ScrOpcTxt = "Interactors"
}
if(ScrOpc == "Wlk"){
  ScrOpcTxt = "Walker"
}

# Construct file path and load results
results_file = paste("Data/Generated/Permutations/Network_Associations/PermutationList_",
                     toupper(ArcArdAss), "_", ScrOpcTxt, "_", NumPer, ".rds", sep="")
results = readRDS(results_file)

# ==============================================================================
# PLOT GENERATION
# ==============================================================================

# Set plot type and data categories
Typ = "Unbalanced"  # Type of analysis
SetArr = c("AgeGen", "ModGen", "DisGen", "HghGen")  # Dataset categories

# Initialize list to store plots
plots = list()
plot_counter = 1
n_pvalue = 16

# Iterate through each network and dataset combination
for(NtwQry in NtwArr){
  for(SetQry in SetArr){
    
    # Extract base plot from results
    p = results[[NtwQry]][[Typ]][[SetQry]]$plot
    
    # --------------------------------------------------------------------------
    # EXTRACT STATISTICAL VALUES FROM SUBTITLE
    # --------------------------------------------------------------------------
    sub_txt <- p$labels$subtitle
    
    p_bonferroni <- function(p, n) {
      p_adj <- p * n
      if (p_adj > 1) p_adj <- 1
      return(p_adj)
    }
    
    pretty_num <- function(x) {
      x <- as.numeric(x)
      ifelse(x == 1, "1", format(x, scientific = TRUE))
    }
    
    # Extract P-value
    PVal <- sub(".*P-value: ([0-9\\.eE-]+).*", "\\1", sub_txt)
    P_adj = PVal %>% as.numeric() %>% p_bonferroni(n_pvalue) %>% pretty_num()
    #P_adj = #sprintf("%.1e", P_adj)
    
    # Extract Null mean
    NulVal <- sub(".*Null mean: ([0-9\\.eE-]+).*", "\\1", sub_txt)
    
    # Extract Group name (e.g., GenAgeHum)
    TrgGrp <- sub(".*red;'>\\s*([^:]+):.*", "\\1", sub_txt)
    
    # Extract Group value
    TrgVal <- sub(".*red;'>[^:]+: ([0-9\\.eE-]+).*", "\\1", sub_txt)
    
    # --------------------------------------------------------------------------
    # EXTRACT VERTICAL LINE INFORMATION
    # --------------------------------------------------------------------------
    gb <- ggplot_build(p)
    
    # Extract vline information from plot layers
    vline_info <- data.frame(
      colour = sapply(gb$plot$layers, 
                      function(l) l$aes_params$colour %||% NA),
      xintercept = sapply(gb$data, 
                          function(d) if (!is.null(d$xintercept)) d$xintercept[1] else NA)
    )
    
    # Filter only valid rows
    vline_info <- vline_info[!is.na(vline_info$xintercept), ]
    
    # Create named list by color
    vline_list <- setNames(as.list(vline_info$xintercept), vline_info$colour)
    
   
    
    # --------------------------------------------------------------------------
    # CUSTOMIZE PLOT
    # --------------------------------------------------------------------------
    plots[[plot_counter]] = p +
      # Remove original titles
      ggtitle(NULL) + 
      xlab(NULL) +
      ylab(NULL) + 
      
      
      
      # Create new subtitle with formatted statistics
      labs(subtitle = paste0(
        #"p-value: ", PVal, "<br>",
        "P_adj: ", P_adj, "<br>",
        "<span style='color:blue;'>Null mean: ", NulVal, "</span><br>",
        "<span style='color:red;'>", TrgGrp, ": ", TrgVal, "</span>"
      )) +
      
      # Theme adjustments
      theme(
        plot.subtitle = ggtext::element_markdown(),
        axis.text.x = element_text(size = 7, angle = 30, hjust = 0.5),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(fill = NA, color = "black", size = 0.1)
      ) +
      
      # Custom histogram styling
      geom_histogram(color = '#C6E4EE', size = 1, fill = NA) + 
      
      # Add vertical lines for statistical values
      geom_vline(xintercept = vline_list[["blue"]], 
                 colour = "blue", linewidth = 1.2, linetype = "dashed") +
      geom_vline(xintercept = vline_list[["red"]], 
                 colour = "red", linewidth = 1.2, linetype = "dashed") 
      
    if(ScrOpc == "Dst"){
      plots[[plot_counter]] = plots[[plot_counter]] + scale_x_continuous(labels = function(x) sprintf("%.1f", x))
    } 
    
    if(ScrOpc == "Int"){
      plots[[plot_counter]] = plots[[plot_counter]] + scale_x_continuous(labels = function(x) sprintf("%.1f", x))
    }
    
    if(ScrOpc == "Wlk"){
      plots[[plot_counter]] = plots[[plot_counter]] + scale_x_continuous(labels = function(x) sprintf("%.0e", x))
      }
      
    plot_counter = plot_counter + 1
  }
}

# ==============================================================================
# CREATE GRID LAYOUT
# ==============================================================================

# Define row and column titles
row_titles <- c("PPI", "COX90", "COX95", "KEGG")
col_titles <- c("GenAge.Hum", "GenAge.Mod", "Diseases", "High ARC-Pleiotropy")

# Create main grid of plots with labels
grid_plots <- plot_grid(
  plotlist = plots,
  ncol = 4,
  align = "hv",
  labels = letters[1:16],    # Labels a through p
  label_size = 8,
  label_fontface = "bold",
  hjust = -0.2,
  vjust = 1.2
)

# Add column titles
grid_with_cols <- plot_grid(
  NULL,  # Spacer
  # Create column title plots
  plot_grid(
    plotlist = lapply(col_titles, function(x) 
      ggdraw() + draw_label(x, fontface = "bold", size = 9)), 
    ncol = 4
  ),
  grid_plots,
  ncol = 1,
  rel_heights = c(0.02, 0.05, 1)  # Relative heights
)

# Add row titles (rotated 90 degrees)
row_labels <- plot_grid(
  plotlist = lapply(row_titles, function(x) 
    ggdraw() + draw_label(x, angle = 90, fontface = "bold", size = 9)), 
  ncol = 1
)

# Add "count" caption for y-axis
row_caption <- plot_grid(
  ggdraw() + draw_label("count", angle = 90, fontface = "plain", size = 8),
  ncol = 1
)

# Combine row labels, caption, and main grid
grid_final <- plot_grid(
  row_labels,          # Row titles
  row_caption,         # "count" label
  grid_with_cols,      # Main plot grid with column titles
  ncol = 3,
  rel_widths = c(0.05, 0.05, 1)  # Relative widths
)

# ==============================================================================
# ADD MAIN TITLE AND CAPTION
# ==============================================================================

if(ScrOpc == "Int"){
  TitTxt = "Permutation results: Number of ARC-Interactions"
  xTxt = "Number of ARC-Interactions"
  SavTxt = "Last_Figures/Supp_Permutation_ARC_Interactions.svg"
} 

if(ScrOpc == "Wlk"){
  TitTxt = "Permutation results: Mean RWR score from ARCs "
  xTxt = "Mean RWR score from ARCs"
  SavTxt = "Last_Figures/Supp_Permutation_ARC_Walker.svg"
}

if(ScrOpc == "Dst"){
  TitTxt = "Permutation results: Shortest path distance to ARCs"
  xTxt = "Mean distance to ARCs"
  SavTxt = "Last_Figures/Supp_Permutation_ARC_Distance.svg"
}


# Create main title
title <- ggdraw() + 
  draw_label(
    #"Permutation results: Number of ARC-Interactions",
    TitTxt,
    fontface = "bold",
    size = 12,
    hjust = 0.5,
    vjust = 1
  )

# Create bottom caption (x-axis label)
caption_bottom <- ggdraw() + 
  draw_label(
    #"Number of ARC-Interactions",
    xTxt,
    fontface = "plain",
    size = 8,
    hjust = 0.5
  )

# Combine all elements
grid_final_with_title_caption <- plot_grid(
  title,
  grid_final,
  caption_bottom,
  ncol = 1,
  rel_heights = c(0.05, 1, 0.05)
)

# Display final plot
print(grid_final_with_title_caption)

# ==============================================================================
# SAVE PLOT
# ==============================================================================

# Load svglite only when needed for saving
library(svglite)

svglite(
  SavTxt,
  width = 6.5,   # Width in inches
  height = 8     # Height in inches
)
print(grid_final_with_title_caption)
dev.off()

