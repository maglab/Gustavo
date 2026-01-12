# ==============================================================================
# GENE ONTOLOGY ENRICHMENT ANALYSIS PIPELINE
# Author: Gustavo Daniel Vega Magdaleno
# Date: 30/10/2025
# Description: Comprehensive GO enrichment analysis for disease, aging, and 
#              network-derived gene sets with redundant term simplification
# ==============================================================================

# Load required libraries ------------------------------------------------------
library(clusterProfiler)    # For GO enrichment analysis
library(org.Hs.eg.db)       # Human genome annotation database
library(dplyr)              # Data manipulation
library(stringr)            # String operations
library(readr)              # File reading/writing
library(GOSemSim)           # GO semantic similarity computation
library(openxlsx)           # Excel file operations
library(EnsDb.Hsapiens.v86) # Ensembl database for human genome annotations (protein-coding genes)

# Note: GOSemSim is loaded twice in original - keeping once is sufficient

# ==============================================================================
# SECTION 1: DATA LOADING AND PREPROCESSING
# ==============================================================================

# Define base data directory path
DATA_DIR_RETRIEVED <- "Data/Retrieved"
DATA_DIR_GENERATED <- "Data/Generated"


# Load all protein-coding genes ------------------------------------------------
AllArr <- genes(EnsDb.Hsapiens.v86, 
                filter = ~ gene_biotype == "protein_coding") %>%
                as.data.frame() %>%
                pull(symbol) %>%
                unique()

# Load disease-associated genes -------------------------------------------------
DisArr <- readRDS(file.path(DATA_DIR_RETRIEVED, "Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds")) %>%
  pull(Gene) %>% 
  unique()

# Extract immune-related genes from disease data -------------------------------
ImmArr <- readRDS(file.path(DATA_DIR_RETRIEVED, "Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds")) %>%
  dplyr::filter(ARC_Meaning == "immunological/systemic disorders") %>%
  pull(Gene) %>% 
  unique()

# HIGH ARC-Pleiotropy ----------------------------------------------------------
ComGenFrm <- readRDS(file.path(DATA_DIR_RETRIEVED, "Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds")) %>%
  dplyr::select(Gene,ARC_Meaning) %>%
  unique() %>% 
  pull(Gene) %>%
  table() %>% 
  as.data.frame() %>%
  dplyr::filter(Freq >= 4) 
HghArr = ComGenFrm$. %>% as.character()

# Load aging-related gene sets -------------------------------------------------
GenAgeHumArr <- readRDS(file.path(DATA_DIR_RETRIEVED, "Genes_and_diseases/HAGR/GenAgeHum_Genes.rds"))
GenAgeModArr <- readRDS(file.path(DATA_DIR_RETRIEVED, "Genes_and_diseases/HAGR/GenAgeMod_Genes.rds"))

# Compute set operations for aging genes ---------------------------------------
GenAgeHumExc <- setdiff(GenAgeHumArr, GenAgeModArr)  # Human-specific aging genes
GenAgeModExc <- setdiff(GenAgeModArr, GenAgeHumArr)  # Model-specific aging genes
GenAgeInt <- intersect(GenAgeModArr, GenAgeHumArr)   # Common aging genes
GenAgeCmb <- c(GenAgeModArr, GenAgeHumArr) %>% unique()  # Combined aging genes

# Load network gene sets -------------------------------------------------------
# PPI network genes
GrpPin <- readRDS(file.path(DATA_DIR_GENERATED, "Networks_and_predictions/Networks/PPI/Lists/GraphsList.rds"))
GenPin <- igraph::V(GrpPin$GenGenNtw)$name

# KEGG network genes
GrpKeg <- readRDS(file.path(DATA_DIR_GENERATED, "Networks_and_predictions/Networks/KEGG/Lists/GraphsList.rds"))
GenKeg <- igraph::V(GrpKeg$GenGenNtw)$name

# COX90 network genes
GrpC90 <- readRDS(file.path(DATA_DIR_GENERATED, "Networks_and_predictions/Networks/COX90/Lists/GraphsList.rds"))
GenC90 <- igraph::V(GrpC90$GenGenNtw)$name

# COX95 network genes
GrpC95 <- readRDS(file.path(DATA_DIR_GENERATED, "Networks_and_predictions/Networks/COX95/Lists/GraphsList.rds"))
GenC95 <- igraph::V(GrpC95$GenGenNtw)$name

# Compute network gene set operations ------------------------------------------
NtwCmb <- c(GenPin, GenKeg, GenC90, GenC95) %>% unique()  # All network genes
NtwInt <- GenPin %>% 
  intersect(GenKeg) %>% 
  intersect(GenC90) %>% 
  intersect(GenC95)  # Genes in all four networks

# Save gene lists to files -----------------------------------------------------
GENES_DIR <- file.path("Data/Retrieved/Genes_and_diseases/Genes_List")

# Create output directory if it doesn't exist
if (!dir.exists(GENES_DIR)) {
  dir.create(GENES_DIR, recursive = TRUE)
}




# Save each gene list with descriptive names
write_gene_list <- function(genes, filename) {
  write.table(genes, 
              file = file.path(GENES_DIR, filename),
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
}

gene_list_path = "Data/Retrieved/Genes_and_diseases/Genes_List/"

write_gene_txt <- function(x, out_file) {
  x <- unique(na.omit(trimws(as.character(x))))
  x <- x[nzchar(x)]
  writeLines(x, out_file)
}


write_gene_txt(AllArr,       file.path(gene_list_path, "Protein_Coding_All.txt"))
write_gene_txt(DisArr,       file.path(gene_list_path, "Diseases_All.txt"))
write_gene_txt(GenAgeHumArr, file.path(gene_list_path, "GenAge_Hum.txt"))
write_gene_txt(GenAgeModArr, file.path(gene_list_path, "GenAge_Mod.txt"))
write_gene_txt(NtwCmb,       file.path(gene_list_path, "Networks_Combination.txt"))
write_gene_txt(NtwInt,       file.path(gene_list_path, "Networks_Intersection.txt"))
write_gene_txt(ImmArr,       file.path(gene_list_path, "Diseases_Immunological.txt"))
write_gene_txt(GenAgeHumExc, file.path(gene_list_path, "GenAge_Hum_Exclusive.txt"))
write_gene_txt(GenAgeModExc, file.path(gene_list_path, "GenAge_Mod_Exclusive.txt"))
write_gene_txt(GenAgeInt,    file.path(gene_list_path, "GenAge_Intersection.txt"))
write_gene_txt(GenAgeCmb,    file.path(gene_list_path, "GenAge_Combination.txt"))
write_gene_txt(HghArr,       file.path(gene_list_path, "High_ARC_Pleiotropy.txt"))


# ==============================================================================
# SECTION 2: PRE-COMPUTE SEMANTIC SIMILARITY DATA
# ==============================================================================

# Precompute semantic data for all GO ontologies to accelerate analysis
# This needs to be done only once per session
sem_list <- list(
  BP = godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE),  # Biological Process
  MF = godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE),  # Molecular Function
  CC = godata('org.Hs.eg.db', ont = "CC", computeIC = FALSE)   # Cellular Component
)

# ==============================================================================
# SECTION 3: GENE ONTOLOGY ENRICHMENT ANALYSES
# ==============================================================================


perform_go_analysis(
  query_genes = GenAgeHumExc,
  background_genes = GenAgeCmb,
  output_filename = paste(Enrichmnt_path,"GenAgeHumExc_AllGenAgeBackground.xlsx",sep=""),
  analysis_name = "Human-Exclusive Aging Genes",
  use_all_background = FALSE
)


# Define enrichment analysis wrapper function ----------------------------------
perform_go_analysis <- function(query_genes, background_genes, 
                                output_filename, analysis_name,
                                use_all_background = FALSE) {
  
  cat("\n=====================================================================\n")
  cat("Performing GO enrichment analysis:", analysis_name, "\n")
  cat("Query genes:", length(query_genes), "\n")
  cat("Background genes:", length(background_genes), "\n")
  cat("=====================================================================\n")
  
  # Execute GO enrichment with redundant term removal
  res <- run_go_enrichment_multi_fast(
    query_genes = query_genes,
    background_genes = background_genes,
    use_all_background = use_all_background,
    keyType = "SYMBOL",
    onts = "BP",  # Focus on Biological Process ontology
    top_n = 80,
    remove_redundant = TRUE,
    similarity_cutoff = 0.8,
    sem_list = sem_list
  )
  
  # Format and sort results
  if (nrow(res) > 0) {
    res <- res %>% 
      relocate(`Adjusted p_value_num`, .after = last_col()) %>%
      arrange(`Adjusted p_value_num`)
    
    # Create output directory if needed
    #output_dir <- file.path(DATA_DIR, "Genes_Enrichment")
    #output_dir <- file.path("Data/Generated/Genes_Enrichment")
    #if (!dir.exists(output_dir)) {
    #  dir.create(output_dir, recursive = TRUE)
    #}
    
    # Save to Excel
    output_path <- output_filename#file.path(output_dir, output_filename)
    write.xlsx(res, file = output_path, rowNames = FALSE)
    cat("Results saved to:", output_path, "\n")
  } else {
    cat("No significant GO terms found for", analysis_name, "\n")
  }
  
  return(res)
}

Enrichmnt_path = "Data/Generated/Genes_Enrichment/"

# Analysis 1: Human-specific aging genes ---------------------------------------
perform_go_analysis(
  query_genes = GenAgeHumExc,
  background_genes = GenAgeCmb,
  output_filename = paste(Enrichmnt_path,"GenAgeHumExc_AllGenAgeBackground.xlsx",sep=""),
  analysis_name = "Human-Exclusive Aging Genes",
  use_all_background = FALSE
)

# Analysis 2: Model-specific aging genes ---------------------------------------
perform_go_analysis(
  query_genes = GenAgeModExc,
  background_genes = GenAgeCmb,
  output_filename = paste(Enrichmnt_path,"GenAgeModExc_AllGenAgeBackground.xlsx",sep=""),
  analysis_name = "Model-Exclusive Aging Genes",
  use_all_background = FALSE
)

# Analysis 3: Common aging genes -----------------------------------------------
perform_go_analysis(
  query_genes = GenAgeInt,
  background_genes = GenAgeCmb,
  output_filename = paste(Enrichmnt_path,"GenAgeIntersection_AllGenAgeBackground.xlsx",sep=""),
  analysis_name = "Common Aging Genes (Intersection)",
  use_all_background = FALSE
)

# Analysis 4: Immune-related disease genes -------------------------------------
perform_go_analysis(
  query_genes = ImmArr,
  background_genes = DisArr,
  output_filename = paste(Enrichmnt_path,"Immunological_AllDiseaseBackground.xlsx",sep=""),
  analysis_name = "Immune-Related Disease Genes",
  use_all_background = TRUE
)


# Analysis 5: High ARC-Pleiotropy ----------------------------------------------
perform_go_analysis(
  query_genes = HghArr,
  background_genes = DisArr,
  output_filename = paste(Enrichmnt_path,"High_ARC_Pleiotropy_AllDiseaseBackground.xlsx",sep=""),
  analysis_name = "High ARC Pleiotropy Genes",
  use_all_background = TRUE
)


# ==============================================================================
# SECTION 4: MACHINE LEARNING PREDICTED GENE ANALYSIS
# ==============================================================================

# Load machine learning predictions --------------------------------------------
#PREDICTIONS_DIR <- file.path(DATA_DIR_GENERATED, "Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions/MachineLearning-Based")
#PrdFrm <- read.csv(file.path(PREDICTIONS_DIR, "All_Int_Mux_Auc88.csv"))
PrdFrm <- read.csv("Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions/Top_Genes_Best_Algorithm/TopGenes_All_Int_Mux_Auc88.csv")


# Extract top 30 predicted non-aging genes
#PrdGen <- PrdFrm %>% 
#  filter(Class == "NotAge") %>% 
#  pull(Label)

PrdGen = PrdFrm$Gen

TopPrdGen <- PrdGen[1:30]

# Save predicted gene list
write_gene_list(TopPrdGen, "Predicted30_All_Int_Mux_Auc88.txt")

# Analysis 5: Predicted genes in multiplex networks ----------------------------
perform_go_analysis(
  query_genes = TopPrdGen,
  background_genes = NtwCmb,
  output_filename = paste(Enrichmnt_path,"Predicted30_All_Int_Mux_Auc88.xlsx",sep=""),
  analysis_name = "Top 30 Predicted Genes (Multiplex Networks)",
  use_all_background = TRUE
)

# Analysis 6: Predicted genes in PPI network -----------------------------------
# Note: Using same genes but different background would require separate network
# For demonstration, using same analysis with focused output
res_ppi <- perform_go_analysis(
  query_genes = TopPrdGen,
  background_genes = GenPin,  # PPI-specific background
  output_filename = paste(Enrichmnt_path,"Predicted30_PPI.xlsx",sep=""),
  analysis_name = "Top 30 Predicted Genes (PPI Network)",
  use_all_background = TRUE
)

# Select key columns for concise reporting
if (nrow(res_ppi) > 0) {
  res_ppi <- res_ppi %>% 
    select(`Term id`, `Term Name`, `Intersection Size`, `Adjusted p_value`)
}

# Analysis 7: Predicted genes in KEGG network ----------------------------------
res_kegg <- perform_go_analysis(
  query_genes = TopPrdGen,
  background_genes = GenKeg,  # KEGG-specific background
  output_filename = paste(Enrichmnt_path,"Predicted30_KEGG.xlsx",sep=""),
  analysis_name = "Top 30 Predicted Genes (KEGG Network)",
  use_all_background = TRUE
)

# Select key columns for concise reporting
if (nrow(res_kegg) > 0) {
  res_kegg <- res_kegg %>% 
    select(`Term id`, `Term Name`, `Intersection Size`, `Adjusted p_value`)
}

# ==============================================================================
# SECTION 5: UTILITY FUNCTIONS
# ==============================================================================

#' Map gene identifiers to ENTREZID
#' 
#' Converts gene symbols or ENSEMBL IDs to ENTREZID format required by clusterProfiler
#' 
#' @param genes Character vector of gene identifiers
#' @param keyType Type of input identifiers: "SYMBOL", "ENSEMBL", or "ENTREZID"
#' @return Character vector of ENTREZID identifiers
.map_to_entrez <- function(genes, keyType = c("SYMBOL", "ENSEMBL", "ENTREZID")) {
  keyType <- match.arg(keyType)
  genes <- unique(genes)
  
  # Return if already ENTREZID
  if (keyType == "ENTREZID") return(genes)
  
  # Perform ID mapping
  suppressMessages({
    mapped <- bitr(genes, 
                   fromType = keyType, 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
  })
  
  return(unique(mapped$ENTREZID))
}

#' Format GO enrichment results
#' 
#' Converts clusterProfiler output to a tidy tibble with formatted columns
#' 
#' @param ego enrichGO result object
#' @param q_size Number of query genes
#' @return Formatted tibble with enrichment results
.format_out <- function(ego, q_size) {
  df <- as.data.frame(ego)
  
  # Handle empty results
  if (NROW(df) == 0) {
    return(tibble(
      `Term Name` = character(),
      `Term id` = character(),
      `Adjusted p_value` = character(),
      `Adjusted p_value_num` = numeric(),
      `Term Size` = integer(),
      `Query Size` = integer(),
      `Intersection Size` = integer(),
      `Ontology` = character()
    ))
  }
  
  # Extract term size from BgRatio column
  term_size <- as.integer(vapply(strsplit(df$BgRatio, "/"), 
                                 function(x) x[[1]], ""))
  
  # Create formatted tibble
  result <- tibble(
    `Term Name` = df$Description,
    `Term id` = df$ID,
    `Adjusted p_value` = toupper(format(df$p.adjust, 
                                        scientific = TRUE, 
                                        digits = 3)),
    `Adjusted p_value_num` = df$p.adjust,  # Numeric for sorting
    `Term Size` = term_size,
    `Query Size` = q_size,
    `Intersection Size` = df$Count,
    `Ontology` = df$ONTOLOGY
  ) %>% 
    arrange(`Adjusted p_value_num`)
  
  return(result)
}

#' Fast GO enrichment analysis with redundancy removal
#' 
#' Performs GO enrichment with precomputed semantic similarity for speed
#' 
#' @param query_genes Query gene set
#' @param background_genes Background gene set (NULL for all annotated genes)
#' @param use_all_background Use all annotated genes as background
#' @param keyType Gene identifier type
#' @param onts GO ontologies to analyze
#' @param p_adjust_method P-value adjustment method
#' @param pvalue_cutoff P-value cutoff
#' @param qvalue_cutoff Q-value cutoff
#' @param min_gs_size Minimum gene set size
#' @param max_gs_size Maximum gene set size
#' @param remove_redundant Remove redundant GO terms
#' @param similarity_cutoff Semantic similarity cutoff for redundancy removal
#' @param simplify_by Metric for term simplification
#' @param top_n Return top N terms
#' @param sem_list Precomputed semantic similarity data
#' @return Data frame with enrichment results
run_go_enrichment_multi_fast <- function(
    query_genes,
    background_genes = NULL,
    use_all_background = TRUE,
    keyType = c("SYMBOL", "ENSEMBL", "ENTREZID"),
    onts = c("BP", "MF", "CC"),
    p_adjust_method = "BH",
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 1,
    min_gs_size = 15,
    max_gs_size = 5000,
    remove_redundant = TRUE,
    similarity_cutoff = 0.8,
    simplify_by = c("p.adjust", "pvalue", "qvalue"),
    top_n = 100,
    sem_list = NULL
) {
  
  # Parameter validation
  keyType <- match.arg(keyType)
  simplify_by <- match.arg(simplify_by)
  onts <- match.arg(onts, several.ok = TRUE)
  
  # Map gene identifiers to ENTREZID
  q_entrez <- .map_to_entrez(query_genes, keyType)

  # Define background gene set
  if (!use_all_background && !is.null(background_genes)) {
    bg_entrez <- .map_to_entrez(background_genes, keyType)
    q_entrez <- intersect(q_entrez, bg_entrez)
    
    if (length(q_entrez) == 0) {
      stop("Query genes share no identifiers with background after mapping.")
    }
  } else {
    bg_entrez <- NULL
  }
  
  # Perform enrichment for each ontology
  results <- list()
  
  for (ont in onts) {
    cat("  Analyzing", ont, "ontology...\n")
    
    # Run GO enrichment
    ego <- enrichGO(
      gene = q_entrez,
      universe = bg_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = p_adjust_method,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = FALSE
    )
    
    # Process results if any significant terms found
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      
      # 1. Filter to top N terms before simplification
      if (!is.null(top_n)) {
        ego@result <- ego@result %>% 
          arrange(p.adjust) %>% 
          dplyr::slice(1:min(n(), top_n))
      }
      
      # 2. Remove redundant terms using semantic similarity
      if (remove_redundant) {
        if (is.null(sem_list) || is.null(sem_list[[ont]])) {
          # Fallback: compute semantic data if not provided
          semData <- godata('org.Hs.eg.db', ont = ont, computeIC = FALSE)
        } else {
          semData <- sem_list[[ont]]
        }
        
        ego <- simplify(
          x = ego,
          cutoff = similarity_cutoff,
          by = simplify_by,
          select_fun = min,
          measure = "Wang",
          semData = semData
        )
      }
    }
    
    # Format and store results
    results[[ont]] <- .format_out(ego, q_size = length(q_entrez))
  }
  
  # Combine results from all ontologies
  final_results <- bind_rows(results[onts])
  
  return(final_results)
}

