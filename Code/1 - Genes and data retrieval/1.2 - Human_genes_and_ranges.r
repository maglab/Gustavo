suppressPackageStartupMessages({
  library(biomaRt)
  library(GenomicRanges)
  library(stringr)
})

out_dir <- "Data/Retrieved/Genes_and_diseases/Ranges"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Settings (paper-safe)
ENSEMBL_VERSION <- 115
GENOME_BUILD    <- "hg38"  # GRCh38
OUTPUT_TAG      <- paste0("hg38_v", ENSEMBL_VERSION)

# ------------------------------------------------------------------------------
# Helper: connect to Ensembl with mirror fallback, and STOP on failure
connect_ensembl <- function(version, mirrors = c("www", "useast", "asia")) {
  # Avoid accidentally reusing an old object
  if (exists("mart", envir = .GlobalEnv)) rm(mart, envir = .GlobalEnv)
  
  last_err <- NULL
  for (m in mirrors) {
    message("Connecting to Ensembl v", version, " (mirror = ", m, ") ...")
    mart_try <- tryCatch(
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        version = version,
        mirror  = m
      ),
      error = function(e) e
    )
    if (!inherits(mart_try, "error")) {
      message("Connected: ", mart_try@host)
      return(mart_try)
    } else {
      last_err <- mart_try
      message("  Failed: ", conditionMessage(mart_try))
    }
  }
  stop("Could not connect to Ensembl (version ", version, "). Last error: ",
       conditionMessage(last_err))
}

# ------------------------------------------------------------------------------
# 1) Connect (hard fail if SSL breaks)
mart <- connect_ensembl(version = ENSEMBL_VERSION)

# Save connection metadata for reproducibility
meta <- list(
  genome_build     = GENOME_BUILD,
  ensembl_version  = ENSEMBL_VERSION,
  ensembl_host     = mart@host,
  retrieved_utc    = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(meta, file.path(out_dir, paste0("ensembl_meta_", OUTPUT_TAG, ".rds")))

# ------------------------------------------------------------------------------
# 2) Retrieve gene table (HGNC, Ensembl, Entrez, coords, strand, biotype, description)
attrs_main <- c(
  "hgnc_symbol",
  "ensembl_gene_id",
  "entrezgene_id",
  "chromosome_name",
  "transcription_start_site",
  "start_position",
  "end_position",
  "strand",
  "gene_biotype",
  "description"
)

genes_raw <- getBM(attributes = attrs_main, mart = mart)

# Filter to autosomes and non-empty symbols
genes_raw <- genes_raw[genes_raw$chromosome_name %in% as.character(1:22), ]
genes_raw <- genes_raw[genes_raw$hgnc_symbol != "" & !is.na(genes_raw$hgnc_symbol), ]
genes_raw <- unique(genes_raw)

# ------------------------------------------------------------------------------
# 3) Annotated gene frame (same column names as your hg19 pipeline)
human_genes_annotated <- genes_raw
colnames(human_genes_annotated) <- c(
  #"Gene",
  #"Ensemble",
  #"Entrez",
  #"chr",
  #"tss",
  #"start_transcript",
  #"end_transcript",
  #"strand",
  #"Function",
  #"description"
  
 "hgnc_symbol",
 "ensembl_gene_id",
 "entrezgene_id",
 "chromosome_name",
 "transcription_start_site",
 "start_position",
 "end_position",
 "strand",
 "gene_biotype",
 "description"
)

human_genes_annotated$description = NULL

saveRDS(
  human_genes_annotated,
  file.path(out_dir, paste0("human_genes_annotated_", OUTPUT_TAG, ".rds"))
)

# ------------------------------------------------------------------------------
# 4) genemap (Entrez, HGNC, Ensembl, description)
genemap <- human_genes_annotated[, c("Entrez", "Gene", "Ensemble", "description")]
colnames(genemap) <- c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "description")
genemap <- unique(genemap)

saveRDS(
  genemap,
  file.path(out_dir, paste0("genemap_", OUTPUT_TAG, ".rds"))
)

# ------------------------------------------------------------------------------
# 5) Build GRanges
genes_gr_df <- data.frame(
  hgnc_symbol      = human_genes_annotated$Gene,
  Ensemble         = human_genes_annotated$Ensemble,
  Entrez           = human_genes_annotated$Entrez,
  chr              = human_genes_annotated$chr,
  start_transcript = human_genes_annotated$start_transcript,
  end_transcript   = human_genes_annotated$end_transcript,
  strand           = human_genes_annotated$strand,
  Function         = human_genes_annotated$Function,
  stringsAsFactors = FALSE
)
genes_gr_df <- unique(genes_gr_df)

genes_gr_df$start  <- genes_gr_df$start_transcript
genes_gr_df$end    <- genes_gr_df$end_transcript
genes_gr_df$strand <- str_replace(genes_gr_df$strand, "-1", "-")
genes_gr_df$strand <- str_replace(genes_gr_df$strand, "1", "+")
genes_gr_df <- genes_gr_df[order(genes_gr_df$chr), ]

GeneCols <- genes_gr_df[, c("hgnc_symbol", "Ensemble", "Entrez", "start_transcript", "end_transcript")]

geneinf <- makeGRangesFromDataFrame(genes_gr_df)
genome(geneinf) <- GENOME_BUILD

for (colnm in colnames(GeneCols)) {
  values(geneinf)[[colnm]] <- GeneCols[[colnm]]
}
names(geneinf) <- geneinf$hgnc_symbol

saveRDS(
  geneinf,
  file.path(out_dir, paste0("GenomicRanges_", OUTPUT_TAG, ".rds"))
)

# ------------------------------------------------------------------------------
# 6) MHC gene list (GRCh38 coordinates)
# chr6:28,510,120-33,480,577
mhc_start <- 28510120
mhc_end   <- 33480577

mhc_idx <- human_genes_annotated$chr == "6" &
  human_genes_annotated$start_transcript >= mhc_start &
  human_genes_annotated$end_transcript   <= mhc_end

MhcGeneArray <- unique(human_genes_annotated$Gene[mhc_idx])

saveRDS(
  MhcGeneArray,
  file.path(out_dir, paste0("MhcGeneArray_", OUTPUT_TAG, ".rds"))
)

# ------------------------------------------------------------------------------
# 7) Ensembl -> external gene name mapping
# (NO filters needed; if you want to restrict, pass values=unique(human_genes_annotated$Ensemble))
HgncEnsembl <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart       = mart
)
HgncEnsembl <- unique(HgncEnsembl)

saveRDS(
  HgncEnsembl,
  file.path(out_dir, paste0("HgncEnsembl_", OUTPUT_TAG, ".rds"))
)

message("Done. Outputs written with tag: ", OUTPUT_TAG)

