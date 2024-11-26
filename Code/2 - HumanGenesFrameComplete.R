library(biomaRt)
library(dplyr)
library("GenomicRanges")

# GENOMIC RANGES OBJECT ########################################################

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c("hgnc_symbol", "ensembl_gene_id" ,"entrezgene_id","chromosome_name", 
               "transcription_start_site", "start_position","end_position", "strand", "gene_biotype"),
  mart = mart)


genes = genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                           "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                           "20", "21", "22"),]

colnames(genes) = c("hgnc_symbol", "Ensemble", "Entrez", "chr", "transcription_start_site", "start_transcript", "end_transcript", "strand", "Function") 

genes$transcription_start_site = NULL
genes = unique(genes)
genes = genes %>% filter(hgnc_symbol != "")
genes$start = genes$start_transcript = genes$start_transcript
genes$end = genes$end_transcript = genes$end_transcript
genes$strand = str_replace(genes$strand, "-1", "-")
genes$strand = str_replace(genes$strand, "1", "+")
genes = genes[order(genes$chr),]

GeneCols = genes[,c("hgnc_symbol", "Ensemble", "Entrez", "start_transcript", "end_transcript")]

# GENE
geneinf = makeGRangesFromDataFrame(genes)
genome(geneinf)='hg19'
for( colnm in colnames(GeneCols)){
  values(geneinf)[[colnm]]=GeneCols[[colnm]]
}

names(geneinf) = geneinf$hgnc_symbol

saveRDS(geneinf, "maglab/Gustavo/Data/Retrieved/Ranges/GenomicRanges.rds")



### HUMAN GENES FRAME COMPLETE (WITH TRANSCRIPTION FACTOR START POINT) #########
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c("hgnc_symbol", "ensembl_gene_id" ,"entrezgene_id","chromosome_name", 
               "transcription_start_site", "start_position","end_position", "strand", "gene_biotype"),
  mart = mart)


genes = genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                           "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                           "20", "21", "22"),]
colnames(genes) = c("Gene", "Ensemble", "Entrez", "chr", "tss", "start_transcript", "end_transcript", "strand", "Function") 

saveRDS(genes, "maglab/Gustavo/Data/Retrieved/Ranges/HumanGenesFrameComplete.rds")


# HUMAN GENES FRAME ############################################################

# Connect to the Ensembl BioMart
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define the chromosome names (1 to 22, X, Y, MT)
chromosome_names <- c(as.character(1:22), "X", "Y", "MT")

# Retrieve HGNC symbols for all chromosomes in a single loop
chromosome_data <- lapply(chromosome_names, function(chr) {
  getBM(
    filters = "chromosome_name",
    attributes = c("chromosome_name", "hgnc_symbol"),
    values = chr,
    mart = ensembl_mart
  )
})

# Combine the data into a single data frame
chromosome_data <- do.call(rbind, chromosome_data)

# Display the resulting data frame
head(chromosome_data)

GenesFrame =  getBM(attributes =  c('hgnc_symbol', "start_position", "end_position"),
                    filters = c('hgnc_symbol'),
                    values = 'ATP6V1G3',#HumanGenes$hgnc_symbol,
                    mart = ensembl_mart)
NumericIndices = !is.na(as.numeric(GenesFrame$chromosome_name))
GenesFrame[NumericIndices,]


HumanGenes['start_position'] = GenesFrame$start_position
HumanGenes['end_position'] = GenesFrame$end_position

HG = HumanGenes
GF = GenesFrame

sHG = HG[match(GF$hgnc_symbol,HG$hgnc_symbol),]
sHG$start_position = GenesFrame$start_position
sHG$end_position = GenesFrame$end_position

colnames(sHG) = c("hgnc_symbol",     "entrezgene_id",   "ensembl_gene_id" ,"gene_biotype"   , "chromosome_name"       ,      "start_position"  ,"end_position")

saveRDS(sHG, 'maglab/Gustavo/Data/Retrieved/Ranges/HumanGenesFrame.rds')


### Query for genemap ##########################################################

genemap <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "description"),
  mart = ensembl
)

saveRDS(genemap, "maglab/Gustavo/Data/Retrieved/Ranges/genemap.rds")


### EXCLUDING MHC ##############################################################

GeneRangeFrm = readRDS('maglab/Gustavo/Data/Retrieved/Ranges/HumanGenesFrame.rds')

HlaGenFrm = GeneRangeFrm %>% filter(chromosome_name == 6 & start_position >= 28477797 & end_position <= 33448354)

MhcGenArr = HlaGenFrm$hgnc_symbol %>% unique()

saveRDS(MhcGenArr,"maglab/Gustavo/Data/Retrieved/Ranges/MhcGeneArray.rds")


### HGNC AND ENSEMBLE RELATIONSHIP #############################################

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name"),
                mart = mart)

saveRDS(G_list, "maglab/Gustavo/Data/Retrieved/Ranges/HgncEnsembl.rds")


