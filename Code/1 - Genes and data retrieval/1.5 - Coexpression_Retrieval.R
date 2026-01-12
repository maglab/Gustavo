# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Col - Column
# Cor - Correlation (e.g., Cor_long, Cor)
# Cox - Coexpression (e.g., CoxPrwFrm90, CoxPrwFrm95)
# Ens - Ensemble (e.g., EnsHgnFrm, DupEnsNam90)
# Nam - Name (e.g., DupEnsNam90, DupEnsNam95)
# Dup - Duplicated (e.g., DupEnsNam90, DupEnsNam95)
# Frm - Frame (e.g., EnsHgnFrm, FltFltEnsHgnFrm, CoxPrwFrm90)
# Hgn - Hgnc (e.g., EnsHgnFrm, FltEnsHgnFrm)
# Flt - Filtered (i.e., genes removed) (e.g., FltEnsHgnFrm, FltFltEnsHgnFrm)
# Gen - Gene (e.g., AllGen90, AllGen95)
# All - All (e.g., AllGen90, AllGen95)
# Abs - Absolute value (e.g., abs_Cor_long)
# Pre - Previous (e.g., PreEnsGenCor, inferred from the commented code)
# Prw - Pairwise (e.g., CoxPrwFrm90, CoxPrwFrm95)
# Phn - Phenotype

### LIBRARIES ##################################################################

library(tidyverse)
library(data.table)
library(stringr)
library(biomaRt)

################################################################################
# COX NETWORKS 
################################################################################

# LOAD GENE CORRELATION/COEXPRESION MATRIX 
Cor = read_tsv("Data/Retrieved/Network_Sources/Coexpression/GeneFriends/human_genes_correlation_matrix.tsv")

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# LONG FORMAT OF FRAME
Cor_long <- pivot_longer(Cor, cols = -Ensemble, names_to = "column", values_to = "value")
abs_Cor_long = Cor_long %>% mutate(value = abs(value))

# Coexpression 90
filtered_Cor_long90 = abs_Cor_long %>% filter(value > 0.90) %>% dplyr::select(Ensemble, column)
colnames(filtered_Cor_long90) = c("row","col")
write.table(filtered_Cor_long90, "Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/MatrixCox90_Ensembl.csv")

EnsCox90 = c(filtered_Cor_long90$row, filtered_Cor_long90$col) %>% unique()

GenCox90 = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                  filters = 'ensembl_gene_id', 
                  values = EnsCox90, 
                  mart = ensembl)

GenCox90 %>% saveRDS("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/GenesCox90.rds")



# Coexpression 95
filtered_Cor_long95 = abs_Cor_long %>% filter(value > 0.95) %>% Dupyr::select(Ensemble, column)
colnames(filtered_Cor_long95) = c("row","col")
write.table(filtered_Cor_long95, "Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/MatrixCox95_Ensembl.csv")

EnsCox95 = c(filtered_Cor_long95$row, filtered_Cor_long95$col) %>% unique()

GenCox95 = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                    filters = 'ensembl_gene_id', 
                    values = EnsCox95, 
                    mart = ensembl)

DupEnsNam95 %>% saveRDS("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/GenesCox95.rds")


################################################################################
# COEXPRESSION MATRIX FOR AGEING AND DISEASE-RELATED GENES 
################################################################################

# Extract Ensembl gene IDs from the row and column names of PreEnsGenCor matrix
gene_ids <- unique(c(rownames(Cor), colnames(Cor)))

EnsHgnFrm = readRDS("Data/Retrieved/Genes_and_diseases/Ranges/HgncEnsembl.rds") %>% 
  filter(ensembl_gene_id %in% gene_ids)

row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
FltEnsHgnFrm = EnsHgnFrm %>% filter(external_gene_name != "")

NonDupEns = FltEnsHgnFrm$ensembl_gene_id[!duplicated(FltEnsHgnFrm$external_gene_name)]
FltFltEnsHgnFrm = FltEnsHgnFrm[NonDupEns,]

# MATRIX WITH HGNA NAMES -------------------------------------------------------
HgnGenCor = Cor[FltFltEnsHgnFrm$ensembl_gene_id,FltFltEnsHgnFrm$ensembl_gene_id]

EnsRow = row.names(HgnGenCor)
EnsCol = colnames(HgnGenCor)

HgnRow = EnsHgnFrm[EnsRow,"external_gene_name"]
HgnCol = EnsHgnFrm[EnsCol,"external_gene_name"]

# Replace Ensembl gene IDs with gene names in the row and column names of the PreEnsGenCor matrix
rownames(HgnGenCor) <- HgnRow
colnames(HgnGenCor) <- HgnCol

# FILTER FOR DISEASE AND AGEING-RELATED GENES ----------------------------------

GenArr = GenPhnFrm$Gene %>% unique()
ArdAgeGenArr = intersect(GenArr, row.names(HgnGenCor))

HgnGenCor_ArdAge = HgnGenCor[PhnGenArr,PhnGenArr]

saveRDS(HgnGenCor_ArdAge,"Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/MatrixCox_AgeARD_Hgnc.rds") 
