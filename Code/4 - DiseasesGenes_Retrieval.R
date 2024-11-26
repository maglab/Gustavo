### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Arc - Ageing-related Community (e.g., ArdArcFrm, GenArdArcFrm)
# Ard - Ageing-related disease (e.g., GenArdFrm_10K$Ard, GenArdFrm_10K_NoMhc$Ard)
# Frm - Frame (e.g., GenArdFrm_10K, GenArdFrm_10K_NoMhc)
# Gen - Gene (e.g., GenArdFrm_10K$Gen, GenArdFrm_10K_NoMhc$Gen)
# gwas - Genome Wide Analysis Studies (e.g., gwasRes, sgwasRes)
# Hrc - Hierarchy (e.g., HrcFrm)
# Hum - Human (e.g., GenAge.Hum)
# Qry - Query
# Mhc - Major Histocompatibility Complex (e.g., MhcGenArr)
# Mod - Model Organisms (e.g., GenAge.Hum)
# N - Number of (e.g., N_Nodes)
# Res - Results
# sgp - SNP-Gene-Phenotype (e.g., sgpList_d10K, sgpList_10K)
# Thr - Threshold
# Tss - Transcription start site

### LIBRARIES ##################################################################

library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)
library(qqman)

##### GENES FRAME ####################################################################################################################

genemap = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/genemap.rds") 
genes.granges = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/GenomicRanges.rds") 
GenesFrame = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/HumanGenesFrameComplete.rds") 
DiseasesCategories = read.csv("maglab/Gustavo/Data/Retrieved/Showcase/Ukb_DiseaseHierarchy.csv") 
colnames(DiseasesCategories)
colnames(GenesFrame)[colnames(GenesFrame) == "hgnc_symbol"] ="Genes"

##############################################################################################################################
### START ####################################################################################################################
##############################################################################################################################

Tss = FALSE 

Nodes = DiseasesCategories[DiseasesCategories$Ageing %in% c(1,2),]$Node
Names = DiseasesCategories [DiseasesCategories$Ageing %in% c(1,2),]$Meaning
Nodes = as.character(Nodes)

Interval_Gene = list()
sgpList_10K = list()

N_Nodes = length(Nodes)

for(i in 1:N_Nodes){
  
  QryTimeStart <- Sys.time()
  
  paste("### PHENOTYPE: ", i, "/", N_Nodes," #############################################", sep="") %>% print()
  noquote("")
  
  Node = Nodes[i]
  Name = Names[i]
  
  # GET PHENOTYPE DATA
  
  FilePath = paste("maglab/Gustavo/Data/Retrieved/Ukbb/a", Node, ".imp.stats", sep = "") # We did not upload these files to GitHub due to theis size.They can be downloaded at: BioStudies (S- BSST407) 
  PathExists = FilePath %>% file.exists()
  
  if(PathExists){
    
    # GET PHENOTYPE DATA
    gwasRes <- read_tsv( FilePath )
    PgwasRes = gwasRes
    gwasRes$P_BOLT_LMM = as.numeric(gwasRes$P_BOLT_LMM)
    gwasRes = AdaptGwasFrame(gwasRes)
    sgwasRes = gwasRes[gwasRes$pvalue < 5E-8,]
    
    sgpList_10K[[as.character(Node)]] = summary2genes_intervals(sgwasRes, Node, Name, thr = 10000, GenesFrame, genes.granges)
    
  }
}

# SELECT DISEASE CODE AND GENES
GenArdFrm_10K = FrameFromList(sgpList_10K) %>% dplyr::select(PhenoCode, Gene) %>% unique()
colnames(GenArdFrm_10K) = c("Ard","Gen")
GenArdFrm_10K$Ard = paste("D",GenArdFrm_10K$Ard,sep="")

# FILTERING OUT MHC
MhcGenArr = readRDS("maglab/Gustavo/Data/Retrieved/Ranges/MhcGeneArray.rds")
GenArdFrm_10K_NoMhc = GenArdFrm_10K %>% dplyr::filter(!(Gen %in% MhcGenArr))

# INTEGRATING GENAGE ###########################################################

# AGEING-RELATED GENES
GenAgeHum = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeMod = readRDS("maglab/Gustavo/Data/Generated/HAGR/GenAgeMod_Genes.rds")

GenAgeHumFrm = data.frame(Ard='HumAge',Gen=GenAgeHum)
GenAgeModFrm = data.frame(Ard='ModAge',Gen=GenAgeMod)

# BINDING ARDs and AGEING
ArdAgeFrm = GenArdFrm_10K_NoMhc %>% rbind(GenAgeHumFrm) %>% rbind(GenAgeModFrm)
colnames(ArdAgeFrm) = c("Code", "Gene")
saveRDS(ArdAgeFrm,"maglab/Gustavo/Data/Generated/General/ARD_GenAge_GeneFrame.rds")

# CREATE ARD - ARC FRAME BASED ON THE HIERARCHY
HrcFrm = read.csv("maglab/Gustavo/Data/Generated/Showcase/Ukb_DiseaseRoots.csv")
HrcFrm$X = NULL
ArdArcFrm = HrcFrm %>% dplyr::select(Node,Meaning,RootMeaning)
colnames(ArdArcFrm) = c('Code','ARD_Meaning','ARC_Meaning')

GenArdArcFrm = merge(ArdAgeFrm,ArdArcFrm,by='Code') %>% 
  dplyr::select('Code','ARD_Meaning','ARC_Meaning', 'Gene') #%>%

saveRDS(GenArdArcFrm,"maglab/Gustavo/Data/Generated/General/ARD_ARC_GeneFrame.rds")


###################################################################################################################################
##### FUNCTIONS ###################################################################################################################
###################################################################################################################################

gwas2GRanges <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}



gwas2GRanges_chr <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  if (all(substr(seqlevels(snpinf),1,3)!='chr')) {
    seqlevels(snpinf)=paste('chr',seqlevels(snpinf),sep='')
  }
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}


snp2gene_eQTL <- function(gwasRes, eQTLfile, genemap){
  tissue <- strsplit(sapply(strsplit(eQTLfile,'/'),function(x)x[length(x)]),'[.]')[[1]][1]
  eqtl <- read_tsv(eQTLfile) %>%
    inner_join(gwasRes) %>%
    dplyr::select(SNP,CHR, BP, ALLELE1, ALLELE0, gene_id, slope, pval_beta, tissue)%>%
    unique()
  genemap <- genemap %>%
    unique()%>%
    dplyr::mutate(entrezgene=as.character(entrezgene))%>%
    dplyr::rename(gene_id = ensembl_gene_id)
  eqtl <- left_join(eqtl, genemap)
  eqtl$entrezgene[eqtl$entrezgene=='']=NA
  eqtl$hgnc_symbol[eqtl$hgnc_symbol=='']=NA
  eqtl$description[eqtl$description=='']=NA
  eqtl
}



FrameFromList <- function(sgpList){
  sgpFrame = data.frame()
  for (i in 1:length(sgpList)){
    Frame = sgpList[[i]]
    if(dim(Frame)[1] > 0)
      sgpFrame = rbind(sgpFrame, Frame)
  }
  return(sgpFrame)
}



range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = as.numeric(df[,2]), end = as.numeric(df[,3]))
  )
  return(gr)
}


# EXTEND RANGES
extend <- function(x, upstream=0, downstream=0)     
{
  names = names(x)
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
  names(x) = names
  return(x)
}


# ADAPT HEALTHSPAN
AdaptHealthspan <- function(gwasRes)     
{
  names(gwasRes)[names(gwasRes) == 'SNPID'] <- 'SNP'
  names(gwasRes)[names(gwasRes) == 'chr'] <- 'CHR'
  names(gwasRes)[names(gwasRes) == 'pos'] <- 'BP'
  names(gwasRes)[names(gwasRes) == 'EA'] <- 'ALLELE1'
  names(gwasRes)[names(gwasRes) == 'RA'] <- 'ALLELE0'
  names(gwasRes)[names(gwasRes) == 'EAF'] <- 'A1FREQ'
  names(gwasRes)[names(gwasRes) == 'beta'] <- 'BETA'
  names(gwasRes)[names(gwasRes) == 'se'] <- 'SE'
  return(gwasRes)
}



AdaptGwasFrame <- function(gwasRes)     
{
  names = colnames(gwasRes)
  if ("X.log10.p.value." %in% names){
    gwasRes$pvalue = 10^(gwasRes$X.log10.p.value.)
  }
  if ("P_BOLT_LMM" %in% names){
    gwasRes$pvalue = gwasRes$P_BOLT_LMM
    gwasRes$X.log10.p.value. = -log10(gwasRes$pvalue)
  }
  return(gwasRes)
}



GetSnp2GeneDistance <- function(sgpFrame)     
{
  distance = c()
  for(i in 1:dim(sgpFrame[1])){
    if(sgpFrame[i,"BP"] < sgpFrame[i,"start_position"]){
      distance[i] = paste("-",as.character(sgpFrame[i,"start_position"] - sgpFrame[i,"BP"]))
    }
    else if(sgpFrame[i,"BP"] > sgpFrame[i,"end_position"]){
      distance[i] = paste("+",as.character(sgpFrame[i,"BP"] - sgpFrame[i,"end_position"]))
    }
    else{
      distance[i] = "within"
    }
  }
  sgpFrame$Snp2Gene_distance = distance
  return(sgpFrame)
}



GetSnp2GeneDistance_Apply <- function(SnpGeneFrame, Tss){
  
  
  SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(
      Snp2Gene_distance = case_when(
        BP < start_position ~ BP - start_position,
        BP > end_position ~ BP - end_position, 
        TRUE ~ 0
        )
      )
  
  if(Tss){
    SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(Snp2Tss_distance = BP - transcription_start_site)
  }
  
  
  return(SnpGeneFrame)
}  
 


# ADAPT PVALUE APPLY
GetSnp2GeneDistance_Apply_For <- function(SnpGeneFrame, Tss){
  
  Snp2Gene_distance = c()
  
  if(!Tss){
    for(k in 1:dim(SnpGeneFrame)[1]){
      if(SnpGeneFrame$BP[k] < SnpGeneFrame$start_position[k]){
        Snp2Gene_distance[k] = SnpGeneFrame$BP[k] - SnpGeneFrame$start_position[k]
      }else if(SnpGeneFrame$BP[k] > SnpGeneFrame$end_position[k]){
        Snp2Gene_distance[k] = SnpGeneFrame$BP[k] - SnpGeneFrame$end_position[k]
      } else{
        Snp2Gene_distance[k] = 0
      }
    }
  }
  
  SnpGeneFrame$Snp2Gene_distance = Snp2Gene_distance
  
  if(Tss){
    SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(Snp2Tss_distance = BP - transcription_start_site)
  }
  
  return(SnpGeneFrame)
}  



summary2genes_intervals <- function(sgwasRes, Node, Name, thr, GenesFrame, genes.granges){
  significants = dim(sgwasRes)[1]
  if (significants > 0){
    sgwasRes = AdaptGwasFrame(sgwasRes)
    
    gwas_as_GR <- gwas2GRanges(sgwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')
    
    r1 = extend(genes.granges, upstream=thr, downstream=thr) #genes.granges
    r2 = gwas_as_GR
    
    names(r2) = r2$RefSNP_id
    
    overlap <- GenomicRanges::findOverlaps(r1, r2)
    snps <- names(r2)[overlap@to]
    genes <- names(r1)[overlap@from]
    snps_genes_frame = data.frame(snps, genes)
    
    colnames(snps_genes_frame) = c("SNP", "Genes")

    SnpGeneFrame = merge(sgwasRes[,c("CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "X.log10.p.value.", "pvalue", "SNP")], snps_genes_frame,by="SNP",all=TRUE)

    SnpGeneFrame = merge(SnpGeneFrame,GenesFrame,by="Genes",all=FALSE)
    
    if(dim(SnpGeneFrame)[1] > 0){

      SnpGeneFrame$PhenoName = Name
      SnpGeneFrame$PhenoCode = Node
      
      
      ChrCond = SnpGeneFrame$chromosome_name == SnpGeneFrame$CHR
      
      SnpGeneFrame$ChrCond = ChrCond
      
      SnpGeneFrame$Interval_any = (SnpGeneFrame$BP >= SnpGeneFrame$start_position - thr) & (SnpGeneFrame$BP <= SnpGeneFrame$end_position + thr)
      SnpGeneFrame$Interval_gene = (SnpGeneFrame$BP >= SnpGeneFrame$start_position) & (SnpGeneFrame$BP <= SnpGeneFrame$end_position)
      SnpGeneFrame$Interval_near = xor(SnpGeneFrame$Interval_gene, SnpGeneFrame$Interval_any)
      
      SnpGeneFrame = GetSnp2GeneDistance(SnpGeneFrame) 
      
      # -------------------------------------------------------------------------------------------------
      
      gwas_as_GR <- gwas2GRanges_chr(sgwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')
      # FURTHER ASSOCIATIONS
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      varInfo <- unique(granges(gwas_as_GR))
      allvar <- locateVariants(varInfo, txdb,
                               AllVariants(promoter=PromoterVariants(upstream=thr,
                                                                     downstream=thr)))
      overs <- findOverlaps(gwas_as_GR,allvar)
      
      NumData = dim(as.data.frame(allvar, row.names = NULL))[1]
      
      if (NumData > 0 ) {
        
        for(colnm in colnames(mcols(gwas_as_GR))){
          mydat <-(mcols(gwas_as_GR)[[colnm]])[queryHits(overs)]
          mcols(allvar)[[colnm]] = '*'
          mcols(allvar)[[colnm]][subjectHits(overs)] <- mydat
        }
        geneids <- unique(allvar$GENEID)
        geneids <- geneids[complete.cases(geneids)]
        genemap <- genemap %>% unique() %>%
          mutate(entrezgene=as.character(entrezgene))
        
        mcols(allvar)$RefSNP_id
        mcols(allvar)$LOCATION
        
        SnpDef = unique(data.frame("SNP" = mcols(allvar)$RefSNP_id, "Location" = mcols(allvar)$LOCATION))
        
        GeneDef = unique(data.frame("Genes" = genemap$hgnc_symbol, "Description" = genemap$description))

        SnpGeneFrame = merge(SnpGeneFrame, unique(SnpDef), by = "SNP", all.x = TRUE, all.y = FALSE)
        SnpGeneFrame = merge(SnpGeneFrame, GeneDef, by = "Genes", all.x = TRUE, all.y = FALSE)
        
        
        SnpGeneFrame = SnpGeneFrame[,c("PhenoCode", "PhenoName", "SNP", "pvalue", "X.log10.p.value.", "CHR", "BP", "Location",
                                      "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE",
                                      "Genes", "chromosome_name", "start_position", "end_position", "ChrCond", "Interval_any",
                                      "Interval_gene", "Interval_near", "Snp2Gene_distance", "gene_biotype", "Description")] #Z would be penultime if required
        
        SnpGeneFrame <- SnpGeneFrame %>% 
          dplyr::rename(CHR_SNP = CHR,
                        BP_SNP = BP,
                        Location_SNP = Location,
                        Gene = Genes,
                        CHR_Gene = chromosome_name,
                        Gene_Start = start_position,
                        Gene_End = end_position,
                        CHR_Congruency = ChrCond,
                        Gene_Description = Description)

      } else {
        SnpGeneFrame = data.frame()
      }
    } else {
      SnpGeneFrame = data.frame()
    }
  } else {
    SnpGeneFrame = data.frame()
  }
  
  return(SnpGeneFrame)
}



MinP_Gene = function(G, Expanded_Summary){
  min_p = Expanded_Summary %>% filter(Gen == G) %>% pull(pvalue) %>% min()
  Expanded_Summary %>% filter(Gen == G) %>% mutate(min_pvalue = min_p)
}
