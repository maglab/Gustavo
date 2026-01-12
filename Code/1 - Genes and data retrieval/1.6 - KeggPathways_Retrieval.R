library("KEGGlincs")
library("KEGGgraph")
library("StarBioTrek")
library("dplyr")
library("limma")
library("stringr")
library("plyr")
library('org.Hs.eg.db')

# Paquetes de Bioconductor
#BiocManager::install("limma")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("KEGGgraph")
#BiocManager::install("KEGGlincs")

### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Arr - Array (e.g., KegArr)
# Col - Column (e.g., OrgCol, NewCol)
# Cod - Code (e.g., KegCod)
# Dst - Destination (e.g., DstFile)
# Fcn - Function (e.g., SelFcn)
# Frm - Frame (e.g., KegFrm)
# Keg - KEGG (e.g., KegFrm, KegCod)
# Len - Length (e.g., KegLen)
# Nam - Name (e.g., RenNamFcn)
# Org - Original (e.g., OrgCol)
# Pul - Pull (e.g., PulFrm)
# Ren - Rename (e.g., RenNamFcn)
# Sel - Select (e.g., SelFcn)
# Sym - Symbol (e.g., GenCodKegFrm$Symbol)
# Tmp - Temporary (e.g., tmp)

###  DOWNLOAD LIST OF GENES-RELATED KEGG PATHWAYS ##############################                                                                                                                                            
GenCodKegFrm <- getGeneKEGGLinks(species="hsa")                                                                                 
GenCodKegFrm$Symbol <-  as.character(mapIds(org.Hs.eg.db, GenCodKegFrm$GeneID, 'SYMBOL', 'ENTREZID'))#mapIds(org.Hs.eg.db, tab$GeneID, column="SYMBOL", keytype="ENTREZID")
GenCodKegFrm$PathwayID = str_remove(GenCodKegFrm$PathwayID, "path:")
GenCodKegFrm$GeneID =  paste("hsa:", GenCodKegFrm$GeneID ,sep="")
head(GenCodKegFrm)

GenCodKegFrm %>% saveRDS('Data/Retrieved/Network_Sources/KEGG/Other/KeggCodingFrame.rds' )


 
### RETRIEVE KEGG PATHWAYS AND CREATE FRAME ####################################

GenCodKegFrm = readRDS("Data/Retrieved/Network_Sources/KEGG/Other/KeggCodingFrame.rds")

KegArr = GenCodKegFrm$PathwayID %>% unique()

KegFrm = data.frame()
KegLen = length(KegArr)
i=1
tmp<- tempfile()
for(i in 1:KegLen){
  print(i)
  KegCod = KegArr[i]
  KegUrl = KEGGgraph::retrieveKGML(KegCod,"hsa",tmp)
  #DstFile = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/Pathways/",KegCod,".xml",sep="")
  DstFile = paste("Data/Retrieved/Network_Sources/KEGG/Pathways/",KegCod,".xml",sep="")
  download.file(url=KegUrl,destfile=DstFile)
  sKegFrm = parseKGML2DataFrame(file=DstFile,expandGenes=T)
  if(nrow(sKegFrm)>0){
    sKegFrm$pathway = KegCod
    KegFrm = rbind(KegFrm, sKegFrm)
  }     
}

KegFrm %>% saveRDS("Data/Retrieved/Network_Sources/KEGG/Other/KeggInfo.rds")

### CREATE FINAL VERSION OF KEGG FRAME WITH GENE CODES #########################

GenCodKegFrm = readRDS('Data/Retrieved/Network_Sources/KEGG/Other/KeggCodingFrame.rds' )

KegFrm = readRDS("Data/Retrieved/Network_Sources/KEGG/Other/KeggInfo.rds") %>% as.data.table()

Keg_from = GenCodKegFrm %>% SelFcn(c('GeneID', 'Symbol')) %>% RenNamFcn(OrgCol=c('GeneID', 'Symbol'),NewCol=c("from", "Gen_from")) 
Keg_to   = GenCodKegFrm %>% SelFcn(c('GeneID', 'Symbol')) %>% RenNamFcn(OrgCol=c('GeneID', 'Symbol'),NewCol=c("to", "Gen_to")) 

Frm1 = inner_join(KegFrm, Keg_from, by="from") %>% unique()
Frm = inner_join(Frm1, Keg_to, by="to") %>% unique()

GenGenKeg = Frm %>% dplyr::select(Gen_from, Gen_to) %>% unique()

GenGenKeg %>% saveRDS("Data/Retrieved/Network_Sources/KEGG/Other/KEGG_Interactions.rds")

GenGenKeg = readRDS("Data/Retrieved/Network_Sources/KEGG/Other/KEGG_Interactions.rds")

################################################################################
# FUNCTIONS
################################################################################

RenNamFcn = function(OrgFrm, OrgCol, NewCol){
  ColNam = colnames(OrgFrm)
  LenRep = length(OrgCol)
  for(i in 1:LenRep){
    colnames(OrgFrm)[ColNam == OrgCol[i]] = NewCol[i]#TrnFrm
  }
  return(OrgFrm)
}

SelFcn = function(Frm, ColArr){
  return(Frm[ColArr])
}

PulFcn = function(Frm, Col){
  return(Frm[[Col]])
}