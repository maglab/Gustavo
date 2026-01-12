# NOTE: Start by running the functions at the end of this file, if any. 

### LIBRARIES ##################################################################

library("biomaRt")
library("DBI")
library("hom.Hs.inp.db")
library('homologene')
library('VennDiagram')
library('dplyr')


# Install from CRAN
install.packages(c("biomaRt", "DBI", "homologene", "VennDiagram"))

# Install from Bioconductor (if needed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hom.Hs.inp.db")  # If this is a Bioconductor package

################################################################################
# HUMAN GENES 
################################################################################

AgeingHumans = read.csv('Data/Retrieved/Genes_and_diseases/HAGR/genage_human.csv')
AgeingHumans$symbol %>% unique() %>% saveRDS('Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds')

################################################################################
# MODEL GENES 
################################################################################

# RETRIEVE OMA ORTHOLOGS
# Homo.sapiens"            "Mus.Musculus"             "Relations"              "ID"
#"Homo.sapiens"            "Drosophila.Melanogaster" "Relations"               "ID" 
#"Homo.sapiens"           "Caenorhabditis.elegans" "Relations"              "ID"                    
#"Homo.sapiens"             "Saccharomyces.cerevisiae" "Relations"                "ID"

#Mouse
Oma_Mouse = read.csv('Data/Retrieved/Genes_and_diseases/OMA/OMA_Mouse.csv') %>% 
  dplyr::select(EntrezGene_2, EntrezGene_1, RelType, OmaGroup) %>%
  rename(Homo.sapiens=EntrezGene_2,
         Mus.Musculus=EntrezGene_1,
         Relations=RelType,
         ID=OmaGroup)

# Fly
Oma_Fly = read.csv('Data/Retrieved/Genes_and_diseases/OMA/OMA_Fly.csv') %>% 
  dplyr::select(EntrezGene_2, EntrezGene_1, RelType, OmaGroup) %>%
  rename(Homo.sapiens=EntrezGene_2,
         Drosophila.Melanogaster=EntrezGene_1,
         Relations=RelType,
         ID=OmaGroup)

# Worm
Oma_Worm = read.csv('Data/Retrieved/Genes_and_diseases/OMA/OMA_Worm.csv') %>% 
  dplyr::select(EntrezGene_2, EntrezGene_1, RelType, OmaGroup) %>%
  rename(Homo.sapiens=EntrezGene_2,
         Caenorhabditis.elegans=EntrezGene_1,
         Relations=RelType,
         ID=OmaGroup)

# Yeast
Oma_Yeast = read.csv('Data/Retrieved/Genes_and_diseases/OMA/OMA_Yeast.csv') %>% 
  dplyr::select(EntrezGene_2, EntrezGene_1, RelType, OmaGroup) %>%
  rename(Homo.sapiens=EntrezGene_2,
         Saccharomyces.cerevisiae=EntrezGene_1,
         Relations=RelType,
         ID=OmaGroup)

# AGEING AND DIETARY RESTRICTION
AgeingModels = read.csv('Data/Retrieved/Genes_and_diseases/HAGR/genage_models.csv')

##### LOAD GENAGE MODELS ##################################################################################################

aMusculus = AgeingModels[AgeingModels$organism == 'Mus musculus',]$symbol
aMusculusId = AgeingModels[AgeingModels$organism == 'Mus musculus',]$entrez.gene.id

aDrosophila  = AgeingModels[AgeingModels$organism == 'Drosophila melanogaster',]$symbol
aDrosophilaId  = AgeingModels[AgeingModels$organism == 'Drosophila melanogaster',]$entrez.gene.id

aElegans = AgeingModels[AgeingModels$organism == 'Caenorhabditis elegans',]$symbol
aElegansId = AgeingModels[AgeingModels$organism == 'Caenorhabditis elegans',]$entrez.gene.id

aCerevisiae  = AgeingModels[AgeingModels$organism == 'Saccharomyces cerevisiae',]$symbol
aCerevisiaeId  = AgeingModels[AgeingModels$organism == 'Saccharomyces cerevisiae',]$entrez.gene.id


##### MART MODELS ##################################################################################################

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
fly = useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
worm = useMart("ensembl", dataset = "celegans_gene_ensembl")
yeast = useMart("ensembl", dataset = "scerevisiae_gene_ensembl")


# MOUSE
ModelSymbol = 'mgi_symbol'
ModelMart = mouse
ModelName = 'Mus musculus'

eaHumanMouse = OmaFrame(RelationsDataset = Oma_Mouse, ModelSymbol, ModelGenes = aMusculus,
                        ModelIDs = aMusculusId, ModelMart, Database = 'entrezgene_id', ModelName)



# FLY
ModelSymbol = 'flybasename_gene'
ModelMart = fly
ModelName = 'Drosophila melanogaster'

eaHumanFly = OmaFrame(RelationsDataset = Oma_Fly, ModelSymbol, ModelGenes = aDrosophila,
                      ModelIDs = aDrosophilaId, ModelMart, Database = 'entrezgene_id', ModelName)


# NEMATODE
ModelSymbol = 'external_gene_name'
ModelMart = worm
ModelName = 'Caenorhabditis elegans'

eaHumanWorm = OmaFrame(RelationsDataset = Oma_Worm, ModelSymbol, ModelGenes = aElegans,
                       ModelIDs = aElegansId, ModelMart, Database = 'entrezgene_id', ModelName)


# YEAST
ModelSymbol = 'external_gene_name'
ModelMart = yeast
ModelName = 'Saccharomyces cerevisiae'

eaHumanYeast = OmaFrame(RelationsDataset = Oma_Yeast, ModelSymbol, ModelGenes = aCerevisiae,
                        ModelIDs = aCerevisiaeId, ModelMart, Database = 'entrezgene_id', ModelName)


##### COMBINE HUMAN GENES ###############################################################################################################

ModelAgeingGenes = unique(c(eaHumanMouse$Human, eaHumanFly$Human, eaHumanWorm$Human, eaHumanYeast$Human))

saveRDS(ModelAgeingGenes, "Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")


#########################################################################################################################################
##### FUNCTIONS #########################################################################################################################
#########################################################################################################################################


EntrezToName <- function(ModelGenes, ModelIDs, Dataset, ModelName){
  NameEntrez = data.frame(ModelGenes,ModelIDs)
  NameEntrez = NameEntrez[!duplicated(NameEntrez),]
  Dataset$ModelEntrez = as.character(Dataset$ModelEntrez)
  
  NameEntrezList = list()
  for (i in 1:dim(NameEntrez)[1]){
    NameEntrezList[[as.character(NameEntrez$ModelIDs[i])]] = NameEntrez$ModelGenes[i]
  }
  
  Model = c()
  Human = c()
  for (i in 1:(dim(Dataset)[1])){
    ds = Dataset[i,]
    Model[i] = NameEntrezList[[ds$ModelEntrez]]
    Human[i] = ds$HumanSymbol
  }
  Organism = rep(ModelName,length(Model))
  RelationsFrame = data.frame(Organism, Model, Human)
  return(RelationsFrame)
}


OmaFrame <- function(RelationsDataset, ModelSymbol, ModelGenes, ModelIDs, ModelMart, Database, ModelName){
  Oma = OmaConversion(RelationsDataset, ModelSymbol, ModelIDs, ModelMart, Database)
  RelationsFrame = EntrezToName(ModelGenes, ModelIDs, Dataset = Oma, ModelName)
  return(RelationsFrame)
}



# OMA CONVERSION
OmaConversion <- function(RelationsDataset, ModelSymbol, ModelIDs, ModelMart, Database){
  RelationsFrame = data.frame('HumanEns' = RelationsDataset$Homo.sapiens, 'ModelEns' = RelationsDataset[[2]])
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  Hframe <- getBM(filters= Database, attributes= c(Database,"hgnc_symbol","entrezgene_id"), values=RelationsFrame$HumanEns, mart= human)
  
  Oframe <- getBM(filters= Database, attributes= c(Database, ModelSymbol,"entrezgene_id"), values=RelationsFrame$ModelEns, mart= ModelMart)
  
  sRelationsFrame = RelationsFrame[RelationsFrame$ModelEns %in% Oframe[[Database]],]
  ModelHomologuesList = HumanModelGenesM(sRelationsFrame = sRelationsFrame, Hframe = Hframe, Oframe = Oframe, ModelSymbol = ModelSymbol, Database = Database)
  hModelInParanoid = ModelHomologuesList$sHumanModel[ModelHomologuesList$sHumanModel$ModelEntrez %in% ModelIDs,]
  hModel = hModelInParanoid[,c('ModelSymbol','HumanSymbol','ModelEntrez')]
  reps = dim(hModel)[1]
  Organism = rep(ModelName,reps)
  hModel$Organism = Organism
  hmodel = hModel[,c(3,1,2)]
  row.names(hmodel) = 1:reps
  return(hModel)
}



# OMA CONVERSION
OmaConversionProtr <- function(RelationsDataset, ModelSymbol, ModelIDs, ModelMart, Database){
  RelationsFrame = data.frame('HumanEns' = RelationsDataset$Homo.sapiens, 'ModelEns' = RelationsDataset[[2]])
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  Hframe <- getBM(filters= Database, attributes= c(Database,"hgnc_symbol","uniprotswissprot"), values=RelationsFrame$HumanEns, mart= human)
  Oframe <- getBM(filters= Database, attributes= c(Database, ModelSymbol,"uniprotswissprot"), values=RelationsFrame$ModelEns, mart= ModelMart)
  sRelationsFrame = RelationsFrame[RelationsFrame$ModelEns %in% Oframe[[Database]],]
  ModelHomologuesList = HumanModelGenesM(sRelationsFrame = sRelationsFrame, Hframe = Hframe, Oframe = Oframe, ModelSymbol = ModelSymbol, Database = Database)
  hModelInParanoid = ModelHomologuesList$sHumanModel[ModelHomologuesList$sHumanModel$ModelEntrez %in% ModelIDs,]
  hModel = hModelInParanoid[,c('ModelSymbol','HumanSymbol','ModelEntrez')]
  reps = dim(hModel)[1]
  Organism = rep(ModelName,reps)
  hModel$Organism = Organism
  hmodel = hModel[,c(3,1,2)]
  row.names(hmodel) = 1:reps
  return(hModel)
}



# HUMAN MODEL GENES
HumanModelGenesM <- function(sRelationsFrame, Hframe, Oframe, ModelSymbol, Database){
  
  SymbolArray = c()
  EnsArray = c()
  EntrezArray = c()
  dimension = dim(sRelationsFrame)[1]
  for (i in 1:dimension){
    print(paste('1-',paste(i, paste('/',dimension))))
    #print(i)
    Ens = sRelationsFrame[i,'ModelEns']
    SymbolArray[i] = Oframe[Oframe[[Database]] == Ens,][[ModelSymbol]]
    EnsArray[i] = Oframe[Oframe[[Database]] == Ens,][[Database]]
    EntrezArray[i] = Oframe[Oframe[[Database]] == Ens,]$entrezgene_id
  }
  rOframe = data.frame('ModelEns' = EnsArray, 'ModelSymbol' = SymbolArray, 'ModelEntrez' = EntrezArray)
  
  SymbolArray = c()
  EnsArray = c()
  
  for (i in 1:dimension){
    print(paste('2-',paste(i, paste('/',dimension))))
    Ens = sRelationsFrame[i,'HumanEns']
    if (length(Hframe[Hframe[[Database]] == Ens,]$hgnc_symbol) != 0){
      SymbolArray[i] = Hframe[Hframe[[Database]] == Ens,]$hgnc_symbol
      EnsArray[i] = Hframe[Hframe[[Database]] == Ens,][[Database]]
      EntrezArray[i] = Hframe[Hframe[[Database]] == Ens,]$entrezgene_id
    }
    else
    {
      SymbolArray[i] = '-'
      EnsArray[i] = '-'
      EntrezArray[i] = '-'
    }
  }
  rHframe = data.frame('HumanEntrez' = EntrezArray, 'HumanSymbol' = SymbolArray, 'HumanEns' = EnsArray)
  
  HumanModel = data.frame(rHframe, sRelationsFrame, rOframe)
  OverlapRows = rowSums(HumanModel[,] != '-') == dim(HumanModel)[2]
  sHumanModel = HumanModel[OverlapRows,]
  NonNArows = !grepl("NA", row.names(sHumanModel), fixed=TRUE)
  sHumanModel =  sHumanModel[NonNArows,]
  HumanModelHomologues = sHumanModel[c('HumanSymbol', 'ModelSymbol')]
  HomologuesList = list()
  HomologuesList$sHumanModel = sHumanModel
  HomologuesList$HumanModelHomologues = HumanModelHomologues
  
  return(HomologuesList)
}



# HUMAN MODEL RELATIONS
HumanModelRelations <- function(RelationsList){
  HumanEns = c()
  ModelEns = c()
  RelationsFrame = data.frame()
  for (i in 1:length(RelationsList)){
    RelationsElement = RelationsList[[i]]
    for (j in 1:length(RelationsElement)){
      HumanEns = c(HumanEns, names(RelationsList)[i])
      ModelEns = c(ModelEns, RelationsElement[j])
    }
  }
  Relations = list()
  Relations[['HumanEns']] = HumanEns
  Relations[['ModelEns']] = ModelEns
  return(Relations)
}



# HUMAN MODEL RELATIONS
ComparisonTable <- function(n,h){
  u = unique(c(n,h))
  A = u %in% n
  B = u %in% h
  DF = data.frame(A,B)
  row.names(DF) = u
  return(DF)
}