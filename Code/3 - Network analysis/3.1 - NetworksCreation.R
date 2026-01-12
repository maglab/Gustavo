##########################################################################################################
# GRAPH ##################################################################################################
##########################################################################################################

# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Edg - Edge (e.g., QryEdg, OrgQryEdg)
# Ens - Ensemble (e.g., EnsNam)
# Frm - Frame (e.g., GenPhnFrm, EdgFrm)
# Gen - Gene (e.g., GenPhnFrm)
# Grp - Graph (e.g., QryGrp)
# Hgn - HGNC (e.g., hgnc_symbol)
# Nam - Name (e.g., EnsNam)
# Nod - Node (e.g., Nod_A, Nod_B)
# Ntw - Network (e.g., NtwArr)
# Org - Original (e.g., OrgQry)
# Phn - Phenotype (i.e., GenAge.Hum, GenAge.Mod or ARDS) (e.g., GenPhnFrm)
# Qry - Query (e.g., QryGrp, QryEdg, NtwQry)

################################################################################

library(rlang)
library(dplyr)
library(igraph)
library(ggrepel)
library(data.table)

####################################################################################################
### START ##########################################################################################
####################################################################################################

# OVERALL PARAMETERS ############################################################################### 

GenPhnFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARDcode_GenAge.rds") %>%
  dplyr::rename(Phn=Code, Gen=Gene)

GenPhnFrm$Phn %>%table()

NtwArr = c('PPI','COX90','COX95', 'KEGG')
#NtwArr = c('KEGG')#,'COX90','COX95', 'KEGG')
NtwLen = length(NtwArr)

i=1
for(i in 1:NtwLen){
  
  print(i)

  NtwQry = NtwArr[i]
  
  ###########################################################################################################
  #### CREATE GRAPH FROM NETWORK DATASETS 
  ###########################################################################################################
  
  if(NtwQry == "PPI"){

    #OldEdgFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Networks/BIOGRID-MV-Physical-4.4.204.tab3.txt", sep="\t") # 04/12 - 164MB
    
    EdgFrm = read.csv("Data/Retrieved/Network_Sources/PPI/BIOGRID-MV-Physical-4.4.249.tab3.txt", sep="\t") # 04/12 - 164MB
    
    
    
    HumEdgFrm = EdgFrm%>%dplyr::filter(Organism.ID.Interactor.A==9606,Organism.ID.Interactor.B==9606)
    
    QryGrp = graph.data.frame(data.frame("Nod_A" = HumEdgFrm$Official.Symbol.Interactor.A, 
                                         "Nod_B" = HumEdgFrm$Official.Symbol.Interactor.B))
    QryGrp <- delete_vertices(QryGrp, V(QryGrp)[name == "NA"])
    QryGrp = as.undirected(QryGrp)
    QryGrp = igraph::simplify(QryGrp)
    
    QryEdg = as.data.frame(get.edgelist(QryGrp))
    colnames(QryEdg) = c("Nod_A", "Nod_B")
    
  } 
  if(NtwQry == "COX90"){
    
    EdgFrm = read.table("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/MatrixCox90_Ensembl.csv")
    EdgFrm = EdgFrm %>% dplyr::filter(row!=col)
    
    OrgEnsNam = readRDS("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/GenesCox90.rds")
    
    EnsNam = OrgEnsNam[!duplicated(OrgEnsNam$ensembl_gene_id),]
    row.names(EnsNam) = EnsNam$ensembl_gene_id
    
    QryGen_A = EnsNam[EdgFrm$row,]$hgnc_symbol
    QryGen_B = EnsNam[EdgFrm$col,]$hgnc_symbol
    
    OrgQryEdg = data.frame(QryGen_A, QryGen_B)
    
    OrgQryEdg = OrgQryEdg[OrgQryEdg$QryGen_A != '', ]
    OrgQryEdg = OrgQryEdg[OrgQryEdg$QryGen_B != '', ]
    OrgQryEdg = OrgQryEdg[!is.na(OrgQryEdg$QryGen_A), ]
    OrgQryEdg = OrgQryEdg[!is.na(OrgQryEdg$QryGen_B), ]

    
    QryGrp = graph.data.frame(data.frame("Nod_A" = OrgQryEdg$QryGen_A, 
                                         "Nod_B" = OrgQryEdg$QryGen_B))
    QryGrp <- delete_vertices(QryGrp, V(QryGrp)[name == "NA"])
    QryGrp = as.undirected(QryGrp)
    QryGrp = igraph::simplify(QryGrp)
    

    QryEdg = as.data.frame(get.edgelist(QryGrp))
    colnames(QryEdg) = c("Nod_A", "Nod_B")
    
    dim(QryEdg)
    
    
  }
  if(NtwQry == "COX95"){
    
    EdgFrm = read.table("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/MatrixCox95_Ensembl.csv")
    EdgFrm = EdgFrm %>% dplyr::filter(row!=col)
    
    OrgEnsNam = readRDS("Data/Retrieved/Network_Sources/Coexpression/Genes_and_edges/GenesCox95.rds")
    
    EnsNam = OrgEnsNam[!duplicated(OrgEnsNam$ensembl_gene_id),]
    row.names(EnsNam) = EnsNam$ensembl_gene_id
    
    QryGen_A = EnsNam[EdgFrm$row,]$hgnc_symbol
    QryGen_B = EnsNam[EdgFrm$col,]$hgnc_symbol
    
    OrgQryEdg = data.frame(QryGen_A, QryGen_B)
    
    OrgQryEdg = OrgQryEdg[OrgQryEdg$QryGen_A != '', ]
    OrgQryEdg = OrgQryEdg[OrgQryEdg$QryGen_B != '', ]
    OrgQryEdg = OrgQryEdg[!is.na(OrgQryEdg$QryGen_A), ]
    OrgQryEdg = OrgQryEdg[!is.na(OrgQryEdg$QryGen_B), ]
    
    
    QryGrp = graph.data.frame(data.frame("Nod_A" = OrgQryEdg$QryGen_A, 
                                         "Nod_B" = OrgQryEdg$QryGen_B))
    QryGrp <- delete_vertices(QryGrp, V(QryGrp)[name == "NA"])
    QryGrp = as.undirected(QryGrp)
    QryGrp = igraph::simplify(QryGrp)
    
    QryEdg = as.data.frame(get.edgelist(QryGrp))
    colnames(QryEdg) = c("Nod_A", "Nod_B")
  }
  
  
  if(NtwQry == "KEGG"){
    
    EdgFrm = readRDS("Data/Retrieved/Network_Sources/KEGG/Other/KEGG_Interactions.rds") 
    colnames(EdgFrm) = c("Nod_A", "Nod_B")
    EdgFrm = EdgFrm %>% dplyr::filter(Nod_A!=Nod_B)
    
    QryGrp = graph.data.frame(data.frame("Nod_A" = EdgFrm$Nod_A, 
                                         "Nod_B" = EdgFrm$Nod_B))
    QryGrp <- delete_vertices(QryGrp, V(QryGrp)[name == "NA"])
    
    QryGrp = as.undirected(QryGrp)
    QryGrp = igraph::simplify(QryGrp)
    
    
    QryEdg = as.data.frame(get.edgelist(QryGrp))
    colnames(QryEdg) = c("Nod_A", "Nod_B")
    
    
  } 
  
  
  ##############################################################################
  # INTEGRATING NETWORK AND GENE-PHENOTYPE
  ##############################################################################
  
  Nod = list()
  
  Nod$Phn = GenPhnFrm$Phn   # Phenotypes
  
  ### EDGES ##################################################################################################
  
  Edg = list()
  
  AllEdgFrm = GenPhnFrm
  colnames(AllEdgFrm) = c("Nod_A", "Nod_B")
  
  Edg$All = rbind(AllEdgFrm,QryEdg)  # Gene-gene and gene-phenotype edges
  Edg$GenPhn = GenPhnFrm  # Gene-Phenotype edges
  Edg$GenGenNtw = QryEdg  # Gene-Gene edges at the network
  
  dim(Edg$GenGenNtw)
  
  unique(c(Edg$GenGenNtw$Nod_A, Edg$GenGenNtw$Nod_B)) %>% length()
  
  ### MORE ON NODES ##########################################################################################
  
  Nod$GenPhn = Edg$GenPhn$Gen %>% unique()  # Phenotype-related genes
  
  Nod$GenNtw = c(Edg$GenGenNtw$Nod_A, Edg$GenGen$Nod_B) %>% unique()  # Network-related genes

  Nod$GenAll = c(Nod$GenPhn, Nod$GenNtw) %>% unique()  # Phenotype- and network-related genes
  
  Nod$All = c(Nod$Phn,Nod$GenAll) %>% unique()   # Phenotypes and Phenotype- and network-related genes
  
  length(Nod$GenNtw)
  length(Nod$GenAll)
  
  ### GRAPHS ##################################################################################################
  
  Grp = list()

  Grp$All = Edg$All %>% graph_from_data_frame(directed = FALSE)  # Network- and phenotype-related genes + phenotype nodes
  Grp$GenGenAll = Grp$All %>% delete_vertices(Nod$Phn) # Genes of the network and phenotype-related genes
  
  Grp$GenGenNtw = Edg$GenGenNtw %>% graph_from_data_frame(directed = FALSE) # Network-related genes
  
  ##############################################################################
  # SAVING 
  ##############################################################################
  
  saveRDS(Nod,paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/Lists/NodesList.rds", sep=''))
  saveRDS(Edg,paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/Lists/EdgesList.rds", sep=''))
  saveRDS(Grp,paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/Lists/GraphsList.rds", sep=''))
 
}


