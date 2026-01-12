### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgePhnHrc_Frm)
# All - All (e.g., GenAll)
# Arc - Ageing-related Community (e.g., GenArdArc)
# Ard - Ageing-related Disease (e.g., GenArdArc)
# Arr - Array (e.g., NtwArr)
# Cat - Category (e.g. GenCatFltCol)
# Cod - Code (e.g., Cod=Code)
# Col - Column (e.g. AllColNam) or Colour (e.g. DstColArr, AllFltCol)
# Dif - Difference (e.g. MaxDifFcn)
# Dis - Disease (e.g. DisDst)
# Dst - Distance (e.g., SrtDstNumDif)
# Edg - Edges (e.g., Edg=EdgesList)
# Fcn - Function (e.g. MaxDifFcn)
# Fct - Factor (e.g. AllColFct)
# Flt - Flattened (i.e., Matrix to vector)
# Frm - Frame (e.g., PrxQryFrm)
# Gen - Gene (e.g., GenAll)
# Grp - Graph(e.g., QryGrp)
# Hrc - Hierarchy (e.g., AgePhnHrc_Frm)
# Htm - Heatmap (e.g., HtmLst)
# Leg - Legend (e.g. SrtLegFcn)
# Len - Length (e.g., NtwLen)
# Lst - List (e.g., HtmLst)
# Mod - Model_Organisms (e.g., SepFct=Mod)
# Ngb - Neighbours (e.g. DisNgb)
# Nod - Node (e.g., Nod=NodesList)
# Ntw - Network (e.g., NtwArr)
# Oth - Others (e.g. OthGen)
# Phn - Phenotype (e.g., GenPhn)
# Prx - Proximity (e.g., PrxQryFrm)
# Qry - Query (e.g., NtwQry)
# Sep - Separison (e.g., SepFct)
# Srt - Sorted (e.g. SrtDstNumDif)
# Tmp - Temporal (e.g. TmpGenTypCol)
# Txt - Text (e.g., TxtPth)
# Typ - Type (e.g., GenTyp) 

### LIBRARIES ##################################################################

library(rje)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(igraph)

### START ######################################################################
  
# LOADING ARD-HIERARCHY
AgePhnHrc_Frm = read.csv("Data/Retrieved/Genes_and_diseases/Diseases/Ukb_DiseaseRoots.csv")


# - PLEIOTROPY -----------------------------------------------------------------

GenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)

ArcPltFrm = GenArdArcFrm %>% 
  dplyr::select(Arc,Gen) %>% 
  unique() %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
colnames(ArcPltFrm) = c('Gen','Frq')


# Outer pleitropy
GenHghPlt = ArcPltFrm %>% filter(Frq>=4) %>% pull(Gen) %>% as.character()
GenLowPlt = ArcPltFrm %>% filter(Frq<4) %>% pull(Gen) %>% as.character()


# ------------------------------------------------------------------------------

NtwArr = c("PPI", "COX90", "COX95", "KEGG")
NtwLen = length(NtwArr)
HtmLst =list()
x=4
for(x in 1:NtwLen){
  print(x)
  NtwQry = NtwArr[x]
  TxtPth = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/",sep="")
  
  Nod = readRDS(paste(TxtPth,"Lists/NodesList",".rds",sep=""))
  Edg = readRDS(paste(TxtPth,"Lists/EdgesList",".rds",sep=""))
  Grp = readRDS(paste(TxtPth,"Lists/GraphsList",".rds",sep=""))
  QryGrp = Grp$GenGenAll

  PrxQryFrm = readRDS(paste(TxtPth, "Proximity_and_Distances/Proximities_To_ARC_ARD.rds",sep=""))[["Closest.Proximity2ARC"]] 
  
  PrxQryFrm = PrxQryFrm[Nod$GenNtw,]
  
  Nod$GenAll %>% length()
  Nod$GenPhn %>% length()
  Nod$GenNtw %>% length()
  
  #  PLEIOTROPIES IN NETWORK
  HghGen = intersect(GenHghPlt,Nod$GenNtw)
  LowGen = intersect(GenLowPlt,Nod$GenNtw)
  
  # GET DISTANCE PLOTTING FRAME
  DstQryFrm = (1/PrxQryFrm)-1    
  QryFrm = DstQryFrm 
  MaxVal = QryFrm[QryFrm!=Inf] %>% max()
  QryFrm = abs(QryFrm-MaxVal)
  QryFrm[QryFrm==Inf] = -1
  
  AllFrm = QryFrm %>% as.matrix() %>% round() %>% t()   
  AllGen = colnames(AllFrm)
  
  # Genes of Diseases and GenAge
  DisGen = Edg$GenPhn %>% filter(!(Phn %in% c("HumAge","ModAge"))) %>% pull(Gen) %>% unique() %>% intersect(AllGen) 
  AgeGen = Edg$GenPhn %>% filter(Phn == "HumAge") %>% pull(Gen) %>% unique() %>% intersect(AllGen) 
  ModGen = Edg$GenPhn %>% filter(Phn == "ModAge") %>% pull(Gen) %>% unique() %>% intersect(AllGen)
  
  # GENE-GENE DISTANCES
  GenGenDstFrm = readRDS(paste(TxtPth, "Proximity_and_Distances/GenGen_DistanceFrame.rds",sep=""))
  GenGenDstFrm = GenGenDstFrm[AllGen, AllGen]
  
  # NEIGBHOURS
  DisNgb = QryGrp %>% ego(order = 1, nodes = DisGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  AgeNgb = QryGrp %>% ego(order = 1, nodes = AgeGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  ModNgb = QryGrp %>% ego(order = 1, nodes = ModGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  
  # GENE DISTANCES TO AGEING AND DISEASES
  AgeDst = GenGenDstFrm[AllGen,AgeGen] %>% rowMins()
  ModDst = GenGenDstFrm[AllGen,ModGen] %>% rowMins()
  DisDst = GenGenDstFrm[AllGen,DisGen] %>% rowMins()
  
  # GENE MEAN DISTANCE
  MeanPrx2Arc = rowMeans(PrxQryFrm, na.rm=TRUE)
  HtmMeanPrx2Arc = abs(MeanPrx2Arc-max(MeanPrx2Arc))
  MeanDst2Arc = rowMeans(DstQryFrm, na.rm=TRUE)
  
  # GENE ARC-INTERACTORS COUNT
  NumArcNgb = DstQryFrm %>% apply(1, function(x){sum(x==1)})
  
  ### ROW SPLITTING ################################################################
  
  # ALL AND EXLCUISE NEIGHBOURS
  AllNgb = c(DisNgb, AgeNgb) %>% unique()
  DisAgeNgbGen = DisNgb %>% c(AgeGen) %>% c(DisGen) %>% unique()
  OthGen = setdiff(AllGen, DisAgeNgbGen)

  
  # DIFFERENTIAL DATA
  AllColNam = AllFrm %>% colnames()
  DisGen_NonAge = setdiff(DisGen, AgeGen)
  DisGen_NonMod = setdiff(DisGen, ModGen)
  AgeGen_NonDis = setdiff(AgeGen, DisGen)
  ModGen_NonDis = setdiff(ModGen, DisGen)
  
  # FACTOR
  AllColFct = AllColNam
  # FACTOR FOR GenAgeHum
  AllColFct = ifelse(AllColFct %in% AgeGen, "Age", AllColFct)
  AllColFct = ifelse(AllColFct %in% DisGen_NonAge, "Dis", AllColFct)
  AllColFct = ifelse((AllColFct %in% AllNgb), "Ngb", AllColFct)
  AllColFct = ifelse(!(AllColFct %in% c("Dis", "Age", "Ngb")), "Oth", AllColFct)
  AllColFct = as.factor(AllColFct)

  # FACTOR FOR GenAgeMod
  AllColFct = ifelse(AllColFct %in% ModGen, "Mod", AllColFct)
  AllColFct = ifelse(AllColFct %in% DisGen_NonMod, "Dis", AllColFct)
  AllColFct = ifelse((AllColFct %in% AllNgb), "Ngb", AllColFct)
  AllColFct = ifelse(!(AllColFct %in% c("Dis", "Mod", "Ngb")), "Oth", AllColFct)
  AllColFct = as.factor(AllColFct)

  
  # GENE TYPE FRAME
  GenTypFrm = data.frame(Gen=AgeGen, GenTyp="AgeGen") %>% 
    rbind( data.frame(Gen=ModGen, GenTyp="ModGen") ) %>%
    rbind( data.frame(Gen=DisGen, GenTyp="DisGen") ) %>%
    rbind( data.frame(Gen=DisGen, GenTyp="PhnGen") ) %>%
    rbind( data.frame(Gen=AgeNgb, GenTyp="AgeNgb") ) %>%
    rbind( data.frame(Gen=ModNgb, GenTyp="ModNgb") ) %>%
    rbind( data.frame(Gen=DisNgb, GenTyp="DisNgb") ) %>%
    rbind( data.frame(Gen=DisNgb, GenTyp="PhnNgb") ) %>%
    rbind( data.frame(Gen=HghGen, GenTyp="HghGen") ) %>% 
    rbind( data.frame(Gen=LowGen, GenTyp="LowGen") ) %>% 
    rbind( data.frame(Gen=AllGen, GenTyp="AllGen") ) %>% 
    rbind( data.frame(Gen=OthGen, GenTyp="OthGen") ) 
    
  # CATEGORICAL FRAME
  GenCatFrm = data.frame("Gen"=AllColNam, "Cat"=AllColFct)
  
  ###  COLUMNS SPLITTING ###########################################################
  
  PhnNam = AllFrm %>% row.names()
  
  AllRowFct = PhnNam %>% as.factor
  AllPhnFct = AgePhnHrc_Frm %>% filter(RootMeaning %in% PhnNam) %>% pull(Meaning) %>% as.factor() 
  
  
  ## ANNOTATION SETS ##########################################################
  
  AgeFrm = AllFrm[,AgeGen] %>% as.matrix() %>% round() %>% t() 
  ModFrm = AllFrm[,ModGen] %>% as.matrix() %>% round() %>% t()  
  DisFrm = AllFrm[,DisGen] %>% as.matrix() %>% round() %>% t()
  NgbFrm = AllFrm[,DisNgb] %>% as.matrix() %>% round() %>% t()
  
  setdiff(DisNgb,colnames(AllFrm))
  
  ### GENE-GENE_LABEL FRAMES FOR PLOTTING #####################################
  
  # FRAMES FOR PLOTTING
  PltAgeFrm = AgeFrm
  PltModFrm = ModFrm
  PltDisFrm = DisFrm
  PltNgbFrm = NgbFrm
  
  PltAgeGenOrg = row.names(PltAgeFrm)
  PltModGenOrg = row.names(PltModFrm)
  PltDisGenOrg = row.names(PltDisFrm)
  PltNgbGenOrg = row.names(PltNgbFrm)
  
  row.names(PltAgeFrm) = paste("Age",1:nrow(PltAgeFrm),sep="")
  row.names(PltModFrm) = paste("Mod",1:nrow(PltModFrm),sep="")
  row.names(PltDisFrm) = paste("Dis",1:nrow(PltDisFrm),sep="")
  row.names(PltNgbFrm) = paste("Ngb",1:nrow(PltNgbFrm),sep="")
  
  PltAllFrm = PltAgeFrm %>% rbind(PltModFrm) %>% rbind(PltDisFrm) %>% rbind(PltNgbFrm) %>% t()
  colnames(PltAllFrm)[1:(nrow(PltAgeFrm)+nrow(PltModFrm)+1)]
  
  PltAllColFct = c(rep("Age",nrow(PltAgeFrm)), rep("Mod",nrow(PltModFrm)), rep("Dis",nrow(PltDisFrm)), rep("Ngb",nrow(PltNgbFrm))) %>% as.factor()
  
  PltAgeGenNew = row.names(PltAgeFrm)
  PltModGenNew = row.names(PltModFrm)
  PltDisGenNew = row.names(PltDisFrm)
  PltNgbGenNew = row.names(PltNgbFrm)
  
  # GENE AND LABEL RELATIONSHIP FRAMES
  PltAgeGenFrm = data.frame(Gen = PltAgeGenOrg, New = PltAgeGenNew)
  PltModGenFrm = data.frame(Gen = PltModGenOrg, New = PltModGenNew)
  PltDisGenFrm = data.frame(Gen = PltDisGenOrg, New = PltDisGenNew)
  PltNgbGenFrm = data.frame(Gen = PltNgbGenOrg, New = PltNgbGenNew)
  row.names(PltAgeGenFrm) = PltAgeGenOrg
  row.names(PltModGenFrm) = PltModGenOrg
  row.names(PltDisGenFrm) = PltDisGenOrg
  row.names(PltNgbGenFrm) = PltNgbGenOrg
  
  # BIND FRAMES FOR PLOTTING
  PltAllGenFrm = PltAgeGenFrm %>% rbind(PltModGenFrm) %>% rbind(PltDisGenFrm) %>% rbind(PltNgbGenFrm)
  
  ### GENE_LABEL FRAME WITH DISTANCES ###########################################
  
  PltAgeAllFrm = PltAllFrm[,PltAllColFct%in%"Age"]
  PltModAllFrm = PltAllFrm[,PltAllColFct%in%"Mod"]
  PltDisAllFrm = PltAllFrm[,PltAllColFct%in%"Dis"]
  PltNgbAllFrm = PltAllFrm[,PltAllColFct%in%"Ngb"]
  
  NumPltGen = ncol(PltAgeAllFrm)
  
  PltNewDisAllFrm = PltDisAllFrm[, sample(1:ncol(PltDisAllFrm), NumPltGen, replace = TRUE) ]
  PltNewAgeAllFrm = PltAgeAllFrm
  PltNewModAllFrm = PltModAllFrm[, sample(1:ncol(PltModAllFrm), NumPltGen, replace = TRUE) ]
  PltNewNgbAllFrm = PltNgbAllFrm[, sample(1:ncol(PltNgbAllFrm), NumPltGen, replace = TRUE) ]
  
  # COLOUR BY GENE TYPE
  GenTypCol = CatFrmColFcn(GenCatFrm, Pal="Pastel 1") 
  TmpGenTypCol = GenTypCol %>% unique()
  
  PltAllFrm = PltNewDisAllFrm %>% cbind(PltNewAgeAllFrm) %>% cbind(PltNewModAllFrm) %>% cbind(PltNewNgbAllFrm)
  PltAllColFct = c(rep("Dis",NumPltGen), rep("Age",NumPltGen), rep("Mod",NumPltGen), rep("Ngb",NumPltGen) )
  PltGenTypCol = c( rep(TmpGenTypCol[1],NumPltGen) , rep(TmpGenTypCol[2],NumPltGen) , rep(TmpGenTypCol[3],NumPltGen) )
  PltGen = colnames(PltAllFrm)
  
  ### LABELLING GENES BY GENE TYPE ###############################################
  
  PltAllGenFrm = PltAgeGenFrm %>% rbind(PltModGenFrm) %>% rbind(PltDisGenFrm) %>% rbind(PltNgbGenFrm)
  
  PltPrxDisMean = HtmMeanPrx2Arc[PltDisGenFrm$Gen]
  PltPrxAgeMean = HtmMeanPrx2Arc[PltAgeGenFrm$Gen]
  PltPrxModMean = HtmMeanPrx2Arc[PltModGenFrm$Gen]
  PltPrxNgbMean = HtmMeanPrx2Arc[PltNgbGenFrm$Gen]
  
  PltDstDisMean = MeanDst2Arc[PltDisGenFrm$Gen]
  PltDstAgeMean = MeanDst2Arc[PltAgeGenFrm$Gen]
  PltDstModMean = MeanDst2Arc[PltModGenFrm$Gen]
  PltDstNgbMean = MeanDst2Arc[PltNgbGenFrm$Gen]
  
  PltDisNgb = NumArcNgb[PltDisGenFrm$Gen]
  PltAgeNgb = NumArcNgb[PltAgeGenFrm$Gen]
  PltModNgb = NumArcNgb[PltModGenFrm$Gen]
  PltNgbNgb = NumArcNgb[PltNgbGenFrm$Gen]
  
  names(PltPrxDisMean) = PltDisGenFrm$New
  names(PltPrxAgeMean) = PltAgeGenFrm$New
  names(PltPrxModMean) = PltModGenFrm$New
  names(PltPrxNgbMean) = PltNgbGenFrm$New
  
  names(PltDstDisMean) = PltDisGenFrm$New
  names(PltDstAgeMean) = PltAgeGenFrm$New
  names(PltDstModMean) = PltModGenFrm$New
  names(PltDstNgbMean) = PltNgbGenFrm$New
  
  names(PltDisNgb) = PltDisGenFrm$New
  names(PltAgeNgb) = PltAgeGenFrm$New
  names(PltModNgb) = PltModGenFrm$New
  names(PltNgbNgb) = PltNgbGenFrm$New
  
  PltMeanPrx2Arc = c(PltPrxDisMean, PltPrxAgeMean, PltPrxModMean, PltPrxNgbMean) #PltMeanPrx2Arc
  PltMeanDst2Arc = c(PltDstDisMean, PltDstAgeMean, PltDstModMean, PltDstNgbMean) #PltMeanDst2Arc
  PltNumArcNgb = c(PltDisNgb, PltAgeNgb, PltModNgb, PltNgbNgb)
  
  PltNumArcNgb = PltNumArcNgb[PltGen]
  PltMeanPrx2Arc = PltMeanPrx2Arc[PltGen]
  PltMeanDst2Arc = PltMeanDst2Arc[PltGen]
  UseAllFrm = PltAllFrm
  
  ### COLOURS ####################################################################
  
  NseFrm = (PltAllFrm + runif(length(PltAllFrm), min = 0, max = 1e-10)) %>% t()
  AllDstFltCol = CatHtmColFcn(UseAllFrm, Pal1 = c("blue", "white"), Pal2 = c("white", "red"), MidNum = 2)
  AllDstFltCol[AllDstFltCol=="#000000"] = "#464646"
  
  # CATEGORICAL
  GenTypCol = CatFrmColFcn(GenCatFrm, Pal="Pastel 1")
  #ComTypCol = CatArrColFcn(AllRowFct, Pal="Dynamic")
  
  ### CONTINUOUS COLORS ##########################################################
  
  Col_Ngb = colorRamp2(c(min(PltNumArcNgb, na.rm=TRUE), max(PltNumArcNgb, na.rm=TRUE)/2, max(PltNumArcNgb, na.rm=TRUE)), c("blue", "white", "red"))
  Col_Mean = colorRamp2(c(min(PltMeanPrx2Arc, na.rm=TRUE), max(PltMeanPrx2Arc, na.rm=TRUE)/2, max(PltMeanPrx2Arc, na.rm=TRUE)), c("blue", "white", "red"))
  Col_Dst = colorRamp2(c(min(PltMeanPrx2Arc, na.rm=TRUE), max(PltMeanPrx2Arc, na.rm=TRUE)/2, max(PltMeanPrx2Arc, na.rm=TRUE)), c("blue", "white", "red"))
  
  
  ### DISTANCES ##################################################################
  
  PltMeanPrx2ArcDst = PltMeanDst2Arc
  NanPltMeanPrx2ArcDst = PltMeanPrx2ArcDst
  NanPltMeanPrx2ArcDst[NanPltMeanPrx2ArcDst==Inf] = NA
  MaxValMeanPrx2ArcDst = max(NanPltMeanPrx2ArcDst, na.rm=TRUE)
  MaxPltMeanPrx2ArcDst = NanPltMeanPrx2ArcDst
  MaxPltMeanPrx2ArcDst[is.na(MaxPltMeanPrx2ArcDst)] = MaxValMeanPrx2ArcDst
  
  ### HEATMAP TITLES #############################################################
  
  if(NtwQry=="PPI"){
    NtwQryTxt = "PPI-based"
  }
  if(NtwQry=="COX90"){
    NtwQryTxt = "COX.90-based"
  }
  if(NtwQry=="COX95"){
    NtwQryTxt = "COX.95-based"
  }
  if(NtwQry=="KEGG"){
    NtwQryTxt = "KEGG-based"
  }
  
  
  DisTxt = paste("Diseases\n",length(DisGen)," Genes",sep="")
  AgeTxt = paste("GenAge.Hum\n",length(AgeGen)," Genes",sep="")
  ModTxt = paste("GenAge.Mod\n",length(ModGen)," Genes",sep="")
  NgbTxt = paste("Neighbours\n",length(DisNgb)," Genes",sep="")
  
  
  ### HEATMAP ####################################################################
  
  PltAllColFctFct = factor(PltAllColFct, levels = c("Dis", "Mod", "Age", "Ngb"))
  
  cn = row.names(UseAllFrm)
  FntSiz = 8
  TxtDeg = 45
  PltAllFrm = as.matrix(UseAllFrm + runif(length(UseAllFrm), min = 0, max = 1e-10))
  
  Htm = Heatmap(t(UseAllFrm), name = "Distance", 
                column_title = paste("Distance to ARCs (",NtwQryTxt,")",sep=""), 
                
                row_title = c(ModTxt, AgeTxt, NgbTxt, DisTxt),
                show_row_names = FALSE, show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 9),
                row_split = PltAllColFctFct,
                border = TRUE,
                col = AllDstFltCol,
                
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(row.names(UseAllFrm), rot = TxtDeg, gp = gpar(fontsize = FntSiz), offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(AllPhnFct) * sin(TxtDeg*pi/180) * FntSiz/11
                ),
                
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 9),   
                  labels_gp = gpar(fontsize = 8),  
                  at = SrtLegFcn(AllDstFltCol) %>% as.numeric() %>% sort() %>% rev() %>% as.character(),
                  labels = MaxDifFcn(AllDstFltCol, MaxVal) %>% as.numeric() %>% sort() %>% as.character() %>% DirIndFcn() 
                )
                
                
                
                
  ) +
  
  
  Heatmap(PltNumArcNgb,
          width = unit(5, "mm"),
          name = "ARC-Interactions",
          border=TRUE,
          col = circlize::colorRamp2( c( min(PltNumArcNgb), ( max(PltNumArcNgb) + min(PltNumArcNgb) )/3 , max(PltNumArcNgb) ), c("#1E8AC6" , "white", "#F8696B") ),
          
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 9),   
            labels_gp = gpar(fontsize = 8)),  
          
          show_row_names = FALSE, show_column_names = FALSE,
          bottom_annotation = HeatmapAnnotation(
            text = anno_text("ARC-Interactions", rot = TxtDeg, gp = gpar(fontsize = 8), offset = unit(1, "npc"), just = "right"),
            annotation_height = max_text_width(cn) * 0.7
          )
  ) + 
    Heatmap(as.numeric(MaxPltMeanPrx2ArcDst),
            width = unit(5, "mm"),
            name = "Mean distance",
            border=TRUE,
            col = circlize::colorRamp2( c(min(MaxPltMeanPrx2ArcDst) , 2 , max(MaxPltMeanPrx2ArcDst) ), c("#F8696B" , "#FED280", "#1E8AC6") ) ,
            heatmap_legend_param = list(
              title_gp = gpar(fontsize = 9),   
              labels_gp = gpar(fontsize = 8)),  
            show_row_names = FALSE, show_column_names = FALSE,
            bottom_annotation = HeatmapAnnotation(
              text = anno_text("Mean distance", rot = TxtDeg, gp = gpar(fontsize = 8), offset = unit(1, "npc"), just = "right"),
              annotation_height = max_text_width(cn) * 0.7
            )
  )
  
  
  # DRAW HEATMAP
  HtmLst[[NtwQry]] = draw(Htm, padding = unit(c(-15, 2, 1.5, 1), "mm")) 

  
  
  ### TESTS ######################################################################
  
  GenTypNew = AllColFct
  
  GenAgeDst = AgeDst
  GenAgeDst[GenAgeDst==Inf] = NA
  GenAgeDst = GenAgeDst %>% as.numeric()
  
  GenModDst = ModDst
  GenModDst[GenModDst==Inf] = NA
  GenModDst = GenModDst %>% as.numeric()
  
  DisDst[DisDst==Inf] = NA
  DisDst = DisDst %>% as.numeric()
  
  # CREATE GENE COMPARISON FRAME
  sArcPrx_ArcInt_Frm = data.frame(Gen=names(NumArcNgb), 
                         AgeDst=GenAgeDst, ModDst=GenModDst, DisDst, 
                         NumArcNgb, 
                         MeanPrx2Arc, MeanDst2Arc,
                         Ntw=NtwQry) %>% merge(GenTypFrm)
  
  # MERGE GENE COMPARISON FRAME
  if(x == 1){
    ArcPrx_ArcInt_Frm = sArcPrx_ArcInt_Frm
  } else{
    ArcPrx_ArcInt_Frm = rbind(ArcPrx_ArcInt_Frm, sArcPrx_ArcInt_Frm)
  }

}

# RENAME GENE COMPARISON FRAME FOR SAVING
ArcPrx_ArcInt_Frm = ArcPrx_ArcInt_Frm %>% rename(Gene = Gen,
                        GeneType=GenTyp,
                        Distance_to_GenAgeHum=AgeDst,
                        Distance_to_GenAgeMod=ModDst,
                        Distance_to_ARCs=DisDst,
                        Number_Of_Neighbouring_ARCs=NumArcNgb,
                        MeanProximity_to_ARCs=MeanPrx2Arc,
                        MeanDistance_to_ARCs=MeanDst2Arc,
                        Network=Ntw
                        )


write.csv(ArcPrx_ArcInt_Frm,"Data/Generated/Networks_and_predictions/Topological_data/Distance_and_centrality/Proximity_and_Interactors.csv",row.names=FALSE)
saveRDS(HtmLst,"Data/Generated/Networks_and_predictions/Topological_data/Distance_and_centrality/ARC_Distance_Heatmaps_List.rds")

################################################################################
# FUNCTIONS
################################################################################

# COLOR FOR CATEGORICAL ARRAY
CatFrmColFcn = function(GenCatFrm, Pal="Pastel 1"){
  
  # UNIQUE DISTANCES AND LENGTH
  GenCatTbl = GenCatFrm$Cat %>% unique() %>% sort() %>% as.character()
  GenCatLen = length(GenCatTbl)
  
  # DEFINE PALETTE
  GenColArr = hcl.colors(GenCatLen, palette = Pal)
  
  # PALETTE-DST TABLE
  GenColFrm = data.frame('Col'=GenColArr, 'Cat'=GenCatTbl)
  row.names(GenColFrm) = GenColFrm$Cat 
  GenColNamArr = GenColFrm$Col
  names(GenColNamArr) = GenColFrm$Cat %>% as.character()
  
  # COL-DST FRAME
  GenCatColFrm = GenColFrm[as.character(GenCatFrm$Cat),]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  GenCatFltCol = GenCatColFrm$Col
  names(GenCatFltCol) = GenCatColFrm$Cat %>% as.vector() %>% as.character()
  
  return(GenCatFltCol)
  
}




# COLOR FOR CATEGORICAL ARRAY
CatArrColFcn = function(GenCatArr, Pal="Pastel 1"){
  
  # UNIQUE DISTANCES AND LENGTH
  GenCatTbl = GenCatArr %>% unique() %>% sort() %>% as.character()
  GenCatLen = length(GenCatTbl)
  
  # DEFINE PALETTE
  GenColArr = hcl.colors(GenCatLen, palette = Pal)
  
  # PALETTE-DST TABLE
  GenColFrm = data.frame('Col'=GenColArr, 'Cat'=GenCatTbl)
  row.names(GenColFrm) = GenColFrm$Cat
  GenColNamArr = GenColFrm$Col
  names(GenColNamArr) = GenColFrm$Cat %>% as.character()
  
  # COL-DST FRAME
  GenCatColFrm = GenColFrm[as.character(GenCatArr),]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  GenCatFltCol = GenCatColFrm$Col
  names(GenCatFltCol) = GenCatColFrm$Cat %>% as.vector() %>% as.character()
  
  return(GenCatFltCol)
  
}


CatHtmColFcn = function(AllFrm, Pal1 = c("green", "white"), Pal2 = c("white", "red"), MidNum = 3){
  
  MidNum = MidNum + 1
  
  # DEFINE PALETTE
  ColPal1 <- colorRampPalette(Pal1)
  ColPal2 <- colorRampPalette(Pal2)
  
  # CREATE DISTANCE TABLE
  AllDstTbl = AllFrm %>% as.matrix() %>% as.vector() %>% unique() %>% sort()
  AllDstLen = length(AllDstTbl)
  
  # CREATE COLOR FOR TABLE
  if(AllDstTbl[1] == -1){
    DstColArr = c( "#000000", ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1)[2:(MidNum) ])
  } else if(AllDstTbl[VecDstLen] == Inf){
    DstColArr = c( ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1)[2:(MidNum), "#000000" ])
  } else{
    DstColArr = c( ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1+1)[2:(MidNum+1)])
  }
  
  # PALETTE-DST TABLE
  DstColFrm = data.frame('Col'=DstColArr, 'Dst'=AllDstTbl)
  row.names(DstColFrm) = DstColFrm$Dst %>% round()
  DstColNamArr = DstColFrm$Col
  names(DstColNamArr) = DstColFrm$Dst %>% round() %>% as.character()
  
  # FLATTERN DATAFRAME
  AllFltCol = AllFrm %>% as.matrix() %>% as.vector() %>% round() %>% as.character()
  
  # COL-DST FRAME
  AllFltColFrm = DstColFrm[AllFltCol,]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  AllDstFltCol = AllFltColFrm$Col
  names(AllDstFltCol) = AllFltColFrm$Dst %>% round() %>% as.vector() %>% as.character()
  
  return(AllDstFltCol)
  
}


SrtLegFcn = function(DstFltCol){
  SrtNam = DstFltCol %>% names() %>% unique() %>% sort() %>% rev()
  return(SrtNam)
}


DirIndFcn = function(Dst){
  Dst[Dst=="0"] = "0 (Direct association)"
  Dst[Dst=="1"] = "1 (ARC-Interaction)"
  return(Dst) 
}


MaxDifFcn = function(Dst,MaxVal){
  SrtDst = Dst %>% names() %>% unique() %>% sort() %>% rev()
  SrtDstNum = SrtDst %>% as.numeric()
  SrtDstNum[SrtDstNum==-1]=Inf
  SrtDstNumDif = abs(SrtDstNum-MaxVal) %>% round() %>% as.character()
  return(SrtDstNumDif)
}
