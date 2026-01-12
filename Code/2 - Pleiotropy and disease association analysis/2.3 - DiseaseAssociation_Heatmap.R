# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Ageing (e.g., GenAgeHumArr)
# Arc - Ageing-related community (e.g., Arc=ARC_Meaning)
# Ard - Ageing-related disease (e.g., Ard=ARD_Meaning)
# Arr - Array (e.g., GenAgeHumArr)
# Ass - Association (e.g., ArcAssFrm)
# Cod - Code (e.g., Cod=Code)
# Cnt - Count (e.g., ArcArr)
# Dis - Disease (e.g., DisGenArr)
# Dst - Distance (e.g., HumArcDstColArr)
# Frm - Frame (e.g., GenArdArcFrm)
# Gen - Gene (e.g., GenArdArcFrm)
# Hrc - Hierarchy (e.g., ArdHrcFrm)
# Htm - Heatmap (e.g., HtmFrm)
# Hum - Human (e.g., AgeHumGen)
# Inf - Infinite (e.g., GenArcInf)
# Len - Length (e.g., DisGenLen)
# Mng - Meaning (e.g., RootMng)
# Mod - Model_Organisms (e.g., GenAgeModArr)
# Non - Non-Human (e.g., NonHumArc)
# Num - Number (e.g., ArcNumFrm)
# Qry - Query (e.g., GenQry)
# Row - Row (e.g., RowArcAssFrm)
# Yes - Yes (e.g., YesHumArc)

### LIBRARIES ##################################################################

library(dplyr)
library(reshape2)
library(dplyr)
library(igraph)
library(rje)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(devtools)

### RETRIEVE DATA #################################################

GenArdArcFrm = readRDS("Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARD_ARC.rds") %>% 
  dplyr::rename(Cod=Code, Ard=ARD_Meaning, Arc=ARC_Meaning, Gen=Gene)
GenArcAgeFrm = readRDS('Data/Retrieved/Genes_and_diseases/Diseases/GeneFrame_ARDcode_GenAge.rds') %>%
  dplyr::rename(Cod=Code, Gen=Gene)

GenAgeHumArr = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeHum_Genes.rds")
GenAgeModArr = readRDS("Data/Retrieved/Genes_and_diseases/HAGR/GenAgeMod_Genes.rds")

ArcArr = GenArdArcFrm$Arc %>% unique()

# GenAge Humans Association of ARCs
YesHumArc = GenArdArcFrm %>% filter(Gen %in% AgeHumGen) %>% pull(Arc) %>% unique()
NonHumArc = setdiff(ArcArr,YesHumArc)
YesHumArcFrm = data.frame(Dst = 0, Arc = YesHumArc)
NonHumArcFrm = data.frame(Dst = 1, Arc = NonHumArc)
HumArcFrm = rbind(YesHumArcFrm,NonHumArcFrm)
row.names(HumArcFrm) = HumArcFrm$Arc


# GenAge Models Association of ARCs
YesModArc = GenArdArcFrm %>% filter(Gen %in% AgeModGen) %>% pull(Arc) %>% unique()
NonModArc = setdiff(ArcArr,YesModArc)
YesModArcFrm = data.frame(Dst = 0, Arc = YesModArc)
NonModArcFrm = data.frame(Dst = 1, Arc = NonModArc)
ModArcFrm = rbind(YesModArcFrm,NonModArcFrm)
row.names(ModArcFrm) = ModArcFrm$Arc

# DISEASE HERARCHY FRAME
ArdHrcFrm = read.csv("Data/Retrieved/Genes_and_diseases/Diseases/Ukb_DiseaseRoots.csv") %>%
  dplyr::rename(Nod=Node, Cod=Coding, Mng=Meaning,RootMng=RootMeaning)

# FILTERING CLUSTERS 1 AND 2 (AGEING-RELATED DISEASES IN DONERTAS)
ArdHrcFrm = ArdHrcFrm %>% filter(Cluster %in% c(1,2))

# GENES ASSOCIATED SITH ARDs and ARCs
DisGenArr = GenArdArcFrm$Gen %>% unique()
DisGenLen = length(DisGenArr)

# NUMBER OF ASSOCIATED ARCs
ArcNumFrm = data.frame()
ArcArr = ArdHrcFrm %>% pull(RootMng) %>% unique()
ArcAssFrm = data.frame(Cnt=ArcArr)
row.names(ArcAssFrm) = ArcAssFrm$Cnt 
ArcAssFrm$Cnt = 0
ArcAssFrm = t(ArcAssFrm)
sArcAssFrm = ArcAssFrm
RowArcAssFrm = ArcAssFrm
for(i in 1:DisGenLen){
  sArcAssFrm = RowArcAssFrm
  print(i)
  GenQry = DisGenArr[i]
  ArcArr = GenArdArcFrm %>% filter(Gen == GenQry) %>% pull(Arc)
  NumArc = ArcArr %>% length()
  NumArc = ArcArr %>% unique() %>% length()
  ArcTbl = table(ArcArr)
  sArcAssFrm[,names(ArcTbl)] = as.numeric(ArcTbl)
  row.names(sArcAssFrm) = GenQry
  if(i == 1){
    ArcAssFrm = sArcAssFrm
  } else{
    ArcAssFrm = rbind(ArcAssFrm, sArcAssFrm)
  }
  sArcNumFrm = data.frame(Gen=GenQry, NumArc)
  ArcNumFrm = rbind(ArcNumFrm, sArcNumFrm)
  
}

# MATRIX OF GENES (Rows) ASSOCIATEION TO ARCs (Columns)
ArcAssFrm = data.frame(ArcAssFrm)
ArcGenNum = ArcAssFrm %>% colSums()
ArcGenNum = ArcGenNum[ArcGenNum>0] %>% names()
ArcAssFrm = ArcAssFrm[ArcGenNum]

ArcAssFrm$Gen = row.names(ArcAssFrm)
ArcAssNumFrm = merge(ArcAssFrm, ArcNumFrm, by="Gen") %>% arrange(desc(NumArc))

##########################################################################################################
### PLOT #################################################################################################
##########################################################################################################

QryFrm = ArcAssNumFrm[,!(colnames(ArcAssNumFrm) %in% c("NumArc"))]
row.names(QryFrm) = QryFrm$Gen
QryFrm$Gen = NULL
QryFrm = QryFrm[rowSums(QryFrm)>0,]
tQryFrm = t(QryFrm)


# COLUMN-SCORE FRAMES
FrmMng = ArcAssNumFrm$NumArc 
FrmRoot = ArcAssNumFrm$NumArc 

MaxMng = max(FrmMng)
MedMng = max(FrmMng)/2
MinMng = min(FrmMng)

MaxRoot = max(FrmRoot)
MedRoot = max(FrmRoot)/2
MinRoot = min(FrmRoot)

MaxFrm = max(QryFrm)
MinFrm = min(QryFrm)

# SORTING COLOUR PALETTES

WhtRedPal <- colorRampPalette(c("white", "red"))   # Apply colorRampPalette
WhtBluPal <- colorRampPalette(c("white", "blue"))   # Apply colorRampPalette
WhtBlkPal <- colorRampPalette(c("white", "black"))   # Apply colorRampPalette

WhtGrnPal <- colorRampPalette(c("green", "white", "orange"))   # Apply colorRampPalette

FltQryFrm = reshape(cbind(QryFrm, id=1, time=rownames(QryFrm)),direction="wide", sep="_", new.row.names=0)[-1]
FltQryVal = as.numeric(FltQryFrm)


RootColArr = colorRamp2(c(0, MaxRoot), c("white", "blue"))
MngColArr = colorRamp2(c(0, MedMng, MaxMng), c("green", "white", "orange"))
FrmColArr = colorRamp2(c(0, MaxFrm), c("white", "red"))

FrmCol = WhtRedPal(MaxFrm-MinFrm+1)
RootCol = WhtBluPal(MaxRoot-MinRoot+1+1)[(1+1):(MaxRoot-MinRoot+1+1)]

MngCol = WhtGrnPal(MaxMng-MinMng+1)
FrmColArr = FrmCol[FltQryVal-MinFrm+1]
RootColArr = RootCol[FrmRoot-MinRoot+1]
MngColArr = MngCol[FrmMng-MinMng+1]


names(FrmColArr) = FltQryVal
names(RootColArr) = FrmRoot
names(MngColArr) = FrmMng


RootLeg = FrmRoot %>% unique() %>% sort(decreasing = TRUE)
MngLeg = FrmMng %>% unique() %>% sort(decreasing = TRUE)
FrmLeg = FltQryVal %>% unique() %>% sort(decreasing = TRUE)

### COLUMNS ANNOTATION ###########################################################################

NewFrmRoot = FrmRoot
NewRootColArr = RootColArr
NewMngColArr = MngColArr
NewRootLeg = RootLeg
NewMngLeg = MngLeg


ha = HeatmapAnnotation(
  Total_ARCs = FrmRoot,
  Total_ARDs = FrmMng,
  
  col = list(Total_ARCs = RootColArr, Total_ARDs = MngColArr),
  
  annotation_height = unit(4, "mm"), 

  annotation_name_gp= gpar(fontsize = 8, fontface = "bold"),
  
  annotation_legend_param = list(


                                 Total_ARCs = list(title = "Total_ARCs",
                                                          at = RootLeg,
                                                          title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                          border = "black", 
                                                          labels_gp = gpar(col = "black", fontsize = 8), 
                                                          legend_height = unit(6, "cm"),
                                                          direction = "horizontal"
                                                          ),
                                 
                                 Total_ARDs = list(title = "Total_ARDs", 
                                                         at = MngLeg,
                                                         title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                         border = "black", 
                                                         labels_gp = gpar(col = "black", fontsize = 8), 
                                                         legend_height = unit(6, "cm"),
                                                         legend_direction = "horizontal"
                                                         )#,
                                 
                                 
                                 )
)


### ROW ANNOTATION ##############################################################################

# PREPARE GEN_AGE DST DATA
RowNam = row.names(tQryFrm)

RowNam = RowNam %>% str_replace("[.]", "/") %>% str_replace("eye.", "eye/") %>% str_replace("[.]", " ") %>% str_replace("genital/tract", "genital tract") # %>% str_replace("/cancer", " cancer")


### MODEL AGE ###################################################################################

HumArcArr = HumArcFrm[RowNam,]
ModArcArr = ModArcFrm[RowNam,]

# ADJUST GEN_AGE_ DST COLORS
ArcDstColNam = hcl.colors(2, palette = "Dynamic")#c("#000000", hcl.colors(ArcDstLen-1, palette = "Pastel 1"))

# ADJUST GEN_AGE_ DST COLORS
ArcGenSumArr = (tQryFrm > 0) %>% rowSums() %>% as.data.frame()
colnames(ArcGenSumArr) = "NumGen"


MaxScr = max(ArcGenSumArr$NumGen)
MedScr = MaxScr/2

ModPltArcDstArr = ifelse(ModArcArr$Dst == 0, "Ovelapping", "Non.Ovelapping")
ModArcDstColArr = ArcDstColNam[ModArcArr$Dst+1]
names(ModArcDstColArr) = ModArcArr
ModPltArcDstColArr = ModArcDstColArr
names(ModPltArcDstColArr) = ifelse(ModArcArr$Dst=="0","Ovelapping","Non.Ovelapping")


HumPltArcDstArr = ifelse(HumArcArr$Dst == 0, "Ovelapping", "Non.Ovelapping")
HumArcDstColArr = ArcDstColNam[HumArcArr$Dst+1]
names(HumArcDstColArr) = HumArcArr
HumPltArcDstColArr = HumArcDstColArr
names(HumPltArcDstColArr) = ifelse(HumArcArr$Dst=="0","Ovelapping","Non.Ovelapping")



# note how we set the width of this empty annotation
ra = rowAnnotation( Gene_Count = anno_barplot(ArcGenSumArr$NumGen),
                    GenAge.Hum = HumPltArcDstArr,
                    GenAge.Mod = ModPltArcDstArr,

                    col = list(GenAge.Hum = HumPltArcDstColArr, GenAge.Mod=ModPltArcDstColArr),
                    
                    annotation_name_gp= gpar(fontsize = 8, fontface = "bold"),
                    show_legend = c(TRUE, TRUE, FALSE),
                    annotation_legend_param = list(
                    GenAge.Hum = list(title = "GenAge (Hum & Mod)",
                                                              title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                              border = "black", 
                                                              labels_gp = gpar(col = "black", fontsize = 8), 
                                                              legend_height = unit(6, "cm"),
                                                              direction = "horizontal")
          ),
          annotation_name_rot = 90  

)




### HEATHMAP ######################################################################################

ncol(t(QryFrm))
cn = colnames(QryFrm)
rn = row.names(QryFrm)
LenRow = length(rn)
LenCol = length(cn)

length(FrmColArr) / 9
dim(tQryFrm)

Htm = Heatmap(as.matrix(tQryFrm), name = "Number of diseases", top_annotation = ha, 
              right_annotation = ra,

              col = FrmColArr,
              
              column_title = paste("Genes (",LenRow,")",sep=""), row_title = paste("ARCs (",LenCol,")",sep="" ),
              use_raster = FALSE, # FEDAULT IS TRUE

              show_column_names = FALSE, 
              border="black",

              row_title_gp=gpar(fontsize=10,fontface="bold"),
              column_title_gp=gpar(fontsize=10,fontface="bold"),
              row_names_gp = gpar(fontsize = 8),
              show_row_names=TRUE,
              
              heatmap_legend_param = list(
                title = "ARDs", 
                at = FrmLeg,
                title_gp = gpar(fontsize = 8, fontface = "bold"), 
                border = "black", 
                labels_gp = gpar(col = "black", fontsize = 8), 
                legend_height = unit(6, "cm"),
                
                direction = "horizontal",

                legend_direction = "horizontal"
                )
                                            
              

) 

### PRINT HEATMAP ##############################################################

draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)

### GEN AGE ASSOCIATIONS #######################################################

GenArdArcFrm %>% dplyr::select(Gen, Arc) %>% unique() %>% dplyr::filter(Gen %in% GenAgeHumArr) %>% pull(Arc) %>% table() %>% as.data.frame()
GenArdArcFrm %>% dplyr::select(Gen, Arc) %>% unique() %>% dplyr::filter(Gen %in% GenAgeModArr) %>% pull(Arc) %>% table() %>% as.data.frame()

