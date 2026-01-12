### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgeArr)
# Ang - Angle (e.g., TckAng)
# Ann - Annotation (e.g., FntSizAnn)
# Arr - Array (e.g., PhnArr)
# Box - Box (e.g., BoxPhnScrFrm)
# Cat - Category (e.g., CatVecCol)
# Col - Column (e.g., ColNam)
# Deg - Degree (e.g., TxtDeg)
# Fnt - Font (e.g., FntSiz)
# Frm - Frame (e.g., PhnMeanFrm)
# Gen - Gene (e.g., GenAge)
# Htm - Heatmap (e.g., Htm)
# Hum - Human (e.g., GenAge.Hum)
# Len - Length (e.g., FrmPhnLen)
# Lbl - Label (e.g., BoxPhnScrFrm$Lbl)
# Lst - List (e.g., MatPhnPhnLst)
# Mat - Matrix (e.g., MatPhnPhnLst)
# Mod - Moduel_Organisms (e.g., ModAge)
# Nam - Name (e.g., SrtNam)
# Phn - Phenotype (e.g., PhnArr)
# Ple - Pleiotropy (e.g., PleArr)
# Plt - Plot (e.g., PltFrm)
# Qry - Query (e.g., FrmPhnQry)
# Row - Row (e.g., Row_Phenotype)
# Scr - Score (e.g., PhnScrArr)
# Siz - Size (e.g., FntSiz)
# Spt - Split (e.g., SptVec)
# Srt - Sorted (e.g., SrtNam)
# Tau - Tau (e.g., Tau)
# Tck - Tick (e.g., TckAng)
# Vec - Vector (e.g., CatVecRow)

### LIBRARIES ##################################################################

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(circlize)

################################################################################
# RETRIEVE VAIBALES
################################################################################

MatPhnPhnLst = readRDS("Data/Generated/Specificity_and_coexpression/Specificity/Specificity_Matrices_List.rds") 
PhnMeanFrm = readRDS("Data/Generated/Specificity_and_coexpression/Specificity/Specificity_Mean_Frame.rds")

### GENE SETS ##################################################################

PhnArr = c("cardiovascular", "endocrine/diabetes", "gastrointestinal/abdominal", "renal/urology", 
           "haematology/dermatology", "immunological/systemic disorders", "musculoskeletal/trauma", "neurology/eye/psychiatry")

#PleArr = c("ARC-Pleiotropy_1", "ARC-Pleiotropy_2", "ARC-Pleiotropy_3", "ARC-Pleiotropy_4", "ARC-Pleiotropy_5", "ARC-Pleiotropy_6")
PleArr = c("Low ARC-Pleiotropy", "High ARC-Pleiotropy")

AgeArr = c("GenAge", "ModAge")

PhnPleAgeArr = c(PhnArr, PleArr, AgeArr)
PhnAgeArr = c(PhnArr, AgeArr)
PleAgeArr = c(PleArr, AgeArr)

################################################################################
### BOXPLOT ####################################################################
################################################################################

PhnPhnScrFrm = readRDS("Data/Generated/Specificity_and_coexpression/Specificity/Specificity_Pairwise_Frame.rds") %>%
  filter(Column_Phenotype == Row_Phenotype) %>% 
  mutate(Phn = Column_Phenotype) %>% 
  select(Column_Phenotype, Scr) %>% 
  rename(Phn=Column_Phenotype)

# TEXT ANGLE
TckAng = -30

# SORT COLUMNS BY MEDIAN
FrmPhnArr = PhnPhnScrFrm$Phn %>% unique()
FrmPhnLen = length(FrmPhnArr)
PhnScrArr = c()
for(i in 1:FrmPhnLen){
  FrmPhnQry = FrmPhnArr[i]
  PhnScrArr[i] = PhnPhnScrFrm %>% filter(Phn == FrmPhnQry) %>% pull(Scr) %>% median(na.rm=TRUE)
}
names(PhnScrArr) = FrmPhnArr
SrtNam = PhnScrArr %>% sort() %>% names() %>% rev()


### BOXPLOT ####################################################################

BoxPhnScrFrm = PhnPhnScrFrm %>% filter(Phn %in% PhnAgeArr)

BoxPhnScrFrm$Phn = ifelse(BoxPhnScrFrm$Phn == "GenAge", "GenAge.Hum", BoxPhnScrFrm$Phn)
BoxPhnScrFrm$Phn = ifelse(BoxPhnScrFrm$Phn == "ModAge", "GenAge.Mod", BoxPhnScrFrm$Phn)

# Define custom color palette
custom_palette <- c("GenAge" = '#F8776D', "immunological/systemic disorders" = "#5A8AC6", "Other_ARCs" = "#CAE2EE")

# Add a new column for custom coloring
BoxPhnScrFrm$color_group <- ifelse(BoxPhnScrFrm$Phn %in% c("GenAge.Hum", "GenAge.Mod"), "GenAge",
                                   ifelse(BoxPhnScrFrm$Phn == "immunological/systemic disorders", "immunological/systemic disorders", "Other_ARCs"))


# ACTUAL PLOT
BoxPhnScrFrm$Lbl = "Tissue specificity in genes associated with ARCs and GenAge"
# PCompute the counts for each box
counts <- as.data.frame(table(BoxPhnScrFrm$Phn))

SrtNam = ifelse(SrtNam == "GenAge", "GenAge.Hum", SrtNam)
SrtNam = ifelse(SrtNam == "ModAge", "GenAge.Mod", SrtNam)


# Create the plot
p = ggboxplot(BoxPhnScrFrm, x = "Phn", y = "Scr",
              outlier.size = 0.4,
              outlier.shape = NA,
              combine = TRUE,
              fill = "color_group",
              facet.by = "Lbl",
              order = SrtNam) + scale_fill_manual(values = custom_palette)

counts$FreqPhnma = format_with_commas(counts$Freq)

p = p+geom_dotplot(binaxis='y', stackdir='center',
                   position=position_dodge(1), dotsize=0.03, color="#202020")

p = p + geom_text(data = counts, aes(x = Var1, y = -0.0, label = FreqPhnma), size=3)

p = p + theme(axis.text.x = element_text(angle = TckAng, vjust = 0.6, hjust=0),
              strip.text = element_text(face="bold")) 

p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL)

p = ggpar(p, xlab ="ARC and GenAge sets", ylab = "Tau", legend.title = "") +
  font("title", size = 10, color = "black", face = "bold") +
  font("xlab", size = 10, color = "black", face = "bold") +
  font("ylab", size = 10, color = "black", face = "bold") +
  font("xy.text", size = 9, color = "black") +
  theme(legend.position = "none")

# Ensure legend is visible
p = p + theme(legend.position = "top")  

p = p + theme(plot.margin = unit(c(0, 1, 0, 0), "cm"))  # Adjust plot margins

# Modify the x-axis labels
p = p + scale_x_discrete(labels = c(
  "immunological/systemic disorders" = "immunological/systemic disorders\n(High ARC-Pleiotropy genes)", 
  "haematology/dermatology" = "haematology/dermatology", 
  "endocrine/diabetes" = "endocrine/diabetes",
  "renal/urology" = "renal/urology",
  "gastrointestinal/abdominal" = "gastrointestinal/abdominal",
  "cancer" = "cancer",
  "neurology/eye/psychiatry" = "neurology/eye/psychiatry",
  "cardiovascular" = "cardiovascular",
  "musculoskeletal/trauma" = "musculoskeletal/trauma",
  "GenAge.Hum" = "GenAge.Hum",
  "GenAge.Mod" = "GenAge.Mod"
))

p


BoxPhnScrFrm = BoxPhnScrFrm %>%
  rename(Phenotype=Phn,
         Score=Scr,
         Label=Lbl)

write.table(BoxPhnScrFrm,'Data/Generated/Specificity_and_coexpression/Specificity/BoxPlot_Data.csv', row.names = FALSE, sep=",")

################################################################################
### INTRA COMMUNITY DIFFERENCE (HEATMAP) #######################################
################################################################################

### INITIAL FRAMES AND PARAMETERS ##############################################

FntSiz = 8
TxtDeg = 60

PltFrm = MatPhnPhnLst$Dif
TxtFrm = MatPhnPhnLst$TstBnfTxt

RowNam = PltFrm$Phenotype

row.names(PltFrm) = PltFrm$Phenotype
row.names(TxtFrm) = TxtFrm$Phenotype

PltFrm$Phenotype = NULL
TxtFrm$Phenotype = NULL


# SORTING ROWS AND COLUMNS 
TxtFrm = TxtFrm[PhnPleAgeArr,PhnPleAgeArr]
PltFrm = PltFrm[PhnPleAgeArr,PhnPleAgeArr]

# GENE GROUPS LABELLING
RowNam = row.names(PltFrm)
SptVec = rep("D",length(RowNam))
SptVec[RowNam %in% c("GenAge","ModAge")] = "H"
SptVec[RowNam %in% PleArr] = "P"
PhnMeanFrm = PhnMeanFrm[row.names(PltFrm),] %>% data.frame()
colnames(PhnMeanFrm) = "Mean"
row.names(PhnMeanFrm) = row.names(PltFrm)

# TOP ANNOTATION OF HEATMAP ####################################################

TopAnn = HeatmapAnnotation(
  
  Mean_Specificity = anno_simple(PhnMeanFrm$Mean, 
                                  col = colorRamp2(c(min(PhnMeanFrm$Mean, na.rm=TRUE), 
                                                     ( min(PhnMeanFrm$Mean, na.rm=TRUE) + max(PhnMeanFrm$Mean, na.rm=TRUE) )/2, 
                                                     max(PhnMeanFrm$Mean, na.rm=TRUE)), 
                                                   c("#5A8AC6", "white", "#F8696B")),
                                  pch = as.character(PhnMeanFrm$Mean),
                                  gp= gpar(col = "black"),
                                  height = unit(8, "mm"),,
                                  pt_size = unit(1, "pt")*6
  ),
  annotation_name_gp= gpar(fontsize = 8),
  
  text = anno_text(replace_dirplt(colnames(PltFrm)), rot = TxtDeg, gp = gpar(fontsize = FntSiz))
)


# HEATMAP ######################################################################

FntSiz = 8
TxtDeg = 60
Htm = Heatmap(as.matrix(PltFrm), name = "Differential Specificity\nbetween Gene-sets:\nColumn.Set − Row.Set", 
              column_title = "Column Gene-set (ARC, GenAge and ARC-Pleiotropy)", 
              row_title = "Row Gene-set (ARC, GenAge and ARC-Pleiotropy)",
              column_title_gp = gpar(fontsize = 11), 
              row_title_gp = gpar(fontsize = 11),     
              show_row_names = FALSE, 
              show_column_names = FALSE,
              border = "black",
              row_split =  SptVec,
              column_split =  SptVec,
              rect_gp = grid::gpar(col = "grey", lwd = 0.5), 
              
              bottom_annotation = TopAnn,
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(TxtFrm[i, j], x, y, gp = gpar(fontsize = FntSiz))
              }
) +
  Heatmap(PhnMeanFrm$Mean,
          width = unit(9, "mm"),
          border = "black",
          rect_gp = grid::gpar(col = "black", lwd = 0.5), 
          name = "Gene-set mean Specificity",
          
          right_annotation = rowAnnotation(
            text = anno_text(replace_dirplt(row.names(PltFrm)), gp = gpar(fontsize = FntSiz))
          ),
          col = colorRamp2(c(min(PhnMeanFrm$Mean, na.rm=TRUE), 
                             ( min(PhnMeanFrm$Mean, na.rm=TRUE) + max(PhnMeanFrm$Mean, na.rm=TRUE) )/2, 
                             max(PhnMeanFrm$Mean, na.rm=TRUE)), 
                           c("#5A8AC6", "white", "#F8696B")),
          show_row_names = FALSE, show_column_names = FALSE,
          
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text( PhnMeanFrm[i, j], x, y, gp = gpar(fontsize = 7))
          },
          
          
  ) 


draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     show_heatmap_legend = FALSE, show_annotation_legend = FALSE)


################################################################################
### FUNCTION
################################################################################

replace_dirplt <- function(arr) {
  # Use sub() function to replace all occurrences of "DirPlt" with "dPleiotropy"
  arr[arr=="GenAge"] = "GenAge.Hum" 
  arr[arr=="ModAge"] = "GenAge.Mod" 
  #return(sub("Plt", "ARC-Pleiotropy_", arr))
  return(arr)
}

format_with_commas <- function(x) {
  # Check if x is a numeric vector or a matrix
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector or matrix.")
  }
  
  # Format each element in x with commas
  formatted <- format(x, big.mark = ",", scientific = FALSE, justify = "none")
  
  # Remove leading spaces by trimming the formatted strings
  formatted <- trimws(formatted)
  
  # Return the formatted vector
  return(formatted)
}
