# DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# All - All (e.g., AllLen)
# ARC - Ageing-related Community (e.g., ARC_Pleiotropy)
# Arr - Array (e.g., NtwArr)
# Col - Color (e.g., ColFrm)
# Dir - Direct (e.g., DirIndOrg)
# Frm - Frame (e.g., ColFrm)
# Gen - Gene (e.g., NtwGen)
# Ind - Indirect (e.g., iARC_Interactor)
# Len - Length (e.g., NtwLen)
# Log - Logarithm (e.g., AllLog)
# Ntw - Network (e.g., NtwArr)
# Org - Original (e.g., DirIndOrg)
# Per - Percentage (e.g., AllPer)
# Phn - Phenotype (e.g., PhnPhn)
# Pth - Path (e.g., Pth = paste(...))
# Qry -  Query (e.g., NtwQry)
# Reg - Region of the plot (e.g., Reg = PhnPhnQry)
# Scatter -  (e.g., PieScatterFrm)
# Set - Set (e.g., MaxSet)


# LIBRARIES ######################################################################

library(dplyr)
library(ggplot2)
library(viridis)
library(ggtext)
library(cowplot)

# START ########################################################################

# LOGARITMIC PARAMETRS OF COLOURING BASED ON NUMBER OF GENES. 
numbers <- 1:2125
log_values <- log10(ifelse(numbers == 0, 1, numbers))
normalized <- (log_values - min(log_values)) / (max(log_values) - min(log_values))
color_array <- viridis(length(unique(normalized)))[as.integer(normalized * (length(unique(normalized)) - 1)) + 1]
ColFrm <- data.frame(AllLen = numbers, Col = color_array)


# LOOP OVER THE NETWORKS
NtwArr = c("PPI", "COX90", "COX95", "KEGG")
NtwLen = length(NtwArr)
pLst = list()
for(z in 1:NtwLen){
  
  print(z)
  NtwQry = NtwArr[z]
  
  Pth = paste("Data/Generated/Networks_and_predictions/Networks/",NtwQry,"/",sep="")
  
  # GGRAPH NODES
  Nod = paste(Pth,"Lists/NodesList.rds",sep="") %>% readRDS()
  NtwGen = Nod$GenNtw
  
  # LOADING PLEIOTROPY AND INDIRECT INTERACTORS DATA
  DirIndOrg = paste(Pth,"Proximity_and_Distances/Pleiotropy_and_IndirectInteractors.csv",sep="") %>% read.csv() # DirIndFrm
  DirIndOrg = DirIndOrg %>% filter(Gen %in% NtwGen)
  
  # INTEGRATING THE NUMBER OF ARC_Pleiotropy and iARC_Interactor INTO A SINGLE VARIABLE FOR COORDINATES
  DirIndOrg$PhnPhn = paste(DirIndOrg$ARC_Pleiotropy,DirIndOrg$iARC_Interactor,sep="") 
  PhnPhnArr = DirIndOrg$PhnPhn %>% unique()
  PhnPhnLen = length(PhnPhnArr)
  
  # COMPUTING CORRELATION
  cor_value <- cor(DirIndOrg$ARC_Pleiotropy, DirIndOrg$iARC_Interactor, method = "pearson")  # You can change the method if needed
  
  # GENERATING THE SCATTER FRAME
  PieScatterFrm = data.frame()
  GenLen = nrow(DirIndOrg)
  for(x in 1:PhnPhnLen){
    PhnPhnQry = PhnPhnArr[x]
    DirIndOrgQry = DirIndOrg %>% filter(PhnPhn == PhnPhnQry)
    AllLen = nrow(DirIndOrgQry)
    
    sPieScatterFrm = data.frame(Dir = unique(DirIndOrgQry$ARC_Pleiotropy), 
                            Ind = unique(DirIndOrgQry$iARC_Interactor), 
                            AllLen, 
                            AllLog = log10(AllLen),
                            Reg = PhnPhnQry
    )
    PieScatterFrm = rbind(PieScatterFrm,sPieScatterFrm)
    
  }
  
  PieScatterFrm$AllPer = PieScatterFrm$AllLen / GenLen

  
  ### PLOTTING THE SCATTERING FRAME #############################################
  
  MaxSet = max(PieScatterFrm$AllLog)
  PieScatterFrm$Rad = 0.5*PieScatterFrm$AllLog / MaxSet
  
  PieScatterFrm$G = NtwQry
  
  PieScatterFrm = merge(PieScatterFrm, ColFrm, by="AllLen")
  
  
  # PlOT
  p <- ggplot(PieScatterFrm, aes(x = Dir, y = Ind)) +
    geom_point(aes(size = 1, fill = AllLen), shape = 21, color = "black", stroke = 1) +
    scale_fill_gradientn(colors = PieScatterFrm$Col, trans="log10", name="Gene Count", # I CAN END HERE WITH ")" + PhnMENT THE FOLLOWING TWO LINES
                         guide = guide_colourbar(direction = "horizontal")) + theme(legend.position = "bottom") +
    guides(fill = "none")
  
  
  p = p +  scale_x_continuous(breaks = c(0,2,4,6)) +   # Ensure x-axis has ticks at every integer from 1 to 7
    scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
    facet_wrap(~ G) +
    theme(strip.text = element_text(face = "bold", size=10)) # This line makes the facet labels bold
  

  p = p + theme(
    axis.title.x = element_markdown(size = 10), 
    axis.title.y = element_markdown(size = 10)   
  )
  
  
  p = p + labs(x="**ARC-Pleiotropy**<br>(Directly connected ARCs)", y="**ARC-Interactions**<br>(Indirectly connected ARCs)")
  
  p = p + annotate("text", x = 6.5, y = 9, label = sprintf("Corr: %.2f", cor_value), hjust = 1, vjust=0, size=3.5)  # Adjust x and y for positioning
  
  
  if(z %in% c(2,4,6)){
    #p = p + theme(axis.title.y = element_blank())
    #p = p + theme(axis.title.y = "")
  }
  

  if(z %in% c(1,2)){
    #p = p + theme(axis.title.x = element_blank())
    #p = p + theme(axis.title.x = "")
  }
  
  p = p + coord_cartesian(ylim = c(-0.5,9.5))
  
  p = p + guides(size = guide_none())
  
  pLst[[NtwQry]] = p
  
}


p = plot_grid(pLst$PPI, pLst$COX90, pLst$COX95, 
          pLst$KEGG,
          labels = "AUTO", label_size = 12, ncol=2)
  

# Añadir título general
title <- ggdraw() + 
  draw_label(
    paste("ARC-Pleiotropy vs ARC-Interactors"),
    fontface = "bold", 
    size = 14,       # << tamaño de la fuente
    x = 0.5, hjust = 0.5
  )

final_plot <- plot_grid(title, p, ncol = 1, rel_heights = c(0.05, 1))

final_plot
