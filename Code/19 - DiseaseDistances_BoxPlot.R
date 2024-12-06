### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgeDst)
# Alp (e.g., Alp=1)
# Ann (e.g., AnnTxtSiz)
# Arc - Ageing-related Community (e.g., ArcPrx_ArcInt_Frm)
# Dis - Disease (e.g., DisDst)
# Dst - Distance (e.g., AgeDst)
# Gen - Gene (e.g., GenTyp)
# Grp - Graph (e.g., GrpLbl)
# Hum - Human (e.g., GenAgeHum)
# Int - Interactor (e.g., ArcPrx_ArcInt_Frm)
# Jmp - Jump (e.g., YtstJmp)
# Lab - Label (e.g., xTckLab)
# Len - Length (e.g., GenSetLen)
# Max - Max (e.g., yLimMax)
# Min - Min (e.g., yLimMin)
# Mod - Model_Organisms (e.g., GenAgeMod)
# Ngb - Neighbours (e.g., DisNgb)
# Ntw - Network(e.g., NtwArr)
# Oth - Other (e.g., OthGen)
# Prx - Proximity (e.g., MeanPrx2Arc)
# Pwc - Pairwise test
# Scr - Score (e.g., Yscr)
# Sep - Separison (e.g., TstSep)
# Set - Set (e.g., GenSet)
# Siz - Size (e.g., TstTxtSiz)
# Txt - Text (e.g., TstTxtSiz)
# Tck - Tick (e.g., xTckLab)
# Tst - Test (e.g., WilTst)
# Typ - Type (e.g., GenTyp)
# Var - Variable (e.g., Xvar)
# Wil - Wilcoxon (e.g., WilTst)

### LIBRARIES ##################################################################

library(rstatix)
library(ggpubr)
library(dplyr)
library(rje)
library(cowplot)

# START ########################################################################

# LOADING FRAME ON GENE ARC-PROXIMITIES AND ARC-INTERACTORS
ArcPrx_ArcInt_Frm = read.csv("maglab/Gustavo/Data/Generated/Proximity_and_Interactors/Proximity_and_Interactors.csv")
  
ArcPrx_ArcInt_Frm = ArcPrx_ArcInt_Frm %>% 
  rename(Gen = Gene,
  GenTyp=GeneType,
  AgeDst=Distance_to_GenAgeHum,
  ModDst=Distance_to_GenAgeMod,
  DisDst=Distance_to_ARCs,
  NumArcNgb=Number_Of_Neighbouring_ARCs,
  MeanPrx2Arc=MeanProximity_to_ARCs,
  MeanDst2Arc=MeanDistance_to_ARCs,
  Ntw = Network)


### Pairwise comparisons #######################################################

# QUERY VARIABLES
Xvar = "GenTyp" 
Yvar = "MeanPrx2Arc"


# - SORTED GROUPS --------------------------------------------------------------
GenSet = c("DisGen", "AgeGen", "ModGen", "DisNgb", "OthGen") 

xTckLab = c("Diseases", "GenAge.Hum", "GenAge.Mod", "Neighbours", "Others")
GenSetLen = length(GenSet)
GenSetGr1 = data.frame(group1=GenSet, xmin=1:GenSetLen)
row.names(GenSetGr1) = GenSetGr1$group1
GenSetGr2 = data.frame(group2=GenSet, xmax=1:GenSetLen)
row.names(GenSetGr2) = GenSetGr2$group2


### INITIAL VARIABLES -- DisGen, AgeGen, aExcNgb -- GenNgbb #####################

WilTst = FALSE
Ymin = 1 
yLimMin = rep(-0.04,6)
yLimMax = rep(1.17, 6) 

YtstMin = rep(0.80,6) 
YtstJmp = rep(0.062, 6) 
YannMin = rep(-0.155,6)

# Parameters for plotting
TstTxtSiz = 2.5 
AnnTxtSiz = 3
DstGrp = FALSE 
AnnTxtExt = "" 
YtstBia = rep(0,6)
TstSep = 0.1

Alp = 1
TckAng = 30

Yscr = "Closest.Proximity to ARCs"

### PLOTTING LOOP ##############################################################

NtwArr = c("PPI", "COX90", "COX95", "KEGG")
GrpLbl = NtwArr#c("PPI", "COX.95", "COX.90", "KEGG")

MeanMedianFrm = data.frame()
x=1
NtwLen = length(NtwArr)
P = list()
for(x in 1:NtwLen){
  
  print(x)
  NtwQry = NtwArr[x]
  ArcPrxFrm = ArcPrx_ArcInt_Frm %>% filter(GenTyp %in% GenSet) %>% filter(Ntw == NtwQry)
  
  # DISTANCES MEANS
  grouped_means <-  ArcPrxFrm %>%
    group_by(GenTyp) %>%
    summarise(mean_Y = mean(MeanPrx2Arc, na.rm = TRUE))
  
  # DISTANCES MEDIANS
  grouped_medians <-  ArcPrxFrm %>%
    group_by(GenTyp) %>%
    summarise(mean_Y = median(MeanPrx2Arc, na.rm = TRUE))
  
  # GROUPING MEANS AND DISTANCES
  sMeanMedianFrm = grouped_means %>% merge(grouped_medians, by='GenTyp')
  sMeanMedianFrm$Network = NtwQry
  MeanMedianFrm = rbind(MeanMedianFrm,sMeanMedianFrm)
  
  # NUMBER OF ELEMENTS
  AgeNum = ArcPrxFrm %>% filter(GenTyp == "AgeGen") %>% nrow()
  ModNum = ArcPrxFrm %>% filter(GenTyp == "ModGen") %>% nrow()
  DisNum = ArcPrxFrm %>% filter(GenTyp == "DisGen") %>% nrow()
  NgbNum = ArcPrxFrm %>% filter(GenTyp == "DisNgb") %>% nrow()
  BeyNum = ArcPrxFrm %>% filter(GenTyp == "OthGen") %>% nrow()
  NumArr = c(DisNum, AgeNum, ModNum, NgbNum, BeyNum)
  
  
  # FILTER DESIRED FRAME
  TstFrm = ArcPrxFrm %>% select(.data[[Xvar]], .data[[Yvar]], Ntw)
  
  colnames(TstFrm) = c("X", "Y", "G")
  TstFrm$Y = TstFrm$Y 
  TstFrm = TstFrm %>% filter(!is.na(X))
  TstFrm$G = GrpLbl[x]
  
  pwc <- TstFrm %>% pairwise_t_test( Y ~ X , p.adjust.method = "bonferroni")
  
  GrpLen_Pwc = c(pwc$group1, pwc$group2) %>% unique() %>% length()
  GrpLen_Org = TstFrm$X %>% unique() %>% length()
  
  yMaxLoc = max(TstFrm$Y,na.rm=TRUE)
  
  # IF NOT ALL PWC ROWS ARE PRESENT WHERE THEY SHOULD BE
  if(GrpLen_Pwc < GrpLen_Org){
    
    Xarr = TstFrm$X %>% unique() %>% sort()
    Xprw <- consecutive_pairs(Xarr)
    
    # COMBINATORY
    PrwFrm = abs (combn(Xarr, 2) - max(Xarr) ) %>% as.data.frame()
    SrtPrwFm = PrwFrm[length(PrwFrm):1]
    Xprw = apply(SrtPrwFm, 2, list)
    
    Xlen = length(Xprw)
    pwc = data.frame()
    GrpMaxArr = c()
    
    # COMPUTE ROWWISE
    for(i in 1:Xlen){
      QryGrpArr = Xprw[[i]] %>% unlist() %>% rev()
      sTstFrm = TstFrm %>% filter(X %in% QryGrpArr)
      GrpMaxArr[i] = sTstFrm %>% pull(Y) %>% max()
      spwc = sTstFrm %>% pairwise_t_test( Y ~ X , p.adjust.method = "bonferroni")
      # IF ROW DOESNT CONTAIN ANYTING
      if(nrow(spwc) == 0){
        npwc = pwc[i-1,]
        GrpLenArr = sTstFrm$X %>% table() %>% as.numeric()
        spwc = npwc %>% mutate(group1 = QryGrpArr[1], group2 = QryGrpArr[2], n1 = GrpLenArr[1], n2 = GrpLenArr[2], p=1, p.signif="ns", p.adj=1, p.adj.signif="ns")
      }
      pwc = rbind(pwc, spwc)
    }
    
  }
  
  
  # WILCOX TEST IF NEEDED
  CapTxt = get_pwc_label(pwc)
  if(Yvar %in% c("GenNgb0", "GenNgb1", "GenNgb10") & WilTst==TRUE){
    pwc2 <- pairwise.wilcox.test( TstFrm$Y , TstFrm$X , p.adjust.method = "bonferroni") #%>% tidy() #%>% mutate(p.adj=p.value, P=p.value, `.y.`="Y") 
    if(nrow(pwc)==0){
      pwc = pwc2
    }
    pwc$p = pwc2$p.value
    pwc$p.adj = pwc2$p.value %>% signif(2)
    CapTxt = "pwc: Wilcox test; p.adjust: Bonferroni"
    Yscr = "Indirect pleitropy"
  }
  
  # - X,Y POSITIONS OF STAT SYMBOLS --------------------------------------------
  pwc = pwc %>% add_xy_position(x = "X")
  PrePwc = pwc
  NaPwx = sum(is.na(PrePwc)) > 0
  
  # IF NA DATA IN X,Y POSITIONS
  if(NaPwx){
    PrePwc$y.position = GrpMaxArr
    PwcLen = nrow(PrePwc)
    for(i in 1:PwcLen){
      PrePwc[i,][['groups']][[1]] = Xprw[[i]] %>% unlist() %>% rev() %>% as.character()
    }
    names(PrePwc$groups) = paste("V",1:PwcLen,sep="")
    
  }
  
  
  PltPwc = PrePwc 
    
  PltPwc$xmin = GenSetGr1[PltPwc$group1,"xmin"]
  PltPwc$xmax = GenSetGr2[PltPwc$group2,"xmax"]
  PltPwc <- PltPwc[order(PltPwc$xmin, PltPwc$xmax), ]
  PltPwc$y.position = sort(PltPwc$y.position)
  
  abs(PltPwc$xmin - PltPwc$xmax)
  
  MinArr = PltPwc$xmin
  MaxArr = PltPwc$xmax
  MinMaxFrm = data.frame(MinArr, MaxArr)
  PltPwc$xmax = rowMaxs(MinMaxFrm)
  PltPwc$xmin = rowMins(MinMaxFrm)
  
  PltPwc$y.position %>% unique()
  
  yPosArr = sort(PltPwc$y.position)
  FstPltPwc = PltPwc %>% filter(y.position == yPosArr[1]) 
  yPosMin = yPosArr[1]
  
  if(FstPltPwc$xmax-FstPltPwc$xmin != 1){
    NewPos = PltPwc %>% filter(xmax == 1+xmin) %>% pull(y.position) %>% min()
    PltPwc[PltPwc$y.position==yPosMin, "y.position"] = NewPos
  }
  
  
  PltPwc = PltPwc %>% mutate(y.position = ifelse(abs(xmin-xmax)==1,yPosMin,y.position) )
  
  PltPwc = PltPwc %>% mutate(xmin = xmin+0.1, xmax = xmax-0.1)
  PltPwc = PltPwc[order(PltPwc$y.position), ]
  YposTbl = PltPwc$y.position %>% table() %>% names()
  YposTblLen = length(YposTbl)
  
  MinPnt = YtstMin[x]

  for(j in 1:YposTblLen){
    PltPwc[PltPwc$y.position == YposTbl[j], "y.position"] = MinPnt + ((j-1)*YtstJmp[x])
  }


  # STATS ASTHERISKS
  PsgnArr = c("ns", "ns", "ns")
  PsgnArr[PltPwc$p > 5e-2] = "ns"
  PsgnArr[5e-2 >= PltPwc$p & PltPwc$p > 1e-2] = "*"
  PsgnArr[1e-2 >= PltPwc$p & PltPwc$p > 1e-3] = "**"
  PsgnArr[1e-3 >= PltPwc$p & PltPwc$p > 1e-4] = "***"
  PsgnArr[1e-4 >= PltPwc$p] = "****"
  
  PadjArr = c("ns", "ns", "ns")
  PadjArr[PltPwc$p.adj > 5e-2] = "ns"
  PadjArr[5e-2 >= PltPwc$p.adj & PltPwc$p.adj > 1e-2] = "*"
  PadjArr[1e-2 >= PltPwc$p.adj & PltPwc$p.adj > 1e-3] = "**"
  PadjArr[1e-3 >= PltPwc$p.adj & PltPwc$p.adj > 1e-4] = "***"
  PadjArr[1e-4 >= PltPwc$p.adj] = "****"
  
  PltPwc$p.signif = PsgnArr
  PltPwc$p.adj.signif = PadjArr
    

  if(class(TstFrm$X) == "numeric"){
    ElmOrd = TstFrm$X %>% unique() %>% sort()
    xTckLab = ElmOrd
  } else{
    ElmOrd = GenSet
  }
  
  ElmLen = length(ElmOrd)
  
  # - PLOT BOX -----------------------------------------------------------------
  p <- ggboxplot(TstFrm, x = "X", y = "Y",
                 fill='#F8776D',
                 facet.by = "G",
                 combine = TRUE,
                 order = ElmOrd, #c("AgeGen", "DisGen", "ExcNgb")
                 outlier.size = 0.1
  ) 
  
  
  p = p + coord_cartesian(ylim = c(yLimMin[x],yLimMax[x]))

  p = p+stat_pvalue_manual(PltPwc, label = "p.adj.signif", tip.length = 0, step.increase = 0.00, label.size = TstTxtSiz) 

  p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL) #x, y, xy
  

  p = p +
    annotate("text",
             x = 1:ElmLen,
             y = YannMin[x],
             label = NumArr %>% format_with_commas(), 
             col = "black",
             size = AnnTxtSiz,
             vjust = - 1)
  
  p =  p + scale_x_discrete(labels=xTckLab)
  
  p = p + theme(axis.text.x = element_text(angle = TckAng, vjust = 1.2, hjust=1)) 
  
  p = p + scale_y_continuous(breaks = c(-0.2,0,0.2,0.4,0.6,0.8,1))   
  
  
  # MEAN SYMBOL
  #p = p + stat_summary(fun = mean, geom = "point", shape = "*", 
  #             size = 6, color = "black", position = position_dodge(0.8))
  
  p = ggpar(p, xlab ="Gene Type", ylab = Yscr, legend.title = "") +
    font("title", size = 10, color = "black", face = "bold") +
    font("xlab", size = 9, color = "black", face = "bold")+
    font("ylab", size = 9, color = "black", face = "bold")+
    font("xy.text", size = 8, color = "black") +
    theme(legend.position = "none") 
  
  if(x %in% c(1,2)){
    p = p + xlab(NULL)
  }
  
  
  if(x %in% c(2,4,5,6)){
    p = p + ylab(NULL)
  }
  
  P[[NtwQry]] = p
  
}


plot_grid(P$PPI, P$COX90, P$COX95, P$KEGG, labels = c("a","b","c","d"), label_size = 12, ncol=2)


################################################################################
# FUNCTIONS
################################################################################

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
