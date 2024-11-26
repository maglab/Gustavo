### DICTIONARY #################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Alp - Alpha (Transparency) (e.g., LinAlp)
# Arr - Array (e.g., PhnArr)
# Big - Big (e.g., BigSiz)
# Blu - Blue (e.g., BluCol)
# CI - Confidence Interval (e.g., WdtICv)
# Cnf - Confidence (e.g., LowCnf.Cox, UppCnf.Tau)
# Col - Column (e.g., ColFrm) or Color (e.g., BluCol)
# Cox - Coexpression (e.g., Mean.Cox, StsCox)
# Frm - Frame (e.g., CoxFrm, TauFrm)
# Gry - Grey (e.g., GryCol)
# h - Horizontal (e.g., WdtIQRh, WdtICh)
# Hum - Human (e.g., GenAge.Hum)
# IRQ - Interquartile Range (e.g., WdtIQRv)
# Lbl - Label (e.g., Lbl)
# Lin - Line (e.g., LinAlp)
# Len - Length (e.g., PhnLen)
# Low - Lower (e.g., LowCnf.Cox)
# Mod - Model Organisms (e.g., GenAge.Mod)
# Nod - Node (e.g., NodAlp)
# Phn - Phenotype (e.g., PhnArr, PhnQry)
# Q1 - First Quartile (e.g., Q1.Cox)
# Q3 - Third Quartile (e.g., Q3.Cox)
# Qry - Query (e.g., PhnQry)
# Scr - Score (e.g., Scr)
# Sml - Small (e.g., SmlSiz)
# Sts - Stats (e.g., StsCox, StsTau)
# Tau - Tau (e.g., TauFrm)
# Upp - Upper (e.g., UppCnf.Cox)
# v - Vertical (e.g., WdtIQRv, WdtICv)
# Wdt - Width (e.g., WdtIQRv)

### LIBRARIES ##################################################################

library(ggplot2)
library(dplyr)
library(ggrepel)

################################################################################
# LOADING SPECIFICITY AND COEXPRESSION DATA 
################################################################################

# LOADING COEXPRESSION DATA
CoxFrm = read.csv("maglab/Gustavo/Data/Generated/Coexpression/BoxPlot_Data.csv") %>%
  rename(Phn=Phenotype, Scr=Score, Lbl=Label) %>%
  filter(!(Phn %in% c("ARC-Pleiotropy_1","ARC-Pleiotropy_2","ARC-Pleiotropy_3",
                      "ARC-Pleiotropy_4","ARC-Pleiotropy_5","ARC-Pleiotropy_6")))

# LOADING SPECIFICITY DATA
TauFrm = read.csv("maglab/Gustavo/Data/Generated/Specificity/BoxPlot_Data.csv") %>%
  rename(Phn=Phenotype, Scr=Score, Lbl=Label) %>%
  filter(!(Phn %in% c("ARC-Pleiotropy_1","ARC-Pleiotropy_2","ARC-Pleiotropy_3",
                      "ARC-Pleiotropy_4","ARC-Pleiotropy_5","ARC-Pleiotropy_6")))

# COMPUTING STATS 
PhnArr = CoxFrm$Phn %>% unique()
PhnLen = length(PhnArr)
StsCox = data.frame()
StsTau = data.frame()
i=1
for(i in 1:PhnLen){
  print(i)
  PhnQry = PhnArr[i]
  CoxArr = CoxFrm %>% filter(Phn %in% PhnQry) %>% pull(Scr)
  TauArr = TauFrm %>% filter(Phn %in% PhnQry) %>% pull(Scr)
  sStsCox = calculateStats(CoxArr, PhnQry)
  sStsTau = calculateStats(TauArr, PhnQry)
  StsCox = rbind(StsCox,sStsCox)
  StsTau = rbind(StsTau,sStsTau)
}

colnames(StsCox) = c("Phn", "Mean.Cox", "Median.Cox", "Q1.Cox", "Q3.Cox", "LowCnf.Cox", "UppCnf.Cox")
colnames(StsTau) = c("Phn", "Mean.Tau", "Median.Tau", "Q1.Tau", "Q3.Tau", "LowCnf.Tau", "UppCnf.Tau")

# GENE-SET LABELLING
ColFrm = data.frame(Phn = c("immunological/systemic disorders", "cardiovascular", "endocrine/diabetes", "gastrointestinal/abdominal", "musculoskeletal/trauma",       
                            "haematology/dermatology", "neurology/eye/psychiatry", "renal/urology", "GenAge.Hum", "GenAge.Mod"),
                    Col = c("Immunological", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "GenAge", "GenAge"))

# STATS FRAME
StsCoxTau = merge(StsCox, StsTau, by="Phn") %>% merge(ColFrm, by="Phn")


################################################################################
# PLOT 
################################################################################

# SETTING PLOT PARAMETERS

# COLORS
RedCol = "#F8696B"
ModCol = "#FF7F50"
BluCol = "#5A8AC6"
GryCol = "#CAE2EE" 

# ALPHA, TRANSPARENCY
LinAlp = 0.7
NodAlp = 0.8

# SIZE OF BALS
BigSiz = 2.5
SmlSiz = 1


# CI AND IQR LINES WIDTHS
WdtCIv = 0.010
WdtCIh = 0.010

WdtIQRv = 0.040
WdtIQRh = 0.065

# Plot the data using ggplot2 and ggrepel
p <- ggplot(StsCoxTau, aes(x = Median.Tau, y = Median.Cox, label = Phn)) +
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "Other"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, size = SmlSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "Other"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Other"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, size = SmlSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Other"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = GryCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "GenAge"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = RedCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "GenAge"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = RedCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "GenAge"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = RedCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "GenAge"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = RedCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Phn == "GenAge.Mod"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = ModCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Phn == "GenAge.Mod"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = ModCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Phn == "GenAge.Mod"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = ModCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Phn == "GenAge.Mod"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = ModCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "Immunological"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = BluCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "Immunological"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = BluCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Immunological"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = BluCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Immunological"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = BluCol, alpha = LinAlp) +
  
  geom_point(aes(color = Col), size = 4, alpha = NodAlp) +
  
 
  scale_color_manual(values = c("GenAge" = "red", "Immunological" = "blue", "ModAge" = "#FF4500", "Other" = "#67ADCF"), 
                     labels = c("GenAge.Hum\n(High iARC-Interactors\nin PPI and KEGG)\n", "Immunological\nDisorders\n(High ARC-Pleiotropy)\n",
                                "GenAge.Mod\n", "Other ARCs")) + #3E6C48, #398A50"
  
  
  scale_linetype_manual(values = c("IQR" = "solid", "CI" = "solid"), 
                        guide = guide_legend(override.aes = list(size = c(1.5, 1)))) +
  
  theme_minimal() +
  labs(title = "Specificity vs Self-Coexpression in GenAge and ARCs",
       x = "Specificity", y = "Self-Coexpression") +
  theme(panel.grid.major = element_line(linetype = 'dashed', color = "grey"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  labs(color = "Group Category") + #+
  
  guides(linetype = guide_legend(title = "Error Bar Type", 
                                 override.aes = list(color = "black", size = c(1.5, 1))))


print(p)


################################################################################
# FUNCTIONS
################################################################################

calculateStats <- function(array, Phn) {
  if (!is.vector(array) || is.null(array)) {
    stop("Input must be a non-null vector")
  }
  
  # Basic statistics
  stats <- data.frame(
    Phn,
    Mean = mean(array, na.rm = TRUE),
    Median = median(array, na.rm = TRUE),
    Q1 = quantile(array, 0.25, na.rm = TRUE),
    Q3 = quantile(array, 0.75, na.rm = TRUE)
  )
  
  # Confidence values
  stats$Lower_Confidence_Value <- quantile(array, 0.05, na.rm = TRUE)
  stats$Upper_Confidence_Value <- quantile(array, 0.95, na.rm = TRUE)
  
  return(stats)
}

