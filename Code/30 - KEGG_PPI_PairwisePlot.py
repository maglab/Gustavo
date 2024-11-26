### DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# All - All (e.g., AllCorr)
# Age (e.g., AgeCorr)
# Bar - Bar (e.g., ColBar)
# Col - Column (e.g., ColBar)
# Corr - Correlation (e.g., AgeCorr)
# lim - Limit (e.g., lim+.0015)

### LIBRARIES #################################################################################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
from scipy.interpolate import interp1d
import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report, plot_confusion_matrix
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sklearn
import seaborn as sns
from scipy.stats import pearsonr
from numpy.polynomial.polynomial import polyfit

# LOAD DATA ##################################################################################

BigFrame = pd.read_csv('maglab/Gustavo/Data/Generated/Pairwise_Predictions/PairwisePred_All.csv')

Test = BigFrame['Test']
Class = BigFrame['Class']
Label = BigFrame['Gen']

lim = 1
Prob_Mean = BigFrame['Prob.Mean']
Prob_KEGG = BigFrame['Prob.KEGG']
Prob_PPI = BigFrame['Prob.PPI']

ColBar = pd.DataFrame({'Label': Label, 'Test': Test, 'Class': Class, 'Prob_Mean': Prob_Mean, 'Prob_KEGG': Prob_KEGG, 'Prob_PPI': Prob_PPI} )

### FIGURE PARAMETERS ####################################################################################

fig, ax = plt.subplots(1,1)
x, y = np.meshgrid(np.linspace(0, lim+.0015, 100), np.linspace(0, lim+.0015, 1000),)#
class_indices = ColBar['Test'] == 0
regplot = False
# For textbox

grey = 'grey'
ticksidesize = 11
tickmainsize = 11
labelsidesize = 11
labelmainsize = 11
titlesize = 11
titleweight = 'bold'
labeltickdist = 0.08

Pal = sns.color_palette()
Colors = Pal.as_hex()

snblue = Colors[0]
snred = Colors[3]

# GRAPH SELECTION 0:A, 1:G
j = 0

### FIGURE ###############################################################################################

L = {}
#0-10 # lines
if j == 0:
    for i in range(11):
        L[i] = MeanLines(ax,i,lim)
if j == 1:
    for i in range(11):
        L[i] = GmeanLines(ax,i,lim)
i = i+1
#11
if j == 0:
    L[i] = MeanPlot(ax, x, y)
if j == 1:
    L[i] = GmeanPlot(ax, x, y)
hi = i
i = i+1
#12
L[i] = ax.scatter(x = ColBar.loc[class_indices, 'Prob_PPI'], y =  ColBar.loc[class_indices, 'Prob_KEGG'],  c = snblue, edgecolor='k', alpha=1, linewidth=1, s = 40, zorder=2)  
nci = i
i = i+1
#13
L[i] = ax.scatter(x = ColBar.loc[~class_indices, 'Prob_PPI'], y =  ColBar.loc[~class_indices, 'Prob_KEGG'],  c = snred, edgecolor='k', alpha=1, linewidth=1, s = 40, zorder=2) 
ci=i
i = i+1

m = ColBar.loc[class_indices, 'Prob_PPI']
n =  ColBar.loc[class_indices, 'Prob_KEGG']
NotAgeCorr = round(pearsonr(m, n)[0] , 2)
m = ColBar.loc[~class_indices, 'Prob_PPI']
n =  ColBar.loc[~class_indices, 'Prob_KEGG']
AgeCorr = round(pearsonr(m, n)[0] , 2)
m = ColBar.loc[:, 'Prob_PPI']
n =  ColBar.loc[:, 'Prob_KEGG']
AllCorr = round(pearsonr(m, n)[0] , 2)

if regplot:
    a = 1
    ls = 'solid'
    lw = 1
    # 14
    m = ColBar.loc[class_indices, 'Prob_PPI']
    n =  ColBar.loc[class_indices, 'Prob_KEGG']
    NotAgeCorr = round(pearsonr(m, n)[0] , 2)
    b, m = polyfit(m, n, 1)
    L[i] = ax.plot(m, b + (m * m), '-', c=snblue, alpha=a, linestyle = ls, linewidth=lw, label=str(NotAgeCorr) , zorder=1)
    ilc = i
    i = i + 1
    # 15
    m = ColBar.loc[~class_indices, 'Prob_PPI']
    n =  ColBar.loc[~class_indices, 'Prob_KEGG']
    AgeCorr = round(pearsonr(m, n)[0] , 2)
    b, m =  polyfit(m, n, 1)
    L[i] = ax.plot(m, b + (m * m), '-', c=snred,  alpha=a, linestyle = ls, linewidth=lw, label=str(AgeCorr), zorder=1)
    iln = i
    i = i + 1
    # 16
    m = ColBar.loc[:, 'Prob_PPI']
    n =  ColBar.loc[:, 'Prob_KEGG']
    AllCorr = round(pearsonr(m, n)[0] , 2)
    b, m =  polyfit(m, n, 1)
    L[i] = ax.plot(m, b + (m * m), '-', c='tab:purple',  alpha=a, linestyle = ls, linewidth=lw, label=str(AllCorr), zorder=1)
    ila = i
    # Legend
    axpos = ax.get_position()
    xmiddleaxe = (axpos.x0 + axpos.x1) / 2
    legend2 = ax.legend(framealpha=0.8, fontsize = 9, bbox_to_anchor=[xmiddleaxe, -0.2], loc='center', ncol=3, )
    legend2.set_title('Legend')
    ax.get_legend().get_title().set_color("red")
    ax.get_legend().get_title().set_weight("bold")
    ax.add_artist(legend2)

# TITLE
if j == 0:
    ax.set_title('Arithmetic mean of Ageing-probabilities\n\n', size = titlesize, weight = titleweight)
else:
    ax.set_title('Geometric mean of Ageing-probabilities\n', size = titlesize, weight = titleweight)

# LATERAL
ax.text(lim2 + labeltickdist-0.01, 0.5, 'Ageing-probabilities mean', rotation=270, va='center', c=grey, size = labelsidesize) # Gray color

# Lateral numbers amean
Lm = 0.02 # Lateral margin
if j == 0:
    for i in range(6,11):
        lim2 = lim + Lm
        ii = i/10
        ylevel =  (2*ii) - lim2
        text = str( ii )
        ax.text(lim2, ylevel , text, fontsize = ticksidesize, verticalalignment = 'center', c=grey)#, ha='left', va='left')#,ha='center', va='center')

# Lateral numbers gmean#
if j == 1:
    for i in range(10):
        lim2 = lim + Lm
        ii = i #i - 0.1
        ylevel =  pow( (i+1)/10, 2) / lim2
        text = str( (i+1)/10 )
        ax.text(lim2, ylevel , text, fontsize = ticksidesize, verticalalignment = 'center', c=grey)#, ha='left', va='left')#,ha='center', va='center')
    
    

# LEFT
ax.set_ylabel("KEGG-based Ageing-probabilities", size = labelmainsize)
ax.tick_params(axis='y', labelsize = labelmainsize)

# DOWN'
if j == 0:
    ydown = - Lm - 0.02 
    for i in range(6):
        lim2 = lim + 0.01
        ii = i #i - 0.1
        xlevel = (2*ii/10) #- lim2
        text = str( (i)/10 )
        ax.text(xlevel, ydown , text, fontsize = ticksidesize, verticalalignment = 'center', c=grey)#, ha='left', va='left')#,ha='center', va='center')
    
    ax.text(0.5, ydown - labeltickdist + 0.01 , 'Ageing-probabilities mean', rotation=0, ha='center', c=grey, size = labelsidesize) # Gray color
 
if j == 1:
    ax.set_xlabel("PPI-based Ageing-probabilities", size = labelmainsize)
    ax.tick_params(axis='x', labelsize = labelmainsize)

# TOP

if j == 0:

    ax.text(0.5, lim2 + labeltickdist + 0.01, "PPI-based Ageing-probabilities", rotation=0, ha='center', size = labelmainsize) # Gray color

    ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=False, bottom=False, top=True, left=True, right=False, labelsize = tickmainsize)


ax.set_xlim([0, lim])
ax.set_ylim([0, lim])

ax.grid(linestyle='dotted', linewidth=0.5, color = 'k')
ax.set_axisbelow('line')

legend1 = ax.legend([ L[ci],L[nci] ],['Ageing', 'Not_Ageing'], loc='upper left', framealpha=0.8, fontsize = 9, handletextpad=0.1)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="10%", pad=0.6)
fig.colorbar(L[hi], cax = cax)

trans = ax.get_xaxis_transform() # x in data untis, y in axes fraction
plt.annotate('Pearson correlation:\nAgeing: ' + str(AgeCorr) + '\tNot_Ageing: ' + str(NotAgeCorr) + '\t$All$: ' + str(AllCorr), xy=(0.5, -0.17), size=11, ha='center', va='top', bbox=dict(boxstyle='round', fc='w', ec="k"), xycoords=trans)

plt.show()

#######################################################################################################################
### FUNCTIONS #########################################################################################################
#######################################################################################################################

def MeanPlot(ax, x, y):
    h = (x + y) / 2
    return ax.pcolormesh(x, y, h, cmap="viridis", zorder=0)


def GmeanPlot(ax, x, y):
    h = np.sqrt(x * y) 
    return ax.pcolormesh(x, y, h, cmap="viridis", zorder=0)   


def GmeanLines(ax, score, lim):
    x  = np.arange(0.001, lim+0.01, 0.001)
    amp  = pow( (score/10) ,2) / x
    return ax.plot(x, amp, c = 'w', alpha=1, linestyle = 'dotted', linewidth=1, zorder=2)# dashed, solid, dotted

def MeanLines(ax, score, lim):
    x  = np.arange(0.001, lim+0.01, 0.001)
    amp  = (2*(score/10)) - x
    return ax.plot(x, amp, c = 'w', alpha=1, linestyle = 'dotted', linewidth=1, zorder=2)# dashed, solid, dotted


