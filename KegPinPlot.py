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

#BigFrame = pd.read_csv('F:/The Project/AJGF.csv')
#BigFrame = pd.read_csv('F:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/KegPin/FrmAll.csv')
BigFrame = pd.read_csv('D:/PhD/EndNtw/KegPin/FrmAll.csv')

#Go=Pin
#dip=Keg

normalised = False

Test = BigFrame['Test']
BigFrame['Class'] =  BigFrame['Test'] == 1
BigFrame.loc[BigFrame['Test'] == 1, 'Class'] = 'CR'
BigFrame.loc[BigFrame['Test'] == 0, 'Class'] = 'NotCR'
Class = BigFrame['Class']
Label = BigFrame['Gen']
#Label = BigFrame['Label']



# Mapped, common
if normalised:
    lim = 1.02
    amProb = BigFrame['Aam_jsrProb']
    gmProb = BigFrame['Agm_jsrProb']
    dProb = BigFrame['AjsrProbDIP']
    gProb = BigFrame['AjsrProbGO']
# Original
else:
    lim = 1
    amProb = BigFrame['Prb.Men']
    gmProb = BigFrame['Prb.Gmn']
    dProb = BigFrame['Prb.Keg']
    gProb = BigFrame['Prb.Pin']
    #amProb = BigFrame['Aam_joProb']
    #gmProb = BigFrame['Agm_joProb']
    #dProb = BigFrame['AjoProbDIP']
    #gProb = BigFrame['AjoProbGO']

ColBar = pd.DataFrame({'Label': Label, 'Test': Test, 'Class': Class, 'nsMean': amProb, 'nsGeom': gmProb, 'dns': dProb, 'gns': gProb} )


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#top=0.844,
#bottom=0.18,
#left=0.11,
#right=0.93,
#hspace=0.2,
#wspace=0.

fig, ax = plt.subplots(1,1)
x, y = np.meshgrid(np.linspace(0, lim+.0015, 100), np.linspace(0, lim+.0015, 1000),)#
class_indices = ColBar['Test'] == 0
regplot = False
# For textbox

# For secondary labels
#grey = '#474B4E'
#grey = '#FF0000'
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


L = {}
#0-10 # lines
if j == 0:
    for i in range(11):
        L[i] = aLines(ax,i,lim)
if j == 1:
    for i in range(11):
        L[i] = gLines(ax,i,lim)
i = i+1
#11
if j == 0:
    L[i] = AmPlot(ax, x, y)
if j == 1:
    L[i] = GmPlot(ax, x, y)
hi = i
i = i+1
#12
L[i] = ax.scatter(x = ColBar.loc[class_indices, 'gns'], y =  ColBar.loc[class_indices, 'dns'],  c = snblue, edgecolor='k', alpha=1, linewidth=1, s = 40, zorder=2)  
nci = i
i = i+1
#13
L[i] = ax.scatter(x = ColBar.loc[~class_indices, 'gns'], y =  ColBar.loc[~class_indices, 'dns'],  c = snred, edgecolor='k', alpha=1, linewidth=1, s = 40, zorder=2) 
ci=i
i = i+1

xx = ColBar.loc[class_indices, 'gns']
yy =  ColBar.loc[class_indices, 'dns']
Ncorr = round(pearsonr(xx, yy)[0] , 2)
xx = ColBar.loc[~class_indices, 'gns']
yy =  ColBar.loc[~class_indices, 'dns']
Ccorr = round(pearsonr(xx, yy)[0] , 2)
xx = ColBar.loc[:, 'gns']
yy =  ColBar.loc[:, 'dns']
Acorr = round(pearsonr(xx, yy)[0] , 2)

if regplot:
    a = 1
    ls = 'solid'
    lw = 1
    # 14
    xx = ColBar.loc[class_indices, 'gns']
    yy =  ColBar.loc[class_indices, 'dns']
    Ncorr = round(pearsonr(xx, yy)[0] , 2)
    b, m = polyfit(xx, yy, 1)
    L[i] = ax.plot(xx, b + (m * xx), '-', c=snblue, alpha=a, linestyle = ls, linewidth=lw, label=str(Ncorr) , zorder=1)
    ilc = i
    i = i + 1
    # 15
    xx = ColBar.loc[~class_indices, 'gns']
    yy =  ColBar.loc[~class_indices, 'dns']
    Ccorr = round(pearsonr(xx, yy)[0] , 2)
    b, m =  polyfit(xx, yy, 1)
    L[i] = ax.plot(xx, b + (m * xx), '-', c=snred,  alpha=a, linestyle = ls, linewidth=lw, label=str(Ccorr), zorder=1)
    iln = i
    i = i + 1
    # 16
    xx = ColBar.loc[:, 'gns']
    yy =  ColBar.loc[:, 'dns']
    Acorr = round(pearsonr(xx, yy)[0] , 2)
    b, m =  polyfit(xx, yy, 1)
    L[i] = ax.plot(xx, b + (m * xx), '-', c='tab:purple',  alpha=a, linestyle = ls, linewidth=lw, label=str(Acorr), zorder=1)
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

# COLBAR

#cbar_ax = fig.add_axes([0.90, 0.15, 0.04, 0.7]) #ok
#fig.colorbar(L[hi], cax = cbar_ax)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="10%", pad=0.6)
fig.colorbar(L[hi], cax = cax)

trans = ax.get_xaxis_transform() # x in data untis, y in axes fraction
plt.annotate('Pearson correlation:\nAgeing: ' + str(Ccorr) + '\tNot_Ageing: ' + str(Ncorr) + '\t$All$: ' + str(Acorr), xy=(0.5, -0.17), size=11, ha='center', va='top', bbox=dict(boxstyle='round', fc='w', ec="k"), xycoords=trans)

plt.show()

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

def go2mean(x):
    return x#2*x


def mean2go(x):
    return x#x/2



def AmPlot(ax, x, y):
    h = (x + y) / 2
    return ax.pcolormesh(x, y, h, cmap="viridis", zorder=0)


def GmPlot(ax, x, y):
    h = np.sqrt(x * y) 
    return ax.pcolormesh(x, y, h, cmap="viridis", zorder=0)   


def gLines(ax, score, lim):
    x  = np.arange(0.001, lim+0.01, 0.001)
    amp  = pow( (score/10) ,2) / x
    return ax.plot(x, amp, c = 'w', alpha=1, linestyle = 'dotted', linewidth=1, zorder=2)# dashed, solid, dotted

def aLines(ax, score, lim):
    x  = np.arange(0.001, lim+0.01, 0.001)
    amp  = (2*(score/10)) - x
    return ax.plot(x, amp, c = 'w', alpha=1, linestyle = 'dotted', linewidth=1, zorder=2)# dashed, solid, dotted


