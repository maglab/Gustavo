### DICTIONARY ###################################################################

# NOTE: Start by running the functions at the end of this file, if any. 

# Age - Age (e.g., AgeGen)
# Alg - Algorithm (e.g., EndAlgFrm)
# Arr - Array (e.g., GrpSrcArr)
# Auc - Area Under the Roc Curve (e.g., AucTrain)
# b - Binary (e.g., b_cvLen )
# BRF - Balanced Random Forest (e.g., )
# c - Continuous (e.g., c_cvLen )
# C - Correlation-based Filtering (e.g., VCPS)
# cv - Cross validation (e.g., c_cvLen)
# Dat - Dataset (e.g., DatArr)
# End - End (e.g., EndAlgFrm)
# Frm - Frame (e.g., EndAlgFrm)
# Grp - Graph (e.g., GrpSrcArr)
# Len - Length (e.g., c_cvLen)
# Mdl - Model (e.g., Mdl) 
# o - Original (e.g., oDataset)
# P - Pvalue-based Filtering (e.g., VCPS)
# S - Scaling (e.g., VCPS)
# Scr - Source (e.g., GrpSrcArr)
# Tst - Test (e.g., TstAuc )
# V - Variance-based Filterin (e.g., VCPS)

### LIBRARIES ###################################################################

import pickle
import scipy
import numpy as np 
scipy.interp = np.interp # temporal parch
import math
import pandas as pd
from numpy import interp
import seaborn as sns
import scikitplot as skplt
from itertools import chain
import statsmodels.api as sm
import inspect, math, itertools
from sklearn import preprocessing
from varname import nameof
from matplotlib import pyplot as plt
from imblearn.combine import SMOTEENN, SMOTETomek
from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import RandomUnderSampler
from sklearn.impute import KNNImputer, SimpleImputer

#from plot_metric.functions import BinaryClassification
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import OneHotEncoder, label_binarize, StandardScaler
from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold, cross_validate
from imblearn.over_sampling import BorderlineSMOTE, SMOTENC, SMOTE, SVMSMOTE, ADASYN, RandomOverSampler
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report
import copy
from functools import reduce
import os

import sys
# Add Python's project path
proyecto_path = r"D:/Respaldo_PHD/Code"
sys.path.append(proyecto_path)
# main.py
from ML_Pipeline import (
    run_ml_pipeline, 
    load_and_format_dataset, 
    get_dataset_array,
    impute_knn,
    BRFparamGridFunc
)



#############################################################################################################
# INDIVIDUAL (SINGLE ALGORITHM-DATASETS COMBINATIONS)
#############################################################################################################

GrpSrcArr = ['PPI','COX90','COX95','KEGG']

# Initialize EndAlgFrm as an empty DataFrame with the specified columns
EndAlgFrm = pd.DataFrame(columns=['Dataset', 'Auc', 'Network'])

DatArr =  get_dataset_array("All")

Dat = DatArr[1]
GrpSrc = GrpSrcArr[1]

for GrpSrc in GrpSrcArr:

    data = {'Dataset': ['X','Y'], 'Auc': [0,0]}
    AlgFrm = pd.DataFrame(data)

    for Dat in DatArr:
        print(Dat)
        #PathPrediction =  'D:/Respaldo_PHD/Nature_Data/Data/Generated/'+GrpSrc+'/Ageing_Prediction/Predictions/'
        oDataset = pd.read_csv('Data/Generated/Networks_and_predictions/Networks/'+GrpSrc+'/Ageing_Prediction/Datasets/'+Dat+'.csv')
        oDataset = oDataset.replace(np.nan, 0)
        oDataset.rename(columns={'AgeGen': 'Class'}, inplace=True)

       # Run ML
        results = run_ml_pipeline(
        dataset=oDataset,
        dataset_name=Dat,
        param_grid_func=BRFparamGridFunc,
        verbose=2,
        output_path = 'Data/Generated/Networks_and_predictions/Networks/'+GrpSrc+'/Ageing_Prediction/Predictions/'
        )
    i=1
i=1

#############################################################################################################
# PER-ALGORITHM (ONE ALGORITHM, ALL GRAPHS)
#############################################################################################################
 
GrpSrcArr = ['PPI','COX90','COX95','KEGG']

# Initialize EndAlgFrm as an empty DataFrame with the specified columns
EndAlgFrm = pd.DataFrame(columns=['Dataset', 'Auc', 'Network'])

DatArr =  get_dataset_array("All")

Dat = DatArr[1]
GrpSrc = GrpSrcArr[1]

NaTyp = "Imp" #["Int", "Imp"]#"Int" #Int, Imp
#DatArr = DatArr[1:]
for Dat in DatArr:

    data = {'Dataset': ['X','Y'], 'Auc': [0,0]}
    AlgFrm = pd.DataFrame(data)

    iDataset_list = []
    for GrpSrc in GrpSrcArr:
        iDataset = load_and_format_dataset(GrpSrc, Dat, base_path="Data/Generated/Networks_and_predictions/Networks")
        iDataset_list.append(iDataset)
        print(iDataset.head())
        i=1
    j=1
    # Progressive merge using "Unnamed: 0" and keeping one single "Class" column
    ooDataset = reduce(
        lambda left, right: pd.merge(left, right, on=["Unnamed: 0", "Class"], how="outer"),
        iDataset_list
    )

    # Handle na values
    oDataset = handle_na_values(ooDataset, NaTyp)

    # Save dataset
    oDataset.to_csv('Data/Generated/Networks_and_predictions/Networks/PER-ALGORITHM/Ageing_Prediction/Datasets/' + Dat + '.csv', index=False)

    # Use ML Function
    results = run_ml_pipeline(
    dataset=oDataset,
    dataset_name=Dat,
    param_grid_func=BRFparamGridFunc,
    verbose=2,
    output_path = 'Data/Generated/Networks_and_predictions/Networks/PER-ALGORITHM/Ageing_Prediction/Predictions'
    )
i=1

#############################################################################################################
# PER-GRAPH (ONE GRAPH, ALL ALGORITHMS)
#############################################################################################################

GrpSrcArr = ['PPI','COX90','COX95','KEGG']

# Initialize EndAlgFrm as an empty DataFrame with the specified columns
EndAlgFrm = pd.DataFrame(columns=['Dataset', 'Auc', 'Network'])

DtsTyp = "Arc" #Arc, Ard, All

DtsTypArr = ["All", "Ard", "Arc"]

#GrpSrc = 'COX95'

# For graph
for GrpSrc in GrpSrcArr:

    data = {'Dataset': ['X','Y'], 'Auc': [0,0]},

    # For array of features
    for DtsTyp in DtsTypArr:

        DatArr = get_dataset_array(DtsTyp)

        #Dat = DatArr[1]
        #GrpSrc = GrpSrcArr[2]

        # Create list of datasets
        AlgFrm = pd.DataFrame(data)
        iDataset_list = []
        for Dat in DatArr:
            print(Dat)
            iDataset = load_and_format_dataset(GrpSrc, Dat, base_path="Data/Generated/Networks_and_predictions/Networks")
            iDataset = iDataset.rename(
                columns=lambda c: c + " (" + DtsTyp + ")" + " (" + Dat + ")" if c not in ["Unnamed: 0", "Class"] else c
            )
            iDataset_list.append(iDataset)
            print(iDataset.head())
            j=1
        i=1

        # Progressive merge using "Unnamed: 0" and keeping one single "Class" column
        oDataset = reduce(
            lambda left, right: pd.merge(left, right, on=["Unnamed: 0", "Class"], how="outer"),
            iDataset_list
        )

        # Save dataset
        oDataset.to_csv('Data/Generated/Networks_and_predictions/Networks/PER-GRAPH/Ageing_Prediction/Datasets/' + GrpSrc + '.csv', index=False)

        # Use ML Function
        results = run_ml_pipeline(
        dataset=oDataset,
        dataset_name=GrpSrc + '_' + DtsTyp,
        param_grid_func=BRFparamGridFunc,
        verbose=2,
        output_path = 'Data/Generated/Networks_and_predictions/Networks/PER-GRAPH/Ageing_Prediction/Predictions'
        )
    i=1
i=1

#############################################################################################################
# MULTIPLEX
#############################################################################################################

GrpSrcArr = ['PPI','COX90','COX95','KEGG']

# Initialize EndAlgFrm as an empty DataFrame with the specified columns
EndAlgFrm = pd.DataFrame(columns=['Dataset', 'Auc', 'Network'])

DtsTyp = "All" #Arc, Ard, All
DtsTypArr = ["All", "Ard", "Arc"]
NaTyp = "Int" #Int, Imp
Mix = True

# For array of features
for DtsTyp in DtsTypArr:

    iDataset_list = []

    # MULTIPLEX -----------------------------------------------------
    if Mix==True:
          if DtsTyp == "Arc" or DtsTyp == "All":
              iDataset = pd.read_csv('Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Datasets/RWR.Proximity2ARC.csv')
              iDataset = iDataset.replace(np.nan, 0)
              iDataset.rename(columns={'AgeGen': 'Class'}, inplace=True)
              new_cols = {col: f"{col} (MUX-RWR.Proximity2ARC)" for col in iDataset.columns if col not in ["Unnamed: 0", "Class"]}
              iDataset.rename(columns=new_cols, inplace=True)
              iDataset_list.append(iDataset)
          if DtsTyp == "Ard" or DtsTyp == "All":
              iDataset = pd.read_csv('Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Datasets/RWR.Proximity2ARD.csv')
              iDataset = iDataset.replace(np.nan, 0)
              iDataset.rename(columns={'AgeGen': 'Class'}, inplace=True)
              new_cols = {col: f"{col} (MUX-RWR.Proximity2ARD)" for col in iDataset.columns if col not in ["Unnamed: 0", "Class"]}
              iDataset.rename(columns=new_cols, inplace=True)
              iDataset_list.append(iDataset)

    #Dat = DatArr[1]
    #GrpSrc = GrpSrcArr[1]

    data = {'Dataset': ['X','Y'],'Auc': [0,0]},

    # Create list of datasets
    AlgFrm = pd.DataFrame(data)

    # Progressive merge using "Unnamed: 0" and keeping one single "Class" column
    ooDataset = reduce(lambda left, right: pd.merge(left, right, on=["Unnamed: 0", "Class"], how="outer"), iDataset_list)
    
    # Handle na values
    oDataset = handle_na_values(ooDataset, NaTyp)

    if DtsTyp=='All':
            oDataset.to_csv('Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Datasets/RWR.Proximity2ARD_ARC.csv', index=False)

    # Use ML Function
    results = run_ml_pipeline(
    dataset=oDataset,
    dataset_name=DtsTyp,
    param_grid_func=BRFparamGridFunc,
    verbose=2,
    output_path = 'Data/Generated/Networks_and_predictions/Networks/MULTIPLEX/Ageing_Prediction/Predictions'
    )
i=1



######################################################################################################
# FUNCTIONS
######################################################################################################

def load_and_format_dataset(GrpSrc, Dat, base_path="Data/Generated"):
    """
    Loads a dataset, replaces NaN with 0, renames AgeGen to Class,
    and appends the suffix (GrpSrc) to all columns except
    'Unnamed: 0' and 'Class'.
    """
    
    file_path = os.path.join(
        base_path,
        GrpSrc,
        "Ageing_Prediction",
        "Datasets",
        f"{Dat}.csv"
    )
    
    iDataset = pd.read_csv(file_path)
    iDataset = iDataset.replace(np.nan, 0)
    
    iDataset.rename(columns={'AgeGen': 'Class'}, inplace=True)
    
    new_cols = {
        col: f"{col} ({GrpSrc})"
        for col in iDataset.columns
        if col not in ["Unnamed: 0", "Class"]
    }
    
    iDataset.rename(columns=new_cols, inplace=True)
    
    return iDataset




def get_dataset_array(DtsTyp):
    """
    Returns the list of datasets according to the dataset type:
    - 'All' : ARD + ARC
    - 'Ard' : ARD only
    - 'Arc' : ARC only
    """
    
    if DtsTyp == "All":
        return [
            'Closest.Proximity2ARD', 'Closest.Proximity2ARC',
            'Average.Proximity2ARD', 'Average.Proximity2ARC',
            'Neighbours2ARD', 'Neighbours2ARC',
            'RWR.Proximity2ARD', 'RWR.Proximity2ARC'
        ]
    
    if DtsTyp == "Ard":
        return [
            'Closest.Proximity2ARD',
            'Average.Proximity2ARD',
            'Neighbours2ARD',
            'RWR.Proximity2ARD'
        ]
    
    if DtsTyp == "Arc":
        return [
            'Closest.Proximity2ARC',
            'Average.Proximity2ARC',
            'Neighbours2ARC',
            'RWR.Proximity2ARC'
        ]
    
    raise ValueError(f"Unknown DtsTyp: {DtsTyp}")




def handle_na_values(ooDataset, NaTyp, k=5, id_cols=('Unnamed: 0', 'Class')):
    """
    NA handling strategies:
    
    - 'Int' : Dataset intersection
              → removes columns containing at least one NA
    - 'Imp' : KNN imputation
              → keeps all rows
    
    Parameters
    ----------
    ooDataset : pd.DataFrame
        Input dataset
    NaTyp : str
        'Int' or 'Imp'
    k : int
        Number of neighbors for KNN (only if NaTyp == 'Imp')
    id_cols : tuple
        Columns that should not be imputed
    """
    
    if NaTyp == "Int":
        # Keep only complete columns (no NA values)
        return ooDataset.loc[:, ooDataset.notna().all()]
    
    if NaTyp == "Imp":
        # Imputation while keeping all rows
        return impute_knn(ooDataset, k=k, id_cols=id_cols)
    
    raise ValueError(f"Unknown NaTyp: {NaTyp}")

