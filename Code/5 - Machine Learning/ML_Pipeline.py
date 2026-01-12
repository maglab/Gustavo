# ml_pipeline.py
"""
Complete Machine Learning Pipeline for biomedical data analysis.
Includes preprocessing, nested cross-validation, and evaluation.
"""

import pandas as pd
import numpy as np
import os
import math
import warnings
warnings.filterwarnings('ignore')

# Scikit-learn imports
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import make_scorer, confusion_matrix, roc_auc_score
from sklearn.impute import KNNImputer, SimpleImputer

# Statsmodels (for backward elimination)
import statsmodels.api as sm

# Imbalanced-learn imports
from imblearn.metrics import geometric_mean_score
from imblearn.ensemble import BalancedRandomForestClassifier
from imblearn.over_sampling import (
    SMOTE, BorderlineSMOTE, SVMSMOTE, ADASYN, RandomOverSampler
)
from imblearn.under_sampling import RandomUnderSampler
from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.metrics import sensitivity_score, specificity_score

# Others
from pandas.api.types import CategoricalDtype

# ============================================================================
# DATA LOADING AND PREPARATION FUNCTIONS (EARLY FUNCTIONS)
# ============================================================================

def load_and_format_dataset(GrpSrc, Dat, base_path="Data/Generated/Networks_and_predictions/Networks"):
    """
    Loads a dataset, replaces NaN with 0, renames AgeGen to Class
    and adds the suffix (GrpSrc) to columns except 'Unnamed: 0' and 'Class'.
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
    Returns the list of datasets according to type:
    - 'All' : ARD + ARC
    - 'Ard' : only ARD
    - 'Arc' : only ARC
    """
    if DtsTyp == "All":
        return [
            'Closest.Proximity2ARD', 'Closest.Proximity2ARC',
            'Average.Proximity2ARD', 'Average.Proximity2ARC',
            'Neighbours2ARD', 'Neighbours2ARC',
            'RWR.Proximity2ARD', 'RWR.Proximity2ARC'
        ]
    elif DtsTyp == "Ard":
        return [
            'Closest.Proximity2ARD',
            'Average.Proximity2ARD',
            'Neighbours2ARD',
            'RWR.Proximity2ARD'
        ]
    elif DtsTyp == "Arc":
        return [
            'Closest.Proximity2ARC',
            'Average.Proximity2ARC',
            'Neighbours2ARC',
            'RWR.Proximity2ARC'
        ]
    else:
        raise ValueError(f"Unknown DtsTyp: {DtsTyp}")

# ============================================================================
# PREPROCESSING FUNCTIONS
# ============================================================================

def backwardElimination(x, Y, sl, columns):
    """
    Backward elimination based on p-values.
    """
    numVars = len(x[0])
    for i in range(0, numVars):
        regressor_OLS = sm.OLS(Y, x).fit()
        maxVar = max(regressor_OLS.pvalues).astype(float)
        if maxVar > sl:
            for j in range(0, numVars - i):
                if regressor_OLS.pvalues[j].astype(float) == maxVar:
                    x = np.delete(x, j, 1)
                    columns = np.delete(columns, j)
    return columns

def PreProcessing(sXtrain, oXtest, Ytrain, MaxClassPercent=1, CorrCoef=1, SL=0.5, 
                  VCPS=[True, False, False, True], Minoc=1, Mds='M'):
    """
    Applies complete preprocessing pipeline.
    """
    Features = sXtrain.columns
    TrainIndices = sXtrain.index
    TestIndices = oXtest.index
    
    Vtext = ''; Ctext = ''; Ptext = ''; Stext = ''
    
    # 1. VARIANCE
    if VCPS[0]:
        Vtext = f'm{math.floor(Minoc)}'
        Variance = MaxClassPercent * (1 - MaxClassPercent)
        selector = VarianceThreshold(Variance)
        selector.fit(sXtrain)
        vSelectedIndices = selector.get_support(indices=True)
        vSelectedColumns = Features[vSelectedIndices]
    else:
        vSelectedColumns = Features
    
    vXtrain = sXtrain[vSelectedColumns]
    vLen = len(vSelectedColumns)
    
    # 2. CORRELATION
    if VCPS[1]:
        Ctext = f'c{round(CorrCoef*100)}'
        corr = vXtrain.corr()
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i, j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cvSelectedColumns = vXtrain.columns[columns]
    else:
        cvSelectedColumns = vSelectedColumns
    
    cvXtrain = vXtrain[cvSelectedColumns]
    cvLen = len(cvSelectedColumns)
    
    # 3. P-VALUES (BACKWARD ELIMINATION)
    if VCPS[2]:
        Ptext = f'p{round(SL)}'
        pcvSelectedColumns = backwardElimination(
            cvXtrain.values, Ytrain, SL, cvSelectedColumns
        )
    else:
        pcvSelectedColumns = cvSelectedColumns
    
    pcvXtrain = cvXtrain[pcvSelectedColumns]
    pcvLen = len(pcvSelectedColumns)
    
    # 4. XTEST
    pcvXtest = oXtest[pcvSelectedColumns]
    
    # 5. SCALING
    if VCPS[3]:
        Stext = 's'
        sc = StandardScaler()
        spcvXtrain = sc.fit_transform(pcvXtrain)
        spcvXtest = sc.transform(pcvXtest)
    else:
        spcvXtrain = pcvXtrain
        spcvXtest = pcvXtest
    
    # Final DataFrames
    Xtrain = pd.DataFrame(data=spcvXtrain, columns=pcvSelectedColumns, index=TrainIndices)
    Xtest = pd.DataFrame(data=spcvXtest, columns=pcvSelectedColumns, index=TestIndices)
    
    VCPtext = f'{Mds}{Vtext}{Ctext}{Ptext}{Stext}'
    
    return Xtrain, Xtest, vLen, cvLen, pcvLen, VCPtext

# ============================================================================
# SAMPLING FUNCTIONS
# ============================================================================

def Sampling(oXtrain, oYtrain, opc="normal"):
    """
    Applies different sampling techniques to balance classes.
    """
    if opc == "normal":
        return oXtrain, oYtrain
    elif opc == "under":
        rus = RandomUnderSampler(random_state=42)
        return rus.fit_resample(oXtrain, oYtrain)
    elif opc == "over":
        ros = RandomOverSampler(random_state=42)
        return ros.fit_resample(oXtrain, oYtrain)
    elif opc == "smote":
        sm = SMOTE(random_state=42)
        return sm.fit_resample(oXtrain, oYtrain)
    elif opc == "blsmote":
        blm = BorderlineSMOTE(random_state=42)
        return blm.fit_resample(oXtrain, oYtrain)
    elif opc == "svmsmote":
        sm = SVMSMOTE(random_state=42)
        return sm.fit_resample(oXtrain, oYtrain)
    elif opc == "smoteenn":
        sme = SMOTEENN(random_state=42)
        return sme.fit_resample(oXtrain, oYtrain)
    elif opc == "smotetomek":
        smt = SMOTETomek(random_state=42)
        return smt.fit_resample(oXtrain, oYtrain)
    elif opc == "adasyn":
        ada = ADASYN(random_state=42)
        return ada.fit_resample(oXtrain, oYtrain)
    else:
        return oXtrain, oYtrain

# ============================================================================
# AUXILIARY FUNCTIONS
# ============================================================================

def unique(list1):
    """Returns list with unique elements."""
    return list(set(list1))

def intersection(List1, List2):
    """Returns intersection of two lists."""
    return set(List1).intersection(set(List2))

def get_binary_columns(feature_array):
    """Identifies binary columns in a DataFrame."""
    binary_columns_indices = []
    for i in range(feature_array.shape[1]):
        column = feature_array.iloc[:, i]
        is_binary = np.all((column == 0) | (column == 1) | (np.isnan(column)))
        if is_binary:
            binary_columns_indices.append(i)
    
    binary_columns = feature_array.iloc[:, binary_columns_indices]
    mask = np.ones(feature_array.shape[1], dtype=bool)
    mask[binary_columns_indices] = 0
    non_binary_columns = feature_array.loc[:, mask]
    
    return binary_columns, non_binary_columns

def impute_knn(df: pd.DataFrame, k: int = 5, id_cols=('Unnamed: 0',), round_ints: bool = True) -> pd.DataFrame:
    """
    Imputes missing values using KNN for numerical and mode for categorical.
    """
    df = df.copy()
    orig_dtypes = df.dtypes.to_dict()
    
    # Present ID columns
    id_cols_present = [c for c in id_cols if c in df.columns]
    
    # Categorical/boolean
    cat_cols = df.select_dtypes(include=['object', 'category', 'bool']).columns.tolist()
    cat_cols = [c for c in cat_cols if c not in id_cols_present]
    
    # Numerical
    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    num_cols = [c for c in num_cols if c not in id_cols_present 
                and not pd.api.types.is_bool_dtype(df[c].dtype)]
    
    # Numerical imputation with KNN
    num_df = pd.DataFrame(index=df.index)
    if num_cols:
        knn = KNNImputer(n_neighbors=k, weights='distance')
        num_imp = knn.fit_transform(df[num_cols])
        num_df = pd.DataFrame(num_imp, columns=num_cols, index=df.index)
        
        for col in num_cols:
            if num_df[col].isna().any():
                col_mean = np.nanmean(num_df[col].values)
                if np.isnan(col_mean):
                    col_mean = 0.0
                num_df[col] = num_df[col].fillna(col_mean)
    
    # Categorical imputation with mode
    cat_df = pd.DataFrame(index=df.index)
    if cat_cols:
        simputer = SimpleImputer(strategy='most_frequent')
        cat_imp = simputer.fit_transform(df[cat_cols])
        cat_df = pd.DataFrame(cat_imp, columns=cat_cols, index=df.index)
    
    # DataFrame reconstruction
    pieces = [df[id_cols_present]]
    if num_cols: 
        pieces.append(num_df)
    if cat_cols: 
        pieces.append(cat_df)
    
    out = pd.concat(pieces, axis=1)
    out = out[[c for c in df.columns if c in out.columns]]
    
    # Restore original types
    for col, dt in orig_dtypes.items():
        if col not in out.columns:
            continue
        try:
            if pd.api.types.is_bool_dtype(dt):
                out[col] = out[col].astype(bool)
            elif pd.api.types.is_integer_dtype(dt) and col in num_cols and round_ints:
                out[col] = pd.to_numeric(np.rint(out[col]), errors='coerce').astype(dt)
            elif isinstance(dt, CategoricalDtype):
                out[col] = out[col].astype('category')
            elif dt == 'object' and col in cat_cols:
                out[col] = out[col].astype('object')
        except Exception:
            if pd.api.types.is_integer_dtype(dt) and col in num_cols and round_ints:
                out[col] = pd.to_numeric(np.rint(out[col]), errors='coerce').fillna(0).astype('int64')
            elif col in cat_cols:
                out[col] = out[col].astype('object')
    
    return out

# ============================================================================
# PARAMETER AND METRIC FUNCTIONS
# ============================================================================

def BRFparamGridFunc(Xtrain):
    """
    Returns parameter grid for Balanced Random Forest.
    """
    NumFeatures = Xtrain.shape[1]
    
    BRFparamGrid = {
        'bootstrap': [True],
        'replacement': [True, False],
        'max_features': ['sqrt', 'log2'],
        'n_estimators': [500],
        'class_weight': ['balanced', 'balanced_subsample', None],
        'sampling_strategy': [1],
    }
    return BRFparamGrid

def GeometricMean(Ytest, Ypred):
    """Calculates Geometric Mean Score."""
    return geometric_mean_score(Ytest, Ypred)

def calculate_metrics(y_true, y_pred, y_prob=None):
    """Calculates evaluation metrics."""
    metrics = {
        'confusion_matrix': confusion_matrix(y_true, y_pred),
        'sensitivity': sensitivity_score(y_true, y_pred),
        'specificity': specificity_score(y_true, y_pred),
        'gmean': geometric_mean_score(y_true, y_pred),
    }
    
    if y_prob is not None:
        metrics['auc'] = roc_auc_score(y_true, y_prob)
    
    return metrics

# ============================================================================
# MAIN PIPELINE FUNCTIONS
# ============================================================================

def apply_sampling(X, y, method='normal'):
    """Wrapper for sampling function."""
    return Sampling(X, y, method)

def apply_preprocessing(X_train, X_test, y_train, **kwargs):
    """
    Applies preprocessing using the PreProcessing function.
    """
    vcps = kwargs.get('vcps', [False, False, False, True])
    corr_coef = kwargs.get('corr_coef', 0.99)
    minoc_val = kwargs.get('minoc_val', 2)
    sl = kwargs.get('sl', 0.5)
    
    return PreProcessing(
        sXtrain=X_train,
        oXtest=X_test,
        Ytrain=y_train,
        CorrCoef=corr_coef,
        Minoc=minoc_val,
        SL=sl,
        VCPS=vcps
    )[:2]  # Returns only X_train, X_test

def process_fold(X_train, y_train, X_test, y_test, **kwargs):
    """
    Processes an individual fold of the outer cross-validation.
    """
    # Extract parameters
    ml_algorithm = kwargs.get('ml_algorithm')
    param_grid_func = kwargs.get('param_grid_func')
    scorer = kwargs.get('scorer')
    sampling_method = kwargs.get('sampling_method', 'normal')
    vcps = kwargs.get('vcps', [False, False, False, True])
    minoc_array = kwargs.get('minoc_array', [2])
    corr_array = kwargs.get('corr_array', [0.99])
    sl = kwargs.get('sl', 0.5)
    inner_skf = kwargs.get('inner_skf')
    n_jobs = kwargs.get('n_jobs')
    verbose = kwargs.get('verbose', 1)
    fold_num = kwargs.get('fold_num', 0)
    
    # 1. Sampling
    X_train_resampled, y_train_resampled = apply_sampling(
        X_train, y_train, method=sampling_method
    )
    
    # 2. Search for best preprocessing parameters
    best_model = None
    best_metrics = {'gmean_train': -1}  # Initialize with negative value
    best_preprocessing_params = None
    
    # Flag to verify if any model was found
    model_found = False
    
    for corr_coef in corr_array:
        for minoc_val in minoc_array:
            if verbose >= 2:
                print(f"  Testing: Corr={corr_coef}, Minoc={minoc_val}")
            
            try:
                # Apply preprocessing
                X_train_processed, X_test_processed = apply_preprocessing(
                    X_train_resampled, X_test, y_train_resampled,
                    vcps=vcps,
                    corr_coef=corr_coef,
                    minoc_val=minoc_val,
                    sl=sl
                )
                
                # Verify we have data after preprocessing
                if X_train_processed.shape[1] == 0:
                    if verbose >= 2:
                        print(f"  Warning: No features after preprocessing")
                    continue
                
                # Optimize hyperparameters
                param_grid = param_grid_func(X_train_processed)
                grid_search = GridSearchCV(
                    estimator=ml_algorithm,
                    param_grid=param_grid,
                    scoring=scorer,
                    cv=inner_skf,
                    n_jobs=n_jobs,
                    verbose=verbose-1 if verbose > 0 else 0
                )
                
                grid_search.fit(X_train_processed, y_train_resampled)
                
                # Evaluate on train
                y_pred_train = grid_search.best_estimator_.predict(X_train_processed)
                gmean_train = geometric_mean_score(y_train_resampled, y_pred_train)
                
                # Save if better
                if gmean_train > best_metrics['gmean_train']:
                    best_model = grid_search.best_estimator_
                    best_metrics = {
                        'gmean_train': gmean_train,
                        'grid_search': grid_search,
                        'param_grid': param_grid
                    }
                    best_preprocessing_params = {
                        'corr_coef': corr_coef,
                        'minoc_val': minoc_val,
                        'X_train_processed': X_train_processed,
                        'X_test_processed': X_test_processed
                    }
                    model_found = True
                    
            except Exception as e:
                if verbose >= 2:
                    print(f"  Error in preprocessing/GridSearch: {e}")
                continue
    
    # Verify that at least one model was found
    if not model_found:
        # Use default values
        if verbose >= 1:
            print(f"  Warning: No optimal model found, using default configuration")
        
        # Use first value of arrays
        corr_coef_default = corr_array[0] if corr_array else 0.99
        minoc_val_default = minoc_array[0] if minoc_array else 2
        
        X_train_processed_default, X_test_processed_default = apply_preprocessing(
            X_train_resampled, X_test, y_train_resampled,
            vcps=vcps,
            corr_coef=corr_coef_default,
            minoc_val=minoc_val_default,
            sl=sl
        )
        
        # Train model with default configuration
        param_grid_default = param_grid_func(X_train_processed_default)
        grid_search_default = GridSearchCV(
            estimator=ml_algorithm,
            param_grid=param_grid_default,
            scoring=scorer,
            cv=inner_skf,
            n_jobs=n_jobs,
            verbose=0
        )
        
        grid_search_default.fit(X_train_processed_default, y_train_resampled)
        
        best_model = grid_search_default.best_estimator_
        best_preprocessing_params = {
            'corr_coef': corr_coef_default,
            'minoc_val': minoc_val_default,
            'X_train_processed': X_train_processed_default,
            'X_test_processed': X_test_processed_default
        }
        
        y_pred_train_default = best_model.predict(X_train_processed_default)
        gmean_train_default = geometric_mean_score(y_train_resampled, y_pred_train_default)
        
        best_metrics = {
            'gmean_train': gmean_train_default,
            'grid_search': grid_search_default,
            'param_grid': param_grid_default
        }
    
    # 3. Final evaluation with best model
    X_test_final = best_preprocessing_params['X_test_processed']
    X_train_final = best_preprocessing_params['X_train_processed']
    
    # Make predictions
    y_pred_test = best_model.predict(X_test_final)
    
    # Get probabilities for BOTH classes
    try:
        y_proba_test = best_model.predict_proba(X_test_final)
        
        # Verify probability shape
        if y_proba_test.shape[1] == 2:
            # Class 0 = Age, Class 1 = NotAge
            y_prob_age = y_proba_test[:, 0]  # Probability of being Age
            y_prob_not_age = y_proba_test[:, 1]  # Probability of being NotAge
        else:
            # Model may have different structure
            y_prob_age = y_proba_test[:, 0] if y_proba_test.shape[1] > 0 else np.zeros(len(y_pred_test))
            y_prob_not_age = 1 - y_prob_age
    except Exception as e:
        if verbose >= 1:
            print(f"  Warning: Error getting probabilities: {e}")
        # If fails, use binary predictions
        y_prob_age = y_pred_test.astype(float)
        y_prob_not_age = 1 - y_prob_age
    
    # For metrics, use y_prob_not_age (as in standard scikit-learn)
    y_prob_test = y_prob_not_age
    
    # Calculate metrics
    y_pred_train_final = best_model.predict(X_train_final)
    
    metrics_train = calculate_metrics(y_train_resampled, y_pred_train_final)
    metrics_test = calculate_metrics(y_test, y_pred_test, y_prob_test)
    
    # Confusion matrices
    conf_matrix_train = confusion_matrix(y_train_resampled, y_pred_train_final)
    conf_matrix_test = confusion_matrix(y_test, y_pred_test)
    
    # Print fold metrics if verbose >= 1
    if verbose >= 1:
        print(f"\n{'='*60}")
        print(f"FOLD {fold_num + 1} - RESULTS")
        print(f"{'='*60}")
        
        print(f"Optimal parameters:")
        print(f"  Corr: {best_preprocessing_params['corr_coef']}, Minoc: {best_preprocessing_params['minoc_val']}")
        
        print(f"\nTrain Metrics:")
        print(f"  G-mean: {metrics_train['gmean']:.4f}")
        print(f"  Sensitivity: {metrics_train['sensitivity']:.4f}")
        print(f"  Specificity: {metrics_train['specificity']:.4f}")
        
        print(f"\nTest Metrics:")
        print(f"  G-mean: {metrics_test['gmean']:.4f}")
        print(f"  Sensitivity: {metrics_test['sensitivity']:.4f}")
        print(f"  Specificity: {metrics_test['specificity']:.4f}")
        if 'auc' in metrics_test:
            print(f"  AUC: {metrics_test['auc']:.4f}")
        
        print(f"\nConfusion Matrix (Test):")
        print(f"  TP: {conf_matrix_test[0,0]}, FP: {conf_matrix_test[0,1]}")
        print(f"  FN: {conf_matrix_test[1,0]}, TN: {conf_matrix_test[1,1]}")
        print(f"{'='*60}")
    
    # Return fold results
    return {
        'fold_num': fold_num,
        'model': best_model,
        'grid_search': best_metrics['grid_search'],
        'preprocessing_params': best_preprocessing_params,
        'y_pred_test': y_pred_test,
        'y_prob_test': y_prob_test,        # Prob NotAge (for metrics)
        'y_prob_age': y_prob_age,          # Prob Age
        'y_prob_not_age': y_prob_not_age,  # Prob NotAge
        'y_test': y_test,
        'metrics_train': metrics_train,
        'metrics_test': metrics_test,
        'conf_matrix_train': conf_matrix_train,
        'conf_matrix_test': conf_matrix_test,
        'best_params': {
            'corr_coef': best_preprocessing_params['corr_coef'],
            'minoc_val': best_preprocessing_params['minoc_val']
        }
    }

def calculate_global_metrics(all_results, y_test, y_pred, y_prob, labels):
    """
    Calculates global metrics by averaging all folds.
    """
    # Calculate global metrics
    global_metrics = calculate_metrics(y_test, y_pred, y_prob)
    
    # Calculate averages and standard deviations
    sensitivities = [r['metrics_test']['sensitivity'] for r in all_results.values()]
    specificities = [r['metrics_test']['specificity'] for r in all_results.values()]
    gmeans = [r['metrics_test']['gmean'] for r in all_results.values()]
    aucs = [r['metrics_test'].get('auc', 0) for r in all_results.values() if 'auc' in r['metrics_test']]
    
    avg_metrics = {
        'avg_sensitivity': np.mean(sensitivities),
        'std_sensitivity': np.std(sensitivities),
        'avg_specificity': np.mean(specificities),
        'std_specificity': np.std(specificities),
        'avg_gmean': np.mean(gmeans),
        'std_gmean': np.std(gmeans),
    }
    
    if aucs:
        avg_metrics['avg_auc'] = np.mean(aucs)
        avg_metrics['std_auc'] = np.std(aucs)
    
    return {**global_metrics, **avg_metrics}

def save_pipeline_results(results, dataset_name, classifier_name, output_path, scoring_method, global_auc):
    """
    Saves results to disk with format similar to original.
    """
    os.makedirs(output_path, exist_ok=True)
    
    global_auc_text = round(global_auc * 100)

    # Save global metrics
    #filename = f"{dataset_name}_{classifier_name}_{scoring_method}.csv"
    filename = f"{dataset_name}_{classifier_name}_AUC{global_auc_text}.csv"
    filepath = os.path.join(output_path,scoring_method, filename)

    # crear el directorio si no existe
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Create DataFrame with results
    results_df = pd.DataFrame([{
        'dataset': dataset_name,
        'classifier': classifier_name,
        'avg_sensitivity': results.get('avg_sensitivity', 'N/A'),
        'std_sensitivity': results.get('std_sensitivity', 'N/A'),
        'avg_specificity': results.get('avg_specificity', 'N/A'),
        'std_specificity': results.get('std_specificity', 'N/A'),
        'avg_gmean': results.get('avg_gmean', 'N/A'),
        'std_gmean': results.get('std_gmean', 'N/A'),
        'avg_auc': results.get('avg_auc', 'N/A'),
        'std_auc': results.get('std_auc', 'N/A'),
        'global_gmean': results.get('gmean', 'N/A'),
        'global_auc': results.get('auc', 'N/A')
    }])
    
    results_df.to_csv(filepath, index=False)
    
    print(f"\nGlobal results saved to: {filepath}")

# ============================================================================
# MAIN RUN_ML_PIPELINE FUNCTION
# ============================================================================

def run_ml_pipeline(
    dataset,  # DataFrame with data
    dataset_name,  # Dataset name
    target_column='Class',
    classifier_name='BRF',
    ml_algorithm=None,
    param_grid_func=None,
    importance_method='Gini',
    scoring_method='Gmean',
    sampling_method='normal',
    vcps=[False, False, False, True],
    minoc_flag=True,
    minoc_array=[2],
    corr_array=[0.99],
    sl=0.5,
    outer_splits=10,
    inner_splits=5,
    n_jobs=None,
    save_results=True,
    save_predictions=True,
    output_path='./results/',
    verbose=1
):
    """
    Executes a complete machine learning pipeline with nested cross-validation.
    """
    # Initialize algorithm if None
    if ml_algorithm is None:
        ml_algorithm = BalancedRandomForestClassifier(random_state=42)
    
    # Initial validations
    if param_grid_func is None:
        raise ValueError("param_grid_func function required")
    
    # Prepare output directory
    if save_results and not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # --------------------------------------------------------------------
    # 1. DATA PREPARATION
    # --------------------------------------------------------------------
    if verbose >= 1:
        print(f"\n{'='*80}")
        print(f"PROCESSING DATASET: {dataset_name}")
        print(f"{'='*80}")
        print(f"Classifier: {classifier_name}")
        print(f"Sampling method: {sampling_method}")
        print(f"Preprocessing (VCPS): {vcps}")
        print(f"Outer splits: {outer_splits}, Inner splits: {inner_splits}")
        print(f"{'='*80}")
    
    # Copy dataset and prepare
    data = dataset.copy()
    
    # Rename index column if necessary
    if data.columns[0] == "Unnamed: 0":
        data = data.rename(columns={"Unnamed: 0": "Genes"})
    elif data.columns[0] == "X":
        data = data.rename(columns={"X": "Genes"})
    
    # Set index
    if 'Genes' in data.columns:
        data.index = data['Genes']
        data = data.drop(columns=['Genes'])
    
    # Save original gene names
    gene_names = data.index.tolist()
    
    # Separate features and target
    features = [col for col in data.columns if col != target_column]
    X = data[features]
    y_label = data[target_column].values
    
    # Encode labels - IMPORTANT: Explicitly define classes
    le = LabelEncoder()
    y = le.fit_transform(y_label)
    
    # Class mapping for reference
    class_mapping = {0: "Age", 1: "NotAge"}
    class_names = le.classes_
    
    # Dataset information
    if verbose >= 1:
        print(f"\nDATASET INFORMATION:")
        print(f"Dataset size: {len(y)} samples (genes)")
        print(f"Number of features: {len(features)}")
        print(f"Class encoding: {dict(zip(le.transform(le.classes_), le.classes_))}")
        print(f"Class 'Age' (0): {sum(y == 0)} genes ({sum(y == 0)/len(y)*100:.1f}%)")
        print(f"Class 'NotAge' (1): {sum(y == 1)} genes ({sum(y == 1)/len(y)*100:.1f}%)")
        print(f"{'='*80}")
    
    # --------------------------------------------------------------------
    # 2. CROSS-VALIDATION CONFIGURATION
    # --------------------------------------------------------------------
    outer_skf = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=42)
    inner_skf = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=42)
    
    # Configure scoring
    if scoring_method == 'Gmean':
        scorer = make_scorer(geometric_mean_score, greater_is_better=True)
    else:
        scorer = scoring_method
    
    # --------------------------------------------------------------------
    # 3. EXTERNAL CROSS-VALIDATION (OUTER CV)
    # --------------------------------------------------------------------
    all_results = {}
    
    # Store all predictions
    all_y_test = []
    all_y_pred = []
    all_y_prob_age = []  # Probability of being Age (class 0)
    all_y_prob_not_age = []  # Probability of being NotAge (class 1)
    all_labels = []  # Gene names
    
    # Store class counts per fold
    train_counts_age = []
    train_counts_not_age = []
    test_counts_age = []
    test_counts_not_age = []
    
    for fold, (train_idx, test_idx) in enumerate(outer_skf.split(X, y)):
        if verbose >= 1:
            print(f"\n{'='*80}")
            print(f"EXTERNAL FOLD {fold + 1}/{outer_splits}")
            print(f"{'='*80}")
        
        # Split current fold
        X_train_fold = X.iloc[train_idx, :]
        y_train_fold = y[train_idx]
        X_test_fold = X.iloc[test_idx, :]
        y_test_fold = y[test_idx]
        labels_fold = X_test_fold.index.tolist()  # Gene names in test
        
        # Save class counts
        train_counts_age.append(sum(y_train_fold == 0))
        train_counts_not_age.append(sum(y_train_fold == 1))
        test_counts_age.append(sum(y_test_fold == 0))
        test_counts_not_age.append(sum(y_test_fold == 1))
        
        if verbose >= 1:
            print(f"\nClass distribution in Fold {fold + 1}:")
            print(f"  Train - Class 'Age' (0): {train_counts_age[-1]} genes")
            print(f"  Train - Class 'NotAge' (1): {train_counts_not_age[-1]} genes")
            print(f"  Test - Class 'Age' (0): {test_counts_age[-1]} genes")
            print(f"  Test - Class 'NotAge' (1): {test_counts_not_age[-1]} genes")
        
        # --------------------------------------------------------------------
        # 4. PREPROCESSING AND OPTIMIZATION PER FOLD
        # --------------------------------------------------------------------
        fold_results = process_fold(
            X_train_fold, y_train_fold, X_test_fold, y_test_fold,
            ml_algorithm=ml_algorithm,
            param_grid_func=param_grid_func,
            scorer=scorer,
            sampling_method=sampling_method,
            vcps=vcps,
            minoc_array=minoc_array,
            corr_array=corr_array,
            sl=sl,
            inner_skf=inner_skf,
            n_jobs=n_jobs,
            fold_num=fold,
            verbose=verbose
        )
        
        # Save fold results
        all_results[fold] = fold_results
        
        # Accumulate predictions
        all_y_test.extend(y_test_fold)
        all_y_pred.extend(fold_results['y_pred_test'])
        all_y_prob_age.extend(fold_results['y_prob_age'])
        all_y_prob_not_age.extend(fold_results['y_prob_not_age'])
        all_labels.extend(labels_fold)
        
        # Print some predictions from this fold
        if verbose >= 2:
            print(f"\nPREDICTIONS - Fold {fold + 1} (first 10 genes):")
            print(f"{'Gene':<20} {'Real':<10} {'Predicted':<10} {'Prob Age':<12} {'Prob NotAge':<12}")
            print(f"{'-'*70}")
            
            # Show first 10 genes of the fold
            for i in range(min(10, len(labels_fold))):
                gene = labels_fold[i]
                real = class_mapping[y_test_fold[i]]
                pred = class_mapping[fold_results['y_pred_test'][i]]
                prob_age = fold_results['y_prob_age'][i]
                prob_not_age = fold_results['y_prob_not_age'][i]
                
                # Highlight incorrect predictions
                marker = "✓" if real == pred else "✗"
                
                print(f"{gene:<20} {real:<10} {pred:<10} {prob_age:.4f}{marker:>2}    {prob_not_age:.4f}")
    
    # --------------------------------------------------------------------
    # 5. SAVE PREDICTIONS PER GENE (similar to original)
    # --------------------------------------------------------------------
    predictions_df = None
    predictions_filepath = None
    
    if save_predictions:
        # Create preprocessing text
        vcps_text = ""
        if vcps[0]: vcps_text += "V"
        if vcps[1]: vcps_text += "C"
        if vcps[2]: vcps_text += "P"
        if vcps[3]: vcps_text += "S"
        if not vcps_text:
            vcps_text = "normal"
        
        # Calculate global AUC to include in filename
        try:
            # Use NotAge probability for AUC calculation (as in scikit-learn)
            global_auc = roc_auc_score(all_y_test, all_y_prob_not_age)
            auc_score_text = f"Auc{int(100 * round(global_auc, 2))}"
        except Exception as e:
            print(f"Warning: Could not calculate AUC: {e}")
            auc_score_text = "AucXX"
        
        # Create DataFrame with predictions (as in original)
        predictions_data = []
        
        for i in range(len(all_labels)):
            gene = all_labels[i]
            real_class_num = all_y_test[i]
            pred_class_num = all_y_pred[i]
            prob_age_val = all_y_prob_age[i]
            prob_not_age_val = all_y_prob_not_age[i]
            
            # Apply transformation as in original
            # In original: Age=0, NotAge=1, but output shows Age=1, NotAge=0
            # That's why it does: Frm['Test'] = 1 - Frm['Test'], Frm['Pred'] = 1 - Frm['Pred']
            
            predictions_data.append({
                'Test': 1 - real_class_num,  # Transform for output
                'Pred': 1 - pred_class_num,  # Transform for output
                'Prob': prob_age_val,  # Probability of Age (already correct)
                'Label': gene,
                'Class': class_mapping[real_class_num],
                'Prob_Age': prob_age_val,
                'Prob_NotAge': prob_not_age_val,
                'Correct': 1 if real_class_num == pred_class_num else 0
            })
        
        predictions_df = pd.DataFrame(predictions_data)
        
        # Sort by Age probability descending (as in original)
        predictions_df = predictions_df.sort_values(by='Prob', ascending=False)
        
        # Save CSV file
        predictions_path = os.path.join(output_path, "Predictions")
        os.makedirs(predictions_path, exist_ok=True)
        
        filename = f"{dataset_name}_{auc_score_text}.csv"
        predictions_filepath = os.path.join(predictions_path, filename)
        
        # Save only necessary columns (as in original)
        output_columns = ['Test', 'Pred', 'Prob', 'Label', 'Class']
        predictions_df[output_columns].to_csv(predictions_filepath, index=False)
        
        if verbose >= 1:
            print(f"\n{'='*80}")
            print(f"PREDICTIONS SAVED - {dataset_name}")
            print(f"{'='*80}")
            print(f"File: {predictions_filepath}")
            print(f"Total genes: {len(predictions_df)}")
            print(f"Column format: {output_columns}")
            
            # Show statistics
            age_pred_count = sum(predictions_df['Pred'] == 1)  # Transformed prediction
            not_age_pred_count = sum(predictions_df['Pred'] == 0)
            
            print(f"\nPrediction distribution:")
            print(f"  Genes predicted as 'Age': {age_pred_count} ({age_pred_count/len(predictions_df)*100:.1f}%)")
            print(f"  Genes predicted as 'NotAge': {not_age_pred_count} ({not_age_pred_count/len(predictions_df)*100:.1f}%)")
            
            # Show top genes
            print(f"\nTop 5 genes with highest probability of being 'Age':")
            print(f"{'Gene':<20} {'Prob Age':<10} {'Real':<10} {'Pred':<10}")
            print(f"{'-'*60}")
            
            top_5 = predictions_df.head(5)
            for _, row in top_5.iterrows():
                gene = row['Label']
                prob = row['Prob']
                real = "Age" if row['Test'] == 1 else "NotAge"
                pred = "Age" if row['Pred'] == 1 else "NotAge"
                
                print(f"{gene:<20} {prob:.4f}     {real:<10} {pred:<10}")
    
    # --------------------------------------------------------------------
    # 6. GLOBAL RESULTS
    # --------------------------------------------------------------------
    if verbose >= 1:
        print(f"\n{'='*80}")
        print(f"GLOBAL RESULTS - {dataset_name}")
        print(f"{'='*80}")
        
        # Print class distribution per fold
        print(f"\nClass distribution across all folds:")
        print(f"Train - Class 'Age' per fold: {train_counts_age}")
        print(f"Train - Class 'NotAge' per fold: {train_counts_not_age}")
        print(f"Test - Class 'Age' per fold: {test_counts_age}")
        print(f"Test - Class 'NotAge' per fold: {test_counts_not_age}")
        
        print(f"\nMETRICS PER FOLD (Test):")
        print(f"{'Fold':<6} {'Sens':<8} {'Spec':<8} {'G-mean':<8} {'AUC':<8} {'Confusion Matrix'}")
        print(f"{'-'*60}")
        
        for fold in range(outer_splits):
            metrics = all_results[fold]['metrics_test']
            conf_matrix = all_results[fold]['conf_matrix_test']
            
            sens = metrics['sensitivity']
            spec = metrics['specificity']
            gmean = metrics['gmean']
            auc_val = metrics.get('auc', 0)
            
            # Compact confusion matrix format
            conf_str = f"[{conf_matrix[0,0]},{conf_matrix[0,1]};{conf_matrix[1,0]},{conf_matrix[1,1]}]"
            
            print(f"{fold+1:<6} {sens:.4f}  {spec:.4f}  {gmean:.4f}  {auc_val:.4f}  {conf_str}")
    
    # Calculate global metrics
    global_results = calculate_global_metrics(
        all_results, all_y_test, all_y_pred, all_y_prob_not_age, all_labels
    )
    
    # Print average metrics
    if verbose >= 1:
        print(f"\n{'='*80}")
        print(f"AVERAGE METRICS:")
        print(f"{'='*80}")
        print(f"Average sensitivity: {global_results['avg_sensitivity']:.4f} ± {global_results['std_sensitivity']:.4f}")
        print(f"Average specificity: {global_results['avg_specificity']:.4f} ± {global_results['std_specificity']:.4f}")
        print(f"Average G-mean: {global_results['avg_gmean']:.4f} ± {global_results['std_gmean']:.4f}")
        if 'avg_auc' in global_results:
            print(f"Average AUC: {global_results['avg_auc']:.4f} ± {global_results['std_auc']:.4f}")
        
        # Show global AUC
        try:
            global_auc = roc_auc_score(all_y_test, all_y_prob_not_age)
            print(f"Global AUC: {global_auc:.4f}")
        except:
            pass
        
        # Show best preprocessing parameters found
        print(f"\nBEST PREPROCESSING PARAMETERS:")
        print(f"{'Fold':<6} {'Corr':<8} {'Minoc':<8}")
        print(f"{'-'*25}")
        for fold in range(outer_splits):
            params = all_results[fold]['best_params']
            print(f"{fold+1:<6} {params['corr_coef']:<8} {params['minoc_val']:<8}")
    
    # --------------------------------------------------------------------
    # 7. SAVE RESULTS
    # --------------------------------------------------------------------
    if save_results:
        save_pipeline_results(
            global_results,
            dataset_name,
            classifier_name,
            output_path,
            scoring_method,
            global_auc
        )
        
        # Save detailed results per fold
        detailed_results = []
        for fold in range(outer_splits):
            metrics = all_results[fold]['metrics_test']
            params = all_results[fold]['best_params']
            
            fold_data = {
                'fold': fold + 1,
                'dataset': dataset_name,
                'classifier': classifier_name,
                'sensitivity': metrics['sensitivity'],
                'specificity': metrics['specificity'],
                'gmean': metrics['gmean'],
                'auc': metrics.get('auc', None),
                'corr_coef': params['corr_coef'],
                'minoc_val': params['minoc_val'],
                'tp': all_results[fold]['conf_matrix_test'][0, 0],
                'fp': all_results[fold]['conf_matrix_test'][0, 1],
                'fn': all_results[fold]['conf_matrix_test'][1, 0],
                'tn': all_results[fold]['conf_matrix_test'][1, 1]
            }
            detailed_results.append(fold_data)
        
        detailed_df = pd.DataFrame(detailed_results)
        results_path = os.path.join(output_path, "Results")
        os.makedirs(results_path, exist_ok=True)
        detailed_path = os.path.join(results_path, f"{dataset_name}_{classifier_name}_{scoring_method}_folds.csv")
        detailed_df.to_csv(detailed_path, index=False)
        
        if verbose >= 1:
            print(f"\nDetailed results per fold saved to: {detailed_path}")
    
    # --------------------------------------------------------------------
    # 8. GENERATE COMPLETE PREDICTIONS REPORT
    # --------------------------------------------------------------------
    if save_predictions and predictions_df is not None:
        # Create complete report
        report_path = os.path.join(output_path, "Reports")
        os.makedirs(report_path, exist_ok=True)
        
        # Report with all metrics
        full_report_df = predictions_df.copy()
        full_report_path = os.path.join(report_path, f"{dataset_name}_full_report.csv")
        full_report_df.to_csv(full_report_path, index=False)
        
        # Summary statistics
        accuracy = sum(predictions_df['Correct']) / len(predictions_df)
        age_correct = len(predictions_df[(predictions_df['Test'] == 1) & (predictions_df['Pred'] == 1)])
        not_age_correct = len(predictions_df[(predictions_df['Test'] == 0) & (predictions_df['Pred'] == 0)])
        
        summary_stats = {
            'Dataset': dataset_name,
            'Total_Genes': len(predictions_df),
            'Accuracy': accuracy,
            'Age_Correct': age_correct,
            'NotAge_Correct': not_age_correct,
            'Age_Predictions': sum(predictions_df['Pred'] == 1),
            'NotAge_Predictions': sum(predictions_df['Pred'] == 0),
            'AUC_Score': auc_score_text.replace('Auc', '')
        }
        
        summary_df = pd.DataFrame([summary_stats])
        summary_path = os.path.join(report_path, f"{dataset_name}_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        
        if verbose >= 1:
            print(f"\nComplete reports saved to: {report_path}")
            print(f"- Complete report: {full_report_path}")
            print(f"- Statistical summary: {summary_path}")
    
    # --------------------------------------------------------------------
    # 9. RETURN RESULTS
    # --------------------------------------------------------------------
    return {
        'dataset_name': dataset_name,
        'classifier_name': classifier_name,
        'global_results': global_results,
        'fold_results': all_results,
        'predictions': predictions_df,
        'predictions_filepath': predictions_filepath,
        'all_predictions_data': {
            'labels': all_labels,
            'y_test': all_y_test,
            'y_pred': all_y_pred,
            'y_prob_age': all_y_prob_age,
            'y_prob_not_age': all_y_prob_not_age
        },
        'class_mapping': class_mapping,
        'config': {
            'outer_splits': outer_splits,
            'inner_splits': inner_splits,
            'sampling_method': sampling_method,
            'scoring_method': scoring_method,
            'vcps': vcps,
            'minoc_array': minoc_array,
            'corr_array': corr_array
        }
    }

def print_fold_summary(fold_results, fold_num, outer_splits):
    """
    Prints a visual summary of fold results.
    """
    print(f"\n{'='*60}")
    print(f"FOLD {fold_num + 1}/{outer_splits} - SUMMARY")
    print(f"{'='*60}")
    
    metrics_test = fold_results['metrics_test']
    conf_matrix = fold_results['conf_matrix_test']
    
    # Create visual table
    print(f"\nConfusion matrix:")
    print(f"{'':<15} {'Predicted 0':<15} {'Predicted 1':<15}")
    print(f"{'-'*45}")
    print(f"{'Actual 0':<15} {conf_matrix[0,0]:<15} {conf_matrix[0,1]:<15}")
    print(f"{'Actual 1':<15} {conf_matrix[1,0]:<15} {conf_matrix[1,1]:<15}")
    
    print(f"\nMetrics:")
    print(f"{'-'*30}")
    for metric, value in metrics_test.items():
        if metric != 'confusion_matrix':
            print(f"{metric.capitalize():<15} {value:.4f}")
    
    print(f"{'='*60}")


def save_predictions_per_gene(all_labels, all_y_test, all_y_pred, all_y_prob_age, 
                              dataset_name, classifier_name, output_path, 
                              sampling_method, vcps_text=""):
    """
    Saves individual predictions for each gene.
    """
    # Create directory for predictions if it doesn't exist
    predictions_path = os.path.join(output_path, "predictions")
    os.makedirs(predictions_path, exist_ok=True)
    
    # Create DataFrame similar to original
    class_mapping = {0: "Age", 1: "NotAge"}
    
    # Create arrays for each column
    labels = []
    real_values = []
    pred_values = []
    prob_age_values = []  # Probability of being Age
    class_names = []
    
    n_samples = len(all_labels)
    
    for i in range(n_samples):
        labels.append(all_labels[i])
        real_values.append(all_y_test[i])
        pred_values.append(all_y_pred[i])
        prob_age_values.append(all_prob_age[i])  # Use prob_age directly
        class_names.append(class_mapping.get(all_y_test[i], "Unknown"))
    
    # Create DataFrame
    data = {
        'Test': real_values,
        'Pred': pred_values,
        'Prob': prob_age_values,  # Now it's Probability of Age
        'Label': labels,
        'Class': class_names
    }
    
    predictions_df = pd.DataFrame(data)
    
    # NOW WE DON'T NEED TO DO 1 - Prob because it's already Age probability
    # We only need to transform Test and Pred values if necessary
    
    # If original encodes Age=0 and NotAge=1, and we want Age=1, NotAge=0 in output:
    predictions_df['Test'] = 1 - predictions_df['Test']
    predictions_df['Pred'] = 1 - predictions_df['Pred']
    
    # Sort by probability descending (higher Age probability first)
    predictions_df = predictions_df.sort_values(by='Prob', ascending=False)
    
    # Generate filename similar to original
    # In original: Df_Nam + "_" + 'Auc' + str(int(100*round(AverageAucTest,2))) + '.csv'
    # We'll create a similar name
    
    # Calculate global AUC to include in name
    try:
        global_auc = roc_auc_score(all_y_test, all_y_prob)
        auc_score_text = f"Auc{int(100 * round(global_auc, 2))}"
    except:
        auc_score_text = "AucXX"
    
    # Filename similar to original
    filename = f"{dataset_name}_{auc_score_text}.csv"
    filepath = os.path.join(predictions_path, filename)
    
    # Save file
    predictions_df.to_csv(filepath, index=False)
    
    print(f"Predictions per gene saved to: {filepath}")
    print(f"Total genes predicted: {len(predictions_df)}")
    
    # Also save a predictions summary
    summary_filename = f"{dataset_name}_predictions_summary.txt"
    summary_filepath = os.path.join(predictions_path, summary_filename)
    
    with open(summary_filepath, 'w') as f:
        f.write(f"Predictions for dataset: {dataset_name}\n")
        f.write(f"Classifier: {classifier_name}\n")
        f.write(f"Sampling method: {sampling_method}\n")
        f.write(f"Preprocessing: {vcps_text}\n")
        f.write(f"Total genes: {len(predictions_df)}\n")
        f.write(f"\nTop 10 genes with highest probability of being Age:\n")
        f.write(f"{'-'*80}\n")
        
        top_10 = predictions_df.head(10)
        for idx, row in top_10.iterrows():
            f.write(f"{row['Label']}: Prob={row['Prob']:.4f}, Real={row['Test']}, Pred={row['Pred']}\n")
    
    return predictions_df, filepath


def generate_predictions_report(all_labels, all_y_test, all_y_pred, all_y_prob, 
                                dataset_name, output_path):
    """
    Generates a detailed report of predictions for each gene.
    """
    report_path = os.path.join(output_path, "predictions_reports")
    os.makedirs(report_path, exist_ok=True)
    
    # Create DataFrame with all predictions
    report_data = []
    
    for i in range(len(all_labels)):
        gene = all_labels[i]
        real_class = "Age" if all_y_test[i] == 0 else "NotAge"
        pred_class = "Age" if all_y_pred[i] == 0 else "NotAge"
        prob_age = 1 - all_y_prob[i]  # Probability of being Age
        
        report_data.append({
            'Gene': gene,
            'Real_Class': real_class,
            'Predicted_Class': pred_class,
            'Probability_Age': prob_age,
            'Probability_NotAge': all_y_prob[i],
            'Correct': 1 if all_y_test[i] == all_y_pred[i] else 0
        })
    
    report_df = pd.DataFrame(report_data)
    
    # Sort by probability of being Age (descending)
    report_df = report_df.sort_values(by='Probability_Age', ascending=False)
    
    # Save complete report
    full_report_path = os.path.join(report_path, f"{dataset_name}_full_predictions_report.csv")
    report_df.to_csv(full_report_path, index=False)
    
    # Save summary report
    summary_stats = {
        'Total_Genes': len(report_df),
        'Correct_Predictions': report_df['Correct'].sum(),
        'Accuracy': report_df['Correct'].sum() / len(report_df),
        'Age_Predictions': len(report_df[report_df['Predicted_Class'] == 'Age']),
        'NotAge_Predictions': len(report_df[report_df['Predicted_Class'] == 'NotAge'])
    }
    
    summary_df = pd.DataFrame([summary_stats])
    summary_path = os.path.join(report_path, f"{dataset_name}_predictions_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    
    # Save top genes by category
    top_age_genes = report_df.head(20)  # Top 20 with highest Age probability
    top_age_path = os.path.join(report_path, f"{dataset_name}_top_age_genes.csv")
    top_age_genes.to_csv(top_age_path, index=False)
    
    top_not_age_genes = report_df[report_df['Predicted_Class'] == 'NotAge'].head(20)
    top_not_age_path = os.path.join(report_path, f"{dataset_name}_top_not_age_genes.csv")
    top_not_age_genes.to_csv(top_not_age_path, index=False)
    
    print(f"\nPrediction reports saved to: {report_path}")
    print(f"- Complete report: {full_report_path}")
    print(f"- Summary report: {summary_path}")
    print(f"- Top Age genes: {top_age_path}")
    print(f"- Top NotAge genes: {top_not_age_path}")
    
    return report_df

# ============================================================================
# MAIN EXECUTION (if run as script)
# ============================================================================

if __name__ == "__main__":
    # Example usage
    print("ML Pipeline - Loading example data...")
    
    # Create example dataset
    n_samples = 200
    n_features = 100
    X_example = np.random.randn(n_samples, n_features)
    y_example = np.zeros(n_samples)
    y_example[:int(0.3 * n_samples)] = 1  # 30% class 1
    
    # Create DataFrame
    data = pd.DataFrame(
        X_example, 
        columns=[f'Gene_{i}' for i in range(n_features)]
    )
    data['Class'] = y_example
    
    # Run pipeline
    try:
        results = run_ml_pipeline(
            dataset=data,
            dataset_name="Example_Dataset",
            classifier_name="BRF",
            param_grid_func=BRFparamGridFunc,
            outer_splits=3,  # Reduced for quick test
            inner_splits=2,
            verbose=1,
            save_results=False
        )
        
        print(f"\nPipeline completed successfully!")
        print(f"Global G-Mean: {results['global_results']['gmean']:.4f}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()