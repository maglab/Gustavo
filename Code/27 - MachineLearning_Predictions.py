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
import math
import numpy as np 
import pandas as pd
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
from plot_metric.functions import BinaryClassification
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import OneHotEncoder, label_binarize, StandardScaler
from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold, cross_validate
from imblearn.over_sampling import BorderlineSMOTE, SMOTENC, SMOTE, SVMSMOTE, ADASYN, RandomOverSampler
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report
import copy

################################################################################################################################################
### START ######################################################################################################################################
################################################################################################################################################
 
GrpSrcArr = ['PPI','COX90','COX95','KEGG']


# Initialize EndAlgFrm as an empty DataFrame with the specified columns
EndAlgFrm = pd.DataFrame(columns=['Dataset', 'Auc', 'Network'])


DatArr =  ['Closest.Proximity2ARD', 'Closest.Proximity2ARC',
           'Average.Proximity2ARD','Average.Proximity2ARC',
           'Neighbours2ARD','Neighbours2ARC']

Dat = DatArr[1]
GrpSrc = GrpSrcArr[1]

for GrpSrc in GrpSrcArr:

    data = {
        'Dataset': ['X','Y'],
        'Auc': [0,0]},


    AlgFrm = pd.DataFrame(data)

    for Dat in DatArr:
        print(Dat)

        PathPrediction =  'maglab/Gustavo/Data/Generated/'+GrpSrc+'/Ageing_Prediction/Predictions/'

        oDataset = pd.read_csv('maglab/Gustavo/Data/Generated/'+GrpSrc+'/Ageing_Prediction/Datasets/'+Dat+'.csv')
        oDataset = oDataset.replace(np.nan, 0)
        oDataset.rename(columns={'AgeGen': 'Class'}, inplace=True)

        ### Continue as normal ############################################################################################################################

        Df_Nam = Dat
        ImportanceMethod = 'Gini' #Gini, Permutation
        nDataset = Dat
        nnDataset = nDataset
        ClassLabel = 'Class'
        Target = ClassLabel

        # CLASSIFIER 
        Classifier = 'BRF'
        MLalgorithm = BalancedRandomForestClassifier(random_state=42)
        ParamGridFunc = BRFparamGridFunc
        nJobs = None

        # SCORE
        MyScore = make_scorer(GeometricMean, greater_is_better=True)
        Score = 'Gmean'

        # SAMPLING
        SamplingMethod = 'normal'

        # PREPROCESSING
        VCPS = [False, False, False, True]   # For preprocessing, V: Variance-based Filtering, C: Correlation-based Filtering, P: Pvalue-based Filtering, S: SCALING
        bMask = [True, False, False, False]
        cMask = [False, False, False, True]
        MinocFlag = True
        if MinocFlag:
            Minoc = 4                      # Minimun allowed ocurrence in minority feature
            MinocArray = [2]
        else:
            NZV = 0.995

        CorrArray = [0.99]
        SL = 0.5

        # CROSS VALIDATIONS
        Method = 'cv'
        OuterSplits = 10
        CvOuterMethod = 'o'
        OuterSkf = StratifiedKFold(n_splits=OuterSplits, shuffle=True, random_state = 42)
        InnerSplits = 5
        CvInnerMethod = 'i'
        InnerSkf = StratifiedKFold(n_splits=InnerSplits, shuffle=True, random_state = 42)

        ImportanceSplits = 3
        ImportanceSkf = StratifiedKFold(n_splits=ImportanceSplits, shuffle=True, random_state = 42)

        # SEPARATOR
        Cntr = '_'


        # ------------------------------------------------------------------------------------------------------------------------------------------------------------

        # DEFINE AND PREPARE DATASET
        if oDataset.columns[0] == "Unnamed: 0":
            oDataset = oDataset.rename(columns={"Unnamed: 0": "Genes"})
        else:
            oDataset = oDataset.rename(columns={"X": "Genes"})

        oDataset.index = oDataset.iloc[:,0]            # Set the first column as row names 
        oDataset = oDataset.drop(columns=['Genes'])    # Drop first column
        Dataset = oDataset
 
        # DEFINE FEATURES AND TARGET
        Features = Dataset.columns.tolist()
        del Features[-1]                             # Remove class from features
        oLen = len(Features)

        # FEATURES AND TARGET ARRAYS
        X = Dataset[Features]

        # Y LABEL
        Ylabel = Dataset[Target].values

        # LABEL ENCODING
        le = preprocessing.LabelEncoder()
        le.fit(Ylabel)
        le.classes_
        y = le.transform(Ylabel)

        ####################################################################################################################################################################
        ##### OUTER CROSS VALIDATION #######################################################################################################################################
        ####################################################################################################################################################################

        ##### OUTER SPLITTING ##############################################################################################################################################

        # PARAMETERS SPLITTING
        kXtrain = []
        kYtrain = []
        kXtest = [] 
        kYtest = []
        kYlabel = []
        kYtrainCountAge = []
        kYtestCountAge = []
        kYtrainCountNotAge = []
        kYtestCountNotAge = []

        cont = 0
        for TrainIndex, TestIndex in OuterSkf.split(X, y):
            kXtrain.append(X.iloc[TrainIndex,:])
            kXtest.append(X.iloc[TestIndex]) 
            kYtrain.append(y[TrainIndex])
            kYtest.append(y[TestIndex])
            kYtrainCountAge.append(sum(kYtrain[cont] == 0))
            kYtestCountAge.append(sum(kYtest[cont] == 0))
            kYtrainCountNotAge.append(sum(kYtrain[cont] == 1))
            kYtestCountNotAge.append(sum(kYtest[cont] == 1))
            kYlabel.append(kXtrain[0].index)
            cont = cont+1

            print("DATASET SET:", len(y), "\n")
            print("MAJORITY CLASS IN DATASET:", sum(y == 1), "\n")
            print("MINORITY CLASS IN DATASET:", sum(y == 0), "\n\n")

            print("MAJORITY CLASS IN TRAINIG SET SPLITS:", kYtrainCountNotAge, "\n")
            print("MINORITY CLASS IN TRAINING SET SPLITS:", kYtrainCountAge, "\n")
            print("MAJORITY CLASS IN TESTING SET SPLITS:", kYtestCountNotAge, "\n")
            print("MINORITY CLASS IN TESTING SET SPLITS:", kYtestCountAge, "\n")

        # Minoc = MinocArray[0]

        AllYtest = []
        AllYpredTest = []
        AllYprobTest = []
        AllYlabelTest = []

        ModelsList = {}
        k = 0
        for k in range(0, OuterSplits):
            oXtrain = kXtrain[k]
            oYtrain = kYtrain[k]
            oXtest = kXtest[k]
            Ytest = kYtest[k]

            print("-------------------------------------------------- ", Dat, " OUTER ", k, " ----------------------------------------------------------\n")

            ##### PREPROCESSING ################################################################################################################################################

            # RESAMPLING  
            sXtrain, Ytrain = Sampling(oXtrain, oYtrain, opc = SamplingMethod) # opc = {normal, under, over, smote, blsmote, csmote, svmsmote, smoteenn, smotetomek, adasyn}

            iXtrain = sXtrain
            iXtest = oXtest

            # --- PREPROCESSING  -----------------------------------------------------------------------------------------------------------------------------------------
            # DEFINE MINIMUN NUMBER OF FEATURE OCURRENCES FOR NEAR ZERO VARIANCE
  
            GmeanTrainMax = 0
            for CorrCoefK in CorrArray:
                for MinocK in MinocArray:
                    MinocK = MinocK - 1
                    print("-------------------------------------------------- ", Dat, " OUTER ", k, ": ", MinocK, " ----------------------------------------------------------\n")

                    Nrows = iXtrain.shape[0]
                    Nrows
                    if MinocFlag:
                        MCP = 1 - (MinocK / Nrows)
                        MCP # Put in MaxClassPercent to get a determined repetition of features across instances
                        print('MinocFlag')
                    if MinocFlag == False:
                        MCP = NZV
                        MinocK = (1-MCP) * Nrows
                        print('NoMinocFlag')

                    print('NZV = ' +  str(MCP))

                    # SPLITTING BINARY AND CONTINUOUS
                    biXtrain, ciXtrain = get_binary_columns(iXtrain)
                    biXtest = iXtest[biXtrain.columns]
                    ciXtest = iXtest[ciXtrain.columns]

                    #biXtest, ciXtest = get_binary_columns(iXtest)

                    b_oLenK = len(biXtrain.columns)
                    c_oLenK = len(ciXtrain.columns)

                    oLenK = b_oLenK + c_oLenK

                    # FEATURE SELECTION (NEAR ZERO VARIANCE, CORRELATION, PVALUE-BASED ELIMINATION, AND SCALING-CENTERING)
                    XtrainK, XtestK, vLenK, cvLenK, pcvLenK, VCPtextK = PreProcessing(iXtrain, iXtest, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = VCPS, Minoc = MinocK+1, Mds = 'M')
                    b_vLenK = vLenK
                    b_cvLenK = cvLenK
                    b_pcvLenK = pcvLenK
                    c_vLenK = vLenK
                    c_cvLenK = cvLenK
                    c_pcvLenK = pcvLenK

                    print("Original Size     :", oLenK , "\n")
                    print("Near Zero Var Size:", vLenK , "\n")
                    print("Nzv + Corr Size   :", cvLenK , "\n")
                    # ------------------------------------------------------------------------------------------------------------------------------------------------------------

                    ##### INNER CV ###############################################################################################################################################

                    # PARAMETERS SPLITTING
                    kXlearn = []
                    kXlearn = []
                    kYlearn = []
                    kXval = [] 
                    kYval = []
                    kYlearnCountAge = []
                    kYvalCountAge = []
                    kYlearnCountNotAge = []
                    kYvalCountNotAge = []

                    cont = 0
                    for LearnIndex, ValIndex in InnerSkf.split(XtrainK, Ytrain):
                        kXlearn.append(XtrainK.iloc[LearnIndex,:])
                        kXval.append(XtrainK.iloc[ValIndex]) 
                        kYlearn.append(Ytrain[LearnIndex])
                        kYval.append(Ytrain[ValIndex])
                        kYlearnCountAge.append(sum(kYlearn[cont] == 0))
                        kYvalCountAge.append(sum(kYval[cont] == 0))
                        kYlearnCountNotAge.append(sum(kYlearn[cont] == 1))
                        kYvalCountNotAge.append(sum(kYval[cont] == 1))
                        cont = cont+1

                    print("---------------------------------------- INNER ------------------------------------------------\n")
                    print("TRAINING SET:", len(Ytrain), "\n")
                    print("TESTING SET:", len(Ytest), "\n\n")

                    print("MAJORITY CLASS IN TRAINING SET:", sum(Ytrain == 1), "\n")
                    print("MINORITY CLASS IN TRAINING SET:", sum(Ytrain == 0), "\n\n")

                    print("MAJORITY CLASS IN TESTING SET:", sum(Ytest == 1), "\n")
                    print("MINORITY CLASS IN TESTING SET:", sum(Ytest == 0), "\n\n")

                    print("MAJORITY CLASS IN LEARNING SET SPLITS:", kYlearnCountNotAge, "\n")
                    print("MINORITY CLASS IN LEARNING SET SPLITS:", kYlearnCountAge, "\n")
                    print("MAJORITY CLASS IN VALIDATION SET SPLITS:", kYvalCountNotAge, "\n")
                    print("MINORITY CLASS IN VALIDATION SET SPLITS:", kYvalCountAge, "\n")

                    #### CREATE CLASSIFIER ####################################################################

                    # INNER CROSS VALIDATION AND FIT
                    ParamGridK = ParamGridFunc(XtrainK)
                    ParamGridK = ParamGridFunc(XtrainK)

                    Ylabel = XtestK.index

                    GridSearchK = GridSearchCV(estimator = MLalgorithm, param_grid = ParamGridK, scoring = MyScore, cv = InnerSkf, verbose = 1, return_train_score = True)

                    GridSearchK.fit(XtrainK, Ytrain)
                    XtrainK.shape
                    len(Ytrain)


                    # BEST PARAMETERS
                    BestGridK = GridSearchK.best_estimator_
                    YpredTestK = BestGridK.predict(XtestK)
                    YprobTestK = BestGridK.predict_proba(XtestK)[:, 1]
                    YpredTrainK = BestGridK.predict(XtrainK)
                    YprobTrainK = BestGridK.predict_proba(XtrainK)[:, 1]

                    ##### PERFORMANCE ###################################################################################################

                    # TRAINING METRICS
                    ConfusionTrainK = confusion_matrix(YpredTrainK, Ytrain)     # It is reversed for easy of interpretation (it is just the transpose)
                    print(ConfusionTrainK)

                    SensitivityTrainK = sensitivity_score(Ytrain, YpredTrainK) 
                    SpecificityTrainK = specificity_score(Ytrain, YpredTrainK) 
                    GmeanTrainK = geometric_mean_score(Ytrain, YpredTrainK) 
                    AucTrainK = roc_auc_score(Ytrain, YprobTrainK)

                    if GmeanTrainK > GmeanTrainMax:
                        GmeanTrainMax = GmeanTrainK
                        BestGrid = BestGridK
                        GridSearch = GridSearchK
                        YpredTest = YpredTestK
                        YprobTest = YprobTestK
                        YpredTrain = YpredTrainK
                        YprobTrain = YprobTrainK
                        ConfusionTrain = ConfusionTrainK
                        SensitivityTrain = SensitivityTrainK
                        SpecificityTrain = SpecificityTrainK
                        GmeanTrain = GmeanTrainK
                        AucTrain = AucTrainK
                        ParamGrid = ParamGridK

                        b_oLen = b_oLenK
                        c_oLen = c_oLenK
                        oLen = oLenK
                        b_vLen = b_vLenK
                        b_cvLen = b_cvLenK
                        b_pcvLen = b_pcvLenK
                        c_vLen = c_vLenK
                        c_cvLen = c_cvLenK
                        c_pcvLen = c_pcvLenK
                        VCPtext = VCPtextK
                        Xtrain = XtrainK
                        Xtest = XtestK
                        vLen = vLenK
                        cvLen = cvLenK
                        pcvLen = pcvLenK
                        Minoc = MinocK
                        CorrCoef = CorrCoefK
        
                    print('SensitivityTrain: ', SensitivityTrainK)
                    print('SpecificityTrain: ', SpecificityTrainK)
                    print('GmeanTrain: ', GmeanTrainK)
                    print('AucTrain: ', AucTrainK)
                MinocArrayEnd= 1


            print("-------------------------------------------------- BEST MINOC: ",  Minoc, " ----------------------------------------------------------\n")
    
            print('SensitivityTrain: ', SensitivityTrain)
            print('SpecificityTrain: ', SpecificityTrain)
            print('GmeanTrain: ', GmeanTrain)
            print('AucTrain: ', AucTrain)

            # TESTING METRICS 
            ConfusionTest = confusion_matrix(YpredTest, Ytest)        # It is reversed for easy of interpretation (it is just the transpose)
            print(ConfusionTest)

            SensitivityTest = sensitivity_score(Ytest, YpredTest) 
            SpecificityTest = specificity_score(Ytest, YpredTest) 
            GmeanTest = geometric_mean_score(Ytest, YpredTest) 
            AucTest = roc_auc_score(Ytest, YprobTest)

            print('SensitivityTest: ', SensitivityTest )
            print('SpecificityTest: ', SpecificityTest)
            print('GmeanTest: ', GmeanTest)
            print('AucTest: ', AucTest)

            # LISTING RESULTS
            Metrics = {'ConfusionTrain': ConfusionTrain, 'SensitivityTrain': SensitivityTrain, 'SpecificityTrain': SpecificityTrain, 'GmeanTrain': GmeanTrain, 'AucTrain': AucTrain,\
                        'ConfusionTest': ConfusionTest, 'SensitivityTest': SensitivityTest, 'SpecificityTest': SpecificityTest, 'GmeanTest': GmeanTest, 'AucTest': AucTest}
            Models = {'BestModel': BestGrid, 'Grid': GridSearch, 'MLalgorithm': MLalgorithm, 'ParamGrid': ParamGrid, 'Score':MyScore}
            Targets = {'Ytest': Ytest, 'YtestPred': YpredTest, 'YtestProb': YprobTest, 'Ylabel': Ylabel}
            Preprocessing = {'VCPS': VCPS, 'Mcp': MCP, 'Minoc': Minoc, 'MinocArray': MinocArray, 'CorrArray': CorrArray, 'CorrCoef': CorrCoef, 'SL': SL, 'oLen': oLen, 'vLen': vLen, 'cvLen': cvLen, 'pcvLen': pcvLen,\
                                'b_oLen': b_oLen, 'b_vLen': b_vLen, 'b_cvLen': b_cvLen, 'b_pcvLen': b_pcvLen,\
                                'c_oLen': c_oLen, 'c_vLen': c_vLen, 'c_cvLen': c_cvLen, 'c_pcvLen': c_pcvLen}
    
            CrossValidation = {'InnerSkf': InnerSkf, 'OuterSkf': OuterSkf}
            ModelsList[k] = {'Metrics':Metrics, 'Models': Models, 'Targets':Targets, 'Preprocessing':Preprocessing, 'CrossValidation':CrossValidation}

            AllYtest = np.concatenate((AllYtest, Ytest), axis=None)
            AllYpredTest = np.concatenate((AllYpredTest, YpredTest), axis=None)
            AllYprobTest = np.concatenate((AllYprobTest, YprobTest), axis=None)
            AllYlabelTest = np.concatenate((AllYlabelTest, Ylabel), axis=None)

            if k == (OuterSplits-1):
                Test = {'Real':AllYtest, 'Pred':AllYpredTest, 'Prob':AllYprobTest, 'Label':AllYlabelTest}

                # AVERAGING RESULTS
                CumulatedConfusionTrain = ModelsList[0]['Metrics']['ConfusionTrain'] - ModelsList[0]['Metrics']['ConfusionTrain']  #Just to initialize a zeros confusion matrix
                CumulatedSensitivityTrain = 0
                CumulatedSpecificityTrain = 0
                CumulatedGmeanTrain = 0
                CumulatedAucTrain = 0
                CumulatedConfusionTest = ModelsList[0]['Metrics']['ConfusionTest'] - ModelsList[0]['Metrics']['ConfusionTest']  #Just to initialize a zeros confusion matrix
                CumulatedSensitivityTest = 0
                CumulatedSpecificityTest = 0
                CumulatedGmeanTest = 0
                CumulatedAucTest = 0


                for z in range(0, OuterSplits):
                    CumulatedConfusionTrain = CumulatedConfusionTrain + ModelsList[z]['Metrics']['ConfusionTrain']
                    CumulatedSensitivityTrain = CumulatedSensitivityTrain + ModelsList[z]['Metrics']['SensitivityTrain']
                    CumulatedSpecificityTrain = CumulatedSpecificityTrain + ModelsList[z]['Metrics']['SpecificityTrain']
                    CumulatedGmeanTrain = CumulatedGmeanTrain + ModelsList[z]['Metrics']['GmeanTrain']
                    CumulatedAucTrain = CumulatedAucTrain + ModelsList[z]['Metrics']['AucTrain']
                    CumulatedConfusionTest = CumulatedConfusionTest + ModelsList[k]['Metrics']['ConfusionTest']
                    CumulatedSensitivityTest = CumulatedSensitivityTest + ModelsList[z]['Metrics']['SensitivityTest']
                    CumulatedSpecificityTest = CumulatedSpecificityTest + ModelsList[z]['Metrics']['SpecificityTest']
                    CumulatedGmeanTest = CumulatedGmeanTest + ModelsList[z]['Metrics']['GmeanTest']
                    CumulatedAucTest = CumulatedAucTest + ModelsList[z]['Metrics']['AucTest']


                AverageConfusionTrain = CumulatedConfusionTrain / OuterSplits
                AverageSensitivityTrain = CumulatedSensitivityTrain / OuterSplits
                AverageSpecificityTrain = CumulatedSpecificityTrain / OuterSplits
                AverageGmeanTrain = CumulatedGmeanTrain / OuterSplits
                AverageAucTrain = CumulatedAucTrain / OuterSplits
                AverageConfusionTest = CumulatedConfusionTest / OuterSplits
                AverageSensitivityTest = CumulatedSensitivityTest / OuterSplits
                AverageSpecificityTest = CumulatedSpecificityTest / OuterSplits
                AverageGmeanTest = CumulatedGmeanTest / OuterSplits
                AverageAucTest = CumulatedAucTest / OuterSplits

                AveragePerformances = {'AverageConfusionTrain': AverageConfusionTrain, \
                                        'AverageSensitivityTrain': AverageSensitivityTrain, \
                                        'AverageSpecificityTrain': AverageSpecificityTrain, \
                                        'AverageGmeanTrain': AverageGmeanTrain, \
                                        'AverageAucTrain': AverageAucTrain, \
                                        'AverageConfusionTest': AverageConfusionTest, \
                                        'AverageSensitivityTest': AverageSensitivityTest, \
                                        'AverageSpecificityTest': AverageSpecificityTest, \
                                        'AverageGmeanTest': AverageGmeanTest, \
                                        'AverageAucTest': AverageAucTest \
                                        }


                Preprocessing = {'VCPS': VCPS, 'Minoc': Minoc, 'CorrCoef': CorrCoef, 'SL': SL, 'MinocArray': MinocArray}
                CrossValidation = {'OuterSkf':OuterSkf}

                DatasetModel = {'AveragePerformances': AveragePerformances, 'ModelsList': ModelsList, 'Preprocessing':Preprocessing, 'CrossValidation':CrossValidation, 'Score':MyScore, 'Test':Test}


                AvSensString = str(int(round(AverageSensitivityTest*100)))
                AvSpecString = str(int(round(AverageSpecificityTest*100)))

                # --- SAVING --------------------------------------------------------------------------------------------------------------------------------------------------------------------

                tOuterSplits = str(OuterSplits)
                tInnerSplits = str(InnerSplits)

                AlgFrm.NewRow = pd.Series([Dat,
                                           int(100*round(AverageAucTest,2))], index=AlgFrm.columns)


                AlgFrm = pd.concat([AlgFrm, pd.DataFrame([AlgFrm.NewRow])], ignore_index=True)

                TstAuc = str(int(100*round(AverageAucTest,2)))


                # TAKING A LOOK AT EACH OF THE TEST PERFORMANCES
                show = 'Confusion'
                for z in range(0,OuterSplits):
                    if(show == 'Confusion'):
                        print('Split ', z, ':\n', DatasetModel['ModelsList'][z]['Metrics']['ConfusionTest'], '\n\n')
                    elif(show == 'Sensitivity'):
                        print('Split ', z, ':\n', DatasetModel['ModelsList'][z]['Metrics']['SensitivityTest'], '\n\n')
                    elif(show == 'Specificity'):
                        print('Split ', z, ':\n', DatasetModel['ModelsList'][z]['Metrics']['SpecificityTest'], '\n\n')
                    elif(show == 'Gmean'):
                        print('Split ', z, ':\n', DatasetModel['ModelsList'][z]['Metrics']['GmeanTest'], '\n\n')
                    elif(show == 'Auc'):
                        print('Split ', z, ':\n', DatasetModel['ModelsList'][z]['Metrics']['AucTest'], '\n\n')

                print(DatasetModel['AveragePerformances']['AverageConfusionTest'], '\n\n')
                print(DatasetModel['AveragePerformances']['AverageSensitivityTest'], '\n\n')
                print(DatasetModel['AveragePerformances']['AverageSpecificityTest'], '\n\n')
                print(DatasetModel['AveragePerformances']['AverageGmeanTest'], '\n\n')
                print(DatasetModel['AveragePerformances']['AverageAucTest'], '\n\n')

                print('Split ', z, ':\n', DatasetModel['ModelsList'][0]['Metrics']['ConfusionTest'], '\n\n')
                print('Split ', z, ':\n', DatasetModel['ModelsList'][0]['Metrics']['SensitivityTest'], '\n\n')
                print('Split ', z, ':\n', DatasetModel['ModelsList'][0]['Metrics']['SpecificityTest'], '\n\n')


                # SAVING FRAME
                PredictionSavingText = (PathPrediction + 
                                   Df_Nam + 
                                   "_" + 
                                   'Auc' + str(int(100*round(AverageAucTest,2))) + 
                                   '.csv')

                Mdl =  DatasetModel['Test']
                Class = np.repeat("NotAge", len(Mdl['Real']), axis=0)
                Class[Mdl['Real'] == 0] = "Age"
                d = {'Test': Mdl['Real'], 'Pred': Mdl['Pred'], 'Prob':Mdl['Prob'], 'Label':Mdl['Label'], 'Class':Class}
                Frm = pd.DataFrame(data=d)

                Frm['Prob'] = 1-Frm['Prob']
                Frm['Test'] = 1-Frm['Test']
                Frm['Pred'] = 1-Frm['Pred']

                Frm = Frm.sort_values(by='Prob', ascending=False)

                Frm.to_csv(PredictionSavingText)
            a=1
        b=1
    AlgFrm = AlgFrm[~AlgFrm['Alg'].isin(["X", "Y"])]
    AlgFrm["Network"] = GrpSrc
    EndAlgFrm = pd.concat([EndAlgFrm,AlgFrm], ignore_index=True)
EndAlgFrm.to_csv('maglab/Gustavo/Data/Generated/Ageing_Prediction_Performances/aML_Based_AUC.csv')

##########################################################################################
# FUNCTIONS 
##########################################################################################


def BRFparamGridFunc(Xtrain):
    MaxFeaturesGridNum = 5
    PowDivisions = 1/MaxFeaturesGridNum
    MaxFeatures = []
    NumFeatures = Xtrain.shape[1] 
    for i in range(1, MaxFeaturesGridNum + 1):
        MaxFeatures.append(round(pow(NumFeatures,PowDivisions*i)))
    MaxFeatures
    BRFparamGrid = {
        'bootstrap': [True],
        'replacement': [True, False],
        'max_features': ['sqrt','log2'], 
        'n_estimators': [500],
        'class_weight': ['balanced', 'balanced_subsample', None],
        'sampling_strategy': [1],
        }
    return(BRFparamGrid)



def GridCreator(Min, Max, Num, Func):
    Num = Num - 1
    Grid = []
    if func == 'lin':
        for i in range(0,Num+1):
            Slope = (Max - Min)/Num
            iValue = Min + Slope*i
            Grid.append(round(iValue))
    elif func == 'exp':
        for i in range(0,Num+1):
            ExpTerm = pow(Max-Min,i)
            iValue = rpund(Min + ExpTerm - 1 + i)
            Grid.append(round(iValue))
    elif func == 'iexp':
        for i in range(0,Num+1):
            ExpTerm = pow(Max-Min,i)
            iValue = Max - Min - rpund(Min + ExpTerm - 1 + i)
            Grid.append(round(iValue))
    elif func == 'gauss':  # FAKE GAUSS PROJECTION MADE WITH POLYNOMIAL
            Xmidl = 0.5
            Xup = 0.8
            Yup = 0.6
            A = np.array([[0, 0, 0, 1], [1, 1, 1, 1], [pow(Xmidl,3), pow(Xmidl,2), Xmidl, 1], [pow(Xup,3), pow(Xup,2), Xup, 1]])
            B = np.array([Min, Max, (Min+Max)/2, (Yup*Min)+((1-Yup)*max)])
            X = np.linalg.inv(A).dot(B)
            iValue = X[0]*pow(i,3) + X[1]*pow(i,2) + X[2]*i + X[3]
            Grid.append(round(iValue))
    return Grid

##### SCORE ###################################################################################################################################################

def GeometricMean(Ytest,Ypred):
   GM = geometric_mean_score(Ytest, Ypred)
   return GM

##### PRE PROCESSING ################################################################################################################################################

def PreProcessing(sXtrain, oXtest, MaxClassPercent = 1, CorrCoef = 1, SL = 0.5, VCPS = [True, False, False, True], Minoc = 1,  Mds = 'M'):

    Mtext = Mds

    # --- INDICES AND COLUMNS ---------------------------------------------------------------------------------------------------------------------------------------
    Features = sXtrain.columns
    TrainIndices = sXtrain.index
    TestIndices = oXtest.index
    Vtext = ''
    Ctext = ''
    Ptext = ''
    Stext = ''

    # --- FEATURES SELECTION (NEAR ZERO VARIANCE) -------------------------------------------------------------------------------------------------------------------
    if (VCPS[0]):
        Vtext = 'm' + str(math.floor(Minoc))# + 'v' + str(round(MaxClassPercent*100))
        Variance = MaxClassPercent * (1 - MaxClassPercent)
        selector = VarianceThreshold(Variance)
        selector.fit_transform(sXtrain) # Originally oXtrain
        vSelectedIndices = selector.get_support(indices=True)
        vSelectedColumns = Features[vSelectedIndices]
    else:
        vSelectedColumns= Features
    vXtrain = sXtrain[vSelectedColumns]
    vLen = len(vSelectedColumns)

    # --- FEATURES SELECTION (CORRELATION) ---------------------------------------------------------------------------------------------------------------------------
    # CORRELATION SELECTION
    if(VCPS[1]):
        Ctext = 'c' + str(round(CorrCoef*100))
        corr = vXtrain.corr()
        #sns.heatmap(corr)
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i,j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cvSelectedColumns = vXtrain.columns[columns]
    else:
        cvSelectedColumns = vSelectedColumns
    cvXtrain = vXtrain[cvSelectedColumns]
    cvLen = len(cvSelectedColumns)
    
    # --- PVALUES SELECTION (BACKWARD ELIMINATION) -------------------------------------------------------------------------------------------------------------------
    if (VCPS[2]):   
        Ptext = 'p' + str(round(SL))
        pcvSelectedColumns = backwardElimination(cvXtrain.values, Ytrain, SL, cvSelectedColumns)
    else:
        pcvSelectedColumns = cvSelectedColumns
    pcvXtrain = cvXtrain[pcvSelectedColumns]
    pcvLen = len(pcvSelectedColumns)

    # --- XTEST ------------------------------------------------------------------------------------------------------------------------------------------------------
    pcvXtest = oXtest[pcvSelectedColumns]

    # --- SCALER -----------------------------------------------------------------------------------------------------------------------------------------------------

    if (VCPS[3]):
        Stext = 's' 
        sc = StandardScaler()
        spcvXtrain = sc.fit_transform(pcvXtrain)
        spcvXtest = sc.transform(pcvXtest)
    else:
        spcvXtrain = pcvXtrain
        spcvXtest = pcvXtest

    # --- DATA FRAMES ------------------------------------------------------------------------------------------------------------------------------------------------
    Xtrain = pd.DataFrame(data=spcvXtrain, columns = pcvSelectedColumns, index = TrainIndices)
    Xtest = pd.DataFrame(data=spcvXtest, columns = pcvSelectedColumns, index = TestIndices)
    VCPtext = Mtext + Vtext + Ctext + Ptext + Stext

    return Xtrain, Xtest, vLen, cvLen, pcvLen, VCPtext




# BACKWARD ELIMINATION FUCTION
def backwardElimination(x, Y, sl, columns):
    numVars = len(x[0])
    for i in range(0, numVars):
        regressor_OLS = sm.OLS(Y, x).fit()
        maxVar = max(regressor_OLS.pvalues).astype(float)
        if maxVar > sl:
            for j in range(0, numVars - i):
                if (regressor_OLS.pvalues[j].astype(float) == maxVar):
                    x = np.delete(x, j, 1)
                    columns = np.delete(columns, j)
                    
    regressor_OLS.summary()
    return columns


##### SAMPLING #####################################################################################################################################################

def Sampling(oXtrain, oYtrain, opc = "normal"):

    # NORMAL
    if (opc == "normal"):
        sXtrain = oXtrain
        Ytrain = oYtrain

    # RANDOM UNDER SAMPLING
    elif(opc == "under"):
        rus = RandomUnderSampler(random_state=42)
        sXtrain, Ytrain = rus.fit_resample(oXtrain, oYtrain)  

    # RANDOM OVER SAMPLER
    elif(opc == "over"):
        ros = RandomOverSampler(random_state=42)
        sXtrain, Ytrain = ros.fit_resample(oXtrain, oYtrain)               

    # SMOTE
    elif(opc == "smote"):
        sm = SMOTE(random_state=42)
        sXtrain, Ytrain = sm.fit_resample(oXtrain, oYtrain)

    # BORDER LINE SMOTE
    elif(opc == "blsmote"):
         blm = BorderlineSMOTE(random_state=42)                            
         sXtrain, Ytrain = blm.fit_resample(oXtrain, oYtrain)

    #CATEGORIAL SMOTE
    elif(opc == "csmote"):
        smote_nc = SMOTENC(categorical_features=[range(0,len(Features))], random_state=42)
        sXtrain, Ytrain = smote_nc.fit_resample(oXtrain, oYtrain)

    # SVM SMOTE
    elif(opc == "svmsmote"):
        sm = SVMSMOTE(random_state=42)
        sXtrain, Ytrain = sm.fit_resample(oXtrain, oYtrain)

    # SMOTE ENN
    elif(opc == "smoteenn"):
        sme = SMOTEENN(random_state=42)
        uXtrain, Ytrain = sme.fit_resample(oXtrain, oYtrain)

    # SMOTE TOMEK
    elif(opc == "smotetomek"):
        smt = SMOTETomek(random_state=42)
        sXtrain, Ytrain = smt.fit_resample(oXtrain, oYtrain)

    # ADASYN
    elif(opc == "adasyn"):
        ada = ADASYN(random_state=42)
        sXtrain, Ytrain = ada.fit_resample(oXtrain, oYtrain)

    else:
        sXtrain = oXtrain
        Ytrain = oYtrain

    return sXtrain, Ytrain

##### UNIQUE ######################################################################################################################################################


def unique(list1): 
    list_set = set(list1) 
    unique_list = (list(list_set)) 
    return(unique_list)

##### NANINFS #####################################################################################################################################################

def NanInfs(Dataset):
    Instances = Dataset.index
    Features = Dataset.columns
    Nrows = Dataset.shape[0]
    Ncols = Dataset.shape[1]
    NanLoc = []
    InfLoc = []
    NanLocNames = []
    InfLocNames = []
    for row in range(0,Nrows):
        for col in range(0,Ncols-1):
            if math.isnan(Dataset.iloc[row,col]):
                NanLoc.append([row,col])
                NanLocNames.append([Instances[row],Features[col]])
            if np.isinf(Dataset.iloc[row,col]):
                InfLoc.append([row,col])
                InfLocNames.append([Instances[row],Features[col]])
    Nnan = len(NanLoc)
    NanIndex = []
    NanColumn = []
    for i in range(0,Nnan):
        NanIndex.append(NanLoc[i][0])
        NanColumn.append(NanLoc[i][1])
    NanIndex = unique(NanIndex)
    NanColumn = unique(NanColumn)
    NanIndexName = Instances[NanIndex]
    NanColumnName = Features[NanColumn]

    Ninf = len(InfLoc)
    InfIndex = []
    InfColumn = []
    for i in range(0,Ninf):
        InfIndex.append(InfLoc[i][0])
        InfColumn.append(InfLoc[i][1])
    InfIndex = unique(InfIndex)
    InfColumn = unique(InfColumn)
    InfIndexName = Instances[InfIndex]
    InfColumnName = Features[InfColumn]

    NanInfColumn = unique(NanColumn + InfColumn)
    NanInfIndex = unique(NanIndex + InfIndex)

    NanInfColumnName = Features[NanInfColumn]
    NanInfIndexName = Instances[NanInfIndex]
    
    NanInfLoc = NanLoc + InfLoc
    NanInfLocNames = NanLocNames + InfLocNames

    return NanIndexName, NanColumnName, InfIndexName, InfColumnName, NanInfColumnName, NanInfIndexName, NanInfLocNames

def get_binary_columns(feature_array):
    binary_columns_indices = []
    for i in range(feature_array.shape[1]):
        column = feature_array.iloc[:, i]
        is_binary = np.all((column == 0) | (column == 1) | (np.isnan(column)))
        if is_binary:
            binary_columns_indices.append(i)
    binary_columns = feature_array.iloc[:, binary_columns_indices]
    mask = np.ones(feature_array.shape[1], np.bool)
    mask[binary_columns_indices] = 0
    non_binary_columns = feature_array.loc[:, mask]
    return binary_columns, non_binary_columns


# UNIQUE
def unique(list1): 
    list_set = set(list1) 
    unique_list = (list(list_set)) 
    return(unique_list)


# INTERSECTION
def intersection(List1, List2):
    Set1 = set(List1)
    Set2 = set(List2)
    return Set1.intersection(Set2)