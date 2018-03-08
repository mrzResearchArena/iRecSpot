# Avoiding warning
import warnings
def warn(*args, **kwargs): pass
warnings.warn = warn
# _______________________________


# Essential Library
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# _____________________________



# Step 01 : Load the dataset :
iRec = 'selectedDataset.csv' # 425
D = pd.read_csv(iRec, header=None)  # Using pandas
# ___________________________________________________________________________


# Step 02 : Divide features (X) and classes (y) :
# ___________________________________________________________________________
X = D.iloc[:, :-1].values
y = D.iloc[:, -1].values


### Remove columns there is all zero values.
v = []
for i in range(X.shape[1]):
    if not np.all(X[:, i] == 0):
        v.append(i)

X = X[:, v]


from sklearn.utils import shuffle
X, y = shuffle(X, y)  # Avoiding bias


# Step 06 : Scaling the feature
# ______________________________________________________________________________
from sklearn.preprocessing import StandardScaler
scale = StandardScaler()
X = scale.fit_transform(X)

# ______________________________________________________________________________
# Step 04 : Encoding y :
# ______________________________________________________________________________
from sklearn.preprocessing import LabelEncoder
y = LabelEncoder().fit_transform(y)
# ______________________________________________________________________________


# scikit-learn :
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import BaggingClassifier, \
    RandomForestClassifier, \
    AdaBoostClassifier, \
    GradientBoostingClassifier

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

Names = ['LR', 'KNN', 'DT', 'NB', 'Bagging', 'RF', 'AdaBoost', 'GBoost', 'SVM', 'LDA']

Classifiers = [

    LogisticRegression(penalty='l2', C=0.10, max_iter=500, solver='sag'),  #1 (95.14% ***)
    KNeighborsClassifier(n_neighbors=7), #2
    DecisionTreeClassifier(), #3
    GaussianNB(), #4
    BaggingClassifier(), #5
    RandomForestClassifier(), #6
    AdaBoostClassifier(), #7
    GradientBoostingClassifier(), #8
    SVC(kernel='rbf', probability=True), #9
    LinearDiscriminantAnalysis(), #10
    ###
]


def runClassifiers():
    Results = []  # compare algorithms

    from sklearn.metrics import accuracy_score, \
        confusion_matrix, \
        roc_auc_score,\
        average_precision_score,\
        roc_curve, f1_score, recall_score, matthews_corrcoef, auc

    # Step 05 : Spliting with 10-FCV :
    from sklearn.model_selection import StratifiedKFold

    cv = StratifiedKFold(n_splits=10, shuffle=True)

    for classifier, name in zip(Classifiers, Names):

        accuray = []
        auROC = []
        avePrecision = []
        F1_Score = []
        AUC = []
        MCC = []
        Recall = []

        mean_TPR = 0.0
        mean_FPR = np.linspace(0, 1, 100)

        CM = np.array([
            [0, 0],
            [0, 0],
        ], dtype=int)

        print(classifier.__class__.__name__)

        model = classifier
        for (train_index, test_index) in cv.split(X, y):

            X_train = X[train_index]
            X_test = X[test_index]

            y_train = y[train_index]
            y_test = y[test_index]

            model.fit(X_train, y_train)


            # Calculate ROC Curve and Area the Curve
            y_proba = model.predict_proba(X_test)[:, 1]
            FPR, TPR, _ = roc_curve(y_test, y_proba)
            mean_TPR += np.interp(mean_FPR, FPR, TPR)
            mean_TPR[0] = 0.0
            roc_auc = auc(FPR, TPR)
            ##########################################
            # print(FPR)
            # print(TPR)
            ##########################################

            y_artificial = model.predict(X_test)

            auROC.append(roc_auc_score(y_test, y_proba))

            accuray.append(accuracy_score(y_pred=y_artificial, y_true=y_test))
            avePrecision.append(average_precision_score(y_test, y_proba)) # auPR
            F1_Score.append(f1_score(y_true=y_test, y_pred=y_artificial))
            MCC.append(matthews_corrcoef(y_true=y_test, y_pred=y_artificial))
            Recall.append(recall_score(y_true=y_test, y_pred=y_artificial))
            AUC.append(roc_auc)

            CM += confusion_matrix(y_pred=y_artificial, y_true=y_test)

        accuray = [_*100.0 for _ in accuray]
        Results.append(accuray)

        mean_TPR /= cv.get_n_splits(X, y)
        mean_TPR[-1] = 1.0
        mean_auc = auc(mean_FPR, mean_TPR)
        plt.plot(
            mean_FPR,
            mean_TPR,
            linestyle='-',
            label='{} ({:0.3f})'.format(name, mean_auc), lw=2.0)

        print('Accuracy: {0:.4f} %'.format(np.mean(accuray)))
        # print('auROC: {0:.6f}'.format(np.mean(auROC)))
        print('auROC: {0:.6f}'.format(mean_auc))
        print('auPR: {0:.4f}'.format(np.mean(avePrecision))) # average_Precision
        print('F1-score: {0:.4f}'.format(np.mean(F1_Score)))
        print('MCC: {0:.4f}'.format(np.mean(MCC)))
        # print('average_AUC:', np.mean(AUC))
        # tn, fp, fn, tp = CM.ravel()
        TN, FP, FN, TP = CM.ravel()
        print('Recall: {0:.4f}'.format( np.mean(Recall)) )
        # print('AUC: {0:.4f}'.format( np.mean(AUC)) )
        print('Sensitivity (+): {0:.4f} %'.format( float( (TP) / (TP + FN) )*100.0 ))
        print('Specificity (-): {0:.4f} %'.format( float( (TN) / (TN + FP) )*100.0 ))
        print('Confusion Matrix:')
        print(CM)

        print('_______________________________________')

    ### auROC Curve ###
    auROCplot()

    # boxplot algorithm comparison
    boxPlot(Results, Names)
    ### --- ###



def boxPlot(Results, Names):
    ### Algoritms Comparison ###
    # boxplot algorithm comparison
    fig = plt.figure()
    # fig.suptitle('Classifier Comparison')
    ax = fig.add_subplot(111)
    ax.yaxis.grid(True)
    plt.boxplot(Results, patch_artist=True, vert = True, whis=True, showbox=True)
    ax.set_xticklabels(Names)
    plt.xlabel('\nName of Classifiers')
    plt.ylabel('\nAccuracy (%)')

    plt.savefig('AccuracyBoxPlot.png', dpi=100)
    plt.show()
    ### --- ###


def auROCplot():
    ### auROC ###
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k', label='Random')
    plt.xlim([0.0, 1.00])
    plt.ylim([0.0, 1.02])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    # plt.title('Receiver Operating Characteristic (ROC)')
    plt.legend(loc='lower right')

    plt.savefig('auROC.png', dpi=100)
    plt.show()
    ### --- ###


if __name__ == '__main__':
    runClassifiers()


