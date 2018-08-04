# Avoiding warning
import warnings
def warn(*args, **kwargs): pass
warnings.warn = warn
# _______________________________


# Essential Library
import pandas as pd
import numpy as np
# _____________________________


# Step 01 : Load the dataset :
iRec = 'dataset.csv' #
D = pd.read_csv(iRec, header=None)  # Using pandas
# ___________________________________________________________________________


# Step 02 : Divide features (X) and classes (y) :
# ___________________________________________________________________________
X = D.iloc[:, :-1].values
y = D.iloc[:, -1].values


# Remove columns there is all zero values.
def selectKImportance(model, X):
    Values = np.sort(model.feature_importances_)[::-1]
    return X[:, model.feature_importances_.argsort()[::-1][:len(Values[Values>0.00])]]

## Uncomment when need to know the best features index.
# from sklearn.ensemble import AdaBoostClassifier
# model = AdaBoostClassifier(algorithm='SAMME.R', n_estimators=500, learning_rate=1.0) #n_estimators=500
# model.fit(X,y)
# X = selectKImportance(model, X) #425

### iRec (425)
X = X[:, [9954, 30769, 30710, 12413, 13005, 40130, 35526, 4949, 1610, 18664, 20716, 813, 33333, 38464, 22625, 20755, 40214, 35206, 7711, 32559, 2535, 38322, 9325, 21984, 2985, 7903, 6863, 32069, 32588, 18839, 30459, 5848, 8421, 13508, 37409, 34394, 24020, 33987, 27838, 11912, 37378, 1486, 21097, 12169, 1993, 31502, 32006, 32638, 4697, 23060, 23479, 20219, 20483, 11249, 24933, 22411, 37332, 4021, 18452, 20211, 6387, 25548, 13650, 7367, 11380, 16047, 18128, 13069, 23433, 29403, 40324, 10603, 9261, 14066, 7131, 1712, 17305, 35091, 31996, 7139, 3231, 13060, 2414, 3659, 25992, 1701, 17295, 2405, 8612, 23847, 18410, 22611, 33263, 29993, 15424, 36977, 25151, 19953, 7868, 10503, 3273, 31944, 27362, 24247, 1266, 4053, 24520, 24757, 10901, 38219, 9012, 15383, 37173, 34539, 12218, 32396, 32760, 8414, 4042, 2498, 29431, 2477, 36947, 36948, 12550, 1228, 15439, 32703, 19031, 14059, 27281, 27279, 32711, 33722, 5586, 27544, 40354, 6337, 14282, 21480, 14527, 26020, 1252, 29682, 11792, 26262, 21203, 30633, 31328, 2360, 21414, 14606, 9672, 38379, 28653, 9671, 19731, 8484, 11910, 22009, 9223, 9666, 3123, 13154, 16136, 22442, 6219, 40070, 40167, 13547, 25846, 5673, 14700, 26959, 31432, 1039, 9767, 38867, 12362, 6199, 37097, 28228, 33411, 1813, 14636, 33791, 1069, 4239, 24897, 1072, 24900, 4696, 8076, 10732, 12349, 1076, 10745, 23764, 18025, 33785, 5741, 18473, 38859, 26343, 4733, 22528, 23798, 29736, 10802, 40233, 8566, 3067, 30767, 23420, 19098, 13095, 11039, 3214, 31478, 10630, 2930, 14582, 20657, 6270, 33860, 34838, 31486, 20649, 12280, 34833, 26443, 7964, 32645, 8584, 20680, 6255, 9125, 24848, 36377, 18435, 4788, 38338, 14719, 20709, 19777, 16163, 34985, 11207, 25819, 19122, 11966, 9138, 25945, 452, 4166, 18428, 8557, 32135, 14235, 7998, 14130, 19628, 14508, 4704, 25365, 18690, 5269, 2584, 2585, 12905, 37277, 19426, 13955, 31089, 10184, 24607, 25497, 28480, 15042, 20403, 12652, 16393, 6051, 3425, 19366, 35374, 26674, 22175, 28421, 27090, 7664, 5134, 22785, 22155, 24710, 27721, 31719, 23959, 23076, 11498, 6512, 30304, 30386, 21575, 13340, 2022, 17703, 18213, 5104, 16504, 6581, 14332, 30407, 30597, 778, 22177, 31679, 32236, 36690, 35516, 15191, 11694, 20308, 5331, 4385, 4490, 18897, 16372, 24162, 27904, 39092, 20296, 17698, 10090, 5177, 12661, 28442, 35473, 30487, 1568, 21795, 15966, 31202, 27671, 32258, 7786, 8708, 7487, 12196, 5442, 22213, 40435, 14301, 30121, 8681, 14930, 7701, 16313, 3727, 6997, 4427, 23555, 33156, 40660, 34745, 3298, 22737, 23669, 20229, 6675, 7437, 11622, 13406, 36304, 36607, 9313, 3924, 3856, 7831, 27979, 8886, 22730, 26695, 14005, 21851, 32781, 32969, 6552, 20976, 31159, 8769, 2010, 32955, 23097, 36106, 24003, 4533, 7796, 11050, 33012, 6648, 21559, 22209, 38571, 24011, 3888, 6596, 31631, 20325, 8233, 14395, 3849, 18661, 2543, 39207, 5210, 16933, 39209, 27834, 31828, 1451, 33693, 25587, 12616, 34628]]


def saveCSV(X, Y):
    F = open('selectedDataset.csv', 'w')

    for x, y in zip(X, Y):
        for each in x:
            F.write(str(each) + ',')
        F.write(str(y) + '\n')
    F.close()

saveCSV(X, y)


