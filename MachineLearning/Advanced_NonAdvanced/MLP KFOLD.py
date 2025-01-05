import pandas as pd
import pylab as pl
import numpy as np
import scipy.optimize as opt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt

from sklearn.metrics import classification_report, confusion_matrix
import itertools


# cell_df = pd.read_csv("cell_samples.csv")
df = pd.read_csv('Advanced_or_NonAdvanced_NAFLD_data.csv')

df = df.sample(frac = 1 , random_state=50)

df.head()
# df['Simplified_class'].value_counts()
# df['Simplified_class'].hist()


# ax = df[df['Simplified_class'] == 4][0:50].plot(kind='scatter', x='Age', y='Diabet', color='DarkBlue', label='malignant');
# df[df['Simplified_class'] == 2][0:50].plot(kind='scatter', x='Age', y='Diabet', color='Yellow', label='benign', ax=ax);
# plt.show()

from sklearn.preprocessing import LabelEncoder
lblSex = LabelEncoder()
df.SEX = lblSex.fit_transform(df.SEX)


lblDiabet = LabelEncoder()
df.Diabet = lblDiabet.fit_transform(df.Diabet)

lblClass = LabelEncoder()
df.Simplified_class = lblClass.fit_transform(df.Simplified_class)
print(df.head(10))

X = df.drop(['Simplified_class', 'Patient_ID','Run'], axis=1) 

y = df['Simplified_class'].values
print(X)

# feature_df = cell_df[['Clump', 'UnifSize', 'UnifShape', 'MargAdh', 'SingEpiSize', 'BareNuc', 'BlandChrom', 'NormNucl', 'Mit']]
# X = np.asarray(feature_df)

# ****** Feature Selection *******

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import chi2
from sklearn.feature_selection import mutual_info_classif



print( "X before feature selection: " ,X.shape)

# edit top_best for selected best features

top_best = 5

X_new = SelectKBest( f_classif, k=top_best)


X_new.fit_transform( X, y)

cols_idxs = X_new.get_support(indices=True)
features_df_new = X.iloc[:,cols_idxs]
print( "X After feature selection: " ,features_df_new.shape)
print(features_df_new)

save_features_name = features_df_new
#features_df_new.to_csv('fclassif_features'+str(top_best)+'.csv')


print(X)
print('...')
print(features_df_new)
# Data Normalization
scaler = preprocessing.StandardScaler().fit(features_df_new)
features_df_new = scaler.transform(features_df_new.astype(float))
print(features_df_new[0:5])
print(y)



df['Simplified_class'] = df['Simplified_class'].astype('int')
y = np.asarray(df['Simplified_class'])


X_train, X_test, y_train, y_test = train_test_split( features_df_new, y, test_size=0.2, random_state=4)
print ('Train set:', X_train.shape,  y_train.shape)
print ('Test set:', X_test.shape,  y_test.shape)


from sklearn import svm
from sklearn.neural_network import MLPClassifier

clf = MLPClassifier(solver='sgd', alpha=1e-5,hidden_layer_sizes=(300), random_state=4)

mychart = []


from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from numpy import mean
from numpy import std

for x in range(0,100):

    cv = KFold(n_splits=10,  shuffle=True)
    # create model

    # evaluate model
    scores = cross_val_score(clf, features_df_new, y, scoring='accuracy', cv=cv, n_jobs=-1)
    # report performance
    print('Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
    mychart.append(mean(scores))

import matplotlib.pyplot as plt


x = np.arange(start=1, stop=101, step=1)
y = np.array(mychart)


plt.title("Line graph")
plt.plot(x, y, color="red")

plt.show()

input()
