import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score

from sklearn.metrics import classification_report, confusion_matrix
import itertools

def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')


df = pd.read_csv(r'C:\Users\Mohamadreza\Desktop\NAFLD\New_Version_240405\Data\Normal_or_NAFLD_data.csv')

#df = df.sample(frac = 1 , random_state=50)


print("Full DATA is :",df.describe)

print(df['SEX'].value_counts())
print(df['Diabet'].value_counts())
print(df['Simplified_class'].value_counts())





# Preprocess to change Type Values to Numerical


from sklearn.preprocessing import LabelEncoder
lblSex = LabelEncoder()
df.SEX = lblSex.fit_transform(df.SEX)


lblDiabet = LabelEncoder()
df.Diabet = lblDiabet.fit_transform(df.Diabet)

lblClass = LabelEncoder()
df.Simplified_class = lblClass.fit_transform(df.Simplified_class)
print(df.head(10))



#X = df[['SEX', 'BMI_surg','Age', 'Diabet']] .values  #.astype(float)

X = df.drop(['Simplified_class', 'Patient_ID','Run'], axis=1) 

y = df['Simplified_class'].values
print(X)

# ****** Feature Selection *******

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import chi2
from sklearn.feature_selection import mutual_info_classif



print( "X before feature selection: " ,X.shape)

# CHI2 method

# X_new = SelectKBest( chi2, k=5)
# X[X< 0] = 0

# mutual info
top_best = int(input("please enter number of features that you want to be selected:  "))
X_new = SelectKBest( mutual_info_classif, k=top_best)


X_new.fit_transform( X, y)

cols_idxs = X_new.get_support(indices=True)
features_df_new = X.iloc[:,cols_idxs]
print( "X After feature selection: " ,features_df_new.shape)
print(features_df_new)

save_features_name = features_df_new


#features_df_new.to_csv('NORMALorNAFLD_mutual_info_features'+str(top_best)+'.csv')


# Data Normalization
scaler = preprocessing.StandardScaler().fit(features_df_new)
features_df_new = scaler.transform(features_df_new.astype(float))
print(features_df_new[0:5])
print(y)






from sklearn.model_selection import GridSearchCV, train_test_split
X_train, X_test, y_train, y_test = train_test_split( features_df_new, y, test_size=0.2, random_state=4)
print ('Train set:', X_train.shape,  y_train.shape)
print ('Test set:', X_test.shape,  y_test.shape)

from sklearn.neighbors import KNeighborsClassifier
# k = 4

# #Train Model and Predict  
# neigh = KNeighborsClassifier(n_neighbors = k).fit(X_train,y_train)


# HYPER PARAMETER

kf=KFold(n_splits=5,shuffle=True,random_state=42)
parameter={'n_neighbors': np.arange(2, 30, 1)}
knn=KNeighborsClassifier()
knn_cv=GridSearchCV(knn, param_grid=parameter, cv=kf, verbose=1)
knn_cv.fit(X_train, y_train)
print(knn_cv.best_params_)

best_K = knn_cv.best_params_

knn=KNeighborsClassifier(n_neighbors=best_K["n_neighbors"])
knn.fit(X_train, y_train)
y_pred=knn.predict(X_test)
accuracy_score=accuracy_score(y_test, y_pred)*100
print("Accuracy for testing dataset after tuning : {:.8f}%".format(accuracy_score) + " features: " + str(top_best) + " - best K: "+ str(best_K["n_neighbors"]) )

# END OF HYPER PARAMETER

# Compute confusion matrix
neigh = KNeighborsClassifier(n_neighbors = best_K).fit(X_train,y_train)

yhat = neigh.predict(X_test)
cm = confusion_matrix(y_test, yhat)


cm_df = pd.DataFrame(cm,
                     index = ['NAFLD(0)','Normal(1)'], 
                     columns = ['NAFLD(0)','Normal(1)'])

#Plotting the confusion matrix
import seaborn as sns
plt.figure(figsize=(5,4))
sns.heatmap(cm_df, annot=True)
plt.title('Confusion Matrix KNN - feaures: '+ str(top_best))
plt.ylabel('Actal Values')
plt.xlabel('Predicted Values')
plt.show()


mychart = []


from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from numpy import mean
from numpy import std
for x in range(0,100):
    cv = KFold(n_splits=10,  shuffle=True)
    # create model
    neigh = KNeighborsClassifier(n_neighbors = k).fit(X_train,y_train)
    # evaluate model
    scores = cross_val_score(neigh, features_df_new, y, scoring='accuracy', cv=cv, n_jobs=-1)
    # report performance
    print('Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
    mychart.append(mean(scores))

print("Mean of 100 K-fold iteration: " + str(mean(mychart)))


x = np.arange(start=1, stop=101, step=1)
y = np.array(mychart)


plt.title("Line graph")
plt.plot(x, y, color="green")

plt.show()

input()