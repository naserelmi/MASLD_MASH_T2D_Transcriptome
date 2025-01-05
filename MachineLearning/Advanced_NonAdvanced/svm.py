import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, KFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from scipy.stats import ttest_ind
import seaborn as sns
import os

os.chdir('F:/OneDrive/Desktop/NAFLD_Project/article/scripts/Machine_Learning/Advanced or Non Advanced/Advanced or Non Advanced')

# Load the dataset
df = pd.read_csv('Advanced_or_NonAdvanced_NAFLD_data.csv')
df = df.sample(frac=1, random_state=50)

print("Full DATA is:", df.describe())
print(df['SEX'].value_counts())
print(df['Diabet'].value_counts())
print(df['Simplified_class'].value_counts())

# Preprocess to change type values to numerical
lblSex = LabelEncoder()
df['SEX'] = lblSex.fit_transform(df['SEX'])

lblDiabet = LabelEncoder()
df['Diabet'] = lblDiabet.fit_transform(df['Diabet'])

lblClass = LabelEncoder()
df['Simplified_class'] = lblClass.fit_transform(df['Simplified_class'])

# Feature selection
X = df.drop(['Simplified_class', 'Patient_ID', 'Run', 'id'], axis=1)
y = df['Simplified_class'].values

print("X before feature selection:", X.shape)

# SelectKBest with Mutual_info
X_new = SelectKBest(mutual_info_classif, k=20)
X_new.fit_transform(X, y)

cols_idxs = X_new.get_support(indices=True)
features_df_new = X.iloc[:, cols_idxs]
selected_feature_names = features_df_new.columns.tolist()  # Store selected feature names

print("X After feature selection:", features_df_new.shape)
print(features_df_new)

# Save selected features to CSV
features_df_new.to_csv('Advanced_NonAdvanced_features.csv', index=False)

# Data normalization
scaler = preprocessing.StandardScaler().fit(features_df_new)
features_df_new = scaler.transform(features_df_new.astype(float))
print(features_df_new[0:5])
print(y)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(features_df_new, y, test_size=0.2, random_state=4)
print('Train set:', X_train.shape, y_train.shape)
print('Test set:', X_test.shape, y_test.shape)

# Find the optimal svm model
# Define the SVM model with different kernels
from sklearn import svm
from sklearn.model_selection import GridSearchCV

# Define the parameter grid for SVM with different kernels
param_grid = {
    'C': [0.1, 1, 10, 100],
    'degree': [2, 3, 4],  # Used only for polynomial kernel
    'kernel': ['linear', 'poly', 'rbf', 'sigmoid']
}

# Define the SVM model
model = svm.SVC()

# Grid search for the best parameters
grid_search = GridSearchCV(model, param_grid, cv=10, scoring='accuracy')
grid_search.fit(X_train, y_train)
print(f"Best parameters for SVM: {grid_search.best_params_}")
print(f"Best cross-validated accuracy: {grid_search.best_score_:.2f}")

# Train the SVM model with the best parameters
clf = svm.SVC(kernel=grid_search.best_params_['kernel'], 
              C=grid_search.best_params_['C'], 
              degree=grid_search.best_params_.get('degree', 3))
clf.fit(X_train, y_train)

# Cross-validation for accuracy
cv = KFold(n_splits=10, shuffle=True)
scores = cross_val_score(clf, features_df_new, y, scoring='accuracy', cv=cv, n_jobs=-1)
mean_accuracy = np.mean(scores)
std_accuracy = np.std(scores)
print('Accuracy: %.3f (%.3f)' % (mean_accuracy, std_accuracy))
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, cross_val_predict, KFold
# Plot mean confusion matrix
y_pred = cross_val_predict(clf, features_df_new, y, cv=cv)
conf_matrix = confusion_matrix(y, y_pred)
disp = ConfusionMatrixDisplay( confusion_matrix=conf_matrix)
disp.plot(cmap=plt.cm.Blues)
plt.title('Mean Confusion Matrix')
plt.savefig('mean_confusion_matrix.png')
plt.show()
df['Simplified_class'] = df['Simplified_class'].astype(str)
# Plot bar chart for selected features
print(df['Simplified_class'].unique())
for feature in selected_feature_names:
    feature='ENSG00000164114'
    plt.figure(figsize=(10, 6))
    # Create box plot with custom colors
    sns.boxplot(x='Simplified_class', y=feature, data=df, palette={'0': 'blue', '1': 'red'})    # Perform t-test for p-value
    groups = df.groupby('Simplified_class')[feature]
    t_stat, p_val = ttest_ind(groups.get_group('0'), groups.get_group('1'))
    np.mean(groups.get_group('1'))/np.mean(groups.get_group('0'))
    # Set custom titles and labels
    plt.title(f'MAP9 gene expression comparison (p-value: {p_val:.4f})', fontsize=16)
    #plt.title(f'PAPLN gene expression comparison', fontsize=16)
    plt.xlabel('Group', fontsize=16)  # x-axis title
    plt.ylabel('MAP9 gene expression', fontsize=16)  # y-axis title
    # Customize x-axis tick labels to show group names
    plt.xticks(ticks=[0, 1], labels=['Non_Advanced Fibrosis', 'Advanced Fibrosis'],fontsize=14)
    # Save and show the plot
    plt.savefig(f'MAP9_box_plot.png')
    plt.show()