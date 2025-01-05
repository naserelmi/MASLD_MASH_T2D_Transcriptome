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

# SelectKBest with f_classif
X_new = SelectKBest(mutual_info_classif, k=5)
X_new.fit_transform(X, y)

cols_idxs = X_new.get_support(indices=True)
features_df_new = X.iloc[:, cols_idxs]
features_df_new = X.iloc[:, cols_idxs]
selected_feature_names = features_df_new.columns.tolist()  # Store selected feature names

print("X After feature selection:", features_df_new.shape)
print(features_df_new)

# Save selected features to CSV
features_df_new.to_csv('selected_features.csv', index=False)

# Data normalization
scaler = preprocessing.StandardScaler().fit(features_df_new)
features_df_new = scaler.transform(features_df_new.astype(float))
print(features_df_new[0:5])
print(y)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(features_df_new, y, test_size=0.2, random_state=4)
print('Train set:', X_train.shape, y_train.shape)
print('Test set:', X_test.shape, y_test.shape)

# Find the optimal K for KNN
param_grid = {'n_neighbors': np.arange(1, 31)}
knn = KNeighborsClassifier()
grid_search = GridSearchCV(knn, param_grid, cv=10, scoring='accuracy')
grid_search.fit(X_train, y_train)
print(f"Optimal number of neighbors: {grid_search.best_params_['n_neighbors']}")
print(f"Best cross-validated accuracy: {grid_search.best_score_:.2f}")

# Evaluate KNN model
neigh = KNeighborsClassifier(n_neighbors=grid_search.best_params_['n_neighbors']).fit(X_train, y_train)

# Cross-validation for accuracy
cv = KFold(n_splits=10, shuffle=True)
scores = cross_val_score(neigh, features_df_new, y, scoring='accuracy', cv=cv, n_jobs=-1)
mean_accuracy = mean(scores)
std_accuracy = std(scores)
print('Accuracy: %.3f (%.3f)' % (mean_accuracy, std_accuracy))
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, cross_val_predict, KFold
# Plot mean confusion matrix
y_pred = cross_val_predict(neigh, features_df_new, y, cv=cv)
conf_matrix = confusion_matrix(y, y_pred)
disp = ConfusionMatrixDisplay( confusion_matrix=conf_matrix)
disp.plot(cmap=plt.cm.Blues)
plt.title('Mean Confusion Matrix')
plt.savefig('mean_confusion_matrix.png')
plt.show()

# Plot bar chart for selected features
for feature in selected_feature_names:  # Use stored feature names
    plt.figure(figsize=(10, 6))
    
    # Create box plot with custom colors
    sns.boxplot(x='Simplified_class', y=feature, data=df, palette={0: 'blue', 1: 'red'})
    
    # Perform t-test for p-value
    groups = df.groupby('Simplified_class')[feature]
    t_stat, p_val = ttest_ind(groups.get_group(0), groups.get_group(1))
    
    # Set custom titles and labels
    plt.title(f'{feature} by Group (p-value: {p_val:.3f})')
    plt.xlabel('Group')  # x-axis title
    plt.ylabel(feature)  # y-axis title
    
    # Customize x-axis tick labels to show group names
    plt.xticks(ticks=[0, 1], labels=['Group 1', 'Group 2'])
    
    # Save and show the plot
    plt.savefig(f'{feature}_box_plot.png')
    plt.show()