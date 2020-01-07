# QSAR_classification

### Introduction

The data set source is [UCI machine learning repository.](https://archive.ics.uci.edu/ml/datasets/QSAR+biodegradation)

The data set has 1,055 instances (molecules) and their 41 features (chemical and physical properties).
Each molecule is either readily biodegradable (RB) or not readily biodegradable (NRB).
Logistic regression, logistic with lasso, logistic with ridge, radial SVM, 
and random forests are used here to classify each molecule as RB or NRB based on 
its 41 features. 

### Train and Test Error Rates

The following was repeated 100 times to capture the spread of train error, test error, 
and minimum CV error rates for all 5 classification methods. 

Through random sampling, the data set was separated into 90% *train set* and the 
rest was used as *test set* data. The data was weighted through resampling since there 
was an imbalance of 66% NRB and 34% RB.

Using the *train set*, the hyper parameters of logistic lasso, logistic ridge, and 
radial svm were tuned with 10-fold CV. The minimum CV error rates were extracted. 
The fitted models were used to capture train error, test error, false positive (fp)
train error, false negative (fn) train error, fp test error, and fn test error rates. 



![](https://github.com/asyakhl/QSAR_classification/blob/master/img/Error_Rates.png)
