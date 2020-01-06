# QSAR_classification

### Introduction

The data set source is [UCI machine learning repository.](https://archive.ics.uci.edu/ml/datasets/QSAR+biodegradation)

The data set has 1,055 instances (molecules) and their 41 features (chemical and physical properties).
Each molecule is either readily biodegradable (RB) or not readily biodegradable (NRB).
Logistic regression, logistic with lasso, logistic with ridge, radial SVM, 
and random forests are used here to classify each molecule as RB or NRB based on 
its 41 features. 

The following was repeated 100 times to capture train error, test error, and minimum 
CV error spreads for all 5 classification methods. 

	Through random sampling, the data set was separated into 0.9n *train set* and the 
rest was used as *test set* data. The data was weighted through resampling since there 
was an imbalance of 66% NRB and 34% RB.
	Using the *train set*, the hyper parameters of logistic lasso, logistic ridge, and 
radial svm were tuned with 10-fold CV. For each of these methods train error, test error, 
and minimum CV error were recorded. Random forest was fitted with 300 bootstrapped
trees. Logistic was performed as usual. 
	False positive (fp) and false negative(fn) error rates were also recorded. 

![](https://github.com/asyakhl/QSAR_classification/blob/master/img/Error_Rates.png)
