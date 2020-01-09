# Classification of QSAR Biodegradation data set 

### Introduction
Quantitative Structure Activity Relationship (QSAR) biodegradation data set.
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
was some imbalance of 66% NRB and 34% RB.

Using the *train set*, the hyper parameters of logistic lasso, logistic ridge, and 
radial svm were tuned with 10-fold CV. The minimum CV error rates were extracted. 
The fitted models were used to capture train error, test error, false positive (fp)
train error, false negative (fn) train error, fp test error, and fn test error rates. 

As shown in the figure below, of the 5 classification methods rf and svm models have the 
lowest train error rates with the *train set* data. However, these two models perform much 
worse with *test set* data, while the rest of the models, although have more spread
out error rates, on average have error rate similar to those with their *train set* data.
The box plot of **test fn errors** stands out, here, svm and rf have the worst performance
and should not be used for identification of RB molecules, instead **logistic ridge** 
appears to have the best overall performance for purpose of identifying RB molecules.

![Error Rates](https://github.com/asyakhl/QSAR_classification/blob/master/img/Error_Rates.png)

### Minimum Cross-Validation Error Rates
The following figure shows the spread of minimum CV error rates for logistic lasso, 
logistic ridge, and svm.

<img src="/img/min_cv_errors.png" width="400">


### Cross-Validation Error vs L2-Norm of Beta Coefficients for Logistic Lasso and Ridge
From the figure below, the smallest cross-validation errors for logistic lasso and 
logistic ridge are similar to the cv error of unregularized logistic regression.  

<img src="/img/cv_errors_vs_coefficients.png" width="400">

### Heatmap of radial SVM CV Error Rates
The heatmap is used to determine the SVM parameters (gamma and cost) that give the 
smallest cross-validation error rate. In the case of QSAR biodegradation data set,
the best cost parameter appears to be equal to 100 and the best gamma parameter is 0.1.

<img src="/img/heatmap_svm_parameters.png" width="400">

### Variable Importance of Logistic Lasso, Logistic Ridge, and SVM Methods
V19 was too sparse and was not used with any of the classification methods. 
Feature definitions can be found [here.](https://archive.ics.uci.edu/ml/datasets/QSAR+biodegradation) 
The logistic lasso and logistic ridge have similar patterns for coefficient importance.
However, logistic lasso tends to emphasize some 
coefficients and reduce others to zero, while logistic ridge tends to reduce all 
coefficients in proportion to their importance. Hence, lasso coefficients are either 
large or small, while those of ridge are somewhere inbetween. The variable importance 
pattern of random forests method is completely different from the patterns of logistic 
lasso and logistic ridge, because random forests method is non-linear method and logistic 
regression is a linear method, two different methods produce two different patterns for 
variable importance. 

<img src="/img/lasso_parameters_reformatted.png" width="430">   <img src="/img/ridge_parameters_reformatted.png" width="430">

<img src="/img/rf_parameters_reformatted.png" width="430">