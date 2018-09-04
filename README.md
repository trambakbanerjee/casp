<!-- README.md is generated from README.Rmd. Please edit that file -->
casp
====
<!-- [![Build Status](https://travis-ci.org/trambakbanerjee/asus.svg?branch=master)](https://travis-ci.org/trambakbanerjee/asus)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/asus)](https://cran.r-project.org/package=asus)
![](http://cranlogs.r-pkg.org/badges/grand-total/asus)-->
This is an R package for Coordinate-wise Adaptive Shrinkage Prediction (casp) in a high-dimensional non-exchangeable hierarchical Gaussian model with unknown location as well as unknown spiked covariance structure. CASP [1] uses results on the behavior of eigenvalues and eigenvectors of high-dimensional sample covariance matrix ([2], [3], [4]) to develop a bias-correction principle that leads to an efficient approach for evaluating the Bayes predictors corresponding to popular loss functions such as
quadratic, generalized absolute, and linex losses.
 
How to use this repository?
----------

Soon this repository will hold the R-package casp which is currently under development. Meanwhile you can use the scripts in the folder 'Numerical Experiments' to reproduce the analysis in our paper.  

References
=======
[1.] Improved Shrinkage Prediction under a Spiked Covariance Structure    
Banerjee, T., Mukherjee, G. and Paul, D.

[2.] Baik, J. and J. W. Silverstein (2006). Eigenvalues of large sample covariance matrices of spiked population models. Journal of Multivariate Analysis 97(6), 1382–1408.

[3.] Onatski, A. (2012). Asymptotics of the principal components estimator of large factor models with weakly influential factors. Journal of Econometrics 168(2), 244–258.

[4.] Paul, D. (2007). Asymptotics of sample eigenstructure for a large dimensional spiked covariance model. Statistica Sinica, 1617–1642.

