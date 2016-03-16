#Supervised PCA

When data has more features than sample size, then one can use PCA to reduce the dimension of feature space, followed by classification methods for prediction, but the problem with this approach is that PCA is unsupervised algorithm, it doesn't take into account if the variation in the feature space is related with response variable or not. Usually in genetics dataset, there is high variability in the genes but that is not even related to response variable(prostate cancer status in this case), and it is due to some internal activity. One was to approach such problem is to use Supervised PCA instead. Here, in this repository I have compared both the dimensionality reduction approaches, PCA vs Supervised PCA followed by classification methods like LDA, QDA. 

For more information about Supervised PCA: https://web.stanford.edu/~hastie/Papers/spca_JASA.pdf
