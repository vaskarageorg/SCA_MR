# SCA_MR
Sparse Component Analysis in Mendelian Randomization

This repository includes code on performing sparse component analysis (SCA) [https://arxiv.org/abs/2007.00596], a sparse dimensionality reduction approach shown to perform well in highly correlated data, in summary data of genetic variant-exposure associations, and then use this to investigate potential causal associations with Mendelian Randomization. SCA outperformed other sparse modalities in the authors' simulation studies. 

We provide a function that 
a) receives the SNP-exposure and SNP-outcome effect sizes and corresponding standard errors,
b) performs SCA of the SNP-exposure associations and,
c) in a second step, performs an inverse-variance weighted meta-analysis of the SCA-transformed SNP-exposure data and the SNP-outcome association, in line with the two-sample MR approach [https://academic.oup.com/hmg/article/23/R1/R89/2900899].

Function Information:
a) Input: Harmonised SNP-X and SNP-Y effect estimates should be provided
b) Sparsity: The estimation of the sparsity parameter is based on publicly available code from the Witten and Tibshirani _PMA_ R package [https://cran.r-project.org/web/packages/PMA/PMA.pdf]. It involves an n-fold cross-validation and a choice of the sparsity parameter that minimises the sum of the squared errors in the left-out data.
c) Causal Estimation: The effect sizes are not readily interpretable as the linear transformations skew the magnitude of the results. Interpretation of the direction and significance of the result, rather than the magnitude of the point estimate, are more intuitive.
d) It is helpful to visualise the transformation procedure with heatmaps.
