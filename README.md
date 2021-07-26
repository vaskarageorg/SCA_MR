# SCA_MR
Sparse Component Analysis in Mendelian Randomization

This repository includes code on performing sparse component analysis (SCA) [https://arxiv.org/abs/2007.00596], a sparse dimensionality reduction approach shown to perform well in highly correlated data, in summary data of genetic variant-exposure associations, and then use this to investigate potential causal associations with Mendelian Randomization. SCA outperformed other sparse modalities in the authors' simulation studies. 

We provide a function that 
a) receives the SNP-exposure and SNP-outcome effect sizes and corresponding standard errors,
b) performs SCA of the SNP-exposure associations and,
c) in a second step, performs an inverse-variance weighted meta-analysis of the SCA-transformed SNP-exposure data and the SNP-outcome association, in line with the two-sample MR approach [10.1093/hmg/ddu328].

The estimation of the sparsity parameter 
