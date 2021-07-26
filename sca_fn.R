#Function to perform SCA and IVW-MR in summary-level genetic data
#SCA: https://arxiv.org/abs/2007.00596
#For the estimation of the sparsity parameter: https://cran.r-project.org/web/packages/PMA/PMA.pdf

library(epca)

sca_MR=function(beta_X_mat, beta_Y_vec, SE_Y_vec, nfold, spars_length){
SCA.cv=function (x, N_PC, length.gamma=10,nfolds = 10, 
                 center = TRUE) 
{
  #browser()
  if (nfolds < 2) 
    stop("Please enter more than two folds.")
  percentRemove <- min(0.25, 1/nfolds)
  call <- match.call()
  gamma_seq=seq(sqrt(dim(x)[2]*N_PC)/3, 
                sqrt(dim(x)[2]*N_PC)*2,
                length.out = length.gamma)
  
  xfill <- x
  missing <- is.na(x)
  
  errs <- matrix(NA, nrow = nfolds, ncol = length(gamma_seq))
  nonzerovs <- matrix(NA, nrow = nfolds, ncol = length(gamma_seq))
  rands <- matrix(runif(nrow(x)), ncol = 1)
  for (i in 1:nfolds) {
    print(paste(" Fold ", i, " out of ", nfolds))
    rm <- ((i - 1) * percentRemove < rands) & (rands < i * 
                                                 percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for (j in 1:length(gamma_seq)) {
      out=sca(A = xrm[complete.cases(xrm),],
              k = N_PC,gamma = gamma_seq[j],is.cov = F,
              rotate = 'varimax',shrink = 'soft',center = T,scale = T,
              max.iter = 2000)
      
      #xhat <- as.numeric(out$sdev^2) * out$scores %*% t(out$loadings)
      xhat <- out$scores %*% t(out$loadings)
      
      errs[i, j] <- sum(((xhat - x[!rm,]))^2)
      nonzerovs[i, j] <- sum(out$v != 0)
      
    }
  }
  
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd)/sqrt(nfolds)
  nonzerovs.mean <- apply(nonzerovs, 2, mean)
  
  bestgamma <- gamma_seq[which.min(err.means)]
  bestgamma1se <- gamma_seq[min(which(err.means < min(err.means) + 
                                        err.sds[which.min(err.means)]))]
  object <- (list(cv = err.means, cv.error = err.sds, 
                  bestgamma = bestgamma, 
                  nonzerovs = nonzerovs.mean,
                  nfolds = nfolds, 
                  bestgamma1se = bestgamma1se))
  return(object)
}

scb1=suppressWarnings(suppressMessages(SCA.cv(x = beta_X_mat,
                                              N_PC = N_PCA,nfolds = nfold,
                                              length.gamma = spars_length,center = T)))
scb1_fit=suppressWarnings(suppressMessages(sca(A = beta_X_mat,
                                               k = N_PCA,
                                               gamma = scb1$bestgamma1se)))

#IVW MR
# The obtained scores are used instead of the individual exposures in an inverse variance weighted (SE_{SNP-outcome} ^ (-2)) meta-analysis

scores_subs = scb1_fit$scores
obj_MR=summary(lm(beta_Y_vec ~ 0 + scores_subs,
           weights = (SE_Y_vec)^(-2)))$coef[,c(1,2,4)]
colnames(obj_MR) = c('Sparse Component Estimate', 'SE', 'pvalue')

#A list with the MR estimate, the SCA scores and loadings is returned 
return(list(MR_estimate=obj_MR, Scores = scores_subs, Loadings = scb1_fit$loadings))
}
