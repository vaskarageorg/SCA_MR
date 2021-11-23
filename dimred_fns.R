SCA.cv=function (x, N_PC, length.gamma=10,nfolds = 5, 
                 center = TRUE) 
{
  #browser()
  if (nfolds < 2) 
    stop("Please enter more than two folds.")
  percentRemove <- min(0.25, 1/nfolds)
  call <- match.call()
  gamma_seq=seq(sqrt(dim(x)[2]*N_PC)/2, 
                sqrt(dim(x)[2]*N_PC) * 2,
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
      
      xhat <- as.numeric(out$sdev^2) * out$scores %*% t(out$loadings)
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

apply_spca_fns = function(betaX,N_PCA){
  
  #PCA
  pcb1=prcomp(betaX,scale. = T)
  pc_scores=pcb1$x[,1:N_PCA]
  #SCA
  scb1=suppressWarnings(suppressMessages(SCA.cv(x = matrix(unlist(betaX),nrow = nrow(betaX)),
                                                N_PC = N_PCA,
                                                length.gamma = 5,nfolds = 5,center = T)))
  scb1_fit=suppressWarnings(suppressMessages(sca(A = matrix(unlist(betaX),nrow = nrow(betaX)),
                                                 k = N_PCA,gamma = scb1$bestgamma1se)))
  #EN
  spca_en = spca(x = betaX,K = N_PCA,para = rep(round(ncol(betaX)/N_PCA),ncol(betaX)),sparse = 'varnum')
  spca_en_scores=as.matrix(betaX) %*% spca_en$loadings %*% solve((t(spca_en$loadings) %*% spca_en$loadings))
  
  #robust
  b22=opt.TPO(x = betaX,k.max = N_PCA)
  robs1=sPCAgrid(x = betaX,k = N_PCA,lambda = b22$pc.noord$lambda)
  robs1_scores=robs1$scores
  robs1_load=matrix(as.numeric(robs1$loadings),ncol=N_PCA)
  
  list(pc_scores=pc_scores, pc_loadings=pcb1$rotation,
       sca_scores=scb1_fit$scores, sca_loadings=scb1_fit$loadings,
       spca_en_scores=spca_en_scores, spca_en_loadings=spca_en$loadings,
       rob_scores=robs1_scores,rob_loadings=robs1_load)
}
