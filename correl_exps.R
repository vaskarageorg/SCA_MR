correl_expos=function(nX,nG,nbl1,n_ind,sum_stats=T){
  maf= runif(n = nG,min = 0.01,max = 0.5)
  
  G = matrix(nrow=n_ind,ncol = nG)
  for (nh in 1:ncol(G)){ G[,nh] = matrix(rbinom(n_ind, size = 2, prob = maf[nh]),ncol = 1)}
  
  U = matrix(rnorm(n = n_ind, 0, 5))
  
  epsilon_X=matrix(rnorm(n = n_ind*nX,0,4),ncol = nX)
  epsilon_y=matrix(rnorm(n = n_ind*1,0,3),ncol = 1)
  
  n_blocks=nbl1;bre1 = seq(1, nG, nG %/% n_blocks);bre1 = append(x = bre1,values = nG)
  
  list_causal=list()
  for (jl in 2:length(bre1)){ list_causal[[jl-1]] = seq(bre1[jl-1],bre1[jl]) }
  onelist1=which(unlist(lapply(list_causal,length))==1)
  if (length(onelist1)!=0){list_causal = list_causal[-onelist1]}
  
  gamma_SNP_X=matrix(nrow = nG,ncol = nX)
  gamma_SNP_X_mu=seq(4,5,length.out = length(list_causal))
  
  X_per_block1=seq(1, nX, nX %/% n_blocks)
  X_per_block1 = append(x = X_per_block1,values = nX)
  list_causalX=list()
  for (jl in 2:length(X_per_block1)){ list_causalX[[jl-1]] = seq(X_per_block1[jl-1],X_per_block1[jl]) }
  onelist1=which(unlist(lapply(list_causalX,length))==1)
  if (length(onelist1)!=0){list_causalX = list_causalX[-onelist1]}
  
  while (length(list_causalX) > length(list_causal)){
    list_causalX[[length(list_causalX)-1]] = append(list_causalX[[length(list_causalX)-1]],
                                                    list_causalX[[length(list_causalX)]])
    list_causalX[[length(list_causalX)]] =NULL
  }
  
  for (jm in 1:length(list_causalX)){ 
    gamma_SNP_X[,list_causalX[[jm]]] = matrix(rnorm(n = nG*length(list_causalX[[jm]]), 
                                                    mean = gamma_SNP_X_mu[[jm]], sd = 1),
                                              ncol = length(list_causalX[[jm]]))}
  
  X = matrix(nrow = n_ind,ncol = nX)
  for (lk in 1:length(list_causalX)){
    X[,list_causalX[[lk]]] = G[,list_causal[[lk]]] %*% gamma_SNP_X[list_causal[[lk]],list_causalX[[lk]]] + epsilon_X[,list_causalX[[lk]]]}
  
  X=apply(X,2, function(x) x+U)
  
  if (sum_stats==T){
    #Summary Statistics
    beta_matrix = SE_matrix = pval_matrix = data.frame(matrix(nrow = ncol(G), ncol = ncol(X)))
    pv_df=data.frame(matrix(nrow = ncol(G),ncol = ncol(X)))
    pvg=paste('p_val_GX',1:ncol(X),sep='')
    betag=paste('beta_GX',1:ncol(X),sep='')
    SEg=paste('SE_GX',1:ncol(X),sep='')
    
    for (k in 1:ncol(X)){
      lm_obj1=summary(lm(X[,k] ~ .,data = data.frame(G)))$coefficients
      p_val_temp = lm_obj1[-1,4]
      pv_df[,k] =p_val_temp; colnames(pv_df)[k]=pvg[k]
      beta_matrix[,k]=lm_obj1[-1,1]
      colnames(beta_matrix)[k]=betag[k]
      SE_matrix[,k]=lm_obj1[-1,2]
      colnames(SE_matrix)[k]=SEg[k]
    }
    rownames(beta_matrix)=rownames(SE_matrix)=rownames(pv_df)=paste('SNP',1:dim(G)[2],sep = '')
    b1=beta_matrix;se1=SE_matrix
    return(list(X=X,U=U,G=G,betaX=b1,seX=se1))
  }
  return(list(X=X,U=U,G=G))}

