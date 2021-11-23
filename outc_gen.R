outcome_gen = function(X, G,nbl,U){
  per_block=seq(1, ncol(X), ncol(X) %/% nbl) #Breaking the genes in discreet blocks that contribute to blocks of X
  per_block=append(per_block,ncol(X))
  
  list_causal_X=list()
  for (jl in 2:length(per_block)){
    list_causal_X[[jl-1]]=seq(per_block[jl-1], per_block[jl]) }
  
  beta_X_Y = matrix(nrow = ncol(X),ncol = 1)
  #num_causal_blocks=sample(1:nbl,1)
  #causal_exps=unlist(list_causal_X[1:num_causal_blocks])
  #beta_X_Y[causal_exps,]=runif(n = length(causal_exps),min = 1,max = 2)
  #beta_X_Y[-causal_exps,]=0
  beta_X_Y[1:ncol(X),]=0
  causal_1 = sample(1:ncol(X),1)
  causal_index = sample(1:ncol(X),causal_1,replace = F)
  beta_X_Y[causal_index,]=1
  
  epsilon_Y=rnorm(n=nrow(X),0,2)
  Y = X %*% beta_X_Y + U + epsilon_Y
  
  
  by1=summary(lm(Y ~ .,data = data.frame(G)))$coefficients[-1,1]
  byse1=summary(lm(Y ~ .,data = data.frame(G)))$coefficients[-1,2]
  
  return(list(Y=Y,betaY=by1,seY=byse1,causal_vec=which(beta_X_Y!=0)))
} 
