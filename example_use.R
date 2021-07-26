#An example of highly correlated SNP-exposure association estimates and how the SCA MR approach can be used

set.seed(131)

#Data Generation ---------
k_j = sample(5:10,1)
p=100;K=10*k_j; n_ind=20000; maf=0.3
G = matrix(rbinom(n_ind * p,size = 2,prob = maf),ncol = p)
Ux = matrix(rnorm(n = n_ind*K, 0, 10), ncol = K)
epsilon_X=matrix(rnorm(n = n_ind*K,0,14),ncol = K)
epsilon_y=matrix(rnorm(n = n_ind*1,0,20),ncol = 1)

#2. Generate exposure|G,Ux
## To get a block correlation we use only specific SNPs for each block
n_blocks=k_j
#Breaking into causal and non-causal blocks 
bre1 = seq(1, K, K %/% n_blocks)
list_causal=list()
for (jl in 1:length(bre1)){
  list_causal[[jl]]=(seq(bre1[jl], jl*(K %/% n_blocks)))
}

gamma_SNP_X=matrix(nrow = p,ncol = K)
gamma_SNP_X_mu=seq(4,2,length.out = length(list_causal))
for (jm in 1:length(list_causal)){
  gamma_SNP_X[,list_causal[[jm]]] = matrix(rnorm(n = p*length(list_causal[[jm]]),
                                                 mean = gamma_SNP_X_mu[[jm]],
                                                 sd = 1),
                                           ncol = length(list_causal[[jm]]))}

X=X_noise=X_cor=matrix(nrow = n_ind,ncol = K)
G_per_block=seq(1, p, p %/% n_blocks) #Breaking the genes in discreet blocks that contribute to blocks of X
list_causal_G=list()
for (jl in 1:length(G_per_block)){
  list_causal_G[[jl]]=(seq(G_per_block[jl], jl*(p %/% n_blocks)))
}

SNP_X_beta_noise=matrix(2,nrow = p,ncol = K)
for (lk in 1:length(list_causal)){
  #X_noise=X[,-list_causal[[lk]]] = G[,-list_causal_G[[lk]]] %*% gamma_SNP_X[-list_causal_G[[lk]],list_causal[[lk]]] + Ux[,list_causal[[lk]]] + epsilon_X[,list_causal[[lk]]]
  X[,list_causal[[lk]]] = G[,list_causal_G[[lk]]] %*% gamma_SNP_X[list_causal_G[[lk]],list_causal[[lk]]] + Ux[,list_causal[[lk]]] + epsilon_X[,list_causal[[lk]]]
  #X=X_noise+X_cor
}  

#3. Generate outcome|X,G,Uy with different beta's for different blocks----
#We use non-zero coefficients for the first two blocks
beta_X_Y = matrix(nrow = K*1,ncol = 1)
beta_X_Y[-list_causal[[2]],]=0
beta_X_Y[list_causal[[1]],]= runif(n = length(list_causal[[1]]),min = 0.3,max = 4.5)
beta_X_Y[list_causal[[2]],]= runif(n = length(list_causal[[1]]),min = 0.1,max = 2.5)
beta_X_Y[list_causal[[3]],]= runif(n = length(list_causal[[1]]),min = -3.6,max = -0.2)

#Within the causal block, get some to zero. Maybe they will spuriously associate with Y just because of 
#multicoll8inearity
zero_coef = sample(length(list_causal[[1]]),size = sample(1:8,1))
beta_X_Y[list_causal[[2]],][zero_coef]=0
beta_X_Y[list_causal[[3]],][zero_coef]=0

Uy_coef = matrix(runif(n = K*1,min = 10,max = 15),ncol = 1)

Y = X %*% beta_X_Y + Ux %*% Uy_coef + epsilon_y
#Summary Stats

#3b. Extracting summary statistics----
beta_matrix = SE_matrix = pval_matrix = data.frame(matrix(nrow = ncol(G), ncol = ncol(X)))
pv_df=data.frame(matrix(nrow = ncol(G),ncol = ncol(X)))
pvg=paste('p_val_GX',1:ncol(X),sep='')
betag=paste('beta_GX',1:ncol(X),sep='')
SEg=paste('SE_GX',1:ncol(X),sep='')

rsq_t=NULL
for (k in 1:ncol(X)){
  d1=summary(lm(X[,k] ~ .,data = data.frame(G)))
  rsq_t[k]=d1$r.squared
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
by1=summary(lm(Y ~ .,data = data.frame(G)))$coefficients[-1,1]
byse1=summary(lm(Y ~ .,data = data.frame(G)))$coefficients[-1,2]


#Apply Function -----
sc_obj=sca_MR(beta_X_mat = matrix(unlist(beta_matrix),nrow = nrow(beta_matrix)),
       beta_Y_vec = by1,SE_Y_vec = byse1,nfold = 5,spars_length = 5)

library(RColorBrewer); coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(25)
heatmap(sc_obj$Loadings,Colv = NA,Rowv = NA,col=coul)
sc_obj$MR_estimate

# The first three causal blocks of exposures are projected in the first three sparse 
# Components (as can be seen in the heatmap).
# This is in accordance with the known data generating process of six blocks with ten exposures in each block. 
#Additionally, the MR estimates suggest that these blocks of exposures are significantly
#associated with the outcome.
