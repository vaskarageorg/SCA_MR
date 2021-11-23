source("/lustre/home/vk282/simul_NatC/511/allinc/dimred_fns.R")
source("/lustre/home/vk282/simul_NatC/511/allinc/NULL/outc_gen_NULL.R")
source("/lustre/home/vk282/simul_NatC/511/allinc/correl_exps.R")
library(MendelianRandomization);library(ggplot2);library(epca);library(elasticnet);library(mada);library(tidyverse);library(pcaPP);library(MVMR)

res_list= list();nSim = 1200;Fstat=CFstat=vector(length = nSim)

for (l in 1:nSim){
#First Sample
#generate 
nX = sample(seq(41,97),size = 1);nG = sample(seq(153,99),size = 1);nbl1 = sample(5:8,1)

#First Sample, SNP-X
X_gen=correl_expos(nX = nX,nG = nG,nbl1 = nbl1,n_ind = n_ind,sum_stats=T)

#Second Sample
X_gen1=correl_expos(nX = nX,nG = nG,nbl1 = nbl1,n_ind = n_ind,sum_stats=F)
Y_gen1=outcome_gen(X = X_gen1$X, G = X_gen1$G,
                   nbl = nbl1,U = X_gen1$U)

#Methods
pcb1=prcomp(X_gen$betaX,scale. = T)
#pcb_original=prcomp(X,scale. = T)
##Karlis Saporta Spinakis Criterion
kaiser_v=1+2*sqrt((ncol(X_gen$betaX)-1)/(nrow(X_gen$betaX)-1))
N_PCA=which(pcb1$sdev^2 > kaiser_v)

spca_res1=apply_spca_fns(betaX = X_gen$betaX,N_PCA = N_PCA)
scind=which(str_detect(names(spca_res1),'scores'))
load_ind=which(str_detect(names(spca_res1),'loadings'))

spca_res_scores=spca_res1[scind]
spca_res_load=spca_res1[load_ind]


# Results
spca_mv_res = data.frame(matrix(nrow=N_PCA,ncol=length(spca_res_scores)))
rownames(spca_mv_res)=paste('pc',1:N_PCA,sep='')
colnames(spca_mv_res)=str_replace(names(spca_res_scores),'_scores','')
for (ns in 1:length(spca_res_scores)){
  mv_res1=mr_mvivw(mr_mvinput(bx = as.matrix(spca_res_scores[[ns]]),
                      bxse = matrix(0,nrow = nG,ncol = N_PCA),
                      by = as.numeric(Y_gen1$betaY),byse = as.numeric(Y_gen1$seY)))
  spca_mv_res[,ns]=mv_res1@Pvalue }

pval_Bonf = 0.05 / nX
spca_mv_resTF=spca_mv_res<pval_Bonf
methods_res=data.frame(matrix(nrow=length(spca_res_scores)+2,ncol = 4))
rownames(methods_res)=c(str_replace(names(spca_res_scores),'_scores',''),'MVMR','MVMR_B')
colnames(methods_res)=c('TP','TN','FP','FN')
pc_above_cutoff=function(x){which(abs(x)>quantile(x,probs = 0.95))}
noncausal_vec = seq(1,nX)[-Y_gen1$causal_vec]

for (nr1 in 1:length(spca_res_scores)){
  if (nr1!=1){
    load_tf=spca_res_load[[nr1]]!=0
    pcsigt=which(spca_mv_resTF[,nr1])
    positive_spca=which(apply(data.frame(load_tf[,pcsigt]), 1, function(x) sum(x))>0)
    
    if (length(positive_spca)!=0){negative_spca=seq(1,nX)[-positive_spca]}
    if (length(positive_spca)==0){negative_spca=seq(1,nX)}
    
    TP = length(intersect(positive_spca,Y_gen1$causal_vec))
    TN = length(intersect(noncausal_vec,negative_spca))
    FN = length(intersect(negative_spca,Y_gen1$causal_vec))
    FP = length(intersect(positive_spca,noncausal_vec))
    
    methods_res[nr1,]=c(TP,TN,FP,FN) }
  if (nr1==1){
    positive_spca=sort(unique(unlist(as.vector(apply(data.frame(spca_res_load[[nr1]]),2,pc_above_cutoff)))))
    negative_spca=seq(1,nX)[-positive_spca]
    
    TP = length(intersect(positive_spca,Y_gen1$causal_vec))
    TN = length(intersect(noncausal_vec,negative_spca))
    FN = length(intersect(negative_spca,Y_gen1$causal_vec))
    FP = length(intersect(positive_spca,noncausal_vec))
    
    methods_res[nr1,]=c(TP,TN,FP,FN)
    }}

MV11=mr_mvivw(mr_mvinput(bx = as.matrix(X_gen$betaX),bxse = as.matrix(X_gen$seX),
                    by = as.numeric(Y_gen1$betaY),byse = as.numeric(Y_gen1$seY)))
sigMV=which(MV11@Pvalue < 0.05)
nonsigMV=seq(1,nX)[-sigMV]
sigMVB=which(MV11@Pvalue < pval_Bonf)
if (length(sigMVB)==0) {
  nonsigMVB=seq(1,nX)}
if (length(sigMVB)!=0) {
  nonsigMVB= seq(1,nX)[-sigMVB]}

TP = length(intersect(sigMV,Y_gen1$causal_vec))
TN = length(intersect(noncausal_vec,nonsigMV))
FN = length(intersect(nonsigMV,Y_gen1$causal_vec))
FP = length(intersect(sigMV,noncausal_vec))
methods_res[5,]  = c(TP,TN,FP,FN)
TP = length(intersect(sigMVB,Y_gen1$causal_vec))
TN = length(intersect(noncausal_vec,nonsigMVB))
FN = length(intersect(nonsigMVB,Y_gen1$causal_vec))
FP = length(intersect(sigMVB,noncausal_vec))
methods_res[6,]  = c(TP,TN,FP,FN)

res_list[[l]] = methods_res
print(l/nSim)
#Plotting intermediate results
list_sns=lapply(res_list, function(x) x[,1]/(x[,1]+x[,4]))
list_spc=lapply(res_list, function(x) x[,2]/(x[,2]+x[,3]))

sns23=t(matrix(unlist(list_sns),nrow=nrow(res_list[[1]])))
spc23=t(matrix(unlist(list_spc),nrow=nrow(res_list[[1]])))
colnames(spc23)=colnames(sns23)=rownames(res_list[[1]])

df3_sns=data.frame(sns23) %>% dplyr::select(pc,sca,spca_en,rob,MVMR,MVMR_B) %>% 
  pivot_longer(pc:MVMR_B, names_to = "Method", values_to = "Sensitivity")
df3_spc=data.frame(spc23) %>% dplyr::select(pc,sca,spca_en,rob,MVMR,MVMR_B) %>% 
  pivot_longer(pc:MVMR_B, names_to = "Method", values_to = "Specificity")

df3=cbind(df3_sns,Specificity=df3_spc$Specificity)
df3$FPR=1-df3$Specificity

plr2=ggplot(df3, aes(FPR, y = Sensitivity, color = Method)) + 
  geom_point(alpha=0.1) + geom_jitter()+ 
  xlim(0, 1) + ylim(0,1) + theme_classic() + geom_abline(slope = 1,intercept = 0) +
  scale_colour_brewer(palette = "Paired")

pdf('Roc_interim.pdf')
print(plr2)
dev.off()
sink('progr.txt')
print(paste(n_ind, 'Sims', l/nSim))
sink()
save(res_list, file="res_list.RData")

Fstat[l]=mean(apply((X_gen$betaX / X_gen$seX)^2,2,mean))
CFstat[l]=mean(unlist(strength_mvmr(r_input = format_mvmr(BXGs = X_gen$betaX,
                                                     BYG = Y_gen1$betaY,
                                                     seBXGs = X_gen$seX,
                                                     seBYG = Y_gen1$seY))))

}


