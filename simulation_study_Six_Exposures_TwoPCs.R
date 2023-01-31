#generate exposure
source('C:/Users/vk282/OneDrive - University of Exeter/Desktop/natcommsubm/review_elife/Kettunen/sim2010/dimred_fns_1smr.R')
library(epca)
n=400
sd_vec=round(seq(100,4000,length.out = 5))
#31.1.23: range of N 100-800, trying out outher N's to see if it converges
#sd_vec=seq(2,0.5,length.out = 8)

nSim=300
simil_measure=simil_measure2=NULL; simil_measure_l=simil_measure_l2=list()
fstat_l=NULL;fstat_s=NULL
g53=list()

sns_df=data.frame(matrix(nrow = nSim, ncol = ))

mvmr_list=list(); list_per_N=list()
for (sdr in 1:length(sd_vec)){
  n=sd_vec[sdr]
  #n=500
  nG=50; nX=6
  SD_i=sd_vec[sdr]; g52=data.frame(matrix(nrow = nSim,ncol = nX))
  for (l in 1:nSim){
gamma_GX=matrix(rnorm(nG*nX,0,0.3),ncol = nX)

#set to zero those that do not contribute in the block
dim(gamma_GX);gamma_GX[26:50,1:3]=0;gamma_GX[1:25,4:6]=0
#---------------------------------------------------------
G=matrix(rbinom(n*nG,2,0.2), nrow = n); G=scale(G)
U=rnorm(n,0,1);epsilonX=rnorm(n,0,1)

X=G%*%gamma_GX + U + epsilonX
X[,2]=X[,2]+0.2*X[,1]; X[,3]=X[,3]+0.3*X[,2]

X[,5]=X[,5]+0.2*X[,4]; X[,6]=X[,6]+0.3*X[,5]

beta_caus=matrix(nrow = nX,ncol = 1)
beta_caus[1:3,]=0.2; beta_caus[4:6,]=0
epsilonY=rnorm(n,0,1)
Y=X%*%beta_caus+U+epsilonY
#heatmap(cor(X))
#get xhat
X =scale(X )
Xhat=X
for (bv in 1:ncol(X) ){ df_XG=data.frame(X[,bv],G); Xhat[,bv]=predict(lm(df_XG)) }
#Xhat=scale(Xhat)

pcx1=prcomp(Xhat); N_PCA=2
gamma_hat= (matrix(nrow = nG,ncol = nX))
SE_hat= (matrix(nrow = nG,ncol = nX))
for (nvn in 1:nX){
  ss_1=summary(lm(X[,nvn]~G))$coef;gamma_hat[,nvn]=ss_1[-1,1]
  SE_hat[,nvn]=ss_1[-1,2]}

#'''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''
mvmr1=summary(lm(Y~Xhat))$coef[-1,4]
mvmr_PCA=summary(lm(Y~pcx1$x[,1:2]))$coef[-1,4]

scacv1=apply_spca_fns(X = Xhat,N_PCA = 2)
sca_1_loadings=scacv1$sca_loadings;sca_1_scores=scacv1$sca_scores

mvmr_SCA=summary(lm(Y~sca_1_scores))$coef[-1,4]

if (sum(sca_1_loadings[1:3,1]) > sum(sca_1_loadings[4:6,1])){ mvmr_SCA2=mvmr_SCA}
if (sum(sca_1_loadings[1:3,1]) < sum(sca_1_loadings[4:6,1])){ mvmr_SCA2=mvmr_SCA[c(2,1)]}
if (sum(pcx1$rotation[1:3,1]) > sum(pcx1$rotation[4:6,1])){ mvmr_PCA2=mvmr_PCA}
if (sum(pcx1$rotation[1:3,1]) < sum(pcx1$rotation[4:6,1])){ mvmr_PCA2=mvmr_PCA[c(2,1)]}

mvmr_list[[l]]=list(pca=mvmr_PCA2,sca=mvmr_SCA2,mv=mvmr1)

pcg1=prcomp(gamma_hat)


# t(Xhat) %*% Xhat; t(gamma_hat) %*% gamma_hat ; t(gamma_hat) %*% t(G) %*% G %*% gamma_hat
# t(gamma_hat/SE_hat) %*% (gamma_hat/SE_hat)

fstat_l[l]=mean(apply((gamma_hat/SE_hat)^2,2,mean))

for (gV in 1:ncol(pcx1$rotation)){
g51=summary(lm(pcx1$rotation[,gV]~0+pcg1$rotation[,gV]))
g52[l,gV]=g51$r.squared}
  
simil_measure[l]=sum((abs(pcx1$rotation[,1:N_PCA])-abs(pcg1$rotation[,1:N_PCA]))^2)
simil_measure2[l]=sum(( (pcx1$rotation[,1:N_PCA])- (pcg1$rotation[,1:N_PCA]))^2)

print(l/nSim)}
  pcsig1=t(matrix(unlist(lapply(mvmr_list, function(x) x$pca<0.05)),ncol= nSim))
  scsig1=t(matrix(unlist(lapply(mvmr_list, function(x) x$sca<0.05)),ncol= nSim))
  mvsig1=t(matrix(unlist(lapply(mvmr_list, function(x)  x$mv<0.05)),ncol= nSim))
  list_per_N[[sdr]]=list(pca=apply(pcsig1,2,mean),sca=apply(scsig1,2,mean),mv=apply(mvsig1,2,mean))
  
  simil_measure_l[[sdr]]=simil_measure
  simil_measure_l2[[sdr]]=simil_measure2
  fstat_s[sdr]=mean(fstat_l)
  g53[[sdr]]=g52
}

PCARES1=t(matrix(unlist(lapply(list_per_N,function(x) x$pca)),nrow = 2))
SCARES1=t(matrix(unlist(lapply(list_per_N,function(x) x$sca)),nrow = 2))
MVRES1=t(matrix(unlist(lapply(list_per_N,function(x) x$mv)),nrow = 6))

par(mfrow=c(1,3))
plot(SCARES1[,1],type = 'l',ylim = c(0,1),ylab='SCA Rejection Rate');points(SCARES1[,2],col='red',type = 'l')
abline(h=0.05,lty=2,col='red'); abline(h=1,lty=2,col='black')
#add monte carlo SEs from Tim Morris et al. formulae 10.1002/sim.8086 Table 6
library(tidyverse);library(patchwork)
SCARES1_SE = sqrt(SCARES1*(1-SCARES1)/nSim)
dfr1=data.frame(fstat=fstat_s,sca1=SCARES1[,1],sca2=SCARES1[,2],
           sca1se=SCARES1_SE[,1],sca2se=SCARES1_SE[,2])
dfr1e=dfr1 %>% pivot_longer(cols = sca1:sca2,values_to='Est',names_to='Method')
dfr1s=dfr1 %>% pivot_longer(cols = sca1se:sca2se,values_to='SE',names_to='Method')
dfr2=dfr1e; dfr2$SE=dfr1s$SE

j <- ggplot(dfr2, aes(fstat, Est, ymin = Est-1.96*SE, ymax = Est+1.96*SE,color=Method))
scplot=j + geom_pointrange()+geom_line() +scale_color_brewer(palette = 'Dark2')+ylab('SCA Rejection Rate')+
  geom_hline(yintercept = 0.05,linetype=2,color='chocolate1')+
  geom_hline(yintercept = 1,linetype=2,color='seagreen')

#pca plot
plot(PCARES1[,1],type = 'l',ylim = c(0,1),ylab='PCA Rejection Rate');points(PCARES1[,2],col='red',type = 'l')
abline(h=0.05,lty=2,col='red'); abline(h=1,lty=2,col='black')


PCARES1_SE = sqrt(PCARES1*(1-PCARES1)/nSim)
dfr1=data.frame(fstat=fstat_s,sca1=PCARES1[,1],sca2=PCARES1[,2],
                sca1se=PCARES1_SE[,1],sca2se=PCARES1_SE[,2])
dfr1e=dfr1 %>% pivot_longer(cols = sca1:sca2,values_to='Est',names_to='Method')
dfr1s=dfr1 %>% pivot_longer(cols = sca1se:sca2se,values_to='SE',names_to='Method')
dfr2=dfr1e; dfr2$SE=dfr1s$SE

j <- ggplot(dfr2, aes(fstat, Est, ymin = Est-1.96*SE, ymax = Est+1.96*SE,color=Method))
pcplot=j + geom_pointrange()+geom_line() + geom_hline(yintercept = 0.05,linetype=2,color='coral') +
  geom_hline(yintercept = 0.05,linetype=2,color='chocolate1')+
  scale_color_brewer(palette = 'Dark2')+ylab('PCA Rejection Rate')+
  geom_hline(yintercept = 1,linetype=2,color='seagreen')


(scplot|pcplot)+plot_layout(guides = 'collect')
#mvmr plot
plot(MVRES1[,4],type = 'l',ylim = c(0,1),ylab='MVMR Rejection Rate',col='brown');points(MVRES1[,5],col='coral',type = 'l')
points(MVRES1[,6],col='brown1',type = 'l'); points(MVRES1[,1],col='darkgray',type = 'l')
points(MVRES1[,2],col='gray26',type = 'l'); points(MVRES1[,3],col='gray8',type = 'l')
abline(h=0.05,lty=2,col='red'); abline(h=1,lty=2,col='black')

###############################################################################
library(tidyverse);library(patchwork);library(RColorBrewer)
s2f=unlist(lapply(g53,function(x) apply(x,2,mean)))
s3f=data.frame(t(matrix(s2f,ncol = length(sd_vec))))
colnames(s3f)=paste('PC',1:ncol(s3f),sep='')
s3f$Fstat=fstat_s
s4f=s3f %>% pivot_longer(cols = PC1:PC6)

i <- ggplot(s4f, aes(log(Fstat)/log(10) , value,colour=name))
r2pl=i+geom_line()+ylab('R2 of Loadings')+xlab('Log10 F-statistic')+scale_color_brewer(palette = 'Dark2')
ggsave('corr_R2_stronger410.pdf')
###############################################################################
df1=data.frame(log_F=log(fstat_s)/log(10), simil=unlist(lapply(simil_measure_l, mean)),
           similSD=unlist(lapply(simil_measure_l, function(x) sd(x)/sqrt(length(x)))))
i <- ggplot(df1, aes(log_F, simil,ymin = simil-1.96*similSD, ymax = simil+1.96*similSD))
similpl=i + geom_line()+geom_pointrange()+xlab('Log10 F-statistic')+ylab('Loadings Similarity, 1SMR and 2SMR')
(r2pl/similpl)
ggsave('C:/Users/vk282/OneDrive - University of Exeter/Desktop/natcommsubm/review_elife/Kettunen/similplots.pdf')

plot(log(fstat_s),unlist(lapply(simil_measure_l, mean)),type='l')
plot(log(fstat_s),unlist(lapply(simil_measure_l2, mean)),type='l')


plot(density(simil_measure_l[[1]]),title='Similarity Measure, SCA'
     , ylim=c(0,7),xlim=c(0,2) ,xlab='Similarity Measure for 1SMR and 2SMR Loadings, SCA',main = ''
)
points(density(simil_measure_l[[2]]),type='l',col='red')
points(density(simil_measure_l[[3]]),type='l',col='blue')
points(density(simil_measure_l[[4]]),type='l',col='green')
points(density(simil_measure_l[[5]]),type='l',col='orange')


plot(density(simil_measure_l2[[1]]),title='Similarity Measure, SCA'
     #, ylim=c(0,7),xlim=c(0,2) ,xlab='Similarity Measure for 1SMR and 2SMR Loadings, SCA',main = ''
)
points(density(simil_measure_l2[[2]]),type='l',col='red')
points(density(simil_measure_l2[[3]]),type='l',col='blue')
points(density(simil_measure_l2[[4]]),type='l',col='green')
points(density(simil_measure_l2[[5]]),type='l',col='orange')
