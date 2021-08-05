Description
================
This tutorial contains **step-by-step guidelines and `R` code** to implement the **Gibbs sampler** presented in the article [Durante, Dunson and Vogelstein (2017) *Nonparametric Bayes Modeling of Populations of Networks*](https://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1219260?journalCode=uasa20), with a focus on the **simulation study in Section 5**. For implementation purposes, please execute the code below considering the same order in which is presented.

Import data
------------------
To start the analysis, **set the working directory** where the data `simulated_data.RData` are located. Once this has been done, **clean the workspace**, and **load the data along with useful `R` packages**.

``` r
###############################################################################
###############################################################################
# CLEAR WORKSPACE AND LOAD USEFUL LIBRARIES ###################################
###############################################################################
###############################################################################

rm(list=ls())
library(plotrix)
library(BayesLogit)
library(MCMCpack)
library(MASS)
library(TeachingDemos)
library(network)
library(latentnet)
library(graphics)
library(latticeExtra)
library(qgraph)
library(coda)
library(gdata)
library(gtools)
library(igraph)
library(lattice)

#################################################################################
#################################################################################
# LOAD DATA (a VxVxn array with n the number of networks and V the nodes) #######
#################################################################################
#################################################################################

# Load the simulated data used in the paper and define model dimensions
load("simulated_data.RData")
V <- dim(A)[1]
n <- dim(A)[3]
```

Set the hyperparameters and run the Gibbs sampler
------------------
Execute the code below to **set the hyperparameters** as in Section 5 and **run the Gibbs sampler** described in the article. For reproducibility purposes, a `seed` can be set.

``` r
###############################################################################
###############################################################################
# SET NUMBER OF GIBBS SAMPLES AND TRUNCATION LEVELS ###########################
###############################################################################
###############################################################################

N_sampl <- 5000

#-----------------------------------------------------------------------------#
# Upper bound on latent space dimensions
# Comment: Lambda matrix is by definition always positive definite, however when 
# R is set too high than needed the Lambda matrix may be very close to be 
# positive semidefinite and R-software can find problems in inverting it. 
R <- 10

#-----------------------------------------------------------------------------#
# Upper bound on number of classes
H <- 30

###############################################################################
###############################################################################
# DEFINE HYPERPARAMETERS SETTINGS #############################################
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# Hyperparameters for multiplicative inverse Gamma MIG
a_1 <- 2.5
a_2 <- 3.5
b_1 <- 1
b_2 <- 1

#-----------------------------------------------------------------------------#
# Dirichlet hyperparameters
a_dir <- rep(1/H,H)

#-----------------------------------------------------------------------------#
# Prior variance and mean of Z. inv_Sigma_Z[v,u,]=1/sigma_l^2 (l index pair (v,u))
inv_Sigma_Z <- matrix(0.1,V,V)
mu_Z <- matrix(0,V,V)

################################################################################
################################################################################
# CREATE ALLOCATION MATRICES ###################################################
################################################################################
################################################################################

#-----------------------------------------------------------------------------#
# pi^{(h)} save posterior samples for pi^{(h)}
pi_matr <- array(0,c(V,V,H,N_sampl))

#-----------------------------------------------------------------------------#
# bar_X=X%*%Lambda^(1/2). with bar_X[v,r,h,]=bar_X_{v,r}^{(h)}
bar_X <- array(0,c(V,R,H,N_sampl))

#-----------------------------------------------------------------------------#
# Shared behavior in matrix form
Z <- array(0,c(V,V,N_sampl))

#-----------------------------------------------------------------------------#
# Weights for complexity learning with inv_Lambda[h,r,]=1/lambda_r^{(h)}
inv_Lambda <- array(0,c(H,R,N_sampl))

#-----------------------------------------------------------------------------#
# Gammas elements to construct lambdas via MIG
theta <- array(0,c(H,R,N_sampl))

#-----------------------------------------------------------------------------#
# G (grouping variable) G[i,]: Group of network i at a given iterat. i=1,...,n
Group <- matrix(0,n,N_sampl)

#-----------------------------------------------------------------------------#
# Pr(G_i=h) for h=1,...,H and i=1,...,n
P_Group <- array(0,c(n,H,N_sampl))

#-----------------------------------------------------------------------------#
# nu with nu[h,] marginal probability of class h: Pr(G=h)
nu <- matrix(0,H,N_sampl)

#-----------------------------------------------------------------------------#
# Sufficient statistic for each step Y[,,h]=Y^{(h)}
Y <- array(0,c(V,V,H))

#-----------------------------------------------------------------------------#
# Other useful matrices for implementing the Gibbs
P_Group_Temp <- matrix(0,n,H)
Sel_Group <- rep(0,H)
bar_XX <- array(0,c(V,V,H))
Prob_matr <- array(0,c(V,V,H))

################################################################################
################################################################################
# INITIALIZE QUANTITIES ########################################################
################################################################################
################################################################################

#-----------------------------------------------------------------------------#
inv_Lambda[,,1] <- 5
theta[,1,1] <- 5
theta[,2:R,1] <- 1
#-----------------------------------------------------------------------------#
for (h in 1:H){
for (v in 1:V){
bar_X[v,,h,1] <- mvrnorm(1,rep(0,R),diag(1,R,R))}}
#-----------------------------------------------------------------------------#
Z[,,1] <- 0   
#-----------------------------------------------------------------------------#
nu[,1] <- rep(1/H,H)
#-----------------------------------------------------------------------------#
Group[,1] <- sample(c(1:H),n,replace=TRUE)
#-----------------------------------------------------------------------------#

################################################################################
################################################################################
# GIBBS SAMPLER ################################################################
################################################################################
################################################################################

time <- system.time(
for (t in 2:N_sampl){
	
################################################################################
################################################################################
# COMPUTE SIZE OF EACH GROUP
for (h in 1:H){Sel_Group[h] <- sum(Group[,t-1]==h)}

################################################################################
################################################################################
# CREATE BINOMIAL MATRICES Y
for (h in 1:H){Y[,,h] <- apply(A[,,which(Group[,t-1]==h)],c(1,2),sum)}

################################################################################
################################################################################
# SAMPLE AUGMENTED POLYA-GAMMA DATA (STORED IN MATRIX OMEGA)
Omega <- array(0,c(V,V,H))

for (h in 1:H){
if (Sel_Group[h]>0){
Omega_temp <- rpg.devroye(num=0.5*V*(V-1),h=rep(Sel_Group[h],0.5*V*(V-1)),z=lowerTriangle(Z[,,t-1]+bar_X[,,h,t-1]%*%t(bar_X[,,h,t-1])))	
	
lowerTriangle(Omega[,,h]) <- Omega_temp
Omega[,,h] <- Omega[,,h]+t(Omega[,,h])

} else {Omega[,,h] <- 0}}

################################################################################
################################################################################
# SAMPLE bar_X AND inv_Lambda 

for (h in 1:H){
################################################################################
# Sample from full conditional given the data for non-empty groups
if (Sel_Group[h]>0){

#-----------------------------------------------------------------------------#
# Sample first row of bar_X given the others
Sigma_v <- solve(t(bar_X[-1,,h,t-1])%*%diag(Omega[-1,1,h],V-1,V-1)%*%bar_X[-1,,h,t-1]+diag(inv_Lambda[h,,t-1],R,R))
mu_v <- Sigma_v%*%(t(bar_X[-1,,h,t-1])%*%(Y[-1,1,h]-Sel_Group[h]/2-Omega[-1,1,h]*(Z[,,t-1])[-1,1]))	

bar_X[1,,h,t] <- mvrnorm(1,mu_v,Sigma_v)

#-----------------------------------------------------------------------------#
# Sample other rows of bar_X given the others
for (v in 2:V){
bar_X_sampl <- bar_X[,,h,t-1]
bar_X_sampl[1:(v-1),] <- bar_X[1:(v-1),,h,t]	
	
Sigma_v <- solve(t(bar_X_sampl[-v,])%*%diag(Omega[-v,v,h],V-1,V-1)%*%bar_X_sampl[-v,]+diag(inv_Lambda[h,,t-1],R,R))
mu_v <- Sigma_v%*%(t(bar_X_sampl[-v,])%*%(Y[-v,v,h]-Sel_Group[h]/2-Omega[-v,v,h]*(Z[,,t-1])[-v,v]))	

bar_X[v,,h,t] <- mvrnorm(1,mu_v,Sigma_v)}

#-----------------------------------------------------------------------------#
# Sample inv_Lambda
a <- c(a_1,rep(a_2,R-1))
b <- c(b_1,rep(b_2,R-1))
inv_Lambda_temp <- apply(bar_X[,,h,t]^2,2,sum)

for (rr in 1:R){
inv_L_rr <- exp(cumsum(log(theta[h,,t-1])))*c(rep(0,rr-1),rep(1/theta[h,rr,t-1],R-rr+1))
theta[h,rr,t] <- (rgamma(1,a[rr]+0.5*V*(R-rr+1),1))/(b[rr]+0.5*sum(inv_L_rr*inv_Lambda_temp))}
inv_Lambda[h,,t] <- exp(cumsum(log(theta[h,,t])))	

################################################################################
# Sample from full conditional (but without conditioning on data) for currently empty groups
} else {

bar_X[,,h,t] <- mvrnorm(V,rep(0,R),diag(inv_Lambda[h,,t-1]^(-1),R,R))

a <- c(a_1,rep(a_2,R-1))
b <- c(b_1,rep(b_2,R-1))
inv_Lambda_temp <- apply(bar_X[,,h,t]^2,2,sum)

for (rr in 1:R){
inv_L_rr <- exp(cumsum(log(theta[h,,t-1])))*c(rep(0,rr-1),rep(1/theta[h,rr,t-1],R-rr+1))
theta[h,rr,t] <- (rgamma(1,a[rr]+0.5*V*(R-rr+1),1))/(b[rr]+0.5*sum(inv_L_rr*inv_Lambda_temp))}
inv_Lambda[h,,t] <- exp(cumsum(log(theta[h,,t])))	}}

################################################################################
################################################################################
# SAMPLE Z
sel_g <- c(which(Sel_Group>0))

for (h in 1:H){
bar_XX[,,h] <- bar_X[,,h,t]%*%t(bar_X[,,h,t])}

V_Z <- 1/(lowerTriangle(apply(Omega[,,sel_g],c(1,2),sum))+lowerTriangle(inv_Sigma_Z))
mean_Z <- ((lowerTriangle(apply(Y[,,sel_g]-Omega[,,sel_g]*bar_XX[,,sel_g],c(1,2),sum))-n/2)+lowerTriangle(mu_Z)*lowerTriangle(inv_Sigma_Z))*V_Z

lowerTriangle(Z[,,t]) <- rnorm(0.5*V*(V-1),mean_Z,sqrt(V_Z))
Z[,,t] <- Z[,,t]+t(Z[,,t])

################################################################################
################################################################################
# UPDATE CLASS PROBABILITY MATRICES
for (h in 1:H){
Prob_matr[,,h] <- 1/(1+exp(-Z[,,t]-bar_X[,,h,t]%*%t(bar_X[,,h,t])))}

#Just for Computational porpuses
Prob_matr[which(Prob_matr==1,arr.ind=TRUE)] <- 0.9999999999
Prob_matr[which(Prob_matr==0,arr.ind=TRUE)] <- 0.0000000001

################################################################################
################################################################################
# UPDATE CONDITIONAL CLASS PROBABILITIES
for (h in 1:H){
P_Group_Temp[,h] <- apply(log(lowerTriangle(Prob_matr[,,h])*apply(A,3,lowerTriangle)+lowerTriangle(1-Prob_matr[,,h])*apply(1-A,3,lowerTriangle)),2,sum)+log(nu[h,t-1])}

for (h in 1:H){
P_Group[,h,t] <- exp(-log(1+apply(exp(P_Group_Temp[,-h]-P_Group_Temp[,h]),1,sum)))}

################################################################################
# UPDATE CLASS INDICATOR FOR EACH NETWORK
for (k in 1:n){
Group[k,t] <- sample(c(1:H),1,prob=c(P_Group[k,,t]))}

################################################################################
################################################################################
# SAMPLE CLASS PROBABILITIES
Group_freq <- rep(0,H)
for (h in 1:(H)){
Group_freq[h] <- sum(Group[,t]==h)}

nu[,t] <- rdirichlet(1,a_dir+Group_freq)

################################################################################
################################################################################
# UPDATE PI MATR
for (h in 1:H){
	pi_matr[,,h,t] <- 1/(1+exp(-Z[,,t]-bar_X[,,h,t]%*%t(bar_X[,,h,t])))
}	
	
#-----------------------------------------------------------------------------#
# Print iterations and Group compositions
print(c(t,table(Group[,t])))})

################################################################################
################################################################################
# SAVE SAMPLES OF USEFUL QUANTITIES FOR POSTERIOR INFERENCE ####################
################################################################################
################################################################################

save(pi_matr,nu,file="posterior_samples_Gibbs.RData")
```

Some examples of posterior analyses
------------------
Once the posterior samples for the key quantities are available, **posterior inference** can be performed. The code below reproduces the **posterior predictive checks** (related to the results displayed in **Figure 4** in the article). 

``` r
###############################################################################
###############################################################################
# CLEAR WORKSPACE AND LOAD USEFUL LIBRARIES ###################################
###############################################################################
###############################################################################

rm(list=ls())
library(igraph)
library(coda)
library(gdata)
library(plotrix)
library(BayesLogit)
library(MCMCpack)
library(MASS)
library(TeachingDemos)
library(network)
library(latentnet)
library(graphics)
library(latticeExtra)
library(qgraph)
library(reshape)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

#################################################################################
#################################################################################
# LOAD DATA AND POSTERIOR SAMPLES ###############################################
#################################################################################
#################################################################################

load("simulated_data.RData")
load("posterior_samples_Gibbs.RData")

#################################################################################
#################################################################################
# SET DATA, MODEL AND MCMC DIMENSIONS ###########################################
#################################################################################
#################################################################################

#-----------------------------------------------------------------------------#
n <- dim(A)[3]
V< - dim(A)[1]
#-----------------------------------------------------------------------------#
R <- 10
H <- 30
#-----------------------------------------------------------------------------#
N_samp <- 5000
burn <- 1000
#-----------------------------------------------------------------------------#

#################################################################################
#################################################################################
# POSTERIOR PREDICTIVE CHECKS (RELATED TO FIGURE 4 IN THE PAPER) ################
#################################################################################
#################################################################################

#-------------------------------------------------------------------------------#
# For a set of 8 network features we study their posterior predictive distribution 
# which is obtained via simulation methods as discussed in Sect. 4.1 of the paper

features_pre <- matrix(,8,N_samp)

pi_h <- array(,c(V*(V-1)/2,H,N_samp))
for (r in 1:N_samp){
	for (h in 1:H){
pi_h[,h,r] <- lowerTriangle(pi_matr[,,h,r])		
}
print(r)}

for (r in 1:N_samp){
	
A_temp <- matrix(0,V,V)
G <- sample(c(1:H),1,prob=nu[,r])	
lowerTriangle(A_temp) <- rbinom(V*(V-1)/2,1,pi_h[,G,r])
A_temp <- A_temp+t(A_temp)
net_dat <- graph.adjacency(A_temp, mode=c("undirected"), weighted=NULL, diag=FALSE)
net_triad <- graph.adjacency(A_temp, mode=c("directed"), weighted=NULL, diag=FALSE)

features_pre[1,r] <- graph.density(net_dat)
features_pre[2,r] <- triad.census(net_triad)[16]/choose(V,3)
features_pre[3,r] <- transitivity(net_dat)
features_pre[4,r] <- assortativity(net_dat,c(rep(1,V/2),rep(2,V/2)))
features_pre[5,r] <- mean(evcent(net_dat)$vector)
features_pre[6,r] <- average.path.length(net_dat)
features_pre[7,r] <- mean(degree(net_dat))
features_pre[8,r] <- sd(degree(net_dat))

print(r)}	

#-------------------------------------------------------------------------------#
# Compute the features for each of the observed networks for goodness of fit checks

features_obs <- matrix(,8,n)

for (r in 1:n){
A_temp <- A[,,r]
net_dat <- graph.adjacency(A_temp, mode=c("undirected"), weighted=NULL, diag=FALSE)
net_triad <- graph.adjacency(A_temp, mode=c("directed"), weighted=NULL, diag=FALSE)

features_obs[1,r] <- graph.density(net_dat)
features_obs[2,r] <- triad.census(net_triad)[16]/choose(V,3)
features_obs[3,r] <- transitivity(net_dat)
features_obs[4,r] <- assortativity(net_dat,c(rep(1,V/2),rep(2,V/2)))
features_obs[5,r] <- mean(evcent(net_dat)$vector)
features_obs[6,r] <- average.path.length(net_dat)
features_obs[7,r] <- mean(degree(net_dat))
features_obs[8,r] <- sd(degree(net_dat))

print(r)}	

#-------------------------------------------------------------------------------#
# Compare graphically the observed features with their posterior predictive 
# distribution obtained under our model

lab <- c("density","triangle_freq","transitivity","assortativity","mean_eigencentrality","ave_path_lenght","mean_degree","sd_degree")

par(mfrow=c(2,4))
for (f in 1:8){
plot(density(features_pre[f,(burn+1):N_samp]),main=paste(lab[f]))
points(features_obs[f,],rep(0,n))}
```
