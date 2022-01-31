
require(fields)
require(MASS)
require(lava)
require(mvtnorm)
require(RColorBrewer)
require(classInt)
require(geoR)
require(matrixcalc)
require(MCMCpack)
require(convoSPAT)
require(Rcpp)
require(inline)
require(RcppEigen)
require(spBayes)

# UPDATE PATH #
setwd("C:/Users/mhb001/Downloads/Mixture-M-RA-master/Mixture-M-RA-master")
sourceCpp("materncovcode.cpp")

###########################################################################################################
# THE FOLLOWING PORTION OF THE CODE USES CODE WRITTEN BY MATTHIAS KATZFUSS.                               #
# IT IS REFERENCED IN THE PAPER "A Multi-Resolution Approximation for Massive Spatial Datasets",          #
#  WHICH APPEARD IN JASA IN 2017.  THEIR CODE IS AVAILABLE AT https://github.com/katzfuss-group/MRA_JASA. #
###########################################################################################################

## return numbered index from tree index (i.e., inverse to indices)
num.ind=function(tree.ind,J){
  if(length(tree.ind)==0) num.index=1 else {
    l=seq(length(tree.ind)-1,0,by=-1)
    num.index=sum(J^l)+sum((tree.ind-1)*J^l)+1 #tree.ind[length(tree.ind)]
  }
  return(num.index)
}


## calculate A-B*inv(C)*D' using cholesky of C
woodbury=function(A,B,C.chol,D) {
  A -   if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else 
    t(base::forwardsolve(C.chol,t(B)))%*%base::forwardsolve(C.chol,t(D))
}


## calculate B'*inv(C)*D using cholesky of C
quadform=function(B,C.chol,D) {
  if(length(B)<=1 | length(C.chol)==0 | length(D)<=1)  0 else #matrix(0,nrow=ncol(B),ncol=ncol(D)) else
    t(base::forwardsolve(C.chol,B))%*%forwardsolve(C.chol,D)
}

cholesky=function(A) if(length(A)<2 && (length(A)<=1 || A==0)) 0 else 
  t(tryCatch(chol(A),error=function(e) chol(A+1e-4*mean(diag(A)))*diag(nrow(A))) )


## calculate inv(C.chol)*A, where C.chol is cholesky factor
sol=function(C.chol,A)  if(length(A)<=1 | length(C.chol)<=1) 0 else base::forwardsolve(C.chol,A)

## calculate t(A)*B
tp=function(A,B)  if(length(A)<=1 | length(B)<=1) 0 else t(A)%*%B

## calculate X'*X for all (relevant) combinations
outer.tp=function(X) {
  m=length(X)-1
  outer.list=vector("list",m+1)
  for(l in 0:m) {
    outer.list[[l+1]]=vector("list",m+1)
    for(k in l:m)  outer.list[[l+1]][[k+1]] = tp(X[[l+1]],X[[k+1]])      
  }   
  return(outer.list)
}

# THE FUNCTION partition.2d IS USED TO PARTITION THE DOMAIN AND DEFINE KNOTS.
# AS IT IS WRITTEN, THIS WILL CREATE J=4 PARTITIONS PER LEVEL, WITH A FIXED NUMBER OF KNOTS 
#  DIFINED ON AN EQUIDISTANT GRID.
# NOTE THAT THIS WILL REORDER YOUR DATA.

partition.2d <- function(J=4,M,domain=c(0,0,1,1),locs,rlat,rlon,z,predlocs=FALSE,zp=NULL){
  pred=(!is.null(predlocs)) # do prediction if locs are given (o/w return likelihood)
  
  splitRectangle <- function(rectbounds,Jlon=sqrt(J),Jlat=sqrt(J)){
    blon <- seq(rectbounds[1],rectbounds[3],,Jlon+1)
    blat <- seq(rectbounds[2],rectbounds[4],,Jlat+1)
    out <- cbind(expand.grid(blon[1:(length(blon)-1)],blat[1:length(blat)-1]),expand.grid(blon[2:length(blon)],blat[2:length(blat)]))
    return(out)
  }
  
  split2d <- function(rectbounds,J=4){
    return(splitRectangle(rectbounds,sqrt(J),sqrt(J)))
  }
  
  par.ind=function(ind,J){
    full.j=indices[[ind]]
    par.j=full.j[-length(full.j)]
    num.ind(par.j,J)
  }
  
  # list of indices for all nodes in the tree structure
  indices=list(list())
  ind.list=list()
  if(M>0) {
    for(m in 1:M){
      ind.list=c(ind.list,list(1:J))
      indices=c(indices,as.list(data.frame(t(expand.grid(ind.list)[,seq(m,1,by=-1)]))))
    }
  }
  n.ind=length(indices)
  
  indres <- vector("list",M+1)
  indres[1] <- 1
  for(m in 1:M){
    end <- length(indres[[m]])
    indres[[m+1]] =seq(indres[[m]][end]+1,indres[[m]][end]+J^m) 
  }
  #bounds are a 4xn.ind matrix.
  bounds <- matrix(ncol=4,nrow=n.ind)
  bounds[1,] <- domain
  counter <- 1
  for(m in 1:M){
    for(ind in indres[[m]]){
      pb <- bounds[ind,]
      bounds[counter+(1:J),] <- as.matrix(split2d(pb,J))
      counter = counter + J
    }
  }
  
  # create knot tree structure
  knots=vector("list",n.ind)
  data=vector("list",n.ind)
  
  for(ind in 1:n.ind) {
    l <- bounds[ind,1]
    b <- bounds[ind,2]
    r <- bounds[ind,3]
    t <- bounds[ind,4]
    
    if(length(indices[[ind]])<M) {
      lon.dist=(r-l)/rlon
      lat.dist=(t-b)/rlat
      lon.knots = seq(l+lon.dist/2,r-lon.dist/2,,rlon)
      lat.knots = seq(b+lat.dist/2,t-lat.dist/2,,rlat)
      
      knots[[ind]] <- as.matrix(expand.grid(lon.knots,lat.knots))
      data[[ind]]=NULL
    } 
    else {
      ind.sub=which(locs[,1]>bounds[ind,1] & locs[,1]<=bounds[ind,3] & locs[,2]>bounds[ind,2] & locs[,2]<=bounds[ind,4])
      knots[[ind]]=locs[ind.sub,]
      data[[ind]]=z[ind.sub]
    }
  }
  for(ind in 1:n.ind){
    if(length(indices[[ind]]==M)){
      if(length(data[[ind]])!=0){
        if(!all(!is.na(data[[ind]]))){
          knots[[ind]]=knots[[ind]][-which(is.na(data[[ind]])),]
          data[[ind]]=data[[ind]][-which(is.na(data[[ind]]))]
          print(ind)
        }
      }
    }
  }
  # prediction locations
  if(typeof(predlocs)=='logical') pred.locs=NA  else {
    pred.locs=vector("list",n.ind)
    pred.data=vector("list",n.ind)
    for(ind in 1:n.ind) {
      if(length(indices[[ind]])>=M) {
        ind.sub=which(predlocs[,1]>=bounds[ind,1] & predlocs[,1]<bounds[ind,3] & predlocs[,2]>=bounds[ind,2] & predlocs[,2]<bounds[ind,4])
        pred.locs[[ind]]=predlocs[ind.sub,]
        pred.data[[ind]]=zp[ind.sub]
      } else pred.locs[[ind]]=knots[[ind]]
    }
  }
  if(typeof(predlocs)!='logical') return(list(indices=indices,knots=knots,data=data,bounds=bounds,indres=indres,pred.locs=pred.locs,pred.data=pred.data))
  else return(list(indices=indices,knots=knots,data=data,bounds=bounds,indres=indres))
  
}

# THE FUNCTION get_basis_cpp IS USED TO COMPUTE THE BASIS FUNCTIONS AND THE PRIOR VARIANCE OF THE BASIS FUNCTION WEIGHTS.

# IF YOUR LOCATIONS ARE GIVEN IN  LATITUDE AND LONGITUDE, WE STRONGLY RECOMMEND REPLACING INSTANCES OF THE rdist() FUNCTION
#   WITH THE rdist.earth() FUNCTION.

# CALLING THE get_basis_cpp FUNCTION SHOULD BE PRECEEDED BY THE USE OF THE FUNCTION partition.2d.
# THE INPUTS data, knots, indices, and pred.locs, MUST ALL BE IN THE LIST FORMAT THAT RESULTS FROM CALLING partition.2d.

# PAY ATTENTION TO THE INPUT cov.fun.
# WE HAVE PROVIDED A FUNCTION THAT CALLS COMPILED cpp CODE TO COMPUTE THE MATERN COVARIANCE FUNCTION.
# IF YOU PREFER A DIFFERENT COVARIANCE FUNCTION, MAKE SURE TO RECONCILE THE INPUTS WITH WHATEVER FUNCTION YOU CHOOSE.

get_basis_cpp <- function(theta,cov.fun,data,knots,indices,pred.locs=NULL){
  pred=(!is.null(pred.locs)) # do prediction if locs are given (o/w return likelihood)
  ## create prior quantities
  n.ind=length(indices)
  M=length(indices[[n.ind]])
  J=if(n.ind==1) 1 else indices[[n.ind]][length(indices[[n.ind]])]
  indres=vector("list",M+1)
  for(m in 0:M) indres[[m+1]]=if(m==0) 1 else (indres[[m]][length(indres[[m]])]+1):(indres[[m]][length(indres[[m]])]+J^m)
  
  V.prior=vector("list",n.ind)
  V.prior2=vector("list",n.ind)
  B=vector("list",n.ind)
  if(pred) {V.p=vector("list",n.ind); B.p=vector("list",n.ind); Bp.tilde=vector("list",n.ind);V.op=vector("list",n.ind); L=vector("list",n.ind)}
  R.prior=vector("list",n.ind)
  R.prior2=vector("list",n.ind)
  R.prior.chol=vector("list",n.ind)
  for(ind in 1:n.ind) {
    inds=indices[[ind]] # full (j) index
    m=length(inds)
    V.prior[[ind]]=vector("list",m+1)
    B[[ind]]=vector("list",m+1)
    if(pred) {V.p[[ind]]=vector("list",m+1); B.p[[ind]]=vector("list",m+1);
    Bp.tilde[[ind]]=vector("list",m+1);
    V.op[[ind]]=vector("list",m+1); L[[ind]]=vector("list",m+1)}
    for(l in 0:m){  
      V.prior[[ind]][[l+1]]=vector("list",m+1)
      if(pred) {V.p[[ind]][[l+1]]=vector("list",m+1); V.op[[ind]][[l+1]]=vector("list",m+1)}
      ind.lm1=if(l<2) 1 else num.ind(inds[1:(l-1)],J)
      for(k in l:m){
        ind.k=if(k==0) 1 else num.ind(inds[1:k],J)
        if(length(knots[[ind]])==2){
          knots[[ind]] <- matrix(knots[[ind]],1,2)
        }
        if(length(knots[[ind.k]])==2){
          knots[[ind.k]] <- matrix(knots[[ind.k]],1,2)
        }
        
        V.prior[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(rdist(knots[[ind]],knots[[ind.k]]),theta[1],theta[2],theta[3]) else
          woodbury(V.prior[[ind]][[l]][[k+1]],V.prior[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                   V.prior[[ind.k]][[l]][[l]]) 
        
        if(pred){
          if(length(pred.locs[[ind]])==2){
            pred.locs[[ind]] <- matrix(pred.locs[[ind]],1,2)
          }
          if(length(pred.locs[[ind.k]])==2){
            pred.locs[[ind.k]] <- matrix(pred.locs[[ind.k]],1,2)
          }
          if(!is.matrix(pred.locs[[ind]])){
            pred.locs[[ind]] <- as.matrix(pred.locs[[ind]])
          }
          if(!is.matrix(pred.locs[[ind.k]])){
            pred.locs[[ind.k]] <- as.matrix(pred.locs[[ind.k]])
          }
          V.p[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(rdist(pred.locs[[ind]],pred.locs[[ind.k]]),theta[1],theta[2],theta[3]) else
            woodbury(V.p[[ind]][[l]][[k+1]],V.p[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                     V.p[[ind.k]][[l]][[l]])          
          V.op[[ind]][[l+1]][[k+1]]= if(l==0) cov.fun(rdist(knots[[ind]],pred.locs[[ind.k]]),theta[1],theta[2],theta[3]) else
            woodbury(V.op[[ind]][[l]][[k+1]],V.prior[[ind]][[l]][[l]],R.prior.chol[[ind.lm1]],
                     V.p[[ind.k]][[l]][[l]])
        }
      }
      if(m==M) {
        B[[ind]][[l+1]]=V.prior[[ind]][[l+1]][[l+1]]
        if(pred) {
          B.p[[ind]][[l+1]]=V.p[[ind]][[l+1]][[l+1]]
          L[[ind]][[l+1]]=V.op[[ind]][[l+1]][[l+1]]
        }
      }
    }
    R.prior[[ind]]=V.prior[[ind]][[m+1]][[m+1]]
    R.prior.chol[[ind]]=cholesky(R.prior[[ind]])
  }
  
  
  basis.functions=basis.pred=vector("list",M)
  K=vector("list",M)
  finest.subs=indices[indres[[M+1]]]
  index <- c()
  for(i in indres[[M+1]]){
    if(R.prior[indres[[M+1]]][[i-indres[[M+1]][1]+1]][1]!=0){index <-c(index,indres[[M+1]][i-indres[[M+1]][1]+1])}
  }
  temp=do.call(blockdiag,R.prior[indres[[M+1]]])
  
  # TO USE OBSERVATION LOCATIONS AS KNOTS AT THE LAST LEVEL,
  # REPLACE "for(res in 0:(M-1))" WITH for(res in 0:(M))
  if(M>0) {  for(res in 0:(M-1)) {
    K[[res+1]]=do.call(blockdiag,lapply(R.prior[indres[[res+1]]],solve))
    bf=vector("list",length(indres[[res+1]]))
    if(pred){bp=vector("list",length(indres[[res+1]])) }
    for(sub in 1:length(indres[[res+1]])) {
      tree.ind=indices[[indres[[res+1]][sub]]]
      ind.children=sapply(finest.subs,function(x) all(x[1:res]==tree.ind))
      num.ind.children=sapply(finest.subs[ind.children],num.ind,J=J)
      bf[[sub]]=do.call(rbind,sapply(B[num.ind.children],`[`,res+1))
      if(pred){bp[[sub]]=do.call(rbind,sapply(B.p[num.ind.children],`[`,res+1))}
    }
    basis.functions[[res+1]]=do.call(blockdiag,bf)
    if(pred){basis.pred[[res+1]]=do.call(blockdiag,bp)}
  } }
  
  return(list(basis.functions,K,basis.pred))
  
}

######################################################################
# THIS CONCLUDES THE PORTION OF THE CODE THAT IS WRITTEN BY KATZFUSS #
######################################################################

##################
# MCMC FUNCTIONS #
##################


# IF YOUR MODEL INCLUDES A LARGE SCALE MEAN, MAKE SURE THAT THE INPUT y IS DEFINED 
#  AS THE DIFFERENCE BETWEEN THE DATA AND THE LARGE-SCALE MEAN
sample_theta_cpp <- function(theta_old,B_old,K_inv_old,eta,y,tau_sq,var_sigma_sq,var_phi,var_nu,data,knots,indices,pred_locs,earth=FALSE,cov.fun){
  
  
  # GENERATE PROPOSAL VALUES #
  
  # MARGINAL VARIANCE # 
  bounds_sigma_sq <- c(max(0,theta_old[1]-var_sigma_sq),theta_old[1]+var_sigma_sq)
  length_sigma_sq <- bounds_sigma_sq[2]-bounds_sigma_sq[1]
  J_den_sigma_sq <- 1/length_sigma_sq
  prop_sigma_sq <- runif(1,bounds_sigma_sq[1],bounds_sigma_sq[2])
  bounds_prop_sigma_sq <- c(max(0,prop_sigma_sq-var_sigma_sq),prop_sigma_sq+var_sigma_sq)
  length_prop_sigma_sq <- bounds_prop_sigma_sq[2]-bounds_prop_sigma_sq[1]
  J_num_sigma_sq <- 1/length_prop_sigma_sq
  
  # RANGE # 
  bounds_phi <- c(max(0,theta_old[2]-var_phi),theta_old[2]+var_phi)
  length_phi <- bounds_phi[2]-bounds_phi[1]
  J_den_phi <- 1/length_phi
  prop_phi <- runif(1,bounds_phi[1],bounds_phi[2])
  bounds_prop_phi <- c(max(0,prop_phi-var_phi),prop_phi+var_phi)
  length_prop_phi <- bounds_prop_phi[2]-bounds_prop_phi[1]
  J_num_phi <- 1/length_prop_phi
  
  # SMOOTHNESS # 
  bounds_nu <- c(max(0,theta_old[3]-var_nu),min(2,theta_old[3]+var_nu))
  length_nu <- bounds_nu[2]-bounds_nu[1]
  J_den_nu <- 1/length_nu
  prop_nu <- runif(1,bounds_nu[1],bounds_nu[2])
  bounds_prop_nu <- c(max(0,prop_nu-var_nu),min(2,prop_nu+var_nu))
  length_prop_nu <- bounds_prop_nu[2]-bounds_prop_nu[1]
  J_num_nu <- 1/length_prop_nu
  
  # CREATE PRIOR QUANTS FOR PROPOSAL THETA #
  prop_theta <- c(prop_sigma_sq,prop_phi,prop_nu)
  prop_quants <- get_basis_cpp(prop_theta,cov.fun=cov.fun,data,knots,indices,pred.locs=pred_locs)
  prop_B_list <- prop_quants[[1]]
  prop_B <- prop_B_list[[1]]
  prop_K_inv <- solve(blockdiag(prop_quants[[2]]))
  for(i in 2:length(prop_B_list)){prop_B <- cbind(prop_B,prop_B_list[[i]])}
  
  # COMPUTE ACCEPTANCE RATIO #
  log_num <- -1/(2*tau_sq) * sum((y-prop_B%*%eta)*(y-prop_B%*%eta)) + 
    1/2 * determinant(prop_K_inv)$modulus - 1/2*eta %*% prop_K_inv %*% eta + 
    log(J_num_sigma_sq) + log(J_num_phi) + log(J_num_nu) - log(prop_sigma_sq) - log(prop_phi) 
  
  log_den <- -1/(2*tau_sq) * sum((y-B_old%*%eta)*(y-B_old%*%eta)) + 
    1/2 * determinant(K_inv_old)$modulus - 1/2* eta %*% K_inv_old %*% eta + 
    log(J_den_sigma_sq) + log(J_den_phi) + log(J_den_nu) - log(theta_old[1]) - log(theta_old[2])
  
  log_r <- log_num - log_den
  log_u <- log(runif(1))
  
  
  if(log_r > log_u){
    theta <- prop_theta
    K_list <- prop_quants[[2]]
    K <- blockdiag(K_list)
    Bp.list <- prop_quants[[3]]
    Bp <- Bp.list[[1]]
    for(i in 2:length(Bp.list)){
      Bp <- cbind(Bp,Bp.list[[i]])
    }
    out <- list(theta,prop_B,K,Bp)
    #print("ACCEPT")
  }
  else{
    #print("REJECT")
    out <- list(theta_old)
  }
  return(out)
}

# SAMPLE BASIS FUNCTION WEIGHTS
sample.eta <- function(R,B,BTB,tau2,K.inv){
  V.eta <- solve(BTB/tau2 + as.matrix(K.inv) + 0.01*diag(nrow(K.inv)))
  M.eta <- (V.eta %*% t(B) %*% R)/tau2
  out <- mvrnorm(1,M.eta,V.eta)
  return(out)
}


sample.rho <- function(l.vector,l,rho.old,M,indres,var,a,b){
  bounds <- c(max(0,rho.old-var),min(1,rho.old+var))
  length <- bounds[2]-bounds[1]
  J.den <- 1/length
  prop.rho <- runif(1,bounds[1],bounds[2])
  bounds.prop <- c(max(0,prop.rho-var),min(1,prop.rho+var))
  length.prop <- bounds.prop[2]-bounds.prop[1]
  
  J.num <- 1/length.prop  
  log.num <- log.den <- 0
  start <- 1
  for(m in 1:(M-1)){
    ind <- indres[[m+1]]
    sum.z <- sum(l.vector[ind]!=l)
    n.zero <- if(m==0) 0 else sum(l.vector[indres[[m]]]==l)
    log.num <- log.num + sum.z * log(prop.rho^m) + (4^m - 4*n.zero - sum.z)*log(1-prop.rho^m)
    log.den <- log.den + sum.z * log(rho.old^m) + (4^m - 4*n.zero - sum.z)*log(1-rho.old^m)
  }  
  log.num <- log.num + J.num + (a-1)*log(prop.rho) + (b-1)*log(1-prop.rho)
  log.den <- log.den + J.den + (a-1)*log(rho.old) + (b-1)*log(1-rho.old)
  log.r <- log.num-log.den
  log.u <- log(runif(1))
  if(is.na(log.r)){print("PROBLEM!")}
  if(log.r > log.u){rho <- prop.rho}
  else{rho <- rho.old}
  return(rho)
}

# FOR ALL INTENTS AND PURPOSES, THIS FUNCTION SAMPLES THE LATENT BINARY VARIALBE Z_{m,j}.
# FOR CONVENIENCE, WE SCALE THE OUTPUT IN TERMS OF L, THE CONSTANT WE USE TO SCALE THE VARIANCE OF THE SECOND MIXTURE COMPONENT.
sample.l <- function(eta,rho,l,K,m){
  log.num <- log(rho^m*dmvnorm(eta,rep(0,length(eta)),K))
  log.den <- log(rho^m*dmvnorm(eta,rep(0,length(eta)),K) + (1-rho^m)*dmvnorm(eta,rep(0,length(eta)),K/l))
  log.p <- log.num-log.den
  p <- exp(log.p)
  out <- if (is.na(p)) 0 else rbinom(1,1,1-p)*l
  if(is.na(p)){print("ASDF")}
  return(out)
}

# LARGE SCALE MEAN (IN THE CASE WHERE X IS A VECTOR OF 1'S) OR LINEAR COEFFICIENTS (WHEN X IS A MATRIX OF COVARIATES)
sample_beta <- function(R,X,XTX,tau2){
  prec_beta <- XTX/tau2
  V_beta <- solve(prec_beta + diag(0.01,nrow(prec_beta)))
  M_beta <- V_beta %*% t(X) %*% R / tau2
  return(mvrnorm(1,M_beta,V_beta))
}

# RESIDUAL VARIANCE
sample.tau2 <- function(R,prior.shape,prior.rate,N){
  post.shape <- prior.shape+(N/2)
  post.rate <- 0.5*sum((R)^2)+prior.rate
  new.tau2 <- 1/rgamma(1,post.shape,post.rate)
  return(new.tau2)
}

adjust_acceptance <- function(acc_rate,var,target_rate){
  y <- 1 + 10*(acc_rate-target_rate)^3
  if (y < .8){y <- 0.8}
  if (y > 1.2){y <- 1.2}
  out = var*y
  return(out)
}

