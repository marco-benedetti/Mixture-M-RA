require(fields)
require(MASS)
require(lava)
require(mvtnorm)
require(glmnet)
require(RColorBrewer)
require(classInt)
require(geoR)
require(matrixcalc)
require(MCMCpack)
require(Rcpp)
require(inline)
require(RcppEigen)




#set.seed(42141)

# GENERATE LOCATIONS TO COMPUTE INITIAL QUANTITIES
n <- 500

lat <- runif(n)
lon <- runif(n)
locs <- cbind(lon,lat)
locs <- rbind(locs,as.matrix(expand.grid(seq(1e-4,1-1e-4,,16),seq(1e-4,1-1e-4,,16))))
predlocs <- as.matrix(expand.grid(seq(3e-2,1-3e-2,,16),seq(3e-2,1-3e-2,,16)))
locs.full <- rbind(locs,predlocs)
n.full <- nrow(locs.full)
lat <- locs[,2]
lon <- locs[,1]
predind <- 757:1012
predlocs <- locs.full[predind,]
locs <- locs.full[-predind,]
n <- nrow(locs)
n.pred <- nrow(predlocs)


#  NUMBER OF SIMULATIONS, BURNIN, AND ITERATIONS PER SIMULATION.
numsim <- 1
burnin <- 10
niter <- 20

# NUMBER OF LEVELS PLUS 1, M, AND PARTITIONS PER SUBREGION, J.
J <- 4
M <- 4

# TOTAL NUMBER OF PARTITIONS
n.groups <- sum(J^(0:(M-1)))

# NUMBER OF KNOTS PER PARTITION (IN THIS EXAMPLE, IT IS 16)
rlat <- rlon <- 4

# WE RECOMMEND USING L = 100 FOR THESE DATA BASED ON EXAMINATION OF MIXING BEHAVIOR DURING PRELIMINARY SIMULATIONS.
# MAKE SURE THAT L IS LARGE ENOUGH SO THAT THE SECOND MIXTURE COMPONENT OF OUR PRIOR DISTRIBUTION CAN REASONABLY TAKE THE PLACE OF A POINT-MASS AT ZERO.
l <- 100

X <- rep(1,n)
XTX <- t(X)%*%X
simdata <- preddata <- eta <- alllocs <- allpredlocs <- list()


eta.mean  <- matrix(0,numsim,n.groups*rlat*rlon)
z.mean  <- matrix(0,numsim,n.groups)
rho.mean <- rep(0,numsim)
predmean <- predsd <-  matrix(0, numsim,n.pred)
LLpred <- ULpred <- LLkrig <- ULkrig <- matrix(0,numsim,n.pred)
ptm <- proc.time()
for(s in 1:numsim){

  print(paste("SIMULATION",s))
  print("Generating Data")
  
  # NEW LOCATIONS EVERY SIMULATION
  n <- 500
  lat <- runif(n)
  lon <- runif(n)
  locs <- cbind(lon,lat)
  
  # MAKING SURE THAT THERE IS ONE LOCATION IN EACH PARTITION AT THE HIGHEST LEVEL MAKES IT EASIER TO COMPUTE BASIS FUNCTIONS.
  # IF YOU ARE NOT RANDOMLY GENERATED DATA (IE LOCATIONS IN THE US), IT HELPS TO GENERATE PSEUDO-LOCATIONS IN EACH PARTITION WITH MISSING DATA.
  locs <- rbind(locs,as.matrix(expand.grid(seq(1e-4,1-1e-4,,16),seq(1e-4,1-1e-4,,16))))
  predlocs <- as.matrix(expand.grid(seq(3e-2,1-3e-2,,16),seq(3e-2,1-3e-2,,16)))
  locs.full <- rbind(locs,predlocs)
  n.full <- nrow(locs.full)
  lat <- locs[,2]
  lon <- locs[,1]
  predind <- 757:1012
  predlocs <- locs.full[predind,]
  locs <- locs.full[-predind,]
  n <- nrow(locs)
  n.pred <- nrow(predlocs)
  
  # GENERATE DATA ACCORDING TO A MIXTURE OF MULTIVARIATE NORMAL DISTRIBUTIONS.
  
  w1.s <- mvrnorm(1,rep(0,n.full),maternCovcpp(rdist(locs.full,locs.full),1,1,1))
  w2.s <- mvrnorm(1,rep(0,n.full),maternCovcpp(rdist(locs.full,locs.full),1,0.01,1))
  
  in.region.1 <- as.numeric(locs.full[,1]<0.5)
  w.s <- in.region.1*w1.s +(1-in.region.1)*w2.s
  y.full <- in.region.1*w1.s +(1-in.region.1)*w2.s + rnorm(n.full,0,sqrt(0.05))
  
  predind <- 757:1012
  y.pred <- y.full[predind]
  y <- y.full[-predind]
  
  
  quants <- partition.2d(J=4,M=4,domain=c(0,0,1,1),locs,rlat,rlon,y,predlocs=predlocs,y.pred)
  data <- quants$data
  knots <- quants$knots
  indices <- quants$indices
  indres <- quants$indres
  pred.locs <- quants$pred.locs
  pred.data <- quants$pred.data
  
  # RUNNING partition.2d REORDERS THE DATA BASED ON THEIR LOCATIONS WITHIN PARTITIONS AT THE HIGHEST LEVEL.
  # THEREFORE, WE REORDER THE THE DATA.  HEREAFTER, IF YOU USE THE DATA IN THE ORIGINAL ORDER, THE CODE WILL NOT WORK.
  locs.new <- NULL
  predlocs.new <- NULL
  y.new <- NULL
  yp.new <- NULL
  for(i in indres[[length(indres)]]){
    locs.new <- rbind(locs.new,knots[[i]])
    if(length(pred.locs[[i]]>0)){
      predlocs.new <- rbind(predlocs.new,pred.locs[[i]])
    }
    y.new <- c(y.new,data[[i]])
    yp.new <- c(yp.new,pred.data[[i]])
  }
  
  simdata[[s]] <- y.new
  preddata[[s]] <- yp.new
  alllocs[[s]] <- locs.new
  allpredlocs[[s]] <- predlocs.new

	
  # INITIAL VALUES FOR THETA AND L
  theta <- c(1,0.05,1)
  l.vector <- rep(1,n.groups)
  
  # COMPUTE PRIOR QUANTITIES
  prior.quants <- get_basis_cpp(theta,maternCovcpp,data,knots,indices,pred.locs=predlocs)
  B.list <- prior.quants[[1]]
  K.list <- prior.quants[[2]]
  Bp.list <- prior.quants[[3]]
  K <- blockdiag(K.list)
  
  # STORE PRIOR QUANTITIES IN A WAY THAT IS CONVENIENT FOR COMPUTATION
  B <- B.list[[1]]
  Bp <- Bp.list[[1]]
  for(i in 2:length(B.list)){
    B <- cbind(B,B.list[[i]])
  }
  for(i in 2:length(Bp.list)){
    Bp <- cbind(Bp,Bp.list[[i]])
  }
  
  Blist2 <- list()
  Klist2 <- list()
  Kinvlist <- list()
  BTBlist <- list()
  
  start <- 1
  end <- 0
  for(i in 1:n.groups){
    end <- end+length(knots[[i]][,1])
    Blist2[[i]] <- B[,start:end]
    BTBlist[[i]] <- t(Blist2[[i]])%*%Blist2[[i]]
    Klist2[[i]] <- K[start:end,start:end]
    Kinvlist[[i]] <- solve(Klist2[[i]])
    start <- end+1
  }
  
  npred <- nrow(Bp)
  BTB <- t(B)%*%B
  K.inv <- solve(K)
  
  M <- length(indices[[length(indices)]])
  eta.out <- matrix(0,niter,length(B[1,]))
  par.out <- matrix(0,niter,2)
  w.out <- matrix(0,niter,n)
  theta.out <- matrix(0,niter,3)
  beta.out <- rep(0,niter)
  z.out <- matrix(0,niter,n.groups)
  
  pred.out <- matrix(0,niter-burnin,n.pred)
  init.w <- rep(0,n)
  acc_rho <- den_rho <- 0
  var_rho <- 0.1
  var_phi = 0.01
  var_sigma_sq = 0.1
  var_nu = 0.05
  acc = 0

  print("Starting Gibbs Sampler")
  par(mfrow=c(5,4))
  ptm <- proc.time()
  for(k in 1:niter){
    if(k==1){
      tau2 <- 1
      w <- init.w
      eta <- rep(0,n.groups*rlat*rlon)
      beta <- 0
      rho <- 0.1
    }
    R_tau2 <- y.new-beta-w 
    new.tau2 <- sample.tau2(R_tau2,0.001,0.001,length(y.new))
    
    R_beta <- y.new-w
    beta <- sample_beta(R_beta,X,XTX,new.tau2)
    
    start <- 1
    end <- 0
    for(i in 1:(n.groups)){
      end <- end+length(knots[[i]][,1])
      
      R.temp <- y.new - beta - B %*% as.matrix(eta) + Blist2[[i]] %*% as.matrix(eta[start:end])
      eta[start:end] <- sample.eta(R.temp,Blist2[[i]],BTBlist[[i]],new.tau2,Kinvlist[[i]]*(l.vector[i]+1))
      
      level <- length(indices[[i]])
      prev.ind <- if(level<=1) 1 else num.ind(indices[[i]][1:(level-1)],J)
      l.vector[i] <- if(l.vector[prev.ind]==l) l else sample.l(eta[start:end],rho,l,Klist2[[i]],level)
      
      start <- end+1

    }
    
    z.out[k,] <- 1-l.vector/l
    
    rho.old <- rho
    new.rho <- sample.rho(l.vector,l,rho.old,M,indres,var_rho,1,1)
    if(new.rho!=rho.old){acc_rho <- acc_rho+1}
    rho <- new.rho

    new.theta <- sample_theta_cpp(theta,B,blockdiag(Kinvlist),eta,y.new-beta,new.tau2,var_sigma_sq,var_phi,var_nu,data,knots,indices,pred.locs)
    if(!isTRUE(all.equal(new.theta[[1]],theta))){
      theta <- new.theta[[1]]
      B <- new.theta[[2]]
      K <- new.theta[[3]]
      Bp <- new.theta[[4]]
      start <- 1
      end <- 0
      for(i in 1:n.groups){
        end <- end+length(knots[[i]][,1])
        Blist2[[i]] <- B[,start:end]
        Klist2[[i]] <- K[start:end,start:end]
        BTBlist[[i]] <- t(Blist2[[i]])%*%Blist2[[i]]
        Kinvlist[[i]] <- solve(Klist2[[i]])
        start <- end+1
      }
      acc = acc+1
    }
    
    if((k%%100==0)&(k<burnin)){
      var_phi <- adjust_acceptance(acc/100,var_phi,0.25)
      var_sigma_sq <- adjust_acceptance(acc/100,var_sigma_sq,0.25)
      var_nu <-adjust_acceptance(acc/100,var_nu,0.25)
      var_rho <- adjust_acceptance(acc_rho/100,var_rho,0.30)
      acc = 0
      acc_rho = 0
    }
    if(k%%100==0){
      print(paste("Iteration: ",k))

    }
    new.w <- B %*% as.matrix(eta)  
    w <- new.w
    
    tau2 <- new.tau2

    # PREDICTION
    if(k > burnin){
      pred.out[k-burnin,] <-  beta + Bp %*% as.matrix(eta) + as.matrix(rnorm(npred,0,sqrt(tau2)))
      # if(k<(21+burnin)){
      #   plot(yp.new,col=2,type='l')
      #   lines(pred.out[k-burnin,])
      # }
    }
    
    # STORE OUTPUT
    beta.out[k] <- beta
    par.out[k,] <- c(tau2,rho)
    theta.out[k,] <- theta
    eta.out[k,] <- as.vector(eta)
    w.out[k,] <- w
    
  }
  # TIME PER ITERATION
  print((proc.time()-ptm)/k)
  
  # PLOT
  par(mfrow=c(1,1))
  plot(yp.new,type='l')
  lines(apply(pred.out,2,mean),col=2)
  rho.mean[s] <- mean(par.out[burnin:niter,2])
  z.mean[s,] <- apply(z.out[burnin:niter,],2,mean)
  eta.mean[s,] <- apply(eta.out[burnin:niter,],2,mean) 
  predmean[s,] <- apply(pred.out,2,mean)
  predsd[s,] <- apply(pred.out,2,sd)
  LLpred[s,] <- apply(pred.out,2,function(x) return(quantile(x,0.025)))
  ULpred[s,] <- apply(pred.out,2,function(x) return(quantile(x,0.975)))

  
}

 
plot(y,new,type='l')
lines(predmean[1,],col=2)