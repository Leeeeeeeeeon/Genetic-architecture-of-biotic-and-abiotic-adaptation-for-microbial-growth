jin2 = read.csv("./growdata_mic0.csv", header = TRUE, row.names = 1) # Growth data of Staphylococcus aureus in mic=0
jin3 = read.csv("./growdata_mic2.csv", header = TRUE, row.names = 1) # Growth data of Staphylococcus aureus in mic=2
jin4 = read.csv("./growdata_mic4.csv", header = TRUE, row.names = 1) # Growth data of Staphylococcus aureus in mic=4
jin5 = read.csv("./growdata_mic6.csv", header = TRUE, row.names = 1) # Growth data of Staphylococcus aureus in mic=6
SNP = read.csv("./snpdata.csv", header = TRUE, row.names = 1) # SNP data
k8str = read.csv("./FASTSTRUCTURE.csv", header = FALSE) # FASTSTRUCTURE data

t <- jin2[,1]#times point
times <- jin2[,1]#times point
get_strjin <- function(jin,si){
  
  
  BB <- apply(si,1,function(x){which(x==max(x))})
  
  BBt <- table(BB)
  BBn <- as.numeric(names(BBt))
  rjin2 <- c()
  for (i in 1:nrow(jin)) {
    
    
    pdata <-   jin[i,]
    subpop <- list()
    almea <- c()
    for( j in 1:length(BBt)){
      
      index1 <- which(BB==BBn[j])
      subpop[[j]] <- pdata[index1]
      almea[j] <- mean(as.numeric(pdata[index1]))
      
      
    }
    for (j in 1:length(BBt)) {
      chazhi=almea[j]-mean(almea)
      pdata[which(BB==BBn[j])] <- pdata[which(BB==BBn[j])]-chazhi
      
      
      
      
      
    }
    rjin2 <- rbind(rjin2,pdata)
    
    
  }
  
  
  rjin2 <- cbind(jin2[,1],rjin2)
  
  return(rjin2)}
rjin2 <- get_strjin(jin2[,-1],k8str)
rjin3 <- get_strjin(jin3[,-1],k8str)
rjin4 <- get_strjin(jin4[,-1],k8str)
rjin5 <- get_strjin(jin5[,-1],k8str)
ct0<-t(rjin2[,-1])
ct2<-t(rjin3[,-1])
ct4<-t(rjin4[,-1])
ct6<-t(rjin5[,-1])
Mct0<-colMeans(ct0)
Mct2<-colMeans(ct2)
Mct4<-colMeans(ct4)
Mct6<-colMeans(ct6)
library(deSolve)
library(mvtnorm)
###########################拟合初始参数########################################
com.get_mu1 <- function(par, times, x0, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2],
      a1 = par[3]
    );
  }
  
  state0 <- c(X=x0);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-(X/k1)^a1)/a1
            list(c(dX))
          }
          
    ) # end with(as.list ...)
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


com.get_mu2 <- function(par, times, x0, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      b1 = par[1],
      c1 = par[2]
      
      
    );
  }
  
  state0 <- c(X=x0);
  y <- COMP.f2( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2] ) );
}


COMP.f2 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- (b1*X^2+c1*X)/(1+X)
            list(c(dX))
          }
          
    ) # end with(as.list ...)
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

init.est<-function(ct0){
  spp<-ct0
  ty <- cbind(spp)
  cm <- Mct0
  
  sum.curve<-function(par6){
    y<-com.get_mu1(par6[c(1,2,3)],times,x0=0.09428755)-com.get_mu2(par6[c(4,5)],times,x0=0)
    A <- sum((cm-y)^2)
    A
  }
  #init.par<-c(0.2400414,1.5424620,0.5542492,0.01544294)
  init.par<-mu.g
  #init.par<-c(-0.02360943,1.60739709,5.60750879,-10.38635991)
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=55000))
  return(a)
}

AR1.get_mat <- function(par0, times, traits=1, options=list())
{
  par<-par0;
  if (class(par0)=="list")
    par <- unlist(par0);
  
  t_len <- length( times );
  
  Ar.1 <- array(0, dim=c(t_len*traits,t_len*traits));
  for (i0 in 1:traits)
    for (i1 in 1:traits)
    {
      if (i0==i1)
        for (k0 in 1:t_len)
          for (k1 in 1:t_len)
          {
            Ar.1[(i0-1)*t_len+k0,(i1-1)*t_len+k1] <- par[i0*2]^2 * par[i0*2-1]^abs( k0 - k1 );
          }
    }
  
  return(Ar.1);
}

AR1.get_inv_mat <- function(par, times, traits=1, options=list())
{
  Ar.1 <- AR1.get_mat(par, times, traits, options)
  return(solve(Ar.1));
}

AR1.get_mat_det <- function(par, times, traits=1, options=list())
{
  Ar.1 <- AR1.get_mat(par, times, traits, options)
  return(det(Ar.1));
}
curve.mlefunc<-function( par,y1,times)
{
  len.cov <- 2
  par.covar <- par[1:len.cov]
  n  <- length(y1[,1])
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  AR1 <- AR1.get_inv_mat(par.covar,times)
  sigma <- solve( AR1 )
  
  curve.par <- par[(len.cov+1):(len.cov+ 5)]
  mu <- com.get_mu1(curve.par[c(1,2,3)],times,x0=0.09428755)-com.get_mu2(curve.par[c(4,5)],times,x0=0)
  
  yy <- y1
  for ( i in 1:dim(y1)[2] )
  {
    y1.miss <- which( is.na(y1[,i]) )
    yy[y1.miss, i]<- mu[i]
  }
  fy <- dmvnorm(yy,mu,sigma)
  #fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps
  A <- -sum(log(fy))
  #cat("LL=",A,"\n")
  return(A)
}

S.mlefunc <- function(par,y1,times,snp.index,snp.type)
{
  n  <- length(y1[,1])
  
  len.cov <- 2
  par.covar <- par[1:len.cov]
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  ar1 <- AR1.get_inv_mat(par.covar,times)
  sigma <- solve( ar1)
  len.gen <- 5
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    curve.par <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    yy1 <- y1[snp.index[[i]],]
    mu <- com.get_mu1(curve.par[c(1,2,3)],times,x0=0.09428755)-com.get_mu2(curve.par[c(4,5)],times,x0=0)
    nyy1 <- yy1
    fy1 <- dmvnorm( nyy1, mu, sigma)
    #fy1[which(fy1<=.Machine$double.eps)] <- .Machine$double.eps
    A1 <- c(A1,-sum(log(fy1)))
    len <- len + len.gen
  }
  A <- sum(A1)
  #cat("LL=",A,"\n")
  return (A);
}



H0.est<-function(ct0){
  
  spp<-ct0
  allpheno <- cbind(spp)
  ms <- colMeans(spp)
  mpheno<-as.numeric(ms)
  par<-c(0.289026,1.536441,1.227860,0,0)
  covar.par<-c(0.6,1.74)
  parin<-c(covar.par,par)
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  while(loop_k<max_iter && max_err>epsi){
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[3:7])
      AA <- curve.mlefunc(nnpar,y=ct0,times =times)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "Nelder-Mead",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    
    mle.1 <- function(npar){
      #npar<-parin[1:2]
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mlefunc(nnpar,y=ct0,times =times)
      AA
    }
    r1 <- optim(c(parin[3:7]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    
    loop_k <- loop_k+1; 
  }
  LL <- curve.mlefunc(parin2,y=ct0,times =times)
  return(c(LL,parin2))
}


H0<-H0.est(ct0)
par1<- H0[-1]


H1.est <- function(y11,SNP1,init.par=par1,times){
  index <- table(SNP1)
  snp.type <- as.character(names(index))
  phenos<-as.matrix(y11)
  par.covar <- par1[1:2]
  g.par <- c()
  snp.index <- list()
  
  for(j in 1:length(snp.type)){
    index <- which(SNP1==snp.type[j])
    yy <- as.numeric(colMeans(phenos[index,]))
    sum.curve<-function(par6){
      y<-com.get_mu1(par6[c(1,2,3)],times,x0=0.09428755)-com.get_mu2(par6[c(4,5)],times,x0=0)
      A <- sum((yy-y)^2)
      A
    }
    r1<-optim(par1[3:7],sum.curve,method = "BFGS",control = list(maxit=32000))
    snp.index[[j]] <- index
    g.par <- c(g.par,r1$par)
  }
  parin <- c(par.covar,g.par)
  n.par <- length(parin)
  
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    
    
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[3:n.par])
      AA <- S.mlefunc(nnpar,y=phenos,times =times,snp.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "Nelder-Mead",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- S.mlefunc(nnpar,y=phenos,times =times,snp.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[3:n.par]),mle1.g,method = "BFGS",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}

com.DH.est1 <- function(ct0,interval=c(1,25173)){
  
  y1 <- as.matrix(ct0)
  snp <- SNP[-1]
  nm <- dim(snp)[1]
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nm)
    n2 <- nm
  p.res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)
  for(i in n1:n2){
    #i=1
    S.SNP <- as.character(snp[i,])
    par1<- H0[-1]
    h01 <- H0#没有NA不用重复计算
    h02 <- try(H1.est(y1,S.SNP,init.par=par1,times),TRUE)
    if (class(h02) == "try-error") 
      h02 <- NA
    LR <- 2*(h01[1]-h02[1])
    
    p.tmp <- c(LR,H0[1],h02)
    cat("snp", i, "=", p.tmp, "\n");
    p.res[(i-(n1-1)),1:length(p.tmp)] <- p.tmp
  }
  
  return(p.res)
}
