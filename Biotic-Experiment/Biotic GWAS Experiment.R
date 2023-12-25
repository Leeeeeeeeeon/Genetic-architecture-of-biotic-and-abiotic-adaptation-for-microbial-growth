load("C:/Users/AwesomeLeon/Desktop/Gong微生物/logistic/dat(fsnp).RData")
setwd("C:/Users/AwesomeLeon/Desktop/微分方程")
SE.load <- function(s.file="./data/coculture-S1.csv", s.snp.file="./data/gen_s.txt") {
  
  # 读取 coculture-S1.csv 文件
  SP <- read.csv(s.file)
  rownames(SP) <- as.character(SP[,1])
  SP1 <- SP[,-1]
  colnames(SP1) <- NULL
  SP1L <- log(SP1)  # 对数据进行对数转换
  
  # 读取 gen_s.txt 文件
  ssnp <- read.table(s.snp.file, header=TRUE)
  rownames(ssnp) <- ssnp[,1]
  ssnp1 <- ssnp[,-1]
  colnames(ssnp1) <- as.character(SP[,1])
  
  # 过滤 SNP 数据
  del.i <- c()
  for(i in 1:dim(ssnp1)[1]) {
    index <- min(table(as.character(unlist(c(ssnp1[i,]))))/dim(SP)[1])
    del.i <- c(del.i, index)
  }
  nssnp1 <- ssnp1[-which(del.i<0.1),]
  
  # 返回列表
  list(sp.p=SP1L, sample_N=dim(SP)[1], s.snp=nssnp1)
}

dat <- SE.load(s.file="./data/coculture-Saureus.csv",
               s.snp.file="./data/gen_saureus.txt")

library(mvtnorm)
library(deSolve)

ctA<-dat$sp.p
Mcta<-colMeans(dat$sp.p)
dat$s.snp
st<-dat$sample_times
times<-dat$sample_times
atimes<-dat$sample_times
library(deSolve)
library(mvtnorm)
###########################拟合参数_独立生长###########################################
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
#colMeans(dat$sp.pi)

init.est1<-function(dat){
  spp<-dat$sp.pi
  ty <- cbind(spp)
  cm <- colMeans(ty)
  
  sum.curve<-function(par6){
    y<-com.get_mu1(par6[c(1,2,3)],times,x0=8.517193)
    A <- sum((cm-y)^2)
    A
  }
  #init.par<-c(0.01463287,30.97206591,-4.36923919)
  #init.par<-c(0.01461047,30.97984556,-4.36965472)
  init.par<-mu.g1
  #init.par<-c(0.41552618,0.04515011,3.21155078,3.93858834)
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=55000))
  return(a)
}
init.est1(dat)

mu.g1<-init.est1(dat)$par

plot(times,colMeans(dat$sp.pi),xlim = c(0,36),ylim = c(8,28))
par(new=T)
plot(times,com.get_mu1(mu.g1[c(1,2,3)],times,x0=8.517193),type="l",col="red",xlim = c(0,36),ylim = c(8,28))

k <- length(mu.g1)
n <- length(times)
y <- colMeans(dat$sp.pi)
SSR0<-init.est1(dat)$value

AIC0 <- 2*k+n*log(SSR0/n)
AIC0
BIC0<-k*log(n)+n*log(SSR0/n)
BIC0
mu.g1

SSR0<-init.est1(dat)$value
SSR0<-sum((colMeans(dat$sp.p)-(com.get_mu1(mu.g1,times,x0=8.517193)-com.get_mu2(mu.g2,times,x0=0.001)))^2)
AIC0 <- 2*k+n*log(SSR0/n)
AIC0
BIC0<-k*log(n)+n*log(SSR0/n)
BIC0


###########################拟合参数_环境作用###########################################

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



init.est2<-function(dat){
  
  cm <- colMeans(dat$sp.pi)-colMeans(dat$sp.pi)
  
  sum.curve<-function(par6){
    y<-com.get_mu2(par6,times,x0=0)
    A <- sum((cm-y)^2)
    A
  }
  init.par<-c(0,0)
  #init.par<-mu.g2
  
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=15000))
  return(a)
}
init.est2(dat)


###########################拟合参数_汇总###########################################

init.est<-function(dat){
  
  cm <- colMeans(dat$sp.pi)
  
  sum.curve<-function(par6){
    y<-com.get_mu1(par6[c(1,2,3)],times,x0=8.517193)-com.get_mu2(par6[c(4,5)],times,x0=0)
    A <- sum((cm-y)^2)
    
  }
  #init.par<-c(0.07930295,30,-3.07472075,0.2,0.2)
  #init.par<-c(0.01461503,40.97827292,-2.36956363,0.5258307,-1.8593241,1.6807830)
  #init.par<-c(mu.g1,mu.g2)
  init.par<-mu.g
  #init.par<-c(-0.02360943,1.60739709,5.60750879,-10.38635991)
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=55000))
  return(a)
}


mu.g<-c(0.01461787,30.97731918,-4.36950143,0,0)


SAD1.get_mat <- function (par0, times, traits = 1, options = list()) {
  
  par <- par0
  if (class(par0) == "list") 
    par <- unlist(par0)
  t_len <- length(times)
  SAD.1 <- array(0, dim = c(t_len * traits, t_len * traits))
  for (i0 in 1:traits) for (i1 in 1:traits) {
    if (i0 == i1) 
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        SAD.1[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                k1] <- abs(par[i0 * 2]) * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
        SAD.1[(i0 - 1) * t_len + k1, (i1 - 1) * t_len + 
                k0] <- abs(par[i0 * 2]) * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
      }
  }
  return(SAD.1)
}

SAD1.get_inv_mat <- function (par, times, traits = 1, options = list()) {
  
  SAD.1 <- SAD1.get_mat(par, times, traits, options)
  return(solve(SAD.1))
}

SAD1.get_mat_det <- function (par, times, traits = 1, options = list()) {
  
  SAD.1 <- SAD1.get_mat(par, times, traits, options)
  return(det(SAD.1))
}

curve.mlefunc<-function( par,y1,times)
{
  len.cov <- 2
  par.covar <- par[1:len.cov]
  n  <- length(y1[,1])
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  sigma <- solve( SAD1.get_inv_mat(par.covar, times, 1) )
  
  curve.par <- par[(len.cov+1):(len.cov+ 5)]
  mu <- com.get_mu1(curve.par[c(1,2,3)],times,x0=8.517193)-com.get_mu2(curve.par[c(4,5)],times,x0=0)
  
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
  sigma <- solve( SAD1.get_inv_mat(par.covar, times, 1) )
  len.gen <- 5
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    curve.par <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    yy1 <- y1[snp.index[[i]],]
    mu <- com.get_mu1(curve.par[c(1,2,3)],times,x0=8.517193)-com.get_mu2(curve.par[c(4,5)],times,x0=0)
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


H0.est<-function(dat){
  
  mpheno<-colMeans(dat$sp.pi)
  par<-c(0.05427067,25.51046568,-3.70801206,0.51771292,-1.86380314,1.70714315)
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
      AA <- curve.mlefunc(nnpar,y=dat$sp.pi,times =times)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "Nelder-Mead",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    
    mle.1 <- function(npar){
      #npar<-parin[1:2]
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mlefunc(nnpar,y=dat$sp.pi,times =times)
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
  LL <- curve.mlefunc(parin2,y=dat$sp.pi,times =times)
  return(c(LL,parin2))
}

H0<-H0.est(dat)
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
      y<-com.get_mu1(par6[c(1,2,3)],times,x0=8.517193)-com.get_mu2(par6[c(4,5)],times,x0=0)
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

com.DH.est1 <- function(dat,interval=c(1,3208)){
  
  y1 <- as.matrix(dat$sp.pi)
  snp <- dat$s.snp
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
    h01 <- H0
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

