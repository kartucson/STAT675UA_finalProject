library(boot)
library(CircStats)
library(Hmisc)

data_set_all <- read.csv("train_HMM_Sound.csv") 
#data_set <- data_set_all
data_set <- data_set_all[1:5000,]
data_set$P_ID <- factor(data_set$P_ID)

data_in <- data_set[,c("P_ID", "Soundm","Age","Noise_sensitivity","ToDAfternoon",
"ToDEvening", "GenderFemale", "SDNN")]

data_in_split <- split(data_in,data_in$P_ID)

## K method: K is number of participants and that many normal distributions 
#  to be estimated
## VARIABLE NAMES:
## delta = Markov chain initial state distribution
## intercept = Intercept values of transition probability stochastic function
## betaS = Coefficients of covariates in the logit regression for 
#  transition probabilities
## mu = Mean of the Normal desnity assumed for state distribution probability 
#  function or emission probability
## sigma = Std. deviation of the Normal desnity assumed for state distribution 
#  probability function or emission probability

## Natural parameters to working parameters
pn2pw_mix <- function(delta,intercept,betaS,mu,sigma)
{ 
  tdelta <- logit(delta)     ## Transfrom from [0,1] to (-inf,inf)
  tsigma <- log(sigma)        ## Transfrom from [0,inf) to (-inf,inf)
  ## No need of transformation for mu, intercept, betaS as they 
  #  are already unconstrained
  parvect <- c(tdelta,intercept,betaS,mu,tsigma) 
  return(parvect)
}

## Working parameters back to natural parameters
pw2pn_mix <- function(parvect,K)
{
  delta <- inv.logit(parvect[1])  ##Not transformed
  intercept <- matrix(parvect[2:3],nrow=2)
  betaS <- matrix(parvect[4:15],nrow=2)
  mu <- matrix(parvect[16:(2*K+15)],nrow=2)
  sigma <- matrix(exp(parvect[(2*K+16):(15+4*K)]),nrow=2)
  return(list(delta=delta,betaS=betaS,intercept=intercept,mu=mu,sigma=sigma))
}

## Initialize the parameters
mu01<-c(c((min(data_in$SDNN)+4*mean(data_in$SDNN))/5,
          (max(data_in$SDNN)+4*mean(data_in$SDNN))/5))
sigma01 <- c(sd(data_in$SDNN),sd(data_in$SDNN))
delta0 <- 0.1
betaS0<-matrix(rep(0.1,12),nrow=2)
intercept0 <- matrix(c(0.2,0.9),nrow=2)

K <- length(data_in_split)

mu0 <- matrix(rep(mu01,K),nrow=2)
sigma0 <- matrix(rep(sigma01,K),nrow=2)

## Parameter vector, check if we get the initial values back 
# after inverse transformation
parvect <- pn2pw_mix(delta0,intercept0,betaS0,mu0,sigma0) 
pw2pn_mix(parvect,K)

mllk_mix <-function(parvect,X)
{
  allprobs <- list()
  K <- length(X)
  pn<-pw2pn_mix(parvect,K)
  llk <- rep(NA,K)
  for(part in 1:K)
  {
    allprobs[[part]] <- matrix(NA,nrow=nrow(X[[part]]),ncol=2) 
    
    for (j in 1:2)  ##j is hidden state
    {
      allprobs[[part]][,j] <- dnorm(X[[part]][,"SDNN"],mean=pn$mu[j,part],
      sd=pn$sigma[j,part])  ## Conditional probability for each state  
    }
    
    allprobs[[part]] <- ifelse(!is.na(allprobs[[part]]),allprobs[[part]],1)
    lscale <- 0
    
    n <- nrow(X[[part]])

    cov1<-X[[part]][,"Soundm"]
    cov2<-X[[part]][,"Age"]
    cov3<-X[[part]][,"Noise_sensitivity"]
    cov4<-X[[part]][,"ToDAfternoon"]
    cov5<-X[[part]][,"ToDEvening"]
    cov6<-X[[part]][,"GenderFemale"]

    for (i in 1:n)
    {
      if (i==1){ foo<-c(pn$delta,1-pn$delta) } ## Initialize foo
      
      gamma<-matrix(NA,nrow=2,ncol=2)  ## Initialize gamma matrics
      for (j in 1:2)
      {
        gamma[j,j] <- inv.logit( pn$intercept[j]   + pn$betaS[j,1]*cov1[i] + 
                    pn$betaS[j,2]*cov2[i] + pn$betaS[j,3]*cov3[i] + 
                pn$betaS[j,4]*cov4[i]  + pn$betaS[j,5]*cov5[i] + pn$betaS[j,6]*cov6[i] )
        gamma[j,-j] <- 1 - gamma[j,j]
      }
      
      foo<-foo %*% gamma*allprobs[[part]][i,]                                            
      sumfoo <- sum(foo)                                                        
      lscale<-lscale+log(sumfoo)
      foo<-foo/sumfoo
    }
    llk[part] <- - lscale
  }
  mllk <- sum(llk)
  return (mllk)
}   

mle_mix <- function(data_in,X,delta0,intercept0,betaS0,mu0,sigma0)
{
  parvect0<-pn2pw_mix(delta0,intercept0,betaS0,mu0,sigma0)               
  
  ## NUMERICAL OPTIMIZATION of parameters using packages 'nlm' and methods 
  # in 'optim' package. 
  # We tried different methods in optim, 'BFGS' is shown below. 
  # The control parameters for nlm and optim methods were also tweaked 
  # to check speed of estimation 
  
  #mod <- nlm(mllk_mix,X=X,parvect0,print.level=2,iterlim=200,stepmax=5,
  hessian=F,ndigit=4,steptol=1e-4,gradtol=1e-4)
  #mod <- nlm(mllk_mix,X=X,parvect0,print.level=2,iterlim=1500,
  # stepmax=5,hessian=F,ndigit=4)
  mod <- optim(parvect0,mllk_mix,X=X,method="BFGS",hessian=F)
  pn<-pw2pn_mix(mod$estimate,K=length(X))   
  mllk <- mod$minimum 
  np <- length(parvect0)
  AIC <- 2*(mllk+np)
  n <- nrow(data_in)
  BIC <- 2*mllk+np*log(n)
  list(mu=pn$mu,sigma= pn$sigma, delta = pn$delta, intercept = pn$intercept, 
  betaS=pn$betaS,mllk=mllk,AIC=AIC,BIC=BIC)
}

## Repeat below code for each optimization technique 
# (i.e. using nlm, optim(CG), optim(BFGS), etc)
system.time(
  HMM_mix1 <- mle_mix(data_in,data_in_split,delta0,intercept0,betaS0,mu0,sigma0)
)

HMM_mix1