## makeSafeKanonTables
## Calculations for Li2012; these values are independent of the dataset selected
## Creates a lookup table for the parameters to be tested
## Since this table is shared for all datasets, it is saved in the Rcode directory
## R Cardell-Oliver
## Sep 2021

source("./GlobalParameters.R") #for kanon and betas list

# Defns and Function from Li et al
{
  # F(j; n, β) to denotes the cumulative probability mass function
  # for binomial getting exactly j successes in n trials where 
  # each trial succeeds with probability β. And
  # 
  # Theorem 5. Any strongly-safe k-anonymization algorithm 
  # satisfies (β, ǫ, δ)-DPS for any 0 < β < 1, eps ≥ − ln(1 − β), and δ = d(k, β, ǫ), 
  # where the function d is defined as
  # d(k, β, ǫ) = maxn for ni in n:n>ceil(k/g-1) (sum from j>gn to ni (f(j; n, β))
  # where γ = (E^eps−1+β) / E^eps
  
  
  #input ka kanon minimum, beta (in (0,1) proportion for sampling and eps
  #return delta DPS error bound
  getdelta <- function(ka,beta,eps) 
  {
    if (beta<=0 | beta>=1) {
      print("ERROR getdelta requires 0<beta<1"); return(NULL)
    }
    allsi=c()
    g = (exp(1)^eps - 1 + beta) / exp(1)^eps
    n = ceiling( (ka/g) - 1)
    for (ni in n:(n+10)) {
      gn = ceiling(g*ni)
      si = 0
      for (j in gn:ni) {
        si = si + dbinom(j,ni,beta)
      }
      allsi = c(allsi,si)
    }
    return(max(allsi))
  }
}

## Make a lookup table for (beta,eps,delta)-DPS parameter combinations for given kanon
{
  res=c()
  for (ka in kanon) { #k=5 only met for 0.05 sample
    for (beta in betas) { #sampling percent
      epsmin = (-log(1-beta))
      if (epsmin < epsmax) {
        epstry = seq(from=epsmin,to=epsmax,by=0.01) 
        for (eps in epstry) {
          delta = getdelta(ka,beta,eps)
          res = rbind(res, c(ka,beta,eps,delta))
        }
      }
    }
  }
  
  eallb=matrix(NA,nrow=length(betas),ncol=length(kanon))
  for (b in 1:length(betas)) {
    dok = which(res[,2]==betas[b] & res[,4]<=deltamax)
    #print(c(b,length(dok)))
    if (length(dok)>0) {
      bb = boxplot(res[dok,3] ~ res[dok,1],plot=FALSE) #main=b)
      emin = bb$stats[1,] #min e for d<10^-7
      #print(emin)
      epos = which(is.element(kanon,as.integer(bb$names)))
      eallb[b,epos] = emin
    }
  }
  rownames(eallb)=betas #sampling pct
  colnames(eallb)=kanon
}



## TABLE for paper for our kanon and diffprivmatrix for optimisation search
{
  besteps = matrix(nrow=length(betas),ncol=length(kanon))
  diffprivmatrix = c() #version for optimisation
  rownames(besteps)=betas #sampling pct
  colnames(besteps)=kanon #kanon required
  #betas = c(0.05,0.1,0.25,0.5) #sampling percent from GlobalParameters, exclude 0.999
  for (ki in 1:length(kanon)) { #k=5 only met for 0.05 sample so start from 0.1
    k = kanon[ki]
    for (bi in 1:length(betas)) { 
      beta = betas[bi] 
      epsmin = (-log(1-beta))
      epstry = seq(from=epsmin,to=epsmax,by=0.01) 
      epsres = c()
      for (eps in epstry) {
        delta = getdelta(k,beta,eps)
        #print(c(k,beta,eps,delta))
        epsres = rbind(epsres, c(eps,delta))
      }
      dok = which(epsres[,2]<=deltamax) #all acceptable deltas
      if (length(dok)>0) {
        besteps[bi,ki] = min(epsres[dok,1])  #best (min) eps given deltas
        pos = which.min(epsres[dok,1])
        diffprivmatrix = rbind(diffprivmatrix,c(beta,k,epsres[dok,][pos,]))
      } else {
        diffprivmatrix = rbind(diffprivmatrix,c(beta,k,NA,NA)) #no soln found
      }
    }
  }
  ## Matrix used for optimisation
  diffprivmatrix = as.data.table(diffprivmatrix)
  names(diffprivmatrix) = c("beta","kanon","eps","delta")
  saveRDS(diffprivmatrix,file="./diffprivmatrix.rds")
}

## FIGURE for paper
## based on Figure 6 [Li2012]: A graph showing the value of eps 
##satisfied by a given k if delta ≤ 10−6 with varying sampling probabilities.
pdf(paste(FIGURESDIR,"kbetaepsdelta.pdf",sep=""),width=6,height=6)
{
  someka = 1:6 #leave out 100,300
  krange = as.numeric(colnames(eallb))[someka] 
  plot(krange,eallb[1,someka],
       type="b",lty=1,lwd=2,col=1,
       ylim=c(0,1.75),ylab="Epsilon",
       xlab="k for k-anonymity",xaxt="n")
  lines(krange,eallb[2,someka],type="b",lwd=2,col=2,pch=2)
  lines(krange,eallb[3,someka],type="b",lwd=2,col=3,pch=3)
  lines(krange,eallb[4,someka],type="b",lwd=2,col=4,pch=4)
  axis(1,krange,at=krange)
  axis(4)
  abline(v=krange,lty="dotted",col="gray")
  abline(h=c(0:3)*.5,lty="dotted",col="gray")
  legend("topright",paste("Beta=",rownames(eallb),sep=""),
         lty=c(1,1,1,1),lwd=c(2,2,2),pch=1:4,
         col=1:4,bty="n",bg="white")
}
dev.off()



## Make lookup table for parameter combinations - maybe later
if (FALSE)
{
  res=c()
  for (k in 10:50) { #k=5 only met for 0.05 sample
    for (b in c(5,10,20,99)) { #sampling percent
      beta = b/100 
      epsmin = (-log(1-beta))
      if (epsmin < 2.0) {
        epstry = seq(from=epsmin,to=2.0,by=0.01) 
      } else {
        epstry = seq(from=4.5,to=5,by=0.01) #for beta=0.99
      }
      for (eps in epstry) {
        delta = getdelta(k,beta,eps)
        res = rbind(res, c(k,beta,eps,delta))
      }
    }
  }
  
  ## full dataset (no sampling)
  full = res[which(res[,2]==0.99),]
  full[which(full[,3]==4.5),] #min eps
  range(  full[which(full[,3]==4.5),4])
  ##range of delta 0.6050061 0.9043821
  
  
  ## sampling results
  deltaUB = 10^-6
  eallb=c()
  for (b in c(5,10,20)/100) {
    dok = which(res[,2]==b & res[,4]<=deltaUB)
    bb = boxplot(res[dok,3] ~ res[dok,1],plot=FALSE) #main=b)
    emin = bb$stats[1,] #min e for d<10^-7
    eallb = rbind(eallb,emin)
  }
  rownames(eallb)=c(5,10,20)/100
  colnames(eallb)=10:50
}




## END FILE



