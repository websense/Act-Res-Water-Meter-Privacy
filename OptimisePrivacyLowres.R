## OptimisePrivacy for low resolution data
## Use analysis results from Accuracy, Unicity and Sampling
## Find pareto optimal parameter settings for utility and privacy
## RCO July 2022

## Load anaylsis results from RESULTSDIR
{
  source("./GlobalParameters.R")
  RESULTSLOWRES = "../resultsLowres/"
  
  #ActivityID AggHours       MCC
  lowresaccuracy = as.data.table(readRDS(file=paste(RESULTSLOWRES,"lowresaccuracy.rds",sep=""))) 
  #beta kanon       eps        delta
  diffprivmatrix = readRDS(file="./diffprivmatrix.rds") 
  #Agg Kanon pknown beta      Leaks      Habits   Occupancy    Any Use
  lowresfaithfulnessexperiments = readRDS(file=paste(RESULTSLOWRES,"lowresfaithfulnessexperiments.rds",sep=""))
}

## Merge results for pareto exploration: only 4 Aggs here
{
  A = lowresaccuracy[order(AggHours,ActivityID),c("ActivityID","AggHours","MCC")]
  
  A1 = A[,lapply(.SD, min),by=c("AggHours")] #worst case min accuracy
  A1$ActivityID = 123
  
  B = lowresfaithfulnessexperiments
  B1 = B[Kanon<50,lapply(.SD,max),by=c("Agg,Kanon,beta")] #same k as highres
  B1$maxRMSEpct = B1$BestForAllSamplesRMSE*100  
  
  C = merge(A1,B1,by.x=c("AggHours","ActivityID"),by.y=c("Agg","ActivityID"))
  
  D = merge(diffprivmatrix,C,by.x=c("kanon","beta"),by.y=c("Kanon","beta")) 
  
  E = D[!is.na(rowMeans(D)),] #from 96 to 66
  E = E[order(E$eps,decreasing=TRUE),]
  
  if (runalltests) {
    par(mfrow=c(2,2))
    plot(E$MCC,type="l",col=1,main="Accuracy")
    plot(E$maxRMSEpct,col=2,type="l",main="Faithfulness % RMSE")
    plot(E$kanon,type="l",main="K-anon",col=3)
    plot(E$eps,type="l",main="DP epsilon",col=4)
    par(mfrow=c(1,1))
  }
  

}



## Trial 1: Pareto Front based on MCC, RMSE and EPS
# Already limited to Lalpha=AggHours<=4, kanon<=30, N=3557, beta from .05 to .5
# and assuming A=123, p<=8, 
# gives dim E is 52 by 9
# and pareto front of 11 cases
{
# output variables to optimise on: min is best for all of these (MCC reversed)
# for any inputs k,alpha,beta,N
P1 = E[,.(1-MCC,maxRMSEpct,eps,delta)] #same choice as highres kanon
nP = dim(P1)[1] #number of records
P1front = rep(1,times=nP) #1 if this dominates, set to 0 if any other dominates it
for (r in 1:nP) 
{
  #get set that dominates r
  for (i in 1:nP) {
    if (i != r & dominates(P1[i,],P1[r,])) {
      P1front[r] = 0 #this r is dominated by at least one other (that is i)
    }
  }
}

par(mfrow=c(2,2))
pf = which(P1front==1)
plot(1-P1$V1,type="l",col="gray",ylab="MCC",main="Pareto Front for 3 vars")
lines(pf,(1-P1$V1)[pf],type="p",col="red",pch=20)
plot(P1$maxRMSEpct,type="l",col="gray",ylab="Sampling RMSE")
lines(pf,(P1$maxRMSEpct)[pf],type="p",col="red",pch=20)
plot(P1$eps,type="l",col="gray",ylab="Differential Privacy Epsilon")
lines(pf,(P1$eps)[pf],type="p",col="red",pch=20)
par(mfrow=c(1,1))

E[which(P1front==1),][order(eps,decreasing=FALSE),]
}


## Table for paper
## show optimal points for P1 #maybe cut some for space RETURN HERE
xtable(E[which(P1front==1),.(kanon,round(beta*100),AggHours,3557,
           eps,delta,MCC,maxRMSEpct)][order(kanon,V2),],
       digits=c(0,0,0,0, 0, 3,-1, 1,3))


## Pareto Fronts (dominating points marked)
{
  
  ## use min max norm to map all to 0,1
  minmaxnorm <- function(vec) {
    mi = min(vec)
    mx = max(vec)
    normvec = (vec-mi)/(mx-mi)
    return (normvec)
  }
  
  {
    ## E, P1 front Pareto Front based on MCC, RMSE and EPS, DELTA
    pdf(file=paste(FIGURESDIR,"LowresParetoFront52cases11dom.pdf",sep=""),
        width=7,height=6)
    {
      par(mar=c(2, 6, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P1front==1)

      ## privacy min eps
      plot(minmaxnorm(E$eps),type="l",col="red",pch=1,cex=0.25,
           xlab="", #"Parameter Combinations",
           xaxt="n",
           ylim=c(0,5.5),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(E$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
      #title("Low Resolution Pareto Front")
      
      ## privacy delta
      lines((minmaxnorm(E$delta)+1.5),type="l",col="red")
      lines(pf,(minmaxnorm(E$delta)+1.5)[pf],type="p",col="red",pch=19,cex=1)
      
      
      ## utility max minMCC
      lines(minmaxnorm(E$MCC)+3,type="l",col="blue")
      lines(pf,(minmaxnorm(E$MCC)+3)[pf],type="p",col="blue",
            pch=19,cex=1,ylab="",yaxt="n")
      
      
      ## utility min RMSE
      lines(minmaxnorm(E$maxRMSEpct)+4.5,type="l",col="blue")
      lines(pf,(minmaxnorm(E$maxRMSEpct)+4.5)[pf],type="p",col="blue",
            pch=19,cex=1,ylab="",yaxt="n")
      
      
     
      abline(h=c(0,1,1.5,2.5,3,4, 4.5,5.5),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5),tick=FALSE,las=1,
           c("DP Epsilon","DP Delta","MCC","RMSE %"))
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5),las=1,
           c(round(range(E$eps),digits=2), 
             signif(range(E$delta),digits=2),
             round(range(E$MCC),digits=2),
             round(range(E$maxRMSE),digits=2)))
      abline(v=pf,lty="dotted",col="gray")
      
    }
    dev.off()
  }
  
}



