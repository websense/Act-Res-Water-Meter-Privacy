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
{
#variables to optimise on: min is best for all of these (MCC reversed)
P1 = E[,.(1-MCC,maxRMSEpct,eps)] #same choice as highres kanon
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

E[which(P1front==1),][order(kanon,beta),]
}


## Trial 2: Pareto Front based on RMSE and EPS for any MCC > 0.9 and RMSE < 2% (0.02) #1st Q is 1.5
{
#variables to optimise on: min is best for all of these

E2 = E[MCC>0.9 & BestForAllSamplesRMSE<0.015 & kanon<50,][order(eps),]
P2 = E2[,.(BestForAllSamplesRMSE,eps)]
nP2 = dim(P2)[1] #number of records
P2front = rep(1,times=nP2) #1 if this dominates, set to 0 if any other dominates it
for (r in 1:nP2) 
{
  #get set that dominates r
  for (i in 1:nP2) {
    if (i != r & dominates(P2[i,],P2[r,])) {
      P2front[r] = 0 #this r is dominated by at least one other
    }
  }
}

E2[which(P2front==1),][order(kanon,beta),] #same as E1
}


## Table for paper
## show optimal points for P2
xtable(E2[which(P2front==1),.(kanon,round(beta*100),AggHours,
          eps,delta,MCC,BestForAllSamplesRMSE*100)][order(kanon,V2),],
       digits=c(0,0,0,0, 3,-1, 1,1))


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
    
    ## E, P1 front Pareto Front based on MCC, RMSE and EPS
    pdf(file=paste(FIGURESDIR,"LowresParetoFrontMCCRMSEEPS.pdf",sep=""),
        width=7,height=7)
    {
      par(mar=c(2, 9, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P1front==1)
      
      ## utility min eps
      plot(minmaxnorm(E$eps),type="l",col="red",pch=1,cex=0.25,
           xlab="", #"Parameter Combinations",
           xaxt="n",
           ylim=c(0,5.5),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(E$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
      title("Low Resolution Pareto Front Eps,MCC,RMSE")
      
      ## utility min RMSE
      lines(minmaxnorm(E$maxRMSEpct)+1.5,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E$maxRMSEpct)+1.5)[pf],type="p",col="blue",pch=19,cex=1,ylab="",yaxt="n")
      
      ## utility max minMCC
      lines(E$MCC+3,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(E$MCC+3)[pf],type="p",col="blue",
            pch=17,cex=1,ylab="",yaxt="n")
      
      ## free vars
      lines((minmaxnorm(E$delta)+4.5),type="l",col="red",pch=1,cex=0.25)
      lines(pf,(minmaxnorm(E$delta)+4.5)[pf],type="p",col="red",pch=17,cex=1)
      
      abline(h=c(0,1,1.5,2.5,3,4, 4.5,5.5),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5),tick=FALSE,las=1,
           c("DP Epsilon","RMSE","MCC","DP Delta"))
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5),las=1,
           c(round(range(E$eps),digits=2), 
             round(range(E$SampleRMSE*100),digits=2),
             0,1,
             signif(range(E$delta),digits=2) ))
      abline(v=pf,lty="dotted",col="gray")
      
    }
    dev.off()
    
    ## E2,P2 front
    ## Trial 2: Pareto Front E2, P2 front based on min N and min EPS 
    ## for any MCC > 0.9 and RMSE<1.0pct and delta<1-e04: gives 3 contrasting solns, 
    pdf(file=paste(FIGURESDIR,"LowresParetoFrontEPSRMSE.pdf",sep=""),
        width=6,height=4)
    {
      par(mar=c(2, 9, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P2front==1)
      plot(minmaxnorm(E2$eps),type="l",col="red",pch=1,cex=0.5,
           xlab="", #Parameter Combinations",
           xaxt="n",
           ylim=c(0,5.5),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(E2$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
      
      ## utility RMSE
      lines(minmaxnorm(E2$maxRMSEpct)+1.5,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E2$maxRMSEpct)+1.5)[pf],type="p",col="blue",pch=19,cex=1,ylab="",yaxt="n")
      
      
      #utility MCC
      lines(E2$MCC+3,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(E2$MCC+3)[pf],type="p",col="blue",
            pch=17,cex=1,ylab="",yaxt="n")
      
   
      ## privacy bounded
      lines((minmaxnorm(E2$delta)+4.5),type="l",col="red",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E2$delta)+4.5)[pf],type="p",col="red",pch=17,cex=1)
 
      abline(h=c(0,1, 1.5,2.5, 3,4, 4.5,5.5),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5),
           c("DP Epsilon",
             "Sampling RMSE",
             "Activity MCC",
             "DP Delta < 1e-04"),tick=FALSE,las=1)
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5),las=1,
           c(round(range(E2$eps),digits=2), 
             round(range(E2$maxRMSEpct),digits=2), 
             0,1,
             signif(range(E2$delta),digits=1)) )
      abline(v=pf,lty="dotted",col="gray")
    }
    dev.off()
    
  
  }
  
}



