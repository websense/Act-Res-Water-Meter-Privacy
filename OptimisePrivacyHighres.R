## OptimisePrivacy for high res data
## Use analysis results from Accuracy, Unicity and Sampling
## Find pareto optimal parameter settings for utility and privacy
## RCO Aug 2022

## Load highres anaylsis results from RESULTSDIR
{
  source("./GlobalParameters.R")
  RESULTSHIGHRES = "../resultsHighres/"
  
  highresaccuracy = as.data.table(readRDS(file=paste(RESULTSHIGHRES,"highresaccuracy.rds",sep=""))) 
  #ActivityID AggSeconds      TP FP       TN    FN       MCC
  
  diffprivmatrix = readRDS(file="./diffprivmatrix.rds") 
  # beta kanon        eps        delta
  
  exposedres = readRDS(paste(RESULTSHIGHRES,"largePopulationHighres.rds",sep=""))
  #actID Agg Kanon Inliers    N   NN Beta        RMSE
}

## Merge all experimental data
{
  A = highresaccuracy[order(AggSeconds,ActivityID),c("ActivityID","AggSeconds","MCC")]
  #only for activities 1 to 6
  A1 = A[ActivityID<=6,lapply(.SD, min),by=c("AggSeconds")] #worst case min MCC for all activities
  A1$ActivityID = 123456
  A1[,minMCC:=MCC]
  A1$MCC=NULL
  
  ##So far  only Kanon 1 or 5 feasible but more possible with larger populations [Kanon==1,] #Kanon=5 too hard
  B = exposedres
  B$maxRMSEpct = B$RMSE*100
  B$AggSeconds = B$Agg*10
  
  #join activities, report max RMSEpct
  B1 = c()
  for (ai in unique(B$AggSeconds)) {
    for (ki in unique(B$Kanon)) {
      for (bi in unique(B$Beta)) {
        for (ni in unique(B$N)) {
          bb = B[AggSeconds==ai & Kanon==ki & Beta==bi & N==ni,lapply(.SD,max)] #worst case RMSE is max for some activity
          B1 = rbind(B1,bb)
        }}}}
  B1 = as.data.table(B1)
  B1$ActivityID=123456
  B1 = B1[,.(ActivityID,AggSeconds,Kanon,Beta, Inliers,N,maxRMSEpct)]
  #if any activity is infeasible then maxRMSEpct will be NA
  
  
  C = merge(A1,B1,by=c("ActivityID","AggSeconds"))
  
  D = merge(diffprivmatrix,C,by.x=c("kanon","beta"),by.y=c("Kanon","Beta")) 
  
  if (runalltests) {
    plot(D$minMCC,type="b",ylim=c(0,3),col="red")
    lines(D$maxRMSEpct,type="b",col="blue")
    lines(D$eps,type="b",col="green")
  }
  
  E = D[!is.na(rowMeans(D)),] #remove infeasible combinations, leaves 211 combinations from 960
  E = E[order(E$eps,decreasing=TRUE),]
  
  pdf(file=paste(FIGURESDIR,"HighResParetoData.pdf",sep=""),width=10,height=8)
  {
    par(mfrow=c(3,3))
    plot(E$eps,type="l",main="DP epsilon",col=4)
    plot(E$kanon,type="l",main="K-anon",col=3)
    plot(E$beta,type="l",main="Beta",col=2)
    plot(E$AggSeconds,type="l",main="Temporal Resolution",col=1)
    
    plot(E$minMCC,type="l",col=2,main="Accuracy")
    plot(E$maxRMSEpct,col=2,type="l",main="Faithfulness % RMSE")
    plot(E$N,type="l",col=1,main="Population")
    par(mfrow=c(1,1))
  }
  dev.off()
  
  
}


## Trial 1: Pareto Front based on privacy and utility outputs MCC, RMSE and EPS, DELTA where N is an output
{
  #variables to optimise on: min is best for all of these (MCC reversed)
  P1 = E[,.(1-minMCC,maxRMSEpct,eps,delta)] 
  nP = dim(P1)[1] #number of records
  P1front = rep(1,times=nP) #1 if this dominates, will be set to 0 if any other tuple dominates it
  for (r in 1:nP) 
  {
    #get set that dominates r
    for (i in 1:nP) {
      if (i != r & dominates(P1[i,],P1[r,])) {
        P1front[r] = 0 #this r is dominated by at least one other (that is i)
      }
    }
  }
  
  
  
  par(mfrow=c(2,4))
  pf = which(P1front==1)
  plot(1-P1$V1,type="l",col="gray",ylab="minMCC",main="Priv/Utility outputs")
  lines(pf,(1-P1$V1)[pf],type="p",col="red",pch=20)
  plot(P1$maxRMSEpct,type="l",col="gray",ylab="Sampling RMSE (pct)")
  lines(pf,(P1$maxRMSEpct)[pf],type="p",col="red",pch=20)
  plot(P1$eps,type="l",col="gray",ylab="Differential Privacy Epsilon")
  lines(pf,(P1$eps)[pf],type="p",col="red",pch=20)
  plot(P1$delta,type="l",col="gray",ylab="Differential Privacy Delta")
  lines(pf,(P1$delta)[pf],type="p",col="red",pch=20)
  
  plot(E$kanon,type="l",col="gray",ylab="Kanon",main="Sample Inputs")
  lines(pf,(E$kanon)[pf],type="p",col="red",pch=20)
  
  plot(E$beta,type="l",col="gray",ylab="Beta sample")
  lines(pf,(E$beta)[pf],type="p",col="red",pch=20)
  
  plot(E$AggSeconds,type="l",col="gray",ylab="Alpha resolution")
  lines(pf,(E$AggSeconds)[pf],type="p",col="red",pch=20)
  
  plot(E$N,type="l",col="gray",ylab="Population Size")
  lines(pf,(E$N)[pf],type="p",col="red",pch=20)
  par(mfrow=c(1,1))
  
  E[which(P1front==1),] 
    # 9 solutions from 211 based on 4 outputs: eps,delta and mcc,rmse
    # 28 solutions from 211 when N is included in search
    # 8 solutions when opt only MCC,RMSE,eps,delta all are N 5600 or 6400 and a 10 or 30 sec
  
  #table for paper
  xtable(E[which(P1front==1),.(kanon,round(beta*100),AggSeconds,N,
                               eps,delta,minMCC,maxRMSEpct)][order(kanon,V2),],
         digits=c(0,0,0,0, 0, 3,-1, 1,3))
}


  
  ## Plot Pareto Fronts (dominating points marked)
{

  ## use min max norm to map all to 0,1
  minmaxnorm <- function(vec) {
    mi = min(vec)
    mx = max(vec)
    normvec = (vec-mi)/(mx-mi)
    return (normvec)
  }
  
    ## E, P1 front Pareto Front based on MCC, RMSE and EPS
    pdf(file=paste(FIGURESDIR,"HighresParetoFront211cases9dom.pdf",sep=""),
        width=7,height=6)
    {
      par(mar=c(2, 6, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P1front==1)
      
      ## privacy min eps
      plot(minmaxnorm(E$eps),type="l",col="red",
           xlab="", #"Parameter Combinations",
           xaxt="n",
           ylim=c(0,5.5),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(E$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
    
      ## privacy min delta
      lines((minmaxnorm(E$delta)+1.5),type="l",col="red")
      lines(pf,(minmaxnorm(E$delta)+1.5)[pf],type="p",col="red",pch=19,cex=1)
      
      ## utility max minMCC
      lines(minmaxnorm(E$minMCC)+3,type="l",col="blue")
      lines(pf,(minmaxnorm(E$minMCC)+3)[pf],type="p",col="blue",
            pch=19,cex=1,ylab="",yaxt="n")
      
      ## utility min RMSE
      lines(minmaxnorm(E$maxRMSEpct)+4.5,type="l",col="blue")
      lines(pf,(minmaxnorm(E$maxRMSEpct)+4.5)[pf],type="p",col="blue",
            pch=19,cex=1,ylab="",yaxt="n")
      
      abline(h=c(0,1,1.5,2.5,3,4, 4.5,5.5),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5),c("DP Epsilon","DP Delta","MCC","RMSE %"),tick=FALSE,las=1)
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5),las=1,
           c(round(range(E$eps),digits=2), 
             signif(range(E$delta),digits=1), 
             round(range(E$minMCC),digits=2),
             round(range(E$maxRMSEpct),digits=2)))
      abline(v=pf,lty="dotted",col="gray")
      
    }
    dev.off()
  
}

## END HERE

