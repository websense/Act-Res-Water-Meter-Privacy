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
  #actID Agg Kanon Inliers    N   NN Beta        RMSE   rmsepct
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


## Trial 1: Pareto Front based on MCC, RMSE and EPS
{
  #variables to optimise on: min is best for all of these (MCC reversed)
  P1 = E[,.(1-minMCC,maxRMSEpct,eps,N)]
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
  
  
  
  par(mfrow=c(2,2))
  pf = which(P1front==1)
  plot(1-P1$V1,type="l",col="gray",ylab="minMCC",main="Pareto Front for 4 vars")
  lines(pf,(1-P1$V1)[pf],type="p",col="red",pch=20)
  plot(P1$maxRMSEpct,type="l",col="gray",ylab="Sampling RMSE (pct)")
  lines(pf,(P1$maxRMSEpct)[pf],type="p",col="red",pch=20)
  plot(P1$eps,type="l",col="gray",ylab="Differential Privacy Epsilon")
  lines(pf,(P1$eps)[pf],type="p",col="red",pch=20)
  plot(E$N,type="l",col="gray",ylab="Population Size")
  lines(pf,(E$N)[pf],type="p",col="red",pch=20)
  par(mfrow=c(1,1))
  
  E[which(P1front==1),] #28 solutions
  # kanon beta        eps        delta ActivityID AggSeconds    minMCC Inliers    N maxRMSEpct
  # 1:    20 0.50 1.39314718 7.826090e-05     123456       3600 0.3236606    4800 4800  0.5677702
  # 2:     5 0.10 1.15536052 5.500000e-05     123456         10 0.9686868    5600 5600  0.5298040
  # 3:     5 0.10 1.15536052 5.500000e-05     123456        900 0.6587754    2400 2400  1.6096912
  # 4:     5 0.10 1.15536052 5.500000e-05     123456       3600 0.3236606    1600 1600  1.6586423
  # 5:     5 0.05 0.56129329 6.368983e-05     123456         10 0.9686868    3170 3200  0.7673206
  # 6:     5 0.05 0.56129329 6.368983e-05     123456        900 0.6587754    1600 1600  2.3566354
  # 7:    10 0.10 0.38536052 7.898190e-05     123456         10 0.9686868    5600 5600  0.5305722
  # 8:    10 0.10 0.38536052 7.898190e-05     123456         10 0.9686868    6400 6400  0.5300014
  # 9:    10 0.10 0.38536052 7.898190e-05     123456       3600 0.3236606    2400 2400  1.6470453
  # etc
}

## Trial 2: Pareto Front based on min N and min EPS 
## for any MCC > 0.9 and RMSE<1.0pct : gives 3 contrasting solns, 
## there are 28 for the full data
{
  #variables to optimise on: min is best for all of these
  E2 = E[minMCC>0.9 & maxRMSEpct<1.0,][order(eps,N),] #just take an acceptable range of MCC
  P2 = E2[,.(eps,N)]
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
  
  ## show optimal points
  E2[which(P2front==1),]
  
  ## choose first (larger) set for paper
  xtable(E2[which(P2front==1),.(kanon,beta,AggSeconds,N, eps,delta, minMCC,maxRMSEpct)],
         digits=c(0, 0,2,0,0, 3,-1, 3,1))
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
  
  {
  
    ## E, P1 front Pareto Front based on MCC, RMSE and EPS
    pdf(file=paste(FIGURESDIR,"HighresParetoFrontMCCRMSEEPS.pdf",sep=""),
        width=7,height=7)
    {
      par(mar=c(2, 9, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P1front==1)
      
      ## utility min eps
      plot(minmaxnorm(E$eps),type="l",col="red",pch=1,cex=0.25,
           xlab="", #"Parameter Combinations",
           xaxt="n",
           ylim=c(0,7),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(E$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
      title("Pareto Front Eps,MCC,RMSE")
    
      ## utility min RMSE
      lines(minmaxnorm(E$maxRMSEpct)+1.5,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E$maxRMSEpct)+1.5)[pf],type="p",col="blue",pch=19,cex=1,ylab="",yaxt="n")
      
      ## utility max minMCC
      lines(E$minMCC+3,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(E$minMCC+3)[pf],type="p",col="blue",
            pch=17,cex=1,ylab="",yaxt="n")
      
      ## free vars
      lines((minmaxnorm(E$delta)+4.5),type="l",col="red",pch=1,cex=0.25)
      lines(pf,(minmaxnorm(E$delta)+4.5)[pf],type="p",col="red",pch=17,cex=1)
      
      lines(minmaxnorm(E$N)+6,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E$N)+6)[pf],type="p",col="blue",
            pch=15,cex=1,ylab="",yaxt="n")
      
      abline(h=c(0,1,1.5,2.5,3,4, 4.5,5.5,6,7),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5,6.5),c("DP Epsilon","MCC","RMSE","DP Delta","N"),tick=FALSE,las=1)
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5, 6,7),las=1,
           c(round(range(E$eps),digits=2), 
             round(range(E2$maxRMSEpct),digits=2),
             round(range(E2$minMCC),digits=2),
             signif(range(E2$delta),digits=1), 
             range(E2$N) ))
      abline(v=pf,lty="dotted",col="gray")
      
    }
    dev.off()
    
    ## E2,P2 front
    ## Trial 2: Pareto Front E2, P2 front based on min N and min EPS 
    ## for any MCC > 0.9 and RMSE<1.0pct and delta<1-e04: gives 3 contrasting solns, 
    pdf(file=paste(FIGURESDIR,"HighresParetoFrontEPSN.pdf",sep=""),
        width=6,height=4.5)
    {
      par(mar=c(2, 9, 2, 4) + 0.1, mfrow=c(1,1))
      pf = which(P2front==1)
      plot(minmaxnorm(P2$eps),type="l",col="red",pch=1,cex=0.5,
           xlab="", #Parameter Combinations",
           xaxt="n",
           ylim=c(0,7),ylab="",yaxt="n")
      lines(pf,(minmaxnorm(P2$eps))[pf],type="p",col="red",pch=19,cex=1,ylab="",yaxt="n")
      
      
      lines(minmaxnorm(P2$N)+1.5,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(P2$N)+1.5)[pf],type="p",col="blue",
            pch=15,cex=1,ylab="",yaxt="n")
      
      ## privacy bounded
      lines((minmaxnorm(E2$delta)+3),type="l",col="red",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E2$delta)+3)[pf],type="p",col="red",pch=17,cex=1)
      
      ## utility bounded
      lines(minmaxnorm(E2$maxRMSEpct)+4.5,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E2$maxRMSEpct)+4.5)[pf],type="p",col="blue",pch=19,cex=1,ylab="",yaxt="n")
      
      lines(minmaxnorm(E2$minMCC)+6,type="l",col="blue",pch=1,cex=0.5)
      lines(pf,(minmaxnorm(E2$minMCC)+6)[pf],type="p",col="blue",
            pch=17,cex=1,ylab="",yaxt="n")
 
      
      abline(h=c(0,1, 1.5,2.5, 3,4, 4.5,5.5, 6,7),lty="dotted",col="gray")
      axis(2,at=c(0.5,2,3.5,5,6.5),
           c("DP Epsilon","Population N",
             "DP Delta < 1e-04",
             "Sampling RMSE < 1",
             "Activity MCC > 0.9"),tick=FALSE,las=1)
      axis(4,at=c(0,1, 1.5,2.5, 3,4, 4.5,5.5, 6,7),las=1,
          c(round(range(E2$eps),digits=2), range(E2$N), signif(range(E2$delta),digits=1), 
            round(range(E2$maxRMSEpct),digits=2), round(range(E2$minMCC),digits=2) ))
      abline(v=pf,lty="dotted",col="gray")
    }
    dev.off()
      

  }
  
}

## END HERE

