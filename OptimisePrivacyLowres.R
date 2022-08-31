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
  B1$ActivityID=123
  
  C = merge(A1,B1,by.x=c("AggHours","ActivityID"),by.y=c("Agg","ActivityID"))
  
  D = merge(diffprivmatrix,C,by.x=c("kanon","beta"),by.y=c("Kanon","beta")) 
  
  E = D[!is.na(rowMeans(D)),] #from 96 to 66
  E = E[order(E$eps,decreasing=TRUE),]
  
  if (runalltests) {
    par(mfrow=c(2,2))
    plot(E$MCC,type="l",col=1,main="Accuracy")
    plot(E$SampleRMSE*100,col=2,type="l",main="Faithfulness % RMSE")
    plot(E$kanon,type="l",main="K-anon",col=3)
    plot(E$eps,type="l",main="DP epsilon",col=4)
    par(mfrow=c(1,1))
  }
  

}



## Trial 1: Pareto Front based on MCC, RMSE and EPS
{
#variables to optimise on: min is best for all of these (MCC reversed)
P1 = E[,.(1-MCC,SampleRMSE,eps)] #same choice as highres kanon
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
plot(P1$SampleRMSE,type="l",col="gray",ylab="Sampling RMSE")
lines(pf,(P1$SampleRMSE)[pf],type="p",col="red",pch=20)
plot(P1$eps,type="l",col="gray",ylab="Differential Privacy Epsilon")
lines(pf,(P1$eps)[pf],type="p",col="red",pch=20)
par(mfrow=c(1,1))

E[which(P1front==1),][order(kanon,beta),]
}


## Trial 2: Pareto Front based on RMSE and EPS for any MCC > 0.9 and RMSE < 2% (0.02) #1st Q is 1.5
{
#variables to optimise on: min is best for all of these

E2 = E[MCC>0.9 & SampleRMSE<0.02 & kanon<50,][order(eps),]
P2 = E2[,.(SampleRMSE,eps)]
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


## show optimal points
xtable(E[which(P1front==1),.(kanon,round(beta*100),AggHours,eps,delta,MCC,SampleRMSE*100)],
       digits=c(0,0,0,0, 3,-1, 1,1))

## show pareto front results
pdf(file=paste(FIGURESDIR,"LowresParetoFrontMCCgr85pk6.pdf",sep=""),width=6,height=6)
{
par(mar=c(5, 4, 4, 4) + 0.1, mfrow=c(1,1))
pf = which(P2front==1)
plot(P2$eps,type="b",col="red",
     xlab="Parameter Combinations (K-anon indicated)",xaxt="n",
     ylim=c(0,1.5),ylab="",yaxt="n")
     #main="Pareto Front for 2 vars")
lines(pf,(P2$eps)[pf],type="p",col="red",pch=19,cex=1.5,ylab="",yaxt="n")
#axis(1,at=pf,labels=E2$AggHours[pf]) #all selections are 1 hour
axis(1,at=pf,labels=E2$kanon[pf])
axis(2,at=(0:7)*0.2,labels=(0:7)*0.2,ylab="",col.axis="red",col.ticks ="red")
# add a title for the right axis
mtext("Differential Privacy Epsilon", side=2, line=3, cex.lab=1,las=0, col="red")

#show RMSE as percentage (of HHs population)
rscale = 40
lines(P2$SampleRMSE*rscale,type="b",col="blue",ylab="",yaxt="n")
lines(pf,(P2$SampleRMSE*rscale)[pf],type="p",col="blue",pch=19,cex=1.5,ylab="",yaxt="n")
axis(4,at=(0:6)*0.2,labels=100*(0:6)*0.2/rscale,ylab="",col.axis="blue",col.ticks ="blue")
mtext("Sampling RMSE %", side=4, line=3, cex.lab=1,las=0, col="blue")

#lines(E2$kanon/50,type="b",col="green",ylab="",yaxt="n")
#lines(pf,(E2$kanon/50)[pf],type="p",col="green",pch=19,cex=1.5,ylab="",yaxt="n")

#lines(E2$beta,type="b",col="orange",ylab="",yaxt="n")
#lines(pf,(E2$beta)[pf],type="p",col="orange",pch=19,cex=1.5,ylab="",yaxt="n")


grid(ny=NA,nx=NA)
abline(v=pf,col="gray",lty="dotted")
abline(h=c(0,0.4,0.8,1.2),col="gray",lty="dotted")
legend("topleft",title="Pareto Front (dominating points marked)",bg="white",
       c("Differential Privacy Epsilon","Sampling RMSE"),lty=c(1,1),
       col=c("red","blue"),pch=c(19,19),pt.cex=1.5)
par(mfrow=c(1,1))
}
dev.off()


