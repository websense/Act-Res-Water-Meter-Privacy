## High res populations
## Calculate exposure for pk=1 scenario with various ks

## S1. load libraries and data
{
  source("./FaithfulBetaSamplingFunctions.R")
  RESULTSHIGHRES = "../resultsHighres/"
  FIGURESDIR = "../../figures/"
}


## NEW set N to achieve a target K
{
  exposedres = c()
  for (agi in 1:6) { 
    ag = highresaggselection[agi]
    pdf(file=paste(FIGURESDIR,
                   sprintf("highres%dsecTargetPopulationForK.pdf",highresaggseconds[agi]),
                   sep=""),
        width=12,height=8)
    {
      par(mfrow=c(2,3))
      for (activityID in 1:6) {
        activityGT = readRDS(file=paste(RESULTSHIGHRES,sprintf("highresActivityGTAct%d.rds",activityID),sep=""))
        baseAM = activityGT[[ag]]
        T = dim(baseAM)[1]
        q = T/4
        q1 = 1:q
        q2 = q + 1:q
        q3 = 2*q + 1:q
        q4 = 3*q + 1:q
        shuffle1 = c(q3,q4,q1,q2)
        shuffle2 = c(q2,q1,q4,q3)
        shuffle3 = c(q4,q3,q2,q1)
        #create an even bigger dataset
        thisAM = cbind(baseAM,baseAM[shuffle1,],baseAM[shuffle2,],baseAM[shuffle3,])
        N = dim(thisAM)[2]
        
        thesekanons = kanon[1:5] 
        for (ti in 1:length(thesekanons)) {
          targetk = thesekanons[ti]
          
          set.seed(12345) #repeatable random sample population
          res=c()
          for (n in ((1:8)*(N/8))) {  ## test different population sizes
            hhs = sample(1:N,n)
            rs = rowSums(thisAM[,hhs])
            tex = which(rs>0 & rs<=targetk) #exposed times
            if (length(tex)==1) {
              exposedhh = which(thisAM[tex,hhs]>0)
              hex = length(exposedhh)
              safe = n-hex #not exposed
              safehhs = hhs[which(thisAM[tex,hhs]==0)] ## never in a grp <= k
            } else {
              if (length(tex)>1) {
                exposedhh = which(colSums(thisAM[tex,hhs])>0)
                hex = length(exposedhh)
                safe = n-hex #not exposed
                safehhs = hhs[which(colSums(thisAM[tex,hhs])==0)]
              } else {
                hex = 0  
                safe = n 
                exposedhh = c()
                safehhs = hhs
              }
            }
            
            res = rbind(res, c(n,hex))
            
            for (be in betas) {
              rmse = getLargeSampleAccuracy(thisAM,safehhs,be)
              exposedres = rbind(exposedres, c(activityID, ag, targetk, safe, n, N, be, rmse))
            }
          }
          
          if (ti==1) {
            title = sprintf("Exposed HHs %s %s",
                            highresactivitynames[activityID],highresaggnames[agi])
            plot(res,type="b",col=ti,lwd=2,
                 ylim=c(0,N),xlim=c(0,N),
                 xlab="Population Size",ylab="Exposed HHs")
            abline(a=0,b=0.5,lty="dotted") #50% safe line
            abline(a=0,b=1.0,lty="dotted") #100% safe line
            abline(v=1600,lty="dotted")
            axis(4)
            title(title)
          } else {
            lines(res,type="b",col=ti,lwd=2)
          }
        }
        klen=length(kanons)
        legend("topleft",paste("k =",kanons),title="K-anonymity",
               lty=rep(1,times=klen),lwd=rep(2,times=klen),
               col=1:klen)
      }
      par(mfrow=c(1,1))
    }
    dev.off()
    print(tail(exposedres))

  }
  colnames(exposedres) = c("actID","Agg","Kanon","Inliers","N","NN","Beta","RMSE")
  exposedres = as.data.table(exposedres)
  activityID = 2
  
  
  saveRDS(exposedres,paste(RESULTSHIGHRES,"largePopulationHighres.rds",sep=""))
}


#plot single figures for paper: low res all acts, highres shower 10 sec, taps 30s clothes 5min, 
pdf(file=paste(FIGURESDIR,
               sprintf("SideAttackOutliersFigure.pdf",highresaggseconds[agi]),
               sep=""),
    width=10,height=3.2) #was 16,4.5
{
caption = c()
par(mfrow=c(1,4))
  #low res data
  {
    ## Low res data
    RESULTSLOWRES = "../resultsLowres/"
    lowresmergedoutliers = as.data.table(readRDS(paste(RESULTSLOWRES,"lowresOutlierHHs.rds",sep="")))
    
    lowresoutliers = cbind(lowresmergedoutliers[,1:2], rowSums(lowresmergedoutliers[,4:3560]==1))
    NN=3557
    kas = unique(lowresoutliers$kanon)[1:5] #ignore 50 here for comp with highres
    aggs = unique(lowresoutliers$agg)
    for (ki in 1:length(kas)) {
      if (ki==1) {
        res = lowresoutliers[kanon==kas[ki], ]
        title = "Low Resolution All Activities"
        plot(res$V2,type="b",col=ki,lwd=2,
             ylim=c(0,NN),xaxt="n",
             xlab="Temporal Resolution\n(a)",ylab="Exposed HHs")
        axis(1,1:6,lowresaggnames)
        abline(h=NN/2,lty="dotted") #50% safe line
        abline(h=NN,lty="dotted") #100% safe line
        axis(4)
        #title(title) do this in caption
        caption = c(caption, "(a)", title)
      } else {
        res = lowresoutliers[kanon==kas[ki], ]
        lines(res$V2,type="b",col=ki,lwd=2)
      }
    }
    klen=length(kas)
    legend("topleft",paste("k =",kas),title="K-anonymity",bg="white",
           lty=rep(1,times=klen),lwd=rep(2,times=klen),
           col=1:klen)
  }
  
  ## high res data
  {
    exposedres = readRDS(paste(RESULTSHIGHRES,"largePopulationHighres.rds",sep=""))
    kas=unique(exposedres$Kanon)
    NN = exposedres$NN[1]
    #activityID=2

    for (i in 1:3) 
    { #3 subfigs
      activityID = c(2,3,4)[i]
      agi = c(1,2,4)[i]
      ag = highresaggselection[agi]
      subfig = paste("Population Size\n",c("(b)","(c)","(d)")[i],sep="")
      for (ki in 1:5) {
        if (ki==1) {
          res = exposedres[actID==activityID & Agg==ag & Kanon==kas[ki] & Beta==0.1, ]
          res$exposed = res$N - res$Inliers
          title = sprintf("%s %s",
                          highresactivitynames[activityID],highresaggnames[agi])
          plot(res[,.(N,exposed)],type="b",col=ki,lwd=2,
               ylim=c(0,NN),xlim=c(0,NN),
               xlab=subfig,ylab="Exposed HHs")
          abline(a=0,b=0.5,lty="dotted") #50% safe line
          abline(a=0,b=1.0,lty="dotted") #100% safe line
          axis(4)
          #title(title) save for caption
          caption = rbind(caption, c(c("(b)","(c)","(d)")[i], title))
        } else {
          res = exposedres[actID==activityID & Agg==ag & Kanon==kas[ki] & Beta==0.1, ]
          res$exposed = res$N - res$Inliers
          lines(res[,.(N,exposed)],type="b",col=ki,lwd=2)
        }
      }
        klen=length(kas)
        legend("topleft",paste("k =",kas),title="K-anonymity",bg="white",
               lty=rep(1,times=klen),lwd=rep(2,times=klen),
               col=1:klen)
    } 
  }
par(mfrow=c(1,1))
}
dev.off()


## general figure of all cases
pdf(file=paste(FIGURESDIR,"highresPopulationVSresults10secto1min.pdf",sep=""),
    width=15,height=10)
{
  exposedres = readRDS(paste(RESULTSHIGHRES,"largePopulationHighres.rds",sep=""))
  par(mfrow=c(3,5))
  exposedres$rmsepct = exposedres$RMSE*100
  for (agipos in 1:3) { #show 10 sec, 30sec adn 1 min
    agi = highresaggselection[agipos]
    rr = range(exposedres[Kanon<=50 & Beta<=0.25 & Agg==agi ,rmsepct],na.rm=TRUE)
    
    for (bi in betas[1:3]) {
      for (ki in kanon[5:1]) {
        
        feasible = exposedres[Kanon==ki & Beta==bi & Agg==agi,]
        
        plot(feasible[actID==1,.(N,rmsepct)],type="b",lwd=2,
             col=highresactivitycols[1],ylim=rr,
             xlab="Population size",ylab="Sample RMSE (%)",
             main=sprintf("K=%d B=%d%% %s",ki,round(bi*100),
                          highresaggnames[agipos]))
        for (ai in 2:6) {
          lines(feasible[actID==ai,.(N,rmsepct)],
                type="b",lwd=2,col=highresactivitycols[ai])
        }
      }
    }

    legend("topleft",highresactivitynames[1:6],
           lty=rep(1,times=6),lwd=rep(2,times=6),
           col=highresactivitycols[1:6])
  }
  par(mfrow=c(1,1))
}
dev.off()


## smaller figure for paper, needs to show different k and resolution
{
pdf(file=paste(FIGURESDIR,"highresSampleSelection.pdf",sep=""),
      width=10,height=3)
  {
    exposedres = readRDS(paste(RESULTSHIGHRES,"largePopulationHighres.rds",sep=""))
    caption = c()
    par(mfrow=c(1,4))
    exposedres$rmsepct = exposedres$RMSE*100
    casesKBA = rbind(
      c(30,0.05,5),
      c(30,0.05,2),
      c(10,0.05,2),
      c(1,0.25,1))
    subfig = c("(a)","(b)","(c)","(d)")
    rr = range(exposedres[Kanon<=50 & Beta<=0.25 & Agg<=300 ,rmsepct],na.rm=TRUE)
    for (i in 1:4) {
      ki = casesKBA[i,1]
      bi = casesKBA[i,2]
      agipos = casesKBA[i,3]
      agi = highresaggselection[agipos]
      feasible = exposedres[Kanon==ki & Beta==bi & Agg==agi,]
      plot(feasible[actID==1,.(N,rmsepct)],type="b",lwd=2,
               col=highresactivitycols[1],ylim=rr,
               xlab=paste("Population size\n",subfig[i],sep=""),ylab="Sample RMSE (% HHs)")
      title=sprintf("K=%d B=%d%% %s",ki,round(bi*100),
                           highresaggnames[agipos])
      caption = rbind(caption,c(subfig[i],title))
             
          for (ai in 2:6) {
            lines(feasible[actID==ai,.(N,rmsepct)],
                  type="b",lwd=2,col=highresactivitycols[ai])
          }
        grid()
        if (i==4) { #same caption for all 4
        legend("topleft",highresactivitynames[1:6],bg="white",ncol=1,
             lty=rep(1,times=6),lwd=rep(2,times=6),
             col=highresactivitycols[1:6])
        }
       
    }

    par(mfrow=c(1,1))
  }
dev.off()
}

## show other results
if (FALSE)
{
plot(exposedres$RMSE,type="p",col=exposedres$actID)
boxplot(exposedres$RMSE ~ exposedres$actID)

NN = exposedres$NN[1]

feasible = exposedres[!is.na(exposedres$RMSE) & RMSE<0.01,]

boxplot(feasible$RMSE ~ feasible$N)
boxplot(feasible$RMSE ~ feasible$Kanon)
boxplot(feasible$RMSE ~ feasible$actID)
boxplot(feasible$RMSE ~ feasible$Agg,
        names.arg=highresaggnames)
boxplot(feasible$RMSE ~ feasible$Beta)

## TODO join all acts?
table(feasible[Kanon==20,.(Beta,N)])
table(feasible[Kanon==10,.(Beta,N)])
table(feasible[Kanon==5,.(Beta,N)])
}

