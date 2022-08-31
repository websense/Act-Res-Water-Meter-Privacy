## makeOutlierTable
## Generic functions for identifing exposed households for side channel attacks of PK known events 
## count of activities or activity frequecy at some time
## Input: alloutlierslist[[k]][[agg]] (for pruning) - remember to ignore c(-1) sets, 
## Input: allGT activities that could be known by the adversary 
## Input: allAR activities the adversary will match to (may not be same)
## Output: list of exposed households each pk,agg combination 
## (pk is number of prior knowledge events) 
## Author: Rachel Cardell-Oliver
## Version: Apr 2022

## load libraries and data
{
  source("./GlobalParameters.R") 
}

## Experiments to identify HHs vulnerable to prior knolwedge attacks
## based on approach proposed by
# Culnane et al, Stop the Open Data Bus, We Want to Get Off, http://arxiv.org/abs/1908.05004

#list of kas, choice of pk for activityGT
#activityGT = readRDS(file=paste(RESULTSHIGHRES,sprintf("highresActivityGTAct%d.rds",activityID),sep=""))[[ag]]
#kavalues = kanon[1:5]
#returns a K x N table where each entry pk is the number of known events 
#needed to expose a HH h (1:N) in a group of k (kanon) similar HHs (ie sharing the same exposing events)
getUnicityOutliersCulnane <- function(thisactivityGT,pkvalues,kavalues,NTRIES=1000)
{
  AM = thisactivityGT
  N = dim(AM)[2] #number of HHs
  klen = length(kavalues) #number of kanon values to cover
  
  #count exposures
  exposureInfo = matrix(0,nrow=length(kavalues),ncol=N)
  rownames(exposureInfo) = kavalues
  
  ## pk=1 case
  for (r in 1:klen) {
    thiska = kavalues[r]
    exposedTimes = which(rowSums(AM==1)<=thiska) #fewer than thiska HHs active at these times
    if (length(exposedTimes) > 0) {
      if (length(exposedTimes) > 1) {
        exposedHHs = which(colSums(AM[exposedTimes,]==1)>0) #HH is active at least one of the exposed times
      } else {
        exposedHHs = which(AM[exposedTimes,]>0)
      }
      exposureInfo[r,exposedHHs] = 1
    }
  }
  
  ## pk=2+ case #maybe pks could be done together - select 8 then subsets
  for (pk in pkvalues[2:length(pkvalues)]) 
  { 
    for (h in 1:N) { 
      hdone = (exposureInfo[klen,h]>0) ##if this house already known to be exposed then no more searching needed
      if (!hdone) {
        #individual times may not be exposed, but the combination may be unique
        #Very slow at ag=360:
        #400 to 800 simultaneous activity at any time - will there be any unique combinations?
        #each HH has around 80 events from 672 time stamps
        activetimes = which(AM[,h]==1)
        if (length(activetimes) >= pk) { 
          #print(c(h,pk))
          hexp=FALSE
          for (i in 1:NTRIES) {
            #print(".")
            if (!hexp) { #try more combinations of pk events to see if h can be exposed
              #MAYBE TODO separate pk choices so at least separate events? - seldom a problem so don't worry
              #al = length(activetimes)
              #activestep = c(2,activetimes[2:al]-activetimes[1:(al-1)])
              #eventtimes = activetimes[which(activestep>1)] #just one entry per event
              #activecombos = combinations(length(eventtimes),pk,eventtimes)  
              ## FULL combination of activities, but too close to stack limit for 1000+ choose 2 or choose 4 or 8
              ## what about separations? just first time for each event event then 260 instead of 538 
              ## but still gives 33670 combos for pk=2 and generating them is slow
              ## or 186043585 for pk=4 which is still far too many (>45 needs recursion incr)
              ## SO pick NTRIES at random, aware that 10K is still far fewer than combos
              
              exposedTimes = sample(activetimes,pk) #random pk events of h  TODO spread these out - seems to happen naturally
              if (!hexp) {
                exposedHHs = which(colSums(AM[exposedTimes[1:pk],])==pk) #test pk events each time
                if (length(exposedHHs) <= pk) #this HH is exposed for at least some ka
                { 
                  hexp = TRUE #after this, terminate search for exposing events for this HH h
                  #print(c(h,length(activetimes),length(exposedHHs)))
                  for (i in 1:klen) {
                    k = kavalues[i]
                    if (length(exposedHHs)<=k) { #then h is in a kanon group
                      newex = which(exposureInfo[i,exposedHHs]==0) #HHs that were not exposed before
                      if (length(newex)>0) { 
                        exposureInfo[i,exposedHHs[newex]] = pk  #mark all exposed HHs as covered, not just h
                      }
                    }
                  }
                } 
              }
            }
          } 
        }
      }
    }
  }
  
  return(exposureInfo)
}


## make summary of exposures
#activityGT = readRDS(file=paste(RESULTSHIGHRES,sprintf("highresActivityGTAct%d.rds",actID),sep=""))[[aggi]] a matrix
#allresDT <- as.data.table(readRDS(paste(RESULTSHIGHRES,"unicityCulnane.rds")))
#all pkvalues tested 
summaryofExposure <- function(thisactivityGT, actID, aggi, pkvalues, allresDT)
{
  AM = thisactivityGT
  ss = summary(colSums(AM)) #average num of time periods per HH for this end-use 
  N = dim(AM)[2]
  NTRIES = min(allresDT[activityID==actID & ag==aggi,ntries])
  eventsummary = c(round(ss[1]), round(ss[3]), round(ss[6]) ,N,NTRIES)
  kavalues = unique(allresDT$kanon)
  klen = length(kavalues)
  isexposed = as.matrix(allresDT[activityID==actID & ag==aggi,5:(dim(allresDT)[2])]) #table of HH exposures
  pkvalues0 = c(0, pkvalues)
  
  #make summary list, num HHs exposed per pk and kanon
  #pk0 is number of safe HHs then pk1, pk2, pk4, pk8 is number exposed in group <=ka
  #sum of pki = 1600.  pk1 to 8 increases since pki+1 includes pki
  res=c()
  for (ki in 1:klen) {
    pkexp = c()
    for (pk in pkvalues0) {
      pkexp = c(pkexp,sum(isexposed[ki,]==pk))
    }
    res = rbind(res, c( kavalues[ki], 
                        pkexp, #this is cum exposed - to be added up
                        #pkexp[1], cumsum(pkexp[2:length(pkexp)]), this is total exp
                        eventsummary))
  }
  colnames(res)[1:(1+length(pkvalues0))] = c("kanon",paste("pk",pkvalues0,sep=""))
  nc = dim(res)[2]
  colnames(res)[(nc-1):nc] = c("numHHs","numTries")
  res = as.data.table(res)
  return(res)
}



## function to merge per activity outliers into a single set for each param setting
{
  ## make all outliers list by merging time and count outlier households
  ## c(-1) for settings with no outlier HHs
  ## eg mergeUnicityOutliers(lowresunicityexperiments)
  ## eg unicityexperiments as.data.table(saved unicity results)
  
  #unicityexperiments <- as.data.table(readRDS(paste(RESULTSHIGHRES,"unicity1600hh5000triesMACHOMECulnane.rds",sep="")))
  #unicityexperiments <- as.data.table(readRDS(paste(RESULTSLOWRES,"unicityLowRes1000triesCulnane.rds",sep="")))
  #makes a table combining outliers for each pk,ka,agi setting merges all activities
  mergeUnicityOutliers <- function(unicityexperiments)
  {
    mergedoutliers = c()
    nc = dim(unicityexperiments)[2]
    thiskanon = unique(sort(unicityexperiments$kanon))
    thisagghours = unique(sort(unicityexperiments$aggi))
    thispk = unique(sort(c(as.matrix(unicityexperiments[,4:nc]))))
    thisact = unique(sort(unicityexperiments$activityID))
    theseacts = paste(thisact,collapse="") #put them into one name
    
    for (agg in thisagghours) {
      for (ka in thiskanon) {
        for (pk in thispk) {
          if (pk>0) { #ignore safe
            thisAM = unicityexperiments[aggi==agg & kavalues==ka,4:nc]
            thisres = colSums((thisAM>0 & thisAM<=pk)*1)
            thisexposed = sum(thisres>0) #is exposed (!=0 for 0 is safe)
            mergedoutliers = rbind(mergedoutliers, c(agg,ka,pk,theseacts,thisexposed, thisres))
          }
        }
      }
    }
    colnames(mergedoutliers)[1:5] = c("Agg", "Kanon", "pknown", "activity", "NumberExposed")
    mergedoutliers = as.data.table(mergedoutliers)
    return ( mergedoutliers )
  }
}


## END HERE

