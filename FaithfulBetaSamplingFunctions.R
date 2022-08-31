## Create a faithful sample using Beta, k, and mergedoutliers
## Return epsilon,delta achieved
## Rachel Cardell-Oliver
## August 2022

source("./GlobalParameters.R") #kanon choices


## Return RMSE accuracy HHs (in 0..1) for Leaks, Habits, Occupancy

## eg reads "highresActivityGTAct%d.rds" for a given activityID
## lowresActivityGT = readRDS(paste(RESULTSLOWRES,"lowresActivityGT.rds",sep=""))[[agg]][ai,,]
## outlierIDs list of identified outliers for this agg of acts
## eg mergedoutliers = as.data.table(lowresoutliers or highresoutliers)
## uses SAMPLING_REPEATS=50 (set in GlobalParameters) to select the sample with best RMSE
## returns RMSE and the list of HHs that make up the best sample
getSampleAccuracy <- function(activityGT,mergedoutliers,aggi=1,kai=50,beta=0.1) 
{
  
  set.seed(9876) #use a seed to get repeatable results
  activities = activityGT 
  T = dim(activities)[1]
  N = dim(activities)[2]
  
  #baseline time series of number active HHs for different activities from full population 
  allperactivityTS = rowSums(activities)
  
  #time series of number active HHs for different activities from the beta sample
  if (is.null(mergedoutliers)) {
    inliers = 1:N #all HHs
  } else {
    lasthh = dim(mergedoutliers)[2]
    inliers = which(mergedoutliers[agg==aggi & kanon==kai,4:lasthh]==0)
  }
  
  NB = round(N*beta)
  if (NB > length(inliers)) { #check if sampling is possible
    print(sprintf("ERROR: betasize=%d > inliers=%d", NB, length(inliers)))
    return(list("RMSEs"=NA,"sampleHHs"=NULL))
  }
  
  ## else try multiple cases
  res=c()
  for (r in 1:SAMPLING_REPEATS) {
    betacandidates = sample(inliers,NB)
    sampleperactivityTS  = rowSums(activities[,betacandidates])
    
    ## faithfulness results
    RMSE = rmse( allperactivityTS/N, sampleperactivityTS/NB ) 
    res = rbind(res, c(RMSE, betacandidates))
  }
  best = order(res[,1])[1]  #best for all use case
  return(list("RMSE"=res[best,1],"sampleHHs"=res[best,2:dim(res)[2]]))
}

## as above but for extended populations see PopulationSizeHighres

getLargeSampleAccuracy <- function(activityGT,inlierhhs,beta=0.1) 
{
  set.seed(9876) #use a seed to get repeatable results
  activities = activityGT 
  T = dim(activities)[1]
  N = dim(activities)[2]
  
  #baseline time series of number active HHs for different activities from full population 
  allperactivityTS = rowSums(activities)
  
  
  NB = round(N*beta)
  if (NB > length(inlierhhs)) { #check if sampling is possible
    print(sprintf("ERROR: betasize=%d > inliers=%d", NB, length(inlierhhs)))
    return(NA)
  }
  
  ## else try multiple cases
  res=c()
  for (r in 1:SAMPLING_REPEATS) {
    betacandidates = sample(inlierhhs,NB)
    sampleperactivityTS  = rowSums(activities[,betacandidates])
    
    ## faithfulness results
    RMSE = rmse( allperactivityTS/N, sampleperactivityTS/NB ) 
    res = cbind(res, RMSE)
  }
  bestrmse = min(res)  #best rmse from multiple samplings
  return(bestrmse)
}


## END HERE
