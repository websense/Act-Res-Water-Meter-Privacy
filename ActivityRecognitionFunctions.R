## ACTIVITY RECOGNITION and ACCURACY
## Compare ground truth with the activity detected for all aggregations
## Calculate accuracy (matthews correlation coefficient) for activity recognition 
## Rachel CO 2022

source("./GlobalParameters.R")

# upscale volumes for a given ground truth array (T X N) and aggregate agg
# returns aggregated array (T1 x N) where T1 = T/agg
upscale <- function(volumeGT,agg)
{
  if (agg==1) { return ( volumeGT ) }
  #otherwise aggregate all vols by agg period
  N = dim(volumeGT)[2]
  gtts = as.numeric(rownames(volumeGT))
  pbounds = c(1,which(gtts %% agg == 0)+1) #index of period bounds in GT
  nperiods = length(pbounds)-1
  #make raw matrix for agg GT
  agt = array(0,dim=c(nperiods,N),
              dimnames=list(gtts[pbounds[1:nperiods]],paste("H",1:N,sep=""))) 
  for (i in 1:nperiods) {
    oneperiod = volumeGT[pbounds[i]:(pbounds[i+1]-1),]
    # activity vol at agg scale isum all act in this period
    agt[i,] = colSums(oneperiod)
  }
  return (agt)
}


## Function to calculate accuracy (several metrics including mcc)
#Defns https://en.wikipedia.org/wiki/Precision_and_recall
# input mc = c(TP,FP,TN,FN)
accuracyMetrics <- function(mc) 
{
  TP = mc[1]
  FP = mc[2]
  TN = mc[3]
  FN = mc[4]
  precision = TP/(TP+FP)
  recall  = TP/(TP+FN) 
  TPR = TP/(TP+FN) #true positive rate is recall
  TNR = TN/(TN+FP) #true negative rate
  balancedaccuracy = (TPR+TNR)/2
  mcc = mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  fmeasure = (2*precision*recall)/(precision+recall)
  res = c(TP,FP,FN,TN,precision,recall,balancedaccuracy,fmeasure,mcc)
  names(res)=c("TP","FP","FN","TN","precision","recall",
               "balancedaccuracy","fmeasure","MCC")
  return(res)
}


#compare ground truth and prediction
#both are 0,1 (boolean) vectors for corresponding set of days and meters
aggregatedActivityAccuracy <- function(gt,ar)
{
  TP = sum(gt & ar)
  FP = sum(!gt & ar)
  TN = sum(!gt & !ar) 
  FN = sum(gt & !ar)
  
  res = accuracyMetrics(c(TP,FP,TN,FN))
  return(res)
}


## END activity recognition accuracy functions

  