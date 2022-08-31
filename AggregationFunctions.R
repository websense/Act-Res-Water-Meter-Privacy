## Aggregation Function makeActivityGTlist
## create activity GT list using {high/low}resActivityGT for different aggselections
## works for both low and high res inputs
## return list of A x TA x N matrices, one for each aggregation in aggselection


source("./GlobalParameters.R")


## Aggregation function
#Input: highres or lowres activityMatrixArray T x N array of 0,1 (is active or not) for a particular activity type
#Input: agglist = aggselection, temporal aggregations (in 10s intervals) to be calculated
#aggselection = 1   6  30  90 180 360 in 10 second units
#aggseconds = aggselection*10
#Output: list of arrays Tagg x N per agg of vol upscaled GT 0,1 (is active or not)
#Output is a list because each agg matrix has a diff size
# eg lowresActivityGT = makeActivityGTlist(highresMatrixArrayAct1,highresaggselection,highresactivitynames[1])

makeActivityGTlist <- function(activityMatrixArray, aggselection, thisactivityname)
{
  activityGT = list() #list per agg of activityVolumeArrays to return
  #activitynames1 = names(activityVolumeArray[,1,1]) #double check
  #A = length(activitynames) 
  T = dim(activityMatrixArray)[1]
  N = dim(activityMatrixArray)[2]
  gtts = 1:T #basetimeIDs #IDs 1..T of base period data
  for (ag in aggselection) {
    print(sprintf("Aggregation %d %s",ag,thisactivityname)) #show progress because slow
    if (ag==1) { #take baseline 10 second AR or 1 hour data as ground truth
      activityGT[[ag]] = activityMatrixArray
      #(activityVolumeArray>0 | activityVolumeArray==UNOCCUPIED_MARKER)*1 #activity is occurring is vol>0
    } else {
      # upscale for periods from GT
      pbounds = c(1,which(gtts %% ag == 0)+1) #index of period bounds in GT
      nperiods = length(pbounds)-1 #remove last overstep
      #make raw matrix for agg GT
      agt = Matrix(0,nrow=nperiods,ncol=N,sparse=TRUE) 
      for (i in 1:nperiods) {
        oneperiod = (activityMatrixArray[pbounds[i]:(pbounds[i+1]-1),]!=0)*1 #all to 0,1
        # activity vol at agg scale is 1 if **at least one** water use event occurs within agg period
        # event can be +ve most activities or -ve occupancy
        agt[i,] = (colSums(oneperiod)!=0)*1
      }
      activityGT[[ag]] = agt #0,1 is activity
    }
  }
  return(activityGT)
}

## END HERE

