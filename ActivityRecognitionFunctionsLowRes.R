## Low res (hourly) activity recognition functions
## RCO Jun 2022

## These functions take an T x N array of total use volumes at some aggregation agg hours (maybe upscaled)
## eg lowresVolumeGT = readRDS("../resultsLowres/lowresVolumeGT.rds") or
## upscale(lowresVolumeGT,4)
## and return a T x N array of 0,1 activity occurring of not

#get constants for low res activity recognition
{
  source("./GlobalParameters.R")
  source('./HabitDetectionAlgorithm.R') #habit detection algorithm 
}


## marks all use as habit in the hours habits occur - not quite right but woudl need to change seqHDA2021 to return actual vols
habitDetection <- function( volumegt,  agg)
{
  #choose habit threshold for each agg
  km=kmeans(c(volumegt),2,nstart=5)
  threshold=round(mean(km$centers)) #threshold
  matchratio=0.5 #how close is pattern of days?
  minsup=7  # at least 7 events to be a pattern
  
  habits = volumegt*0 #dummies to be filled
  habitslist = c()
  #learn habits for each household
  mts = as.numeric(row.names(volumegt))
  for (hh in 1:dim(volumegt)[2]) {
    hhres = seqHDA3(hh,volumegt[,hh],mts,agg,threshold,matchratio,minsup,
                    ploton=FALSE)
    habits[,hh] = (hhres$habitevents)*volumegt[,hh]
    habitslist = rbind(habitslist, hhres$habitdetails)
  }                 
  return(habits) #and if needed more details: list("habits"=habits,"habitdetails"=habitslist))
}

#activity 1 leak detection (based on min flow per day)
#input totals array volumegt and aggregation hours
#outputs a smaller array with activity recognised
#return array with minvol (per day) if leak occurs and 0 otherwise
leakDetection <- function( volumegt, agg)
{
  if (agg==24) {  #daily data - no need to check min of all periods
    return( (volumegt > LEAK_LIM*agg)*1 ) 
  }
  #otherwise learn leak from MIN vol from all periods per day
  aggperday = 24/agg
  leaks = volumegt*0 #dummies to be filled
  timestamp = as.numeric(rownames(volumegt)) #hour labels
  dayagg = floor(timestamp/24) #same ID for all periods belonging to one day to get leak
  dha = unique(dayagg)
  for (i in 1:length(dha)) {
    #if at least one activity hour in the period then mark agg period with activity
    pos = which(dayagg==dha[i])
    if (length(pos)==aggperday) { #only label full aggregated periods with the total volume in that period, should always be true
      minvol = apply(volumegt[pos,],2,min) #find min flow/period in one day
      minvol[which(minvol < (LEAK_LIM*agg))] = 0 #ignore contin flows < leak limit, else leak is set to min flow
      for (j in pos) { leaks[j,] = minvol } #copy mins for all hhs for all periods of day
    }
  }
  #leaks contains the value per hour where there are leaks
  #convert min flow = 2L/hour all day into min flow one day is XL/period 
  #convert L/Hr to L/Period = L/Hr * Hr/Period
  return( leaks )
}

#activity 3 occupancy (based on max flow per day)
#takes vol per period matrix, finds leak vol (min flow in a day) 
# returns per day UNOCCUPIED_MARKER if house is unoccupied
# does not allow for auto irrigation - TODO maybe 20 hours a day 0 is enough
occupancyDetection <- function( volumegt, agg)
{
  #perform aggregation from hourly data
  if (agg==24) {
    return( (volumegt <= EMPTY_LIM)*UNOCCUPIED_MARKER ) #no need to find max
  }
  #otherwise learn max vol in readings per day
  aggperday = 24/agg
  occupancy = volumegt*0 #dummies to be filled
  timestamp = as.numeric(rownames(volumegt))
  dayagg = floor(timestamp/24) #same ID for all periods belonging to one day
  dha = unique(dayagg)
  for (i in 1:length(dha)) {
    pos = which(dayagg==dha[i])
    if (length(pos)==aggperday) { 
      maxvol = apply(volumegt[pos,],2,max) 
      for (j in pos) { #copy max vol per day to all periods of day 
        occupancy[j,] = (maxvol<=EMPTY_LIM)*UNOCCUPIED_MARKER
      }
    }
  }
  #occupied if vol>EMPTY_LIM*agg then occupied if vol<=EMPTY_LIM*agg then assume unoccupied
  return( occupancy )
}

