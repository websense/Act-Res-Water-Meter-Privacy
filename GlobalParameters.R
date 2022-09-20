## GlobalParameters
## Colour palettes and labels for figures for the paper for hourly and 10 second data
## Author: Rachel Cardell-Oliver
## Version: Sep 2021, Jul 2022

# libraries and common directories
{
require(data.table)
require(RColorBrewer)
require(mltools) #rmse
require(xtable)
require(Matrix) #sparse matrices
  
FIGURESDIR="../../figures/"

runalltests = FALSE #set this to true in local context if needed
}
  

# colour palettes and names for figures for the paper
## Hourly activities
{
aggcol = brewer.pal(8,"Greens")[8:1] #need seq colours, not finishing too light
lowresaggselection = c(1,2,3,4,12,24)
lowresaggnames =  paste(lowresaggselection,"hr",sep="")
  
lowresactivitycols = c("royalblue","forestgreen","orange","darkgray")
lowresactivitynames = c("Leaks","Habits","Absence","Any Use") 
lowresactivitypch = c(6,5,1) #2,

## constants for lowres activity recognition
LEAK_LIM = 2 #>=2 L/h for last 24 hours defines a leak
EMPTY_LIM = 0 #L/h for last 24 hours defines an unoccupied day
UNOCCUPIED_MARKER = -20 #indicator of activity else 0

}

## highres-grained 10 second activities for STREaM data
## Taps (and Any Use) are most frequent
{
    highresactivitycols = brewer.pal(7,"Paired")
    highresactivitynames = c("Toilet","Shower","Taps",
                             "Clotheswasher","Dishwasher","Bathtub","Any Use") 
    highresactivitypch = 1:7

    # in multiples of 10 seconds (so last is 1 min)
    highresaggselection = c(1,3,6,6*5,6*15,6*60) #drop 30 mins
    
    highresaggseconds = highresaggselection * 10 #convert to seconds
    
    oneday10s = 24*60*6
    onehour10s = 60*6
    
    #labels for graphs
    highresaggnames = paste(as.character(round(highresaggselection*10/60)),"min")
    highresaggnames[1:2] = c("10 sec","30 sec") #10 second labels
}


## Aggregations
## uses agghours = c(1,2,3,4,6,8,12,24) from ColoursandNames



## Privacy control inputs
{
## more settings for selecting best privatisation
MCCTHRESHOLD = 0.6 #lowest bound for acceptable activity recognition accuracy (set to 0 to try everything)
  
## OutlierAnalysis.R
kanon = c(1,5,10,20,30,50,100,300) #min size k to test for anonymity

## for makeUnicityOutlierResults.R
#MINPK = 4 #dataset must resist attacks using external knowledge of 4 or more known events usually 6 max
#MAXPK = 6 #max number of prior knowledge events that could be known for one household
#pks=MINPK:MAXPK  #these are the ones we test

#dataset must resist attacks using external knowledge of pkvalues known events 
#these are the ones we test - because lowres plateaus after pk=4
#unicity tests need to be run 10K times so very slow
pkvalues = c(1,2,4,8)  

# for makeSafeKanonTables.R
betas = c(0.05,0.1,0.25,0.5) #,0.999) #require 0<beta<1
## these eps,delta upper bounds are much too high, 
# but are included to allow metrics for all cases (incl full raw data)
epsmax = 2.0 #upper bound for eps, #or use epsmax=12 to test 0.999 case
deltamax = 10^-4 #was 10^-6 but for high res need more search #upper bound for delta #or use 1.0 is to test 0.999 case
##find eps (nat log = ln)
#eps >= (-log(1-beta))
#delta=getdelta(k,b,eps) see below

SAMPLING_REPEATS = 50 #tests when selecting beta sample
}

## Pareto front functions
{
# v1 dominates v2 if every entry of v1 <= v2 that is v2-v1>=0 AND 
# at least one v1 < v2  that is v1 != v2
# here SMALLEST values are best
dominates <- function(v1,v2) {
  return ((min(v2-v1) >= 0) & (max(v2-v1) > 0))
}

incomparable <- function(v1,v2) {
  return (!dominates(v1,v2) & !dominates(v2,v1))
}
}