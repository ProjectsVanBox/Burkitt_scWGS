########################################
# Define function to track plot number #
########################################

plotCount <- -1
plotCounter <- function(x){
    plotCount <<- sprintf("%02d", as.numeric(x) + 1)
    print(plotCount)
}
