# We know that slice 35 is the least different and 
# slice 57 is the most different. We want confirm that 
# using the bottleneck distance.

source('distance.r')

# Load the 64 dataset.
cdm <- readRDS(paste(source, "cdm", 4, ".rds", sep=""))
wdm <- readRDS(paste(source, "wdm", 4, ".rds", sep=""))

cdm_low <- cdm[[35]]
wdm_low <- wdm[[35]]

cdm_hi <- cdm[[57]]
wdm_hi <- wdm[[57]]

cont_low <- contourDist(cdm_low, wdm_low)
cont_hi <- contourDist(cdm_hi, wdm_hi)

print(paste("Contour Min Dist:", cont_low))
print(paste("Contour Max Dist:", cont_hi))

bot_low <- bottleneckDist(cdm_low, wdm_low)
bot_hi <- bottleneckDist(cdm_hi, wdm_hi)

print(paste("Bottleneck Min Dist:", bot_low))
print(paste("Bottleneck Max Dist:", bot_hi))