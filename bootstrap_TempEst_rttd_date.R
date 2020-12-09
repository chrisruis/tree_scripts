#Calculates the significance of a correlation between root-to-tip distance and collection date using bootstrapping
#Takes a file produced using "Export Data" in TempEst
#To run: RScript bootstrap_TempEst_rttd_date.R distance_date.txt number_bootstraps

args <- commandArgs(TRUE)

#Import the root to tip distance versus sampling date from TempEst
distanceDate <- read.table(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = '!')

#Number bootstraps that will be carried out
bootstraps <- as.numeric(args[2])

#Calculate the correlation between sampling date and root-to-tip distance in the real data
distanceDatelm <- lm(distance ~ date, data = distanceDate)
correlation <- summary(distanceDatelm)$r.squared

#Number of samples in the dataset
numberSamples <- length(distanceDate[,1])

#Used as the base data frame for each bootstrap
bootstrapDistance <- data.frame("date" = 0, "distance" = distanceDate$distance)

#Will be filled with the correlations in the bootstraps
bootstrapCorrelation <- c()

#Iterate through the bootstraps, randomise dates, calculate correlation and add to bootstrapCorrelation
for (bootstrap in 1:bootstraps) {
  
  bootstrapDistanceDate <- bootstrapDistance
	bootstrapDistanceDate[,1] <- sample(distanceDate$date, numberSamples)
	
	bootstraplm <- lm(distance ~ date, data = bootstrapDistanceDate)
	bootstrapCorrelation <- c(bootstrapCorrelation,summary(bootstraplm)$r.squared)}

#Calculate the proportion of bootstraps with R2 at least as high as with the real dates
p.value <- length(bootstrapCorrelation[which(bootstrapCorrelation >= correlation)])/bootstraps

print(paste("Correlation between sampling data and root-to-tip distance with real dates: ", correlation, sep = ""))
print(paste("Proportion of bootstraps with equal or greater correlation (p-value): ", p.value, sep = ""))