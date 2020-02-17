
# Biophony and Technophony Computation ------------------------------------

library(seewave)
library(tuneR)
library(soundecology)
library(stringr)
library(tidyr)

setwd("../data/")

# Preparation of Sound Files ----------------------------------------------

# # Get names of all wav files
# shortNames <- list.files(paste(getwd(), "fullRecordings", sep = "/"),
#                          pattern = "*.wav")
# 
# # Split sound files into 1 minute increments ------------------------------
# 
# setwd(paste(getwd(), "fullRecordings", sep = "/"))
# for(i in seq_along(shortNames)) {
#   system(paste("sox", shortNames[i], paste("MINUTE", shortNames[i], sep = "-"),
#                "trim 0 60 : newfile : restart", sep = " "))
# }


# Prepare 1 minute files for computation ----------------------------------

minuteSongs <- list.files(paste(getwd(), "recordings", sep = "/"), pattern = "MINUTE*")
# Remove the 30th minute and 31st minute for all songs due to microphone error
char1 <- 30
char2 <- 31
char3 <- 29
minuteSongs <- minuteSongs[!(substr(minuteSongs, char1, char2) %in% c("30", "31"))]
# Accounts for recording sites with only two letter abbreviations
minuteSongs <- minuteSongs[!(substr(minuteSongs, char3, char1) %in% c("30", "31"))]


# Compute NDSI of each 1 minute file --------------------------------------

# NOTE: cannot use multiple_sounds() function from soundecology because you need
# to extract the separate psd values, which that function does not provide.

ndsiData <- rep(list(list(ndsi_left = NA, ndsi_right = NA, biophony_left = NA,
                     anthrophony_left = NA, biophony_right = NA,
                     anthrophony_right = NA)), length(minuteSongs))
setwd("../data/recordings/")

# With 0.5 as a lower bound
# for (i in seq_along(minuteSongs)) {
#   songFile <- readWave(minuteSongs[i])
#   ndsiData[[i]] <- ndsi(songFile, anthro_min = 500, anthro_max = 2000, bio_max = 11000)
#   rm(songFile)
# }

# Test with 0 as a lower bound
for (i in seq_along(minuteSongs)) {
  songFile <- readWave(minuteSongs[i])
  ndsiData[[i]] <- ndsi(songFile, anthro_min = 0, anthro_max = 2000, bio_max = 12000)
  rm(songFile)
}


# Format the data and export as csv ---------------------------------------

# Reformat the list of lists into a matrix
ndsiMatrix <- matrix(unlist(ndsiData), ncol = 6, byrow = TRUE)
# Convert to data frame
ndsiDataFrame <- as.data.frame(ndsiMatrix)
names(ndsiDataFrame) <- c("ndsiLeft", "ndsiRight", "biophonyLeft",
                          "anthrophonyLeft", "biophonyRight", "anthrophonyRight")

# Need to do some more manipulation of the data frame here to get location

# Manipulate file name to obtain the time of day
time <- minuteSongs %>%
  str_replace(".*_", "") %>%
  str_replace(".wav", "")
time <- as.integer(time)
timeOfDay <- character(length(time))
timeOfDay[time < 100000000] = "Morning"
timeOfDay[(time > 100000000) & time < 180000000] = "Afternoon"
timeOfDay[(time >= 180000000)] = "Evening"
# Add time of day to data frame
ndsiDataFrame$time <- timeOfDay

# Obtain only location from the recording name
location <- minuteSongs %>%
  str_replace("MINUTE-", "")
ndsiDataFrame$location <- gsub("_.*$", "", location)

# # Export data frame if desired for full NDSI values
# write.csv(ndsiDataFrame, file = "ndsiData.csv")
write.csv(ndsiDataFrame, file = 'ndsiData02.csv')

# After doing above, just read in the data frame
ndsiDataFrame <- read.csv("ndsiData.csv")


# Some further manipulation -----------------------------------------------

# Number of recordings
K <- 29
# Number of time periods
J <- 3
# Number of recording sites
l <- 9 * 2
# Total recordings
numRecordings <- l * J
# Total observations
N <- K * J * l

# Take the average of the left channel and right channel
bio <- rowMeans(ndsiDataFrame[, c("biophonyLeft", "biophonyRight")])
anthro <- rowMeans(ndsiDataFrame[, c("anthrophonyLeft", "anthrophonyRight")])

# Only get the name of the forest plot
forestPlot <- gsub("-.*$", "", ndsiDataFrame$location)

# Create column indicating the part of each recording
recordingNumber <- rep(1:K, times = numRecordings)

# Create new data frame with desired values
# Redclare time as a factor to get proper values for M, A, E
ndsiDat <- data.frame(forestPlot, ndsiDataFrame$location,
                      factor(ndsiDataFrame$time, levels = c("Morning", 
                                                            "Afternoon", 
                                                            "Evening")), bio,
                      anthro, recordingNumber)
names(ndsiDat) <- c("forestPlot", "location", "time", "bio", "anthro",
                    "recordingNumber")
head(ndsiDat)

# Create variable containing only "Center" or "Edge"
ndsiDat$location <- as.character(ndsiDat$location)
for (i in seq_along(ndsiDat$location)) {
  if (substr(ndsiDat$location[i], nchar(ndsiDat$location[i]), nchar(ndsiDat$location[i])) == 'C') {
    ndsiDat$locationInd[i] = "Center"
  } else {
    ndsiDat$locationInd[i] = "Edge"
  }
}
ndsiDat$location <- factor(ndsiDat$location)


# Format the data by number, site, time of day ----------------------------

stackedDat <- ndsiDat[with(ndsiDat, order(time, forestPlot, locationInd)), ]

# Since the PSD bands here are 1.5 and the range of bio is from 2 - 11, scale
# the bio values to be theoretically b/w 0 and 1 by dividing by the total number
# of bands, which is 6
bio <- stackedDat$bio / 6
techno <- stackedDat$anthro

# write.table(bio, "../../beta/data/beta", row.names = FALSE, col.names = FALSE, 
#             sep = "\t")
# write.table(techno, "../../beta/data/alpha", row.names = FALSE, col.names = FALSE,
#             sep = "\t")
