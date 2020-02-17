rm(list = ls())
library(dplyr)
library(tuneR)
library(seewave)
set.seed(1830257)
# 07-psdValidation --------------------------------------------------------

# Code just to get file names in the csv and then exported for manual analysis
# minuteSongs <- list.files("../data/recordings/")
# # Remove the 30th minute and 31st minute for all songs due to microphone error
# char1 <- 30
# char2 <- 31
# char3 <- 29
# minuteSongs <- minuteSongs[!(substr(minuteSongs, char1, char2) %in% c("30", "31"))]
# # Accounts for recording sites with only two letter abbreviations
# minuteSongs <- minuteSongs[!(substr(minuteSongs, char3, char1) %in% c("30", "31"))]
# dat <- data.frame(fileNames = minuteSongs)
# write.csv(dat, "../data/psdValidation.csv", row.names = FALSE)

dat <- read.csv("../data/psdValidation.csv", stringsAsFactors = FALSE)
dat$largeFileName <- substr(dat$fileNames, 1, nchar(dat$fileNames) - 6)

n.full.rec <- length(unique(dat$largeFileName))
n <- 2
recordings <- rep(0, n * n.full.rec)
for (i in 1:n.full.rec) {
  recordings[((i-1) * n + 1):(i * n)] <- sample(dat[dat$largeFileName == unique(dat$largeFileName)[i], "fileNames"], 2)
}

write.csv(recordings, "../data/psdSubsetVal.csv", row.names = FALSE)



# Figure of a spectrogram -------------------------------------------------
rec <- readWave("../data/recordings/MINUTE-IFP-C_20160623_125457019.wav")
png("../figures/spectrogram.png", width = 720, height = 480)
spectro(rec, colbg = "white", flim = c(0.5, 11), palette = spectro.colors, grid = FALSE,
        cexlab = 1, scalecexlab = 1)
abline(h = 2, lty = 2)
text(x = 3, y = 5, "Birds", cex = 1.5)
text(x = 5.5, y = 1.7, "Road Noise", cex = 1.5)
dev.off()
