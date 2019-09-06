library(magrittr)

# Load the data: 
hists <- read.csv("History.csv", header=T)
indivs <- read.csv("Individual.csv", header=T)

# Turn the Dead column into a factor with two levels: 
hists$Dead <- hists$Dead %>% factor

# Convert the dates to Date objects:
hists$capture_date %<>% as.character %>% as.Date(optional = T)
hists$processing_date %<>% as.character %>% as.Date(optional = T)
hists$release_date %<>% as.character %>% as.Date(optional = T)

# Look at the different locations:
hists$location %>% unique %>% length

# The number of individual fish in the data: (88 112)
indivN <- hists$individual_id %>% unique %>% length 

# The total number of observations (captures): (406 907)
captures <- nrow(hists)

# The number of different capture, processing and release dates:
captureT <- hists$capture_date %>% unique %>% length
processT <- hists$processing_date %>% unique %>% length
releaseT <- hists$release_date %>% unique %>% length


max(hists$processing_date)
