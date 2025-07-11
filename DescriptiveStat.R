library(tidyverse)
library(lubridate)

# Specify the file path
file_1 <- "./data_1.rds"
file_2 <- "./data_2.rds"

scalar_covariate_1 <- "./scalar_covariates_1.rds"
scalar_covariate_2 <- "./scalar_covariates_2.rds"

# Read the .RDS file
data_1 <- readRDS(file_1)
data_2 <- readRDS(file_2)

scalar_covariate_1 <- readRDS(scalar_covariate_1)
scalar_covariate_2 <- readRDS(scalar_covariate_2)

extract_min <- function(data) {
  res <- c(time = data[[1]]$arg[which.min(data[[1]]$value)], 
           value = min(data[[1]]$value))
  res
}

extract_min(data_1[1, 2])

which.min(data_1[1, 2][[1]]$value)

min_summary <- data.frame(time = NA, value = NA)
for (i in seq_len(nrow(data_1))) {
  min_summary[i, ] <- extract_min(data_1[i, 2])
}
head(min_summary)

plot(data_1[4, 2][[1]]$arg, data_1[4, 2][[1]]$value)

boxplot(min_summary$time)
boxplot(min_summary$value)



