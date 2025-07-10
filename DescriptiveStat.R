library(tidyverse)
library(lubridate)

# Specify the file path
file_1 <- "G:/Consulting/data/data_1.rds"
file_2 <- "G:/Consulting/data/data_2.rds"

scalar_covariate_1 <- "G:/Consulting/data/scalar_covariates_1.rds"
scalar_covariate_2 <- "G:/Consulting/data/scalar_covariates_2.rds"

# # save as RDS file
# saveRDS(data_1[1, ], file = "C:/Users/Xu Yuan/Documents/file_3")

# Read the .RDS file
data_1 <- readRDS(file_1)
data_2 <- readRDS(file_2)

scalar_covariate_1 <- readRDS(scalar_covariate_1)
scalar_covariate_2 <- readRDS(scalar_covariate_2)

# Now you can use the data object in your R session
# Function to extract time-series data from a list
extract_timeseries <- function(data_list) {
  time_vector <- data_list[1]      # First vector (time)
  measurement_vector <- data_list[2]  # Second vector (measurements)
  
  # Create a data frame
  data.frame(
    time = time_vector,
    measurement = measurement_vector
  )
}

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



# extract_timeseries(data_1[1, 2][[1]])
# 
# extract_timeseries(data_1[, 2][[1]])
# 
# lapply(data_1[, 2], extract_timeseries)
# 
# 
# map_dfr(data_1[, 2], extract_timeseries)
# 
# 
# data_1[, 2]
# 
# apply(data_1[, 2], 1, extract_timeseries)
# 
# data_1[, 2][[1]]
# 
# for (i in seq_len(nrow(data_1))) {
#   
# }
# 
# 
# data_1[1, 2][[1]]
# 
# # Extract all time series from both columns
# column1_data <- map_dfr(your_dataset$column1, extract_timeseries, .id = "series_id")
# column2_data <- map_dfr(your_dataset$column2, extract_timeseries, .id = "series_id")
# 
# # Add column identifiers
# column1_data$column <- "column1"
# column2_data$column <- "column2"
# 
# # Combine if needed
# all_data <- bind_rows(column1_data, column2_data)
# 
# 
# 
# data_1[data_1$id == 1, ] %>% pivot_longer(cols = 2, 
#                                             names_to = "measurement_type", 
#                                             values_to = "measurement_data") %>%
#   ggplot(aes(x = measurement_type, y = measurement_data)) + geom_line()
# 
# 
# 
# 
# 
# 
# 
# summary(data_1)
# str(data_1)
# View(data_1)
# 
# data_1[1:3, 2]
# 
# ggplot()
# 
# str(data_1[, 4])
# # Change character to factor
# data_1_df[, 4] <- as.factor(data_1_df[, 4])
# summary(data_1_df[, 4])
# 
# 
# str(data_1_df[, 2])
# length(data_1_df[, 2])
# 
# str(scalar_covariate_1)
# summary(scalar_covariate_1)
# str(scalar_covariate_2)
# summary(scalar_covariate_2)
# which.min(scalar_covariate_2[, 2])
# scalar_covariate_2[28, ]
# scalar_covariate_1[, ]
# 
# # subset data_1_df to only include rows where the id is 2
# data_1_df_subset <- data_1_df[data_1_df$id == 2, ]
# View(data_1_df_subset)
# 
# data_1_id <- unique(data_1[, 1])
# 
# # Create a new data frame with the id and the number of measurements for each id
# num_measurements <- data.frame(
#   id = data_1_id,
#   num_measurements = 
#     sapply(data_1_id, function(x) nrow(data_1_df[data_1_df$id == x, ]))
# )
# # View the new data frame
# num_measurements
# unique(data_1[data_1$id == 64, 4])
# 
# num_measurements[which.max(num_measurements[,2]), ]
# 
# unique(data_1[data_1$id == 16, 4])
# 
# 
# data_1[1, 2]
# data_2[1, 3]
# 
# 
# unique(data_2[, 1])
# unique(data_1[, 1])
# 
# 
# 
# 
# 
# # Load required libraries
# library(ggplot2)
# 
# # Create sample P-V loop data (cardiac cycle example)
# # This represents one cardiac cycle with phases:
# # 1. Isovolumic contraction
# # 2. Ejection
# # 3. Isovolumic relaxation  
# # 4. Filling
# 
# volume <- c(50, 150, 120, 120, 50)  # mL
# pressure <- c(50, 80, 80, 5, 5)     # mmHg
# 
# 
# t <- data_1[1, 2][[1]]$`1`$arg
# sv <- data_1[1, 2][[1]]$`1`$value
# ecg <- c(rep(0, 20), data_1[1, 3][[1]]$`1`$value)
# 
# 
# 
# # Create data frame
# pv_data <- data.frame(
#   Volume = volume,
#   Pressure = pressure,
#   Phase = c("Start", "Isovolumic", "Ejection", "Relaxation", "Filling")
# )
# 
# 
# contraction_data <- data.frame(
#   time = t,
#   sv = sv
#   )
# 
# # Basic plot
# ggplot(se_data, aes(x = sv, y = time)) +
#   geom_path(size = 1.2, color = "blue") +
#   geom_point(size = 3, color = "red") +
#   labs(title = "Pressure-Volume Loop",
#        x = "signal (ev)",
#        y = "Contraction (cm)") +
#   theme_minimal()
# 
# 
# 
