# Load necessary library
library(dplyr)

# Specify the folder where your CSV files are located
folder_path <- "output"


# List all .csv files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty data frame for your summary table
summary_table <- data.frame(Year = character(), Model = character(), mean = numeric(), sd = numeric(),
                            `2.50%` = numeric(), `25%` = numeric(), `50%` = numeric(), `75%` = numeric(),
                            `97.50%` = numeric(), Rhat = numeric(), n.eff = numeric(), stringsAsFactors = FALSE)

# Iterate over each file
for (file in file_list) {
  
  # Extract year and model name from the file name
  file_name <- basename(file)
  year <- gsub(".*[^0-9](\\d{4})[^0-9].*", "\\1", file_name) # Extracts a 4-digit number
  model <- sub(".*_(sst|ttt).*", "\\1", file_name)
  
  # Read the file
  data <- read.csv(file, header = TRUE)
  
  # Ensure column names are standardized
  colnames(data) <- make.names(colnames(data))
  
  # Find the rows where the first column is equal to "Nsuper"
  rows_with_Nsuper <- which(data[[1]] == "Nsuper")
  
  # If there's a matching row, extract the specific columns for that row
  if (length(rows_with_Nsuper) > 0) {
    for (row in rows_with_Nsuper) {
      summary_table <- rbind(summary_table, data.frame(
        Year = year,
        Model = model,
        mean = if ("mean" %in% colnames(data)) data[row, "mean"] else NA,
        sd = if ("sd" %in% colnames(data)) data[row, "sd"] else NA,
        `2.50%` = if ("X2.50." %in% colnames(data)) data[row, "X2.50."] else NA,
        `25%` = if ("X25." %in% colnames(data)) data[row, "X25."] else NA,
        `50%` = if ("X50." %in% colnames(data)) data[row, "X50."] else NA,
        `75%` = if ("X75." %in% colnames(data)) data[row, "X75."] else NA,
        `97.50%` = if ("X97.50." %in% colnames(data)) data[row, "X97.50."] else NA,
        Rhat = if ("Rhat" %in% colnames(data)) data[row, "Rhat"] else NA,
        n.eff = if ("n.eff" %in% colnames(data)) data[row, "n.eff"] else NA
      ))
    }
  }
}

# Print the summary table
summary_table |> arrange(Model, Year)

# Comparison of Nsuper estimate by model and year
summary_table |> pivot_wider(id_cols = c(Year), names_from = Model, values_from = 'X50.')



