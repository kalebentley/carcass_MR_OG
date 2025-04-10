# Function to collapse third dimension
collapse_third_dimension <- function(data_list) {
  # Initialize an empty list to store long format data frames
  long_dfs <- list()
  
  # Loop over each element in the list
  for (level_name in names(data_list)) {
    long_df <- data_list[[level_name]] %>%
      mutate(type = level_name) %>%
      pivot_longer(-c(dates, type), names_to = 'variable', values_to = 'value') %>%
      mutate(value = as.numeric(as.character(value)))
    long_dfs[[level_name]] <- long_df
  }
  
  # Combine all long data frames
  combined_long <- bind_rows(long_dfs)
  
  # Summarize to collapse the levels
  collapsed_df <- combined_long %>%
    group_by(dates, variable) %>%
    summarize(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = variable, values_from = value) |> 
    mutate(dates = as.Date(dates))
  
  return(collapsed_df)
}