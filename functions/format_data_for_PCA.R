### format_data_for_PCA
### Author: Erik Ländström
### Date: 190124


# Description -------------------------------------------------------------

# Formats data produced by DEP function get_df_wide so it's ready for 
# doing PCA.


# Arguments ---------------------------------------------------------------

# tb = imputed data generated with get_df_wide
# num_samples = numeric, number of samples
# variables = character vector, use c(), vector with group names

# Function ----------------------------------------------------------------

format_data_for_PCA <- function(tb, num_samples, variables) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # Read gene name and intensity data from get_df_wide(data_se)
  temp <- tb[, 1:(1 + num_samples)]
  
  # Create a transposed tibble
  temp <- as_tibble(
    cbind(nms = names(temp), t(temp))
  )
  
  # Use gene names as column names
  colnames(temp) <- temp[1, ]
  
  # Remove name row
  temp <- temp[-1, ]
  
  # Make all intensity columns numeric
  temp[, -1] <- temp[, -1] %>% 
    mutate_all(.funs = as.numeric)
  
  # Separate name column into each variable and join all except sample number
  temp <- temp %>% 
    separate(name, into = variables, sep = "_") %>% 
    unite(col = "group", paste(variables[-1]), remove = FALSE)
  
  # Perform PCA
  temp <- temp %>% 
    nest() %>% # Nest data
    mutate(
      pca = map(data, ~ prcomp(.x %>% dplyr::select(-sample:-litter),
                               scale = TRUE, center = TRUE)), # Perform PCA
      pca_aug = map2(pca, data, ~ augment(.x, data = .y)) # Extract PCA
    )
  
  return(temp)
}