### calculate_variance_explained
### Author: Erik Ländström
### Date: 190124


# Description -------------------------------------------------------------

# Calculates the variance, variance explained and cumulative variance explained
# for each PC and returns it in a tibble.


# Arguments ---------------------------------------------------------------

# tb = pca results


# Function ----------------------------------------------------------------
calculate_variance_explained <- function(tb) {
  
  data <- tb %>% 
    unnest(pca_aug) %>% 
    summarise_at(.vars = vars(contains("fittedPC")), .funs = funs(var)) %>% 
    gather(key = pc, value = variance) %>% 
    mutate(
      Variance                        = variance,
      `Variance Explained`            = Variance / sum(Variance),
      `Cumulative Variance Explained` = cumsum(`Variance Explained`),
      PC                              = str_replace(pc, ".fittedPC", "") 
    ) %>% 
    select(-pc, -variance) %>% 
    select(PC, everything())
  
  return(data)
}