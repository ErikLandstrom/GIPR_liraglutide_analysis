### plot_var_exp
### Author: Erik Ländström
### Date: 190124


# Description -------------------------------------------------------------

# Plots variance explained by the principal components.


# Arguments ---------------------------------------------------------------

# tb = pca results


# Function ----------------------------------------------------------------

plot_var_exp <- function(tb) {
  
  library(tidyverse)
  library(ggplot2)
  
  tb %>% 
    mutate(PC = as.numeric(PC)) %>% 
    gather(key = key, value = value, `Variance Explained`:`Cumulative Variance Explained`) %>% 
    ggplot(aes(PC, value, group = key)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ key, scales = "free_y") +
    lims(y = c(0, 1)) +
    labs(
      y     = "Proportion of variance",
      title = "Proportion of variance explained by each principal component"
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(color = "black")
    )
}