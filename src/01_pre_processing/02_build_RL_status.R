library(taxize)
library(tidyverse)

# Helper functions --------------------------------------------------------

get_rl_status <- function(x) {
  message(paste("Searching ", x, "..."))
  res <- taxize::iucn_summary(x)
  rl_status <- res[[1]]$status
  return(list(species = x, status = rl_status))
}

spp_data <- readr::read_csv("data/external/udr/european_tetrapods/spp_codes.csv")

st <- purrr::map_df(spp_data$species, get_rl_status)

readr::write_csv(st, "data/external/udr/european_tetrapods/spp_rl_statuses.csv")
