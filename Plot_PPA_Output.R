library(tidyverse)


# Directories -------------------------------------------------------------

base_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
# base_directory <- getwd() # Ensure that your base directory is the PPA folder

species_file = "/input/PPA_FG5_filtered.csv"
initComm_file = "/input/PPA_initial_state_fg5_secondary.csv"
output_directory = "/output/"

ppa <- readRDS(file = paste0(base_directory, output_directory, "PPA_output.rds"))



# Functions ---------------------------------------------------------------





# Quick plot --------------------------------------------------------------

ggplot(ppa) +
    geom_line(aes(x = Year,
                  y = BasalArea,
                  color = as.factor(SpeciesID))) +
    labs(color = "SpeciesID") +
    theme_minimal()
