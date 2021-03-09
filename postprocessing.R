# TODO Attach species PCA traits to the processed data.frames

# Libraries ---------------------------------------------------------------

library(tidyverse)


# Directories and data sets -----------------------------------------------

base_directory   <- dirname(rstudioapi::getActiveDocumentContext()$path)
output_directory <- paste0(base_directory, "/output/")
input_directory  <- paste0(base_directory, "/input/")
species_file     <- paste0(base_directory, "/input/PPA_FG5_filtered.csv")

cohorts          <- readRDS(paste0(output_directory, "PPA_output_raw_cohort.rds"))
mortality        <- readRDS(paste0(output_directory, "PPA_output_raw_cohort_mortality.rds"))
spVitals         <- read.table(species_file, sep = ",", header = TRUE)


# Generate size-class data set --------------------------------------------

size_class_designations <- c(0, 5, 20, 60, 300)

size_classes <- cohorts %>%
    mutate(SizeClass = cut(Diameter, breaks = size_class_designations)) %>%
    group_by(Model, Year, SpeciesID, SizeClass) %>%
    summarise(Diameter = sum(Diameter),
              BasalArea = sum(BasalArea),
              Biomass = sum(Biomass))

saveRDS(size_classes, file = paste0(output_directory, "/PPA_output_processed_size_classes.rds"))


# Generate species-level data set -----------------------------------------

species <- cohorts %>%
    group_by(Year, SpeciesID) %>%
    summarise(N = sum(N),
              BasalArea = sum(BasalArea),
              Biomass = sum(Biomass))

species_mortality <- mortality %>%
    group_by(Year, SpeciesID) %>%
    summarise(BiomassLoss = sum(Biomass))

species <- inner_join(species, species_mortality)

species <- species %>%
    group_by(SpeciesID) %>%
    arrange(Year) %>%
    mutate(Productivity = Biomass - lag(Biomass) + BiomassLoss,
           Productivity = ifelse(Productivity < 0, 0, Productivity)) %>%
    select(-BiomassLoss)

saveRDS(species, file = paste0(output_directory, "/PPA_output_processed_species.rds"))
