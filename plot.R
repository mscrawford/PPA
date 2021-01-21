# TODO Create a legend for the cut factor of size_classes
#   (As of right now, the darker the color, the smaller the size class)

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)


# Directories -------------------------------------------------------------

base_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
output_directory <- paste0(base_directory, "/output/")
plot_directory <- paste0(base_directory, "/plots/")


# Plotting size-class cohorts ---------------------------------------------

size_classes <- readRDS(paste0(output_directory, "PPA_output_processed_size_classes.rds"))

p <- ggplot(size_classes) +
    geom_line(aes(x = Year,
                  y = BasalArea,
                  color = as.factor(SpeciesID))) +
    facet_wrap(facets = vars(SizeClass), labeller = label_both) +
    labs(x = "Year",
         y = "Basal Area",
         color = "SpeciesID") +
    theme_minimal()

save_plot(plot = p, filename = paste0(plot_directory, "SizeClasses.pdf"))


# Plotting species-level data ---------------------------------------------

species <- readRDS(paste0(output_directory, "PPA_output_processed_species.rds"))

p <- ggplot(species) +
    geom_line(aes(x = Year,
                  # y = BasalArea,
                  y = N,
                  color = as.factor(SpeciesID))) +
    labs(x = "Year",
         color = "SpeciesID") +
    theme_minimal()

save_plot(plot = p, filename = paste0(plot_directory, "Species.PDF"))
