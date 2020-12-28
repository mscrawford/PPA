# TODO Create a legend for the cut factor of size_classes
#   (As of right now, the darker the color, the smaller the size class)

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(cowplot)


# Directories -------------------------------------------------------------

base_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
output_directory = paste0(base_directory, "/output/")


# Plotting size-class cohorts ---------------------------------------------

size_classes <- readRDS(paste0(output_directory, "PPA_output_processed_size_classes.rds"))

plots <- map(.x = size_classes,
             .f = ~ {
                 ggplot(.x) +
                     geom_line(aes(x = Year,
                                   y = BasalArea,
                                   color = SizeClass),
                               show.legend = FALSE) +
                     facet_grid(cols = vars(SpeciesID)) +
                     scale_color_grey() +
                     labs(x = "Year",
                          y = "Basal Area",
                          title = paste("Size class width:", unique(.x$SizeClass_width))) +
                     theme_minimal()
             }
)

plot_grid(plotlist = plots, nrow = length(plots))


# Plotting species-level data ---------------------------------------------

species <- readRDS(paste0(output_directory, "PPA_output_processed_species.rds"))

ggplot(species) +
    geom_line(aes(x = Year,
                  y = BasalArea,
                  color = as.factor(SpeciesID))) +
    labs(color = "SpeciesID") +
    theme_minimal()
