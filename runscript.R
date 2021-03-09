# Perfect Plasticity Approximation model (Strigul et al. 2008)
# Adapted from Rüger et al. 2020 (code written by Caroline Farrior, cfarrior@gmail.com; https://github.com/cfarrior/Ruger_etal_2020)

# TODO
#   Parameterize based on
#       Initial community files
#       Species files
#       * Parameterization will have to account for instances where different species files have also different initial communities
#   run_parallel

# Libraries ---------------------------------------------------------------

library(tictoc)


# Global mutable parameters -----------------------------------------------

DEBUG <- TRUE

USE_INITIAL_COMMUNITY <- TRUE

CALCULATE_INTERNAL_SEED_RAIN <- TRUE
CALCULATE_EXTERNAL_SEED_RAIN <- FALSE

# Directories
base_directory   <- dirname(rstudioapi::getActiveDocumentContext()$path)
output_directory <- paste0(base_directory, "/output/")

# Files
species_file  <- paste0(base_directory, "/input/PPA_FG5_filtered.csv")
initComm_file <- paste0(base_directory, "/input/PPA_initial_state_fg5_secondary.csv")

# Scripts
PPA_script            <- paste0(base_directory, "/PPA.R")
postprocessing_script <- paste0(base_directory, "/postprocessing.R")
plotting_script       <- paste0(base_directory, "/plot.R")

source(PPA_script)


# run ---------------------------------------------------------------------

run <- function()
{
    parameterization <- parameterize()

    results <- run_serial(parameterization)

    saveRDS(results[[1]], file = paste0(output_directory, "/PPA_output_raw_cohort.rds"))
    saveRDS(results[[2]], file = paste0(output_directory, "/PPA_output_raw_cohort_mortality.rds"))
}


# run_serial --------------------------------------------------------------

run_serial <- function(parameterization)
{
    spVitals_list <- parameterization$spVitals_list
    initComm <- parameterization$initComm

    results <- run_simulation(spVitals_list, initComm)

    return(results)
}


# run_parallel ------------------------------------------------------------

run_parallel <- function(parameterization)
{

}


# parameterize ------------------------------------------------------------

parameterize = function()
{
    # Define the species present in the simulation
    spVitals <- read.table(species_file, sep = ",", header = TRUE)
    spVitals_list <- disaggregateSpeciesVitals(spVitals)

    # Define the initial community composition
    initComm <- NULL
    if (USE_INITIAL_COMMUNITY)
    {
        initComm <- read.table(initComm_file, sep = "\t", header = FALSE)
    }

    return(list(spVitals_list = spVitals_list, initComm = initComm))
}


# disaggregateSpeciesVitals --------------------------------------------------
# Disaggregates the `spVitals` dataframe into its component parts, with each documenting one
# aspect of the species' vital rates. These parameters can then be exported to the
# environment of the calling function. This reduces the amount of code needed to incorporate
# a variable number of layers within the simulation.
#
# This function assumes that there is a growth rate and mortality rate for each layer,
# and that they are ordered decreasing with crown class.
disaggregateSpeciesVitals <- function(spVitals)
{
    N   <- nrow(spVitals)                                    # community size (integer)

    ID  <- spVitals %>% select(SpeciesID) %>% pull()         # species ID (vector)
    G   <- spVitals %>% select(contains("G"))                # growth rates (matrix)
    mu  <- spVitals %>% select(contains("mu"))               # mortality rates (matrix)
    Fec <- spVitals %>% select(contains("F")) %>% pull()     # fecundity rate (vector)
    wd  <- spVitals %>% select(contains("wd")) %>% pull()    # wood density (vector)

    # Set the number of layers defined within the input species data.frame
    nLayers <- ncol(G) + 1

    if (DEBUG)
    {
        assertthat::assert_that(!is.null(nLayers) & nLayers > 0)

        assertthat::assert_that(ncol(G) == ncol(mu))

        columnOrder <- as.numeric(substring(colnames(G), 2))
        assertthat::assert_that(!is.unsorted(columnOrder) & columnOrder[length(columnOrder)] > columnOrder[1])
    }

    return(list(nLayers = nLayers, N = N, ID = ID, G = G, mu = mu, Fec = Fec, wd = wd))
}


# Run ---------------------------------------------------------------------
tic(); run(); toc()
source(postprocessing_script)
source(plotting_script)
