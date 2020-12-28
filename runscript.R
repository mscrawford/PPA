# TODO
#   Parameterize based on
#       Initial community files
#       Species files
#       * Parameterization will have to account for instances where different species files have also different initial communities
#   Run parallel


# Libraries ---------------------------------------------------------------

library(tictoc)


# Global mutable parameters -------------------------------------------------------------------

DEBUG <- TRUE

USE_INITIAL_COMMUNITY <- TRUE

CALCULATE_INTERNAL_SEED_RAIN <- TRUE
CALCULATE_EXTERNAL_SEED_RAIN <- FALSE


# Directories
base_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
# base_directory <- getwd() # Ensure that your base directory is the PPA folder
output_directory = paste0(base_directory, "/output/")

species_file = paste0(base_directory, "/input/PPA_FG5_filtered.csv")
initComm_file = paste0(base_directory, "/input/PPA_initial_state_fg5_secondary.csv")
PPA_script = paste0(base_directory, "/PPA.R")
postprocessing_script = paste0(base_directory, "/postprocessing.R")

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
    spVitals <- parameterization$spVitals
    initComm <- parameterization$initComm

    results <- run_simulation(spVitals, initComm)

    return(results)
}


run_parallel <- function(parameterization)
{

}


# parameterize ------------------------------------------------------------
parameterize = function()
{
    # Define the species present in the simulation
    spVitals <- read.table(species_file, sep = ",", header = TRUE)

    # Define the initial community composition
    initComm <- NULL
    if (USE_INITIAL_COMMUNITY)
    {
        initComm <- read.table(initComm_file, sep = "\t", header = FALSE)
        initComm <- as.matrix(initComm)
        names(initComm) <- NULL
        initComm <- cbind(initComm[,c(2,3)], NA, initComm[,1])
        initComm <- calculateLayers(initComm)
    } else {
        initComm <- generateDefaultCommunity(spVitals)
    }

    # A list of all the precursors to one simulation run
    return(list(spVitals = spVitals, initComm = initComm))
}


# generateDefaultCommunity ------------------------------------------------
# Generates a default, initial community. This data matrix closely resembles that
# of the external seed rain, except that by default all individuals within the
# initial community are in the first light layer (the canopy).
generateDefaultCommunity <- function(spVitals)
{
    initComm <- generateExternalSeedRain(spVitals)
    initComm[,3] <- 1
    return(initComm)
}


# Run ---------------------------------------------------------------------

tic() run() toc()
source(file = postprocessing_script)
