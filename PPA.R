# Perfect Plasticity Approximation model (Strigul et al. 2008)
# Adapted from Rüger et al. 2020 (code written by Caroline Farrior)


# Libraries ---------------------------------------------------------------

library(dplyr)

# Global constants ----------------------------------------------------------------------------

# Note:
#   crown radius = 0.5*dbh^0.62; dbh in cm and crown radius in m
#   crown area = phi*dbh^theta

# Commonly used constants
phi <- round(pi * 0.5^2, 4)
theta <- 1.24

# Global constants
PA <- 10000             # plot size: 1 ha
dnot <- 1               # initial tree diameter: 1 cm
mincohortN <- 0.001     # minimum cohort size
deltaT <- 5             # model timestep. This is a fixed parameter (see generateInternalSeedRain).
cutT <- deltaT
maxT <- 200*deltaT      # simulation terminal year
nLayers <- 5

# Global mutable parameters -------------------------------------------------------------------

set.seed(1)

USE_INITIAL_COMMUNITY <- TRUE

CALCULATE_INTERNAL_SEED_RAIN <- TRUE
CALCULATE_EXTERNAL_SEED_RAIN <- FALSE

# Directories
base_directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
# base_directory <- getwd() # Ensure that your base directory is the PPA folder
species_file = "/input/PPA_FG5_filtered.csv"
initComm_file = "/input/PPA_initial_state_fg5_secondary.csv"
output_directory = "/output/"


# run ---------------------------------------------------------------------
run <- function()
{
    parameterization <- parameterize()

    results <- run_serial(parameterization)

    saveRDS(object = results,
            file = paste0(base_directory, output_directory, "/PPA_output.rds"))
}


# run_serial --------------------------------------------------------------
run_serial <- function(parameterization)
{
    spVitals <- parameterization$spVitals
    initComm <- parameterization$initComm
    externalSeedRain <- parameterization$externalSeedRain

    results <- run_simulation(spVitals, initComm, externalSeedRain)

    return(results)
}


# parameterize ------------------------------------------------------------
parameterize = function()
{
    # Define the species present in the simulation
    spVitals <- read.table(paste0(base_directory, species_file), sep = ",", header = TRUE)

    # Define the initial community composition
    initComm <- NULL
    if (USE_INITIAL_COMMUNITY)
    {
        initComm <- read.table(paste0(base_directory, initComm_file), sep = "\t", header = FALSE)
        initComm <- as.matrix(initComm)
        names(initComm) <- NULL
        initComm <- cbind(initComm[,c(2,3)], NA, initComm[,1])
        initComm <- calculateLayers(initComm)
    } else {
        initComm <- generateDefaultCommunity(spVitals)
    }

    # If applicable, define the external seed rain per time step
    externalSeedRain <- NULL
    if (CALCULATE_EXTERNAL_SEED_RAIN)
    {
        externalSeedRain <- generateExternalSeedRain(spVitals)
    }

    # A list of all the precursors to one simulation run
    return(list(spVitals = spVitals, initComm = initComm, externalSeedRain = externalSeedRain))
}


# run_simulation ------------------------------------------------------------------------------
# Runs one simulation with a given set of species (spVitals) and initial community (initComm).
run_simulation <- function(spVitals, initComm, externalSeedRain)
{
    list2env(x = disaggregateSpeciesVitals(spVitals), envir = environment())

    results <- NULL

    # Data matrix has columns:
    # (1) diameter per individual
    # (2) # of individuals
    # (3) crown class (1: Canopy, 2: understory, 3: 2nd understory, ...)
    # (4) species
    data <- initComm
    for (t in seq(0, maxT, by = deltaT))
    {
        startingData <- data

        # Mortality matrix has columns:
        # (1) diameter per dead individual
        # (2) # of dead individuals
        # (3) crown class of dead individuals
        # (4) species
        mortality <- matrix(0, nrow = 0, ncol = 4)

        for (sp in unique(data[,4]))
        {
            # Step 1: Mortality
            muV <- mu[ID == sp, ]
            for (i in seq(1, dim(data)[1])[data[,4] == sp])
            {
                mortality <- rbind(mortality, c(data[i, 1],
                                                data[i, 2] * (1 - (1 - muV[[data[i, 3]]])^deltaT),
                                                data[i, 3],
                                                sp))

                data[i, 2] <- data[i, 2] * (1 - muV[[data[i, 3]]])^deltaT
            }
            data <- data[data[,2] > mincohortN, , drop = FALSE]

            # Step 2: Growth by layer
            GV <- G[ID == sp, ]
            for (layer in unique(data[data[,4] == sp, , drop = FALSE][,3]))
            {
                data[(data[,4] == sp) & (data[,3] == layer), 1] <-
                    data[(data[,4] == sp) & (data[,3] == layer), 1] + GV[[layer]] * deltaT
            }
        }

        data <- data[data[,2] > mincohortN, , drop = FALSE] # Remove any cohort with too few individuals
        data <- data[data[,1] > 0, , drop = FALSE]          # Remove any cohort with a negative avg. diameter

        # Step 3: Reproduce
        # Step 3a. Internal seed rain
        if (CALCULATE_INTERNAL_SEED_RAIN)
        {
            internalSeedRain <- generateInternalSeedRain(data, spVitals)
            data <- rbind(data, internalSeedRain)
        }

        # Step 3b. External seed rain
        if (CALCULATE_EXTERNAL_SEED_RAIN)
        {
            assertthat::assert_that(!is.null(externalSeedRain))
            data <- rbind(data, externalSeedRain)
        }

        # Step 4: Assign crown class
        data <- calculateLayers(data)

        # Step 5: Record
        if (floor(t/cutT) == t/cutT & t != 0) # Records every cutT timesteps
        {
            results <- rbind(results, calculateOutput(data, startingData, mortality, t, spVitals))
        }
    }

    return(results)
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

# generateInternalSeedRain ------------------------------------------------
# This function generates the internal seed rain generated per species across deltaT
# years. It returns a matrix that can, optionally, be added to the main data matrix on a
# yearly basis.
generateInternalSeedRain <- function(data, spVitals)
{
    list2env(x = disaggregateSpeciesVitals(spVitals), envir = environment())

    tmp_internalSeedRain <- NULL
    for (sp in unique(data[,4]))
    {
        t_sp_cohorts <- data[data[,4] == sp, , drop = FALSE]

        nbaabg <- n_ba_agb(t_sp_cohorts, ID, wd)

        tmp_internalSeedRain <- rbind(tmp_internalSeedRain,
                                      c(dnot,
                                        nbaabg$ba * Fec[ID == sp],
                                        nLayers - 1,
                                        sp,
                                        G[ID == sp, ncol(G)]))
    }

    # Generate deltaT cohorts, each at a different stage of growth
    internalSeedRain <- NULL
    for (i in seq(deltaT - 1, 0, -1))
    {
        baby <- tmp_internalSeedRain
        baby[, 1] <- baby[, 1] + baby[, ncol(baby)] * i; # baby[, ncol(baby)] is to the smallest layer's growth rate
        # Because the census occurs on a five year time step, seedling mortality is implicitly considered
        baby <- as.matrix(baby[, -5]) # remove growth column, matrix is now in the "normal" cohort-level data form

        internalSeedRain <- rbind(internalSeedRain, baby)
    }

    return(internalSeedRain)
}


# generateExternalSeedRain ------------------------------------------------
# This function generates the seed rain per species across deltaT years. It
# returns a matrix that can, optionally, be added to the main data matrix on a
# yearly basis.
generateExternalSeedRain <- function(spVitals)
{
    list2env(x = disaggregateSpeciesVitals(spVitals), envir = environment())

    # Seed rain matrix
    tmp_externalSeedRain <- matrix(1, nrow = N, ncol = 4)

    # with columns:
    tmp_externalSeedRain[,1] <- dnot                # cohort diameter
    tmp_externalSeedRain[,2] <- PA/10000 * Fec      # number of individuals
    tmp_externalSeedRain[,3] <- nLayers - 1         # crown class, by default in the lowest understory
    tmp_externalSeedRain[,4] <- ID                  # speciesID

    externalSeedRain <- NULL
    for (i in seq(deltaT - 1, 0, -1))
    {
        baby <- tmp_externalSeedRain
        baby[,1] <- baby[,1] + G[, ncol(G)] * i # G[, ncol(G)] corresponds to the smallest growth rate
        baby[,2] <- baby[,2] * (1 - mu[, ncol(mu)])^i # mu[, ncol(mu)] corresponds to the smallest mortality rate
        externalSeedRain <- rbind(externalSeedRain, baby)
    }

    return(externalSeedRain)
}


# disaggregateSpeciesVitals --------------------------------------------------
# Disaggregates the `spVitals` dataframe into its component parts, which each document one
# aspect of the species' vital rates. These parameters can then be exported to the
# environment of the calling function. This reduces the amount of code needed to incorporate
# a variable number of layers within the simulation.
#
# This function assumes that there is a growth rate and mortality rate for each layer.
disaggregateSpeciesVitals <- function(spVitals)
{
    N <- nrow(spVitals)                                     # community size (integer)

    ID <- spVitals %>% select(SpeciesID) %>% pull()         # species ID (vector)
    G <- spVitals %>% select(contains("G"))                 # growth rates (matrix)
    mu <- spVitals %>% select(contains("mu"))               # mortality rates (matrix)
    Fec <- spVitals %>% select(contains("F")) %>% pull()    # fecundity rate (vector)
    wd <- spVitals %>% select(contains("wd")) %>% pull()    # wood density (vector)

    assertthat::assert_that(ncol(G) == ncol(mu))
    assertthat::assert_that(ncol(G) == nLayers - 1)

    return(list(N = N, ID = ID, G = G, mu = mu, Fec = Fec, wd = wd))
}


# calculateOutput ---------------------------------------------------------
calculateOutput <- function(data, startingData, mortality, year, spVitals)
{
    list2env(x = disaggregateSpeciesVitals(spVitals), envir = environment())

    out <- NULL

    for (sp in unique(data[,4]))
    {
        ss_t1 <- data[data[,4] == sp, , drop = FALSE]
        ss_t0 <- startingData[startingData[,4] == sp, , drop = FALSE]
        ss_t0_mort <- mortality[mortality[,4] == sp, , drop = FALSE]

        nbaabg_t1 <- n_ba_agb(ss_t1, ID, wd)
        nbaabg_t0 <- n_ba_agb(ss_t0, ID, wd)
        nbaabg_t0_mort <- n_ba_agb(ss_t0_mort, ID, wd)

        productivity <- nbaabg_t1$agb - nbaabg_t0$agb + nbaabg_t0_mort$agb
        productivity <- ifelse(productivity < 0, 0, productivity)

        sp.out <- data.frame(Model = "PPA",
                             SpeciesID = sp,
                             `F` = Fec[ID == sp],
                             Year = year,
                             N = nbaabg_t1$n,
                             BasalArea = nbaabg_t1$ba,
                             Biomass = nbaabg_t1$agb,
                             Productivity = productivity)

        out <- rbind(out, sp.out)
    }

    return(out)
}


# n_ba_agb ----------------------------------------------------------------
# This function calculates the abundance, basal area, and aboveground biomass for each species.
# It assumes that the data matrix supplied to it is composed of only one species.
n_ba_agb <- function(data, ID, wd)
{
    assertthat::assert_that(length(unique(data[,4])) %in% c(0, 1))

    sp <- unique(data[,4])
    n <- sum(data[,2])
    ba <- sum( (data[,1]/200)^2 * pi * data[,2] )

    s_wd <- wd[ID == sp]
    agb <- sum(s_wd *
                   exp(-1.499 + 2.148*log(data[,1]) +
                           0.207*log(data[,1])^2 -
                           0.0281*log(data[,1])^3) /
                   1000 * data[,2])

    return(list(n = n, ba = ba, agb = agb))
}


# calculateLayers ---------------------------------------------------------
# This function reassigns each cohort per species to a new layer, depending on the
# PPA assumption. It first calculates the total crown area (CA). If CA is greater than
# the plot area (PA), this function calls the CCassign function. Lastly, this function
# removes any any plants falling below the minimum considered canopy layer.
calculateLayers <- function(data)
{
    CA <- sum(phi * data[,1]^theta * data[,2]) # calculate total crown area of individuals in the plot
    if (CA <= PA)
    {
        data[,3] = 1 # if less than the ground area, everyone is in the canopy.
    } else {
        data <- CCassign_manylayers_continuous(data) # if greater than the ground area, go through CCassign
    }

    data <- data[data[,3] < nLayers, , drop = FALSE] # kill plants that fall into the smallest layer immediately

    return(data)
}


# Main cohorts function -----------------------------------------------------------------------
# assigns crown class based on the PPA assumption
# assumes all individuals have the same crown area and height allometries
# assumes that CAtot > PA
# works for any number of layers
CCassign_manylayers_continuous = function(data)
{
    CAv <- phi*data[,1]^theta
    data <- data[order(CAv,decreasing = TRUE), ]
    CAv <- CAv[order(CAv,decreasing = TRUE)]
    cohortCAv <- CAv*data[,2]
    cacaV <- cumsum(cohortCAv)

    data[,3] <- 1

    for (layers in seq(1, floor(sum(cohortCAv)/PA)))
    {
        # make a vector cacaV where the ith entry is the crown area of the i'th cohort plus all cohorts
        # with individuals of greater diameter.
        CAv <- phi*data[,1]^theta
        data <- data[order(CAv, decreasing = TRUE), ]
        CAv <- CAv[order(CAv, decreasing = TRUE)]
        cohortCAv <- CAv * data[,2]
        cacaV <- cumsum(cohortCAv)

        # pull out individuals that are definitely in the canopy and understory
        und <- data[cacaV > PA * layers, , drop = FALSE]
        can <- data[cacaV <= PA * layers, , drop = FALSE]

        # split the first cohort in the understory to fill the leftover open canopy space
        canCA <- max(0, sum(phi * can[,1]^theta * can[,2]))
        tosplit <- und[1, , drop = FALSE]
        opencan <- layers * PA - canCA
        splitind_incan <- opencan / (phi * tosplit[1, 1]^theta)
        und[1, 2] <- und[1, 2] - splitind_incan
        tosplit[,2] <- splitind_incan
        tosplit[,3] <- layers
        can <- rbind(can, tosplit)
        if (max(can[,3]) != layers)
        {
            print("error")
        }
        und[,3] <- layers + 1

        # piece the data back together
        data <- rbind(can, und)

        data <- data[data[,2] > mincohortN, , drop = FALSE]
    }

    return(data)
}


# Run ---------------------------------------------------------------------
run()


# Quick plot --------------------------------------------------------------

library(ggplot2)

results <- readRDS(file = paste0(base_directory, output_directory, "PPA_output.rds"))

ggplot(results) +
    geom_line(aes(x = Year,
                  y = BasalArea,
                  color = as.factor(SpeciesID))) +
    labs(color = "SpeciesID") +
    theme_minimal()
