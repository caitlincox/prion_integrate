#ifndef STATE
#define STATE

#include "infecteds.h"
#include "susceptibles.h"

#include <string>

//Helper structs kept here to make it easier to read!

// Integration parameters.  All times are in years.
struct IntegrationParams {
    double deltaTime;
    double totalTime;
    double infectionSizeToAgesSizeRatio;
};

//Store human-readable model parameters
struct ModelParams {
    double maxAge;
    double aveLifespan;
    double aveInfectiousPeriod; //We derive an intrinsic growth rate based on ave infectious period
    double aveInitInfectionLoad; //and ave initial infection load. We don't carry over distributional properties
    //because it's not a terrible assumption that there's an intrinsic prion doubling time. Also, that's super complicated.
    //Note that if using to *fit* a model, you would want to keep pairs together! Please don't input actual averages from data.
    double beta;  // This is a constant used to model the Beta function.
    double kappa;
    double initialSusceptiblePop;
    double initialInfectedPop;
    double gammaShapeParam; //in Stringer et al. this param is called "c"
};

// Parameters computed on the fly.
struct ComputedParams {
    // Change in log(infection) between infection buckets.  For 100 buckets, it
    // should be close to 0.001.
    double deltaLogInfection;
    // Currently computed in death.cc (unimplemented).
    double lambda = 0.0;
    size_t ageSize = 0;
    size_t infectionSize = 0;; 
    // Parameters computed per step.
    double totalInfection = 0.0;
    double infectedPop = 0.0;
    double susceptiblePop = 0.0;
    double popSize = 0.0;  // Sum of infected and susceptible.
    double transferRate = 0.0;
    double aveLoad = 0.0;
    double ageDeaths = 0.0;
    double infectionDeaths = 0.0;
    double maxInfectedsPopDensity = 0.0;
    double maxSusceptiblesPopDensity = 0.0;
    // Computed in initialcontinsion.cc.
    double intrinsicGrowthRate = 0.0; //in Stringer et al. this param is called c1
    double firstBucketLogLoad = 0.0;
    // These are the infections loads for column i in the infecteds table.
    // The first column has all 0's meaning no infecteds have zero load.
    // Also there are no 0-age infecteds, so the 0 row is also 0's.
    std::unique_ptr<std::vector<double>> columnLoads;
};

// This structure contains all state in the system, including parameters
// computed on the fly.  It is effectively the shared database, and all modules
// read from and write to it directly.
struct State {
    State(IntegrationParams integrationParams, ModelParams modelParams);
    // Update the state to reflect the state after a deltaTime step.
    void timeStep();
    // After taking a time step, update computed parameters such as popSize.
    void updateComputedParameters();
    // Write a portable greymap (.pgm) image from the infecteds.
    void writeInfectedsPGM(const std::string& filename) const;
    // Write a portable bitmap (.pbm) graph for the suseptibles.
    void writeSusceptiblesPBM(const std::string& filename, uint32_t width) const;

    std::unique_ptr<Susceptibles> susceptibles;
    std::unique_ptr<Infecteds> infecteds;
    IntegrationParams intParms;
    ModelParams modParms;
    ComputedParams compParms;

private:
    void initializeComputedParameters();
};

#endif
