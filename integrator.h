#ifndef INTEGRATOR
#define INTEGRATOR

#include <memory>

#include "susceptibles.h"
#include "infecteds.h"
#include "birth.h"
#include "death.h"

//Helper structs kept here to make it easier to read!

//Store timestep for integration.  All times are in years.
struct IntegrationParams{
    double deltaTime;
    size_t numInfectionLoadBuckets; 
    double deltaLogInfection;
    double totalTime;
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
};

//Integrate with forward Euler scheme
class Integrator {
public:
    Integrator(
        const IntegrationParams& integrationParams,
        const ModelParams& modelParams,
        std::unique_ptr<BirthScheme> births,
        std::unique_ptr<Death> deaths,
        std::unique_ptr<State> state);
    void run();

private:
    // Update state_ for the current time step.
    void update();
    void findSums();

    IntegrationParams integrationParams_;
    ModelParams modelParams_;
    std::unique_ptr<BirthScheme> births_;
    std::unique_ptr<State> state_;
    // Sums computed per step.
    double totalInfection_;
    double infectedPop_;
    double susceptiblePop_;
    double transferRate_;
    // Computed in initialcontinsion.cc.
    double intrinsicGrowthRate_;
    double firstBucketLogLoad_;
    // These are the infections loads for column i in the infecteds table.
    // The first column has all 0's meaning no infecteds have zero load.
    // Also there are no 0-age infectes, so the 0 row is also 0's.
    std::unique_ptr<std::vector<double>> columnLoads_;
};

#endif
