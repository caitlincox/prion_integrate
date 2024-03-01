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
};

//Integrate with forward Euler scheme
class Integrator {
public:
    Integrator(
        IntegrationParams integrationParams,
        ModelParams modelParams,
        std::unique_ptr<BirthScheme> births,
        std::unique_ptr<BirthScheme> deaths,
        std::unique_ptr<State> state)
    {
        integrationParams_ = integrationParams;
        modelParams_ = modelParams;
        births_ = std::move(births);
        state_ = std::move(state);
    }
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
    // Computed in initialcontinsion.cc.
    double intrinsicGrowthRate_;
};

#endif
