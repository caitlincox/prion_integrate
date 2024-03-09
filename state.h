#ifndef STATE
#define STATE

#include "infecteds.h"
#include "susceptibles.h"

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

struct State {
    std::unique_ptr<Susceptibles> susceptibles;
    std::unique_ptr<Infecteds> infecteds;
    IntegrationParams intParms;
    ModelParams modParms;
};

#endif
