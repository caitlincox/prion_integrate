#ifndef INTEGRATOR
#define INTEGRATOR

#include <memory>

#include "susceptibles.h"
#include "infecteds.h"
#include "birth.h"
#include "death.h"

//Helper structs kept here to make it easier to read!

//Store timestep for integration
struct IntegrationParams{
    double deltaTime;
    double deltaLogInfection;
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
    Integrator(BirthScheme b, State state);
    double* doBirths();
    double* doDeaths();
    double* setUpSusceptibles();
    double* setUpInfecteds();

private:
    std::unique_ptr<BirthScheme> births_;
    std::unique_ptr<State> state_;
    IntegrationParams stepSize_;
    ModelParams params_;
};

#endif
