#ifndef INTEGRATOR
#define INTEGRATOR

#include "susceptibles.h"
#include "infecteds.h"
#include "birth.h"
#include "death.h"

//Integrate with forward Euler scheme
class forwardEuler {
public:
    forwardEuler(BirthScheme b); //constructor
    double * doBirths();
    double * doDeaths();
    double * setUpSusceptibles();
    double * setUpInfecteds();

private:
    BirthScheme births;
    Susceptibles S; //Not proper camalCase, but named thus to follow S-I-R convention
    Infecteds I;
    integrationParams stepSize;
    modelParams params;
};

//Helper structs kept here to make it easier to read!

//Store human-readable model parameters
struct modelParams {
    double maxAge;
    double aveLifespan;
    double aveInfectiousPeriod; //We derive an intrinsic growth rate based on ave infectious period
    double aveInitInfectionLoad; //and ave initial infection load. We don't carry over distributional properties
    //because it's not a terrible assumption that there's an intrinsic prion doubling time. Also, that's super complicated.
    //Note that if using to *fit* a model, you would want to keep pairs together! Please don't input actual averages from data.
};

//Store timestep for integration
struct integrationParams{
    double deltaTime;
    double deltaAge;
    double deltaInfection;
};

#endif
