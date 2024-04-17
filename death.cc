#include "death.h"

#include <cmath>
#include <cstdio>

// Constructor
Death::Death(const State& state) {
    weibullSurvivorshipKappa_ = state.modParms.kappa;
    weibullSurvivorshipLambda_ = state.compParms.lambda;
}

// Given age, calculate survival coefficient
double Death::weibullDeathRate(double age) {
    double deathRate = weibullSurvivorshipKappa_ * weibullSurvivorshipLambda_ * pow(weibullSurvivorshipLambda_ * age, weibullSurvivorshipKappa_ - 1);
    return deathRate;
}

// Kill off population due to natural deaths (e.g. getting eaten).
void Death::kill(State& state) {
  //for now will call weibull... potentially will make polymorphic later
  //will also only kill off susceptibles for now
    double susceptibleDeaths = 0; //this is a counter
    double susTotal = 0;
    for(size_t ageIndex = 0; ageIndex < state.compParms.ageSize; ageIndex++){
        double ageInYears = (double) ageIndex * state.intParms.deltaTime;
        std::vector<double>& susVec = *state.susceptibles->getCurrentState(); //check this referene
        double deathRate = weibullDeathRate(ageInYears) * state.intParms.deltaTime; //multiply by dt because return value is in deaths/year units
        double killed = deathRate * susVec[ageIndex];
        susVec[ageIndex] -= killed; //kill
        susceptibleDeaths += killed; //record deaths
        for(size_t infectionIndex = 0; infectionIndex < state.compParms.infectionSize; infectionIndex++){ //this can be parallel
            double infectedsKilled = state.infecteds->getIndex(ageIndex,infectionIndex) * deathRate;
            //state.infecteds->setIndex(ageIndex, infectionIndex, numLeftInfected); //infected deaths for when implemented
        }
    }
    state.compParms.naturalDeaths = susceptibleDeaths;
    double tot = susceptibleDeaths * state.intParms.deltaTime;

}
