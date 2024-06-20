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
    deathRate += age/10 * 0.0; //PREDATOR MODIFIED
    if(deathRate > 1.0) { deathRate = 1.0; } //PREDATOR MODIFIED PAKJDF;KJ
    return deathRate;
}

// Kill off population due to natural deaths (e.g. getting eaten).
void Death::kill(State& state) {
    //for now will call weibull... potentially will make polymorphic later
    //will also only kill off susceptibles for now
    double susceptibleDeaths = 0; //this is a counter
    double infectedDeaths = 0;
    for (size_t ageIndex = 0; ageIndex < state.compParms.ageSize - 1; ageIndex++){ //age size - 1 to avoid double counting with age deaths
        double ageInYears = (double) ageIndex * state.intParms.deltaTime;
        //removes susceptible density & returns amount subtracted
        susceptibleDeaths += killSusceptibles(ageInYears, state, ageIndex);
        infectedDeaths += killInfecteds(ageInYears, state, ageIndex);
    }
    state.compParms.naturalDeaths = (susceptibleDeaths + infectedDeaths) * state.intParms.deltaTime;
}

//Subtract susceptible density for an age in accordance with deathrate. Returns amount subtracted.
double Death::killSusceptibles(double ageInYears, State& state, size_t index){
    double susDeathRate = weibullDeathRate(ageInYears) * state.intParms.deltaTime; //multiply by dt because return value is in deaths/year units
    std::vector<double>& susVec = *state.susceptibles->getCurrentState(); //check this referene
    double killed = susDeathRate * susVec[index];
    susVec[index] -= killed; //kill
    return killed; //record deaths
}

// Subtract infected density for an age in accordance with deathrate. Returns amount subtracted.
// In order to get infected density for age, integrates out infection variable.
double Death::killInfecteds(double ageInYears, State& state, size_t ageIndex){
    // multiply by dt because return value is in deaths/year units
    double infDeathRate = weibullDeathRate(ageInYears) * state.intParms.deltaTime;
    double infDeaths = 0;
    double deltaInfection;
    //infection size -1 to avoid double counting infection deaths
    for (size_t infectionIndex = 0; infectionIndex < state.compParms.infectionSize - 1; infectionIndex++){ 
        double currentInfecteds = state.infecteds->getIndex(ageIndex, infectionIndex);
        double infectedsKilled = currentInfecteds * infDeathRate;
        deltaInfection = state.compParms.deltaInfectionForLoad->at(infectionIndex);
        state.infecteds->setIndex(ageIndex, infectionIndex, currentInfecteds - infectedsKilled);
        infDeaths += infectedsKilled * deltaInfection; //we integrate out infection to find density at age
    }
    return infDeaths; //height of age vs pop density function for infecteds at age
}
