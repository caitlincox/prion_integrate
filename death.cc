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
    double infectedDeaths = 0;
    for(size_t ageIndex = 0; ageIndex < state.compParms.ageSize; ageIndex++){
        double ageInYears = (double) ageIndex * state.intParms.deltaTime;
        susceptibleDeaths += killSusceptibles(ageInYears, state, ageIndex); //removes susceptible density & returns amount subtracted
        infectedDeaths += killInfecteds(ageInYears, state, ageIndex);
    }
    state.compParms.naturalDeaths = susceptibleDeaths + infectedDeaths;
    double dedSus = susceptibleDeaths * state.intParms.deltaTime;
    double dedInf = infectedDeaths * state.intParms.deltaTime;
}

//Subtract susceptible density for an age in accordance with deathrate. Returns amount subtracted.
double Death::killSusceptibles(double ageInYears, State& state, size_t index){
    double susDeathRate = weibullDeathRate(ageInYears) * state.intParms.deltaTime; //multiply by dt because return value is in deaths/year units
    std::vector<double>& susVec = *state.susceptibles->getCurrentState(); //check this referene
    double killed = susDeathRate * susVec[index];
    susVec[index] -= killed; //kill
    return killed; //record deaths
}

//Subtract infected density for an age in accordance with deathrate. Returns amount subtracted.
//In order to get infected density for age, integrates out infection variable.
double Death::killInfecteds(double ageInYears, State& state, size_t ageIndex){
    double infDeathRate = weibullDeathRate(ageInYears) * state.intParms.deltaTime; //multiply by dt because return value is in deaths/year units
    double infDeaths = 0;
    double deltaInfection;
    for(size_t infectionIndex = 0; infectionIndex < state.compParms.infectionSize; infectionIndex++){ 
        double currentInfecteds = state.infecteds->getIndex(ageIndex, infectionIndex);
        double infectedsKilled = currentInfecteds * infDeathRate;
        state.infecteds->setIndex(ageIndex, infectionIndex, currentInfecteds - infectedsKilled);
        if(infectionIndex == state.compParms.infectionSize - 1) { //to prevent indexing out of vector
            deltaInfection = 1 - (*state.compParms.columnLoads)[infectionIndex]; //Need non log scale
        }else {
            deltaInfection = (*state.compParms.columnLoads)[infectionIndex + 1] - (*state.compParms.columnLoads)[infectionIndex];
        }
        infDeaths += infectedsKilled * deltaInfection; //we integrate out infection to find density at age
    }
    
    return infDeaths; //height of age vs pop density function for infecteds at age
}