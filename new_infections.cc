#include "new_infections.h"
#include "initialconditions.h"
#include "tests.h"

#include <cassert>

NewInfections::NewInfections(const State& state) {
    deltaInfections_ = std::make_unique<Infecteds>(state.compParms.ageSize, state.compParms.infectionSize);
    deltaSusceptibles_ = std::make_unique<Susceptibles>(state.compParms.ageSize);
}

void NewInfections::prepInfecteds(State& state) {
    calculateDeltaSusceptibles(state);
}

// Given S(a*, t*), using forward Euler integration in time, finds the
// component of S(a* + da, t* + dt) - S(a*, t*) that is due to infection of
// susceptibles.  It calculates this for all values of a* between 0 and (maxAge
// - da) In other words, finds the density fcn for susceptibles infected
// between time t* and time t* + dt. 
void NewInfections::calculateDeltaSusceptibles(State& state) {
    betaI_ = calculateInfectionCoefficient(state);
    auto& susceptibles = *state.susceptibles->getCurrentState();
    auto& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double deltaT = state.intParms.deltaTime;
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++){
        //Here, scaling by dt bc dt/1 = dt + da / (1 + 1). It's to adjust step size.
        deltaSusceptibles[xAge] = betaI_ * susceptibles[xAge] * deltaT;
        if(susceptibles[xAge] < deltaSusceptibles[xAge]) {
            deltaSusceptibles[xAge] = susceptibles[xAge];
        }
    }
    double tot = betaI_ * state.compParms.susceptiblePop * deltaT;
}

// Given current state, calculates the integral that corresponds to "Beta * I"
// in an ODE SIR model. More details in Stringer et al. paper.
double NewInfections::calculateInfectionCoefficient(State& state) {
    double coefficient = 0;
    double loadBeta;
    double deltaLoad;
    for (size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        loadBeta = betaForLoad(state, state.compParms.rowLoads->at(xLoad)); // higher load -> more shedding
        deltaLoad = state.compParms.deltaInfectionForLoad->at(xLoad); // Integration delta for infection load. Variable b/c log linear.
        for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            coefficient += state.infecteds->getIndex(xAge, xLoad) * deltaLoad * loadBeta;
        }
    }
    return coefficient * state.intParms.deltaTime; // da doens't change by age, so can multiply by it at end. (dt = da)
}

// Linear scaling based on assumption that quantity of prions in body is
// directly proportional to prion shedding
double NewInfections::betaForLoad(State& state, double infectionLoad) {
    return (state.modParms.beta * infectionLoad);
}

// Assigns infection loads to delta susceptibles given gamma fcn.
// Note that susceptibles that move over max load via this step also vanish instantly
void NewInfections::calculateDeltaInfecteds(State& state) {
    // Shorthand variables to make code less verbose
    std::vector<double>& loadVec = *state.compParms.rowLoads;
    std::vector<double>& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double shapeParam = state.modParms.gammaShapeParam;
    double deltaT = state.intParms.deltaTime;
    double totalInfected = 0.0;
    for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        double gammaVal = gammaDist(loadVec[xLoad], state.modParms.aveInitInfectionLoad, shapeParam);
        double deltaLoad = state.compParms.deltaInfectionForLoad->at(xLoad);
        double deltaArea = deltaT * deltaLoad;
        // Loop over infection levels. We include age 0 bc infecteds haven't moved yet.
        for(size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            double infectedPopDensity = deltaSusceptibles[xAge] * gammaVal;
            totalInfected += infectedPopDensity * deltaArea;
            deltaInfections_->setIndex(xAge, xLoad,  infectedPopDensity);
        }
    }   //assertAproxEqual(totalInfected, state.compParms.susceptiblePop * betaI_ * deltaT);
}

void NewInfections::moveInfecteds(State& state){
    auto& susceptibles = *state.susceptibles->getCurrentState(), deltaSus = *deltaSusceptibles_->getCurrentState();
    Infecteds& infecteds = *state.infecteds;
    Infecteds& deltaInf = *deltaInfections_;
    for(size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        susceptibles[xAge] -= deltaSus[xAge];
        for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
            double newDensity = infecteds.getIndex(xAge, xLoad) + deltaInf.getIndex(xAge, xLoad);
            infecteds.setIndex(xAge, xLoad, newDensity);
        }
    }
}


// TEST: count moving infecteds and susceptibles to show equal approx
