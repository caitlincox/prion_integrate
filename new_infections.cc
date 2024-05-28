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
    calculateDeltaInfecteds(state);
}

// Given S(a*, t*), using forward Euler integration in time, finds the
// component of S(a* + da, t* + dt) - S(a*, t*) that is due to infection of
// susceptibles.  It calculates this for all values of a* between 0 and (maxAge
// - da) In other words, finds the density fcn for susceptibles infected
// between time t* and time t* + dt.  also moves this one da, dt forward in
// age/time
void NewInfections::calculateDeltaSusceptibles(State& state) {
    betaI_ = calculateInfectionCoefficient(state);
    auto& susceptibles = *state.susceptibles->getCurrentState();
    auto& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double deltaT = state.intParms.deltaTime;
    totalSusceptible_ = 0.0;
    for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++){
        //Here, scaling by dt bc dt/1 = dt + da / (1 + 1). It's to adjust step size.
        deltaSusceptibles[xAge] = betaI_ * susceptibles[xAge] * deltaT;
        totalSusceptible_ += deltaSusceptibles[xAge] * deltaT;
    }
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
// Moves everything forward by d-theta
// Note that susceptibles that move over max load via this step also vanish instantly
void NewInfections::calculateDeltaInfecteds(State& state) {
    // Shorthand variables to make code less verbose
    std::vector<double>& loadVec = *state.compParms.rowLoads;
    std::vector<double>& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double shapeParam = state.modParms.gammaShapeParam;
    double deltaT = state.intParms.deltaTime;
    double totalInfected = 0.0;
    for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        // Get proportion that lands in previous infection level & advance load by 1...
        double gammaVal = gammaDist(loadVec[xLoad], state.modParms.aveInitInfectionLoad, shapeParam);
        double newDeltaLoad = state.compParms.deltaInfectionForLoad->at(xLoad);
        double newDeltaArea = deltaT * newDeltaLoad;
        // Loop over infection levels. Skip age 0 -- no infections
        for(size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            double infectedPopDensity = deltaSusceptibles[xAge] * gammaVal;
            totalInfected += infectedPopDensity * newDeltaArea;
            deltaInfections_->setIndex(xAge, xLoad,  infectedPopDensity);
        }
    }
    assertAproxEqual(totalInfected, state.compParms.susceptiblePop * betaI_ * deltaT, 0.1);
    adjustForConstantPop(state, totalInfected);
}

// There are small errors introduced by the non-continuous buckets used to represent infecteds.
// Correct the population so it is constant to within rounding errors.
void NewInfections::adjustForConstantPop(State& state, double totalInfected) {
    if (totalInfected == 0.0) {
        // Nothing to do.
        return;
    }
    double ratio = totalInfected / totalSusceptible_;
    std::vector<double>& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    for(size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        deltaSusceptibles[xAge] *= ratio;
    }
}

void NewInfections::moveInfecteds(State& state){
    auto& susceptibles = *state.susceptibles->getCurrentState(), deltaSus = *deltaSusceptibles_->getCurrentState();
    Infecteds& infecteds = *state.infecteds;
    Infecteds& deltaInf = *deltaInfections_;
    double totalSus = 0.0;
    double totalInf = 0.0;
    double deltaT = state.intParms.deltaTime;
    for(size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
        susceptibles[xAge] -= deltaSus[xAge];
        totalSus += deltaSus[xAge] * deltaT;
        for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
            double deltaLoad = state.compParms.deltaInfectionForLoad->at(xLoad);
            totalInf += deltaInf.getIndex(xAge, xLoad) * deltaT * deltaLoad;
            double newDensity = infecteds.getIndex(xAge, xLoad) + deltaInf.getIndex(xAge, xLoad);
            infecteds.setIndex(xAge, xLoad, newDensity);
        }
    }
    assertAproxEqual(totalInf, totalSus, 0.001);
}


// TEST: count moving infecteds and susceptibles to show equal approx
