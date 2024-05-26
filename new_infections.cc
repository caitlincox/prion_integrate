#include "new_infections.h"
#include "initialconditions.h"

NewInfections::NewInfections(const State& state) {
    deltaInfections_ = std::make_unique<Infecteds>(state.compParms.ageSize, state.compParms.infectionSize);
    deltaSusceptibles_ = std::make_unique<Susceptibles>(state.compParms.ageSize);
}

void NewInfections::prepInfecteds(State& state) {
    calculateDeltaSusceptibles(state);
    calculateDeltaInfecteds(state);
}

// Given S(a*, t*), using forward Euler integration in time, finds the component of
// S(a* + da, t* + dt) - S(a*, t*) that is due to infection of susceptibles.
// It calculates this for all values of a* between 0 and (maxAge - da)
// In other words, finds the density fcn for susceptibles infected
// between time t* and time t* + dt.
// also moves this one da, dt forward in age/time
void NewInfections::calculateDeltaSusceptibles(State& state) {
    double betaI = calculateInfectionCoefficient(state);
    double rate;
    auto& susceptibles = *state.susceptibles->getCurrentState();
    auto& deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double count;
    for (size_t xAge = 0; xAge < state.compParms.ageSize - 1; xAge++){
        rate = betaI * susceptibles[xAge];
        //Here, scaling by dt bc dt/1 = dt + da / (1 + 1). It's to adjust step size.
        deltaSusceptibles[xAge + 1] = rate * state.intParms.deltaTime;
        count += deltaSusceptibles[xAge + 1] * state.intParms.deltaTime;
    }
}


// Given current state, calculates the integral that corresponds to "Beta * I"
// in an ODE SIR model. More details in Stringer et al. paper.
double NewInfections::calculateInfectionCoefficient(State& state) {
    double coefficient = 0;
    double loadBeta;
    double deltaLoad;
    for (size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        loadBeta = betaForLoad(state, state.compParms.columnLoads->at(xLoad)); // higher load -> more shedding
        deltaLoad = state.compParms.deltaInfectionForLoad->at(xLoad); // Integration delta for infection load. Variable b/c log linear.
        for (size_t xAge = 0; xAge < state.compParms.ageSize; xAge++) {
            coefficient += state.infecteds->getIndex(xAge, xLoad) * deltaLoad * loadBeta;
        }
    }
    return (coefficient * state.intParms.deltaTime); // da doens't change by age, so can multiply by it at end. (dt = da)

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
    double deltaArea;
    // Shorthand variables to make code less verbose
    std::vector<double>& loadVec = *state.compParms.columnLoads, deltaSusceptibles = *deltaSusceptibles_->getCurrentState();
    double shapeParam = state.modParms.gammaShapeParam;
    // Counter for debugging
    double totalInfected = 0.0;
    //Zero out min age and min infection since none left after move 
    for (size_t zeros = 0; zeros < state.compParms.infectionSize + state.compParms.ageSize; zeros++){
        if (zeros < state.compParms.infectionSize){
            deltaInfections_->setIndex(0, zeros, 0);
        }else {
            deltaInfections_->setIndex(zeros - state.compParms.infectionSize, 0, 0);
        }
    }
    for(size_t xLoad = 1; xLoad < state.compParms.infectionSize; xLoad++) {
        // Get proportion that lands in previous infection level & advance load by 1...
        double gammaVal = gammaDist(loadVec[xLoad - 1], state.modParms.aveInitInfectionLoad, shapeParam);
        deltaArea = state.intParms.deltaTime * state.compParms.deltaInfectionForLoad->at(xLoad);
        // Loop over infection levels. Skip age 0 -- no infections
        for(size_t xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            double infectedVal = deltaSusceptibles[xAge] * gammaVal;
            totalInfected += infectedVal * deltaArea;
            deltaInfections_->setIndex(xAge, xLoad,  infectedVal);
        }
    }
}

void NewInfections::moveInfecteds(State& state){
    auto& susceptibles = *state.susceptibles->getCurrentState(), deltaSus = *deltaSusceptibles_->getCurrentState();
    Infecteds& infecteds = *state.infecteds;
    Infecteds& deltaInf = *deltaInfections_;
    for(size_t xAge = 1; xAge < state.compParms.ageSize; xAge++) {
        susceptibles[xAge] -= deltaSus[xAge];
        for(size_t xLoad = 1; xLoad < state.compParms.infectionSize; xLoad++) {
            double newDensity = infecteds.getIndex(xAge, xLoad) + deltaInf.getIndex(xAge, xLoad);
            infecteds.setIndex(xAge, xLoad, newDensity);
        }
    }
}


// TEST: count moving infecteds and susceptibles to show equal approx
