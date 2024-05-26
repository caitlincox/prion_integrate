#include "susceptibles.h"
#include "initialconditions.h"

#include <cassert>

void Susceptibles::setBirths(double births) {
  assert((*currentState_)[0] == 0.0);
  (*currentState_)[0] = births;
}


//infect
//oh shoot this is why the normalizer is bad
void infect(const State& state) {
    auto& infecteds = *state.infecteds;
    std::vector<double>& loadVec = *state.compParms.columnLoads;
    double shapeParam = state.modParms.gammaShapeParam;
    std::vector<double>& susceptibles = *state.susceptibles->getCurrentState();
    size_t numSusceptibles = state.modParms.initialSusceptiblePop;
    size_t numInfecteds = state.modParms.initialInfectedPop;
    double totalInfected = 0.0;
    double deltaInfection = state.compParms.deltaLogInfection;
    double deltaArea = state.intParms.deltaTime * deltaInfection;

    for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        // Get proportion that lands in this infection level
        // TODO: The gamma distribution function needs work: It needs more
        // parameters to cause the total infecteds to add to numInfecteds.
        double gammaVal = gammaDist(loadVec[xLoad], state.modParms.aveInitInfectionLoad, shapeParam);

        // Loop over infection levels. Skip age 0 -- no infections
        for(size_t xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            // Get normalized density of susceptibles at this age.
            double weight = susceptibles[xAge] / numSusceptibles;
            // Get density of new infecteds at this bucket & age.
            double infectedVal = weight * gammaVal * numInfecteds;
            totalInfected += infectedVal * deltaArea;
            infecteds.setIndex(xAge, xLoad,  infectedVal);
        }
    }
    // Normalize the initial infections to get the correct total infected population.
    double normalizer = state.modParms.initialInfectedPop / totalInfected;
    printf("Initial infected normalizer = %f\n", normalizer);
    for(size_t xLoad = 0; xLoad < state.compParms.infectionSize; xLoad++) {
        for(size_t xAge = 1; xAge < state.compParms.ageSize; xAge++) {
            infecteds.setIndex(xAge, xLoad,  normalizer * infecteds.getIndex(xAge, xLoad));
        }
    }
}
