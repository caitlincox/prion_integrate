#include "state.h"

void State::findSums() {
    size_t maxAgeStep = modParms.maxAge/intParms.deltaTime;
    compParms.totalInfection = 0.0;
    compParms.infectedPop = 0.0;
    compParms.susceptiblePop = 0.0;
    compParms.transferRate = 0.0;
    auto& columnLoads = *compParms.columnLoads;
    for (size_t ageIndex = 0; ageIndex < maxAgeStep; ageIndex++) {
        compParms.susceptiblePop += (*susceptibles->getCurrentState())[ageIndex];
        for (size_t infectionIndex = 1; infectionIndex < intParms.numInfectionLoadBuckets; infectionIndex++) {
            double popAtLoad = infecteds->getIndex(ageIndex, infectionIndex);
            compParms.infectedPop += popAtLoad;
            compParms.transferRate += modParms.beta * columnLoads[infectionIndex] * popAtLoad;
        }
    }
    compParms.popSize = compParms.susceptiblePop + compParms.infectedPop;
}
