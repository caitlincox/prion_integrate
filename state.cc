#include "state.h"

namespace {

double findLambda(double aveLifespan, double kappa) {
  printf("Not implemented");
  return 0.0;
}

}

// Initialize computed parameters, including those computed up-front.
void State::initializeComputedParameters() {
    compParms.ageSize = modParms.maxAge / intParms.deltaTime;
    compParms.lambda = findLambda(modParms.aveLifespan, modParms.kappa);
}

State::State(IntegrationParams integrationParams, ModelParams modelParams) {
    initializeComputedParameters();
    susceptibles = std::make_unique<Susceptibles>(compParms.ageSize);
    infecteds = std::make_unique<Infecteds>(compParms.ageSize, integrationParams.numInfectionLoadBuckets);
    intParms = integrationParams;
    modParms = modelParams;
}

void State::updateComputedParameters() {
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
