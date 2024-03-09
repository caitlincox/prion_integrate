#include "state.h"

#include <cmath>

#include "initialconditions.h"

namespace {

double findLambda(double aveLifespan, double kappa) {
  printf("Not implemented");
  return 0.0;
}

}

// Initialize computed parameters, including those computed up-front.
void State::initializeComputedParameters() {
    compParms.intrinsicGrowthRate = intrinsicGrowthRate(
        modParms.aveLifespan, modParms.aveInitInfectionLoad);
    compParms.ageSize = modParms.maxAge / intParms.deltaTime;
    compParms.lambda = findLambda(modParms.aveLifespan, modParms.kappa);
    // We represent load buckets on a natural log scale.  This model throws
    // away loads that are very small, which otherwise would show up as animals
    // with longer incubation periods.  To model longer incubation periods,
    // increase the number of infection load buckets.
    //
    // The first bucket we force to 0, and the minimum infection load supported is
    // the value in the next column (at index 1).  We have to compute this
    // bucket's load.  Call 'b' the number of infection load buckets.  When
    // load reaches 1, the animal dies, so the max bucket is the one just
    // before reaching 1.  We can compute 1 = e^w * e^(C*dt*(b-1)) => 0 = w +
    // C*dt*(b-1) => w = -C*dt*(b-1).  The first bucket we force to be 0 load,
    // but the next has e^w.  After that the load is
    // (e^w)*e^(C*dt*(bucketIndex-1)), where bucketIndex is on the infection
    // load axis.  We make bucketIndex == 0 a special case where load == 0, and
    // pop == 0, and at bucketIndex == 1, we have the minimum represented load
    // of e^w.
    compParms.firstBucketLogLoad = -compParms.intrinsicGrowthRate *
        intParms.deltaTime * (intParms.numInfectionLoadBuckets - 1);
    compParms.columnLoads = std::make_unique<std::vector<double>>();
    compParms.columnLoads->resize(intParms.numInfectionLoadBuckets);
    // Compute the infection load for each bucket.
    auto& columnLoads = *compParms.columnLoads;
    // First bucket is always 0.
    columnLoads[0] = 0.0;
    // Second bucket is e^w.
    columnLoads[1] = exp(compParms.firstBucketLogLoad);
    for (size_t i = 2; i < intParms.numInfectionLoadBuckets; i++) {
        columnLoads[i] = exp(compParms.firstBucketLogLoad + (i - 1) * intParms.deltaTime *
            compParms.intrinsicGrowthRate);
    }
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
    compParms.aveLoad = 0.0;
    double totalLoad = 0.0;
    auto& columnLoads = *compParms.columnLoads;
    auto& susceptiblesVec = *susceptibles->getCurrentState();
    for (size_t ageIndex = 0; ageIndex < maxAgeStep; ageIndex++) {
        compParms.susceptiblePop += susceptiblesVec[ageIndex];
        for (size_t infectionIndex = 1; infectionIndex < intParms.numInfectionLoadBuckets; infectionIndex++) {
            double popAtLoad = infecteds->getIndex(ageIndex, infectionIndex);
            compParms.infectedPop += popAtLoad;
            compParms.transferRate += modParms.beta * columnLoads[infectionIndex] * popAtLoad;
            totalLoad += columnLoads[infectionIndex] * popAtLoad;
        }
    }
    compParms.popSize = compParms.susceptiblePop + compParms.infectedPop;
// FIX ME: Should this be popSize instead?
    compParms.aveLoad = totalLoad / compParms.infectedPop;
}
