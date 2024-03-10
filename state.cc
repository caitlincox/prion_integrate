#include "state.h"

#include <cmath>
#include <cstring>

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
    setStartingDistribution(*this);
    // We represent load buckets on a natural log scale.  This model throws
    // away loads that are very small, which otherwise would show up as animals
    // with longer incubation periods.  To model longer incubation periods,
    // increase the number of infection load buckets.
    //
    // The first bucket at load index 0 has load e^w, the minimum infection
    // load modeled.  We have to compute this bucket's load.  Call 'b' the
    // number of infection load buckets.  When load reaches 1, the animal dies,
    // so the max bucket is the one just before reaching 1.  We can compute 1 =
    // e^w * e^(C*dt*b) => 0 = w + C*dt*b => w = -C*dt*b.  The first bucketk
    // has load e^w.  After that the load is (e^w)*e^(C*dt*bucketIndex),
    // where bucketIndex is on the infection load axis.  At bucketIndex == 0,
    // we have the minimum represented load of e^w.  The bucket index with
    // infection load == 1 is not represented, as the animals would be dead.
    compParms.firstBucketLogLoad = -compParms.intrinsicGrowthRate *
        intParms.deltaTime * intParms.numInfectionLoadBuckets;
    compParms.columnLoads = std::make_unique<std::vector<double>>();
    compParms.columnLoads->resize(intParms.numInfectionLoadBuckets);
    // Compute the infection load for each bucket.
    auto& columnLoads = *compParms.columnLoads;
    // First bucket is e^w.
    columnLoads[0] = exp(compParms.firstBucketLogLoad);
    for (size_t i = 1; i < intParms.numInfectionLoadBuckets; i++) {
        columnLoads[i] = exp(compParms.firstBucketLogLoad + i * intParms.deltaTime *
            compParms.intrinsicGrowthRate);
    }
}

State::State(IntegrationParams integrationParams, ModelParams modelParams) {
    initializeComputedParameters();
    susceptibles = std::make_unique<Susceptibles>(compParms.ageSize);
    infecteds = std::make_unique<Infecteds>(compParms.ageSize, integrationParams.numInfectionLoadBuckets);
    intParms = integrationParams;
    modParms = modelParams;
    setInitialInfecteds(*this);
    updateComputedParameters();
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
    for (size_t xAge = 0; xAge < maxAgeStep; xAge++) {
        compParms.susceptiblePop += susceptiblesVec[xAge];
        for (size_t xLoad = 1; xLoad < intParms.numInfectionLoadBuckets; xLoad++) {
            double popAtLoad = infecteds->getIndex(xAge, xLoad);
            compParms.infectedPop += popAtLoad;
            compParms.transferRate += modParms.beta * columnLoads[xLoad] * popAtLoad;
            totalLoad += columnLoads[xLoad] * popAtLoad;
        }
    }
    compParms.popSize = compParms.susceptiblePop + compParms.infectedPop;
// FIX ME: Should this be popSize instead?
    compParms.aveLoad = totalLoad / compParms.infectedPop;
}

void State::timeStep() {
    // Compute deaths before the memmove.
    compParms.ageDeaths = 0.0;
    size_t xMaxAge = compParms.ageSize - 1;
    for (size_t xLoad = 0; xLoad < intParms.numInfectionLoadBuckets; xLoad++) {
        compParms.ageDeaths += infecteds->getIndex(xMaxAge, xLoad);
    }
    // Compute deaths from infection.
    compParms.infectionDeaths = 0.0;
    size_t xMaxLoad = intParms.numInfectionLoadBuckets - 1;
    // The special case of the bucket that would die both from age and
    // infection is added to the age deaths.
    for (size_t xAge = 0; xAge < xMaxAge; xAge++) {
        compParms.infectionDeaths += infecteds->getIndex(xAge, xMaxLoad);
    }
    // Now advance time by deltaTime.
    // Move all susceptiblees to one higher age index.
    std::vector<double>& susVec = *susceptibles->getCurrentState();
    double* susPtr = &susVec[0];
    memmove(susPtr + 1, susPtr, xMaxAge);
    *susPtr = 0.0;  // Set births to 0.
    // Move all the infectes both by a deltaTime and increase infection.
    // This simply moves data by 1 in both age and infection dimensions.
    double* infectedPtr = infecteds->getInfectionRow(0);
    size_t numMoving = xMaxAge * intParms.numInfectionLoadBuckets - 1;
    memmove(infectedPtr + intParms.numInfectionLoadBuckets + 1, infectedPtr, numMoving);
    // Zero out the min infection buckets.
    for (size_t xAge = 0; xAge < compParms.ageSize; xAge++) {
        infecteds->setIndex(xAge, 0, 0.0);
    }
}
