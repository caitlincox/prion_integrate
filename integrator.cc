#include "integrator.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "tests.h"

Integrator::Integrator(
    std::unique_ptr<BirthScheme> births,
    std::unique_ptr<Death> deaths,
    std::unique_ptr<State> state,
    std::unique_ptr<NewInfections> newInfections,
    bool expectConstantPop)
{
    births_ = std::move(births);
    deaths_ = std::move(deaths);
    state_ = std::move(state);
    newInfections_ = std::move(newInfections);
    expectConstantPop_ = expectConstantPop;
}

// Determine if we will pass a time epoch boundary this time step.  This is
// used to do occasional compute-intensive processing such as writing graphs of
// the current state.
uint32_t computeTimeEpoch(float time, float deltaTime, float totalTime, uint32_t numEpochs) {
    return static_cast<uint32_t>(time * numEpochs / totalTime);
}

void Integrator::run() {
    uint32_t epoch = 0;
    for (double time = 0.0; time < state_->intParms.totalTime;
           time += state_->intParms.deltaTime) {
        state_->updateComputedParameters();
        uint32_t nextEpoch = computeTimeEpoch(time, state_->intParms.deltaTime,
                state_->intParms.totalTime, 100);
        if (nextEpoch != epoch) {
            printf("Writing graph %u\n", epoch);
            std::string epochStr = std::to_string(epoch);
            state_->writeSusceptiblesPBM("data/suseptibles_" + epochStr + ".pbm", 2);
            state_->writeInfectedsPGM("data/infecteds_" + epochStr + ".pgm");
            epoch = nextEpoch;
        }
        runTests(*state_, expectConstantPop_);
        timeStep();
// temp
//        testInfection(*state_, *newInfections_);
    }
}

// We kill off animals at a max age.  This has contributions from both susceptibles and infecteds.
void Integrator::computeAgeDeaths() {
    ComputedParams& compParms = state_->compParms;
    Infecteds* infecteds = state_->infecteds.get();
    std::vector<double>& susVec = *state_->susceptibles->getCurrentState();
    // Initialize age deaths to the oldest suseptibles.
    compParms.ageDeaths = susVec[compParms.ageSize - 1];
    susVec[compParms.ageSize - 1] = 0.0;
    size_t xMaxAge = compParms.ageSize - 1;
    for (size_t xLoad = 0; xLoad < compParms.infectionSize; xLoad++) {
        compParms.ageDeaths += infecteds->getIndex(xMaxAge, xLoad) * compParms.deltaInfectionForLoad->at(xLoad);
    }
    compParms.ageDeaths *= state_->intParms.deltaTime; //Scale by delta time
}

// Compute deaths from infection.
void Integrator::computeInfectionDeaths() {
    ComputedParams& compParms = state_->compParms;
    Infecteds* infecteds = state_->infecteds.get();
    size_t xMaxAge = compParms.ageSize - 1;
    compParms.infectionDeaths = 0.0;
    size_t xMaxLoad = compParms.infectionSize - 1;
    // The special case of the bucket that would die both from age and
    // infection is added to the age deaths.
    for (size_t xAge = 0; xAge < xMaxAge; xAge++) {
        compParms.infectionDeaths += infecteds->getIndex(xAge, xMaxLoad) * compParms.deltaInfectionForLoad->at(xMaxLoad);
    }
    // Note that infections in a column are from 0 to 1, so total area is just deltaTime.
    compParms.infectionDeaths *= state_->intParms.deltaTime; //scale by delta time
}

void Integrator::timeStep() {
    computeAgeDeaths();
    computeInfectionDeaths();
    ComputedParams& compParms = state_->compParms;
    size_t xMaxAge = compParms.ageSize - 1;
    Infecteds* infecteds = state_->infecteds.get();
    std::vector<double>& susVec = *state_->susceptibles->getCurrentState();
    // TODO: Insert call to compute delta infecteds here.  We do this before
    // the step because we're using forward Euler.  Add the delta in after the
    // time step below.    
// temp
//     newInfections_->prepInfecteds(*state_);
    // Now kill off population due to natural deaths.
    deaths_->kill(*state_);
    // Compute births.  Add it in after the time step.
    double births = births_->calculateBirth(*state_);
    //deaths_->kill(*state_); goes here if not replacement
    
    // Now advance time by deltaTime.
    // Move all susceptiblees to one higher age index.
    double* susPtr = &susVec[0];
    memmove(susPtr + 1, susPtr, xMaxAge * sizeof(double));
    // Move all the infectes both by a deltaTime and increase infection.
    // This simply moves data by 1 in both age and infection dimensions.
    double* infectedPtr = infecteds->getInfectionRow(0);
    size_t numMoving = xMaxAge * compParms.infectionSize - 1;
    memmove(infectedPtr + compParms.infectionSize + 1, infectedPtr, numMoving * sizeof(double));
    // Zero out the min infection buckets.
    for (size_t xAge = 0; xAge < compParms.ageSize; xAge++) {
        infecteds->setIndex(xAge, 0, 0.0);
    }
    // Add in births to suseptibles.
    susVec[0] = births;
// temp
double deltaInfInv = 1.0 / state_->compParms.deltaInfection;
for (size_t xAge = 1; xAge < compParms.ageSize; xAge++) {
    for (size_t xLoad = 1; xLoad < compParms.infectionSize; xLoad++) {
        infecteds->setIndex(xAge, xLoad, infecteds->getIndex(xAge, xLoad) * deltaInfInv);
    }
}
    // Add newly infecteds
// temp
//     newInfections_->moveInfecteds(*state_);
}
