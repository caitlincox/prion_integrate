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
    bool expectConstantPop)
{
    births_ = std::move(births);
    deaths_ = std::move(deaths);
    state_ = std::move(state);
    expectConstantPop_ = expectConstantPop;
}

void Integrator::run() {
    for (double time = 0.0; time < state_->intParms.totalTime;
           time += state_->intParms.deltaTime) {
// temp
state_->writeSusceptiblesPBM("initial_suseptibles.pbm", 2);
state_->writeInfectedsPGM("initial_infecteds.pgm");
        runTests(*state_, expectConstantPop_);
        update();
    }
}

void Integrator::update() {
    state_->updateComputedParameters();
    timeStep();
}

void Integrator::timeStep() {
    // Compute deaths before the memmove.
    ComputedParams& compParms = state_->compParms;
    Infecteds* infecteds = state_->infecteds.get();
    std::vector<double>& susVec = *state_->susceptibles->getCurrentState();
    // Initialize age deaths to the oldest suseptibles.
    compParms.ageDeaths = susVec[compParms.ageSize - 1];
    size_t xMaxAge = compParms.ageSize - 1;
    for (size_t xLoad = 0; xLoad < compParms.infectionSize; xLoad++) {
        compParms.ageDeaths += infecteds->getIndex(xMaxAge, xLoad); //NEED TO SCALE BY DELTA T - NOT DONE
    }
    // Compute deaths from infection.
    compParms.infectionDeaths = 0.0;
    size_t xMaxLoad = compParms.infectionSize - 1;
    // The special case of the bucket that would die both from age and
    // infection is added to the age deaths.
    for (size_t xAge = 0; xAge < xMaxAge; xAge++) {
        compParms.infectionDeaths += infecteds->getIndex(xAge, xMaxLoad);
    }

    // TODO: Insert call to compute delta infecteds here.  We do this before
    // the step because we're using forward Euler.  Add the delta in after the
    // time step below.

    // Now kill off population due to natural deaths.
    deaths_->kill(*state_);
    // Compute births.  Add it in after the time step.
    double births = births_->calculateBirth(*state_);
    //deaths_->kill(*state_);
    
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
}
