#include "integrator.h"

#include <cmath>
#include <cstdint>
#include <vector>

#include "initialconditions.h"

Integrator::Integrator(
    std::unique_ptr<BirthScheme> births,
    std::unique_ptr<Death> deaths,
    std::unique_ptr<State> state)
{
    births_ = std::move(births);
    deaths_ = std::move(deaths);
    state_ = std::move(state);
    state_->compParms.intrinsicGrowthRate = intrinsicGrowthRate(
        state_->modParms.aveLifespan, state_->modParms.aveInitInfectionLoad);
    // We represent buckets on a natrual log scale.  When load reaches 1,
    // the animal dies, so we can compute 1 = e^w * e^(C*dt * b) => 0 = w +
    // C*dt*b => w = -C*dt*b.  The first bucket has 0 load, but the next
    // has e^w.
    state_->compParms.firstBucketLogLoad = -state_->compParms.intrinsicGrowthRate * state_->intParms.deltaTime *
        state_->intParms.deltaTime * state_->intParms.numInfectionLoadBuckets;
    state_->compParms.columnLoads = std::make_unique<std::vector<double>>();
    state_->compParms.columnLoads->resize(state_->intParms.numInfectionLoadBuckets);
    auto& columnLoads = *state_->compParms.columnLoads;
    columnLoads[0] = 0.0;
    columnLoads[1] = exp(state_->compParms.firstBucketLogLoad);
    for (size_t i = 2; i < state_->intParms.numInfectionLoadBuckets; i++) {
        columnLoads[i] = exp(state_->compParms.firstBucketLogLoad + (i - 1) * state_->intParms.deltaTime *
            state_->compParms.intrinsicGrowthRate);
    }
}

void Integrator::run() {
    for (double time = 0.0; time < state_->intParms.totalTime;
           time += state_->intParms.deltaTime) {
        update();
    }
}

void Integrator::update() {
    findSums();
}

void Integrator::findSums() {
    size_t maxAgeStep = state_->modParms.maxAge/state_->intParms.deltaTime;
    auto& compParms = state_->compParms;
    compParms.totalInfection = 0.0;
    compParms.infectedPop = 0.0;
    compParms.susceptiblePop = 0.0;
    compParms.transferRate = 0.0;
    auto& columnLoads = *compParms.columnLoads;
    for (size_t ageIndex = 0; ageIndex < maxAgeStep; ageIndex++) {
        compParms.susceptiblePop += (*state_->susceptibles->getCurrentState())[ageIndex];
        for (size_t infectionIndex = 1; infectionIndex < state_->intParms.numInfectionLoadBuckets; infectionIndex++) {
            double popAtLoad = state_->infecteds->getIndex(ageIndex, infectionIndex);
            compParms.infectedPop += popAtLoad;
            compParms.transferRate += state_->modParms.beta * columnLoads[infectionIndex] * popAtLoad;
        }
    }
}
