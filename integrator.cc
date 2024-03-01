#include "integrator.h"

#include <cmath>
#include <cstdint>
#include <vector>

#include "initialconditions.h"

Integrator::Integrator(
    const IntegrationParams& integrationParams,
    const ModelParams& modelParams,
    std::unique_ptr<BirthScheme> births,
    std::unique_ptr<Death> deaths,
    std::unique_ptr<State> state)
{
    integrationParams_ = integrationParams;
    modelParams_ = modelParams;
    births_ = std::move(births);
    state_ = std::move(state);
    intrinsicGrowthRate_ = InitialConditions::intrinsicGrowthRate(
        modelParams.aveLifespan, modelParams.aveInitInfectionLoad);
    // We represent buckets on a natrual log scale.  When load reaches 1,
    // the animal dies, so we can compute 1 = e^w * e^(C*dt * b) => 0 = w +
    // C*dt*b => w = -C*dt*b.  The first bucket has 0 load, but the next
    // has e^w.
    firstBucketLogLoad_ = -intrinsicGrowthRate_ * integrationParams.deltaTime *
        integrationParams.deltaTime * integrationParams.numInfectionLoadBuckets;
    columnLoads_ = std::make_unique<std::vector<double>>();
    columnLoads_->resize(integrationParams_.numInfectionLoadBuckets);
    auto& columnLoads = *columnLoads_;
    columnLoads[0] = 0.0;
    columnLoads[1] = exp(firstBucketLogLoad_);
    for (size_t i = 2; i < integrationParams_.numInfectionLoadBuckets; i++) {
        columnLoads[i] = exp(firstBucketLogLoad_ + (i - 1) * integrationParams_.deltaTime *
            intrinsicGrowthRate_);
    }
}

void Integrator::run() {
    for (double time = 0.0; time < integrationParams_.totalTime;
           time += integrationParams_.deltaTime) {
        update();
    }
}

void Integrator::update() {
    findSums();
}

void Integrator::findSums() {
    size_t maxAgeStep = modelParams_.maxAge/integrationParams_.deltaTime;
    totalInfection_ = 0.0;
    infectedPop_ = 0.0;
    susceptiblePop_ = 0.0;
    transferRate_ = 0.0;
    auto& columnLoads = *columnLoads_;
    for (size_t ageIndex = 0; ageIndex < maxAgeStep; ageIndex++) {
        susceptiblePop_ += (*state_->susceptibles->getCurrentState())[ageIndex];
        for (size_t infectionIndex = 1; infectionIndex < integrationParams_.numInfectionLoadBuckets; infectionIndex++) {
            double popAtLoad = state_->infecteds->getIndex(ageIndex, infectionIndex);
            infectedPop_ += popAtLoad;
            transferRate_ += modelParams_.beta * columnLoads[infectionIndex] * popAtLoad;
        }
    }
}
