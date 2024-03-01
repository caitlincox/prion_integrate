#include "integrator.h"

#include <cstdint>
#include <vector>

void Integrator::run() {
    for (double time = 0.0; time < integrationParams_.totalTime;
           time += integrationParams_.deltaTime) {
        update();
    }
}

void Integrator::update() {
    findSums();
}

// We represent buckets on a natrual log scale.  When load reaches 1, the
// animal dies, so we can compute 1 = e^w * e^(C*dt * b) => 0 = w + C*dt*b => w
// = -C*dt*b.  The first bucket has 0 load, but the next has e^w.
void Integrator::findSums() {
    size_t maxAgeStep = modelParams_.maxAge/integrationParams_.deltaTime;
    totalInfection_ = 0.0;
    infectedPop_ = 0.0;
    susceptiblePop_ = 0.0;
    for (size_t ageIndex = 0; ageIndex < maxAgeStep; ageIndex++) {
        susceptiblePop_ += (*state_->susceptibles->getCurrentState())[ageIndex];
        for (size_t infectionIndex = 0; infectionIndex < integrationParams_.numInfectionLoadBuckets; infectionIndex++) {
            infectedPop_ += state_->infecteds->getIndex(ageIndex, infectionIndex);
        }
    }
}
