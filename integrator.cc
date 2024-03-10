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
}

void Integrator::run() {
    for (double time = 0.0; time < state_->intParms.totalTime;
           time += state_->intParms.deltaTime) {
        update();
    }
}

void Integrator::update() {
    state->timeStep();
    state_->updateComputedParameters();
}
