#include "integrator.h"

#include <cmath>
#include <cstdint>
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
// state_->writeInfectedsPGM("initial_infecteds.pgm");
        runTests(*state_, expectConstantPop_);
        update();
    }
}

void Integrator::update() {
    state_->timeStep();
    state_->updateComputedParameters();
    double births = births_->calculateBirth(*state_);
}
