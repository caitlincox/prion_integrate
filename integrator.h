#ifndef INTEGRATOR
#define INTEGRATOR

#include <memory>

#include "birth.h"
#include "death.h"
#include "infecteds.h"
#include "state.h"
#include "susceptibles.h"

//Integrate with forward Euler scheme
class Integrator {
public:
    Integrator(
        std::unique_ptr<BirthScheme> births,
        std::unique_ptr<Death> deaths,
        std::unique_ptr<State> state,
        bool expectConstantPop);
    void run();
    const State& getState() {
        return *state_.get();
    }

private:
    // Update state_ for the current time step.
    void update();

    std::unique_ptr<BirthScheme> births_;
    std::unique_ptr<Death> deaths_;
    std::unique_ptr<State> state_;
    bool expectConstantPop_;
};

#endif
