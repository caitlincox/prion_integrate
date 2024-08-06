#ifndef INTEGRATOR
#define INTEGRATOR

#include <memory>

#include "birth.h"
#include "death.h"
#include "infecteds.h"
#include "state.h"
#include "susceptibles.h"
#include "new_infections.h"

//Integrate with forward Euler scheme
class Integrator {
public:
    Integrator(
        std::unique_ptr<BirthScheme> births,
        std::unique_ptr<Death> deaths,
        std::unique_ptr<State> state,
        std::unique_ptr<NewInfections> newInfections,
        bool expectConstantPop);
    void run();
    void dryRun(double runs);
    const State& getState() {
        return *state_.get();
    }

private:
    // Update the state to reflect the state after a deltaTime step.
    void timeStep();

    // These are specific to how we integrate over time.
    void computeAgeDeaths();
    void computeInfectionDeaths();
    void computeNewInfectionDeaths();
    void calculateDeltas();
    void adjustSusceptibles();
    void adjustInfecteds();


    //integrate over space-like variables
    double integrateInfecteds(Infecteds& infecteds, State& state);
    double integrateSusceptibles(Susceptibles& susceptibles);// NOT IMPLEMENTED


    std::unique_ptr<BirthScheme> births_;
    std::unique_ptr<Death> deaths_;
    std::unique_ptr<State> state_;
    std::unique_ptr<NewInfections> newInfections_;
    bool expectConstantPop_;
};

#endif
