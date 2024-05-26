#include "state.h"
#include "susceptibles.h"
#include "infecteds.h"

#ifndef NEWINFECTIONS
#define NEWINFECTIONS
class NewInfections{
public:
    NewInfections(const State& state);
    void calculateDeltaSusceptibles(State& state);
    void calculateDeltaInfecteds(State& state);
    void moveInfecteds(State& state);
    double calculateInfectionCoefficient(State& state);
    double betaForLoad(State& state, double infectionLoad);
    void prepInfecteds(State& state);
    Infecteds* getInfecteds() {return deltaInfections_.get();}
    Susceptibles* getSusceptibles() {return deltaSusceptibles_.get();}
private:
    std::unique_ptr<Infecteds> deltaInfections_;
    std::unique_ptr<Susceptibles> deltaSusceptibles_;
};
#endif