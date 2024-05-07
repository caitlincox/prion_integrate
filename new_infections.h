#include "state.h"
#include "susceptibles.h"
#include "infecteds.h"

#ifndef NEWINFECTIONS
#define NEWINFECTIONS
class NewInfections{
public:
    NewInfections(State& state);
    void calculateNewInfections(State& state);
    void subtractFromSusceptibles(State& state);
    void addToInfecteds(State& state);
private:
    std::unique_ptr<Infecteds> deltaInfections_;
    std::unique_ptr<Susceptibles> deltaSusceptibles_;
};
#endif