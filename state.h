#ifndef STATE
#define STATE

#include "infecteds.h"
#include "susceptibles.h"

struct State {
    std::unique_ptr<Susceptibles> susceptibles;
    std::unique_ptr<Infecteds> infecteds;
};

#endif
