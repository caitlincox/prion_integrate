#ifndef BIRTH
#define BIRTH

#include "state.h"

#include <cstdio>

//Abstract class since we do population replacement in a few different ways
class BirthScheme {

public:
    virtual double calculateBirth(const State& state) = 0;
};

#endif
