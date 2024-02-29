#include "initialconditions.h"
#include <cmath>
    double InitialConditions::intrinsicGrowthRate(double aveLifespan, double aveInitInfectionLoad) {
        return log(1/aveInitInfectionLoad)/aveLifespan;
    }
