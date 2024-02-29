#include "integrator.h"

#include <cstdint>
#include <vector>

/*
    double aveLifespan = 5.0;
    double .aveInitInfectionLoad = 0.1;
    double growthRae = intrinsicGrowthRate(aveLifespan, aveInitInfectionLoad);
*/

void Integrator::run() {
    for (uint32_t timeStep = 0; timeStep < integrationParams_.totalTime;
           timeStep += integrationParams_.deltaTime) {
        update();
    }
}

void Integrator::update() {

}
