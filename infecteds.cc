#include "infecteds.h"
    void Infecteds::rescaleInfecteds(double scalingFactor) {
        for (int xLoad = 0; xLoad < infectionSize_; xLoad++) {
            for (int xAge = 0; xAge < ageSize_; xAge++) {
                double oldVal = getIndex(xAge, xLoad);
                setIndex(xAge, xLoad, oldVal * scalingFactor);
            }
        }
    }
