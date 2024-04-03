#ifndef DEATH
#define DEATH

#include "state.h"

class Death {
public:
    Death(const State& state);
    double weibullDeathRate(double age);
    // Kill off population according to death rate for all ages.
    void kill(State& state);

    //getters and setters
    double getSurvivalKappa() {return weibullSurvivorshipKappa_;}
    double getSurvivalLambda() {return weibullSurvivorshipLambda_;}
    void setSurvivalKappa(double kappa) {weibullSurvivorshipKappa_ = kappa;}
    void setSurvivalLambda(double lambda) {weibullSurvivorshipLambda_ = lambda;}

private:
    double weibullSurvivorshipLambda_;
    double weibullSurvivorshipKappa_;
};

#endif
