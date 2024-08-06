#ifndef DEATH
#define DEATH

#include "state.h"

class Death {
public:
    Death(const State& state);
    virtual double weibullDeathRate(double age);
    // Kill off population according to death rate for all ages.
    void kill(State& state);

    //getters and setters
    double getSurvivalKappa() {return weibullSurvivorshipKappa_;}
    double getSurvivalLambda() {return weibullSurvivorshipLambda_;}
    void setSurvivalKappa(double kappa) {weibullSurvivorshipKappa_ = kappa;}
    void setSurvivalLambda(double lambda) {weibullSurvivorshipLambda_ = lambda;}
    double killSusceptibles(double ageInYears, State& state, size_t index);


protected:
    double weibullSurvivorshipLambda_;
    double weibullSurvivorshipKappa_;

    double killInfecteds(double ageInYears, State& state, size_t index);
};

class PredatorScheme : public Death {
public:
    PredatorScheme(const State& state, double maxAge, double deathConst) : Death(state){
        maxAge_ = maxAge;
        deathConst_ = deathConst;
    };
    double weibullDeathRate(double age) override;
    double maxAge_;
    double deathConst_;
};

class sexScheme : public Death { 
public:
    sexScheme(const State& state, double femRate, double malRate) : Death(state){
        femRate_ =femRate;
        malRate_ = malRate;
    };
    double weibullDeathRate(double age) override;
private:
    double femRate_;
    double malRate_;
};
#endif
