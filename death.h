#ifndef DEATH
#define DEATH

class Death {
public:
    Death(double kappa, double aveLifespan);
    double weibullDeathRate(double age);
    double findLambda(double aveLifespan, double kappa);

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
