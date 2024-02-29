#ifndef SUSCEPTIBLES
#define SUSCEPTIBLES

//Holds pointer to susceptibles array and allows for state to be updated
class Susceptibles{
public:
    void setInitState(double * state) {currentState = state;}
    double * getCurrentState() {return currentState;}
private:
    double * currentState; //I'm goint to be using memmove a lot, too. What could go wrong?
};
#endif