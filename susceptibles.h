#ifndef SUSCEPTIBLES
#define SUSCEPTIBLES

#include <memory>
#include <vector>

//Holds pointer to susceptibles array and allows for state to be updated
class Susceptibles{
public:
    Susceptibles(std::unique_ptr<std::vector<double>> current_state) {
      currentState_ = std::move(current_state);
    }
    std::vector<double>* getCurrentState() {return currentState_.get();}

private:
    std::unique_ptr<std::vector<double>> currentState_;
};

#endif
