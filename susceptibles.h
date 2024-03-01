#ifndef SUSCEPTIBLES
#define SUSCEPTIBLES

#include <memory>
#include <vector>

//Holds pointer to susceptibles array and allows for state to be updated
class Susceptibles{
public:
    Susceptibles(size_t ageSize) {
        currentState_ = std::make_unique<std::vector<double>>();
        currentState_->resize(ageSize);
    }
    std::vector<double>* getCurrentState() {return currentState_.get();}

private:
    std::unique_ptr<std::vector<double>> currentState_;
};

#endif
