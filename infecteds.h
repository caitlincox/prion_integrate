#ifndef INFECTEDS
#define INFECTEDS

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>

class Infecteds {
public:
    Infecteds(std::unique_ptr<std::vector<double>> currentState, size_t ageSize, size_t infectionSize) {
        if (currentState->size() != ageSize * infectionSize) {
            printf("Invalid table size!");
            exit(1);
        }
        currentState_ = std::move(currentState);
        ageSize_ = ageSize;
        infectionSize_ = infectionSize;
    }
    std::vector<double>* getCurrentState() {return currentState_.get();}
    double getIndex(size_t ageIndex, size_t infectionIndex) {
        return (*currentState_)[ageIndex * infectionSize_ + infectionIndex];
    }
    void setIndex(size_t ageIndex, size_t infectionIndex, double value) {
        (*currentState_)[ageIndex * infectionSize_ + infectionIndex] =  value;
    }
    double* getInfectionRow(size_t ageIndex) {
        return &(*currentState_)[ageIndex * infectionSize_];
    }

private:
    std::unique_ptr<std::vector<double>> currentState_;
    size_t ageSize_;
    size_t infectionSize_;
};

#endif
