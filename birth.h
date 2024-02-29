#ifndef BIRTH
#define BIRTH

#include <cstdio>

//Abstract class since we do population replacement in a few different ways
class BirthScheme {

public:
    virtual double calculateBirth() = 0;
private:
    int l;
};

#endif
