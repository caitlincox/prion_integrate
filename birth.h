#ifndef BIRTH
#define BIRTH

//Abstract class since we do population replacement in a few different ways
class BirthScheme {

public:
    virtual double calculateBirth(){};
private:
    int l;
};

//Stupidest possible thing. Useful for debugging/checking intuition. Constant birth not perportional to population.
class ConstantBirthScheme : public BirthScheme {

public:
    double calculateBirth(){};

    //getters and setters
    double getBirthConstant() {return birthConstant;}
    void setBirthConstant(double value) {birthConstant = value;}

private:
    double birthConstant;
};


//SIR-inspired constant population. Sets birth equal to death. Was used in Stringer et al.
class ReplacementBirthScheme : public BirthScheme {

public:
    double calculateBirth(){};
};


//Biologically more accurate birth, but we don't change the birthrate by age.
class BirthRateBirthScheme : public BirthScheme {

public:
    double calculateBirth(){};

    //getters and setters
    double getBirthRate() {return birthRate;}
    void setBirthRate(double value) {birthRate = value;}
    
private:
    double birthRate;
};
#endif