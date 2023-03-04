#ifndef _ATMOSCONTROL_H
#define _ATMOSCONTROL_H
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Add computation of evapotranspiration here.
class AtmosControl {
private:
    std::vector<std::string> dataHeaders;
    std::vector<std::string> dateTime; // change this to decimal date
    std::vector < std::vector<double> > atmosForcingData;
    double time;
    double precipFact;
    double petFact;
public:
    AtmosControl();
    void set_precip_fact(double precipFactNew){precipFact = precipFactNew;}
    void set_pet_fact(double petFactNew){petFact = petFactNew;}
    void set_atmos_forcing_data(std::vector < std::vector<std::string> > atmosForcingDataStr);
    void print_atmos_forcing_data();
    void advance_time(double timeStep){time += timeStep;}
    double get_precip();
    double get_pet();
    double get_air_temp();
};

#endif

