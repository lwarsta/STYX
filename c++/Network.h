#ifndef _NETWORK_H
#define _NETWORK_H
#include "Settings.h"
#include "JuncGeom.h"
#include "JuncWater.h"
#include "LinkGeom.h"

class Network {
private:
    std::vector<JuncGeom> juncs_geom;
    std::vector<JuncWater> juncs_water;
	std::vector<LinkGeom> links_geom;
public:
    void build_network(std::vector<std::vector<std::string>> &tokens);
    void create_water_juncs();
    void init_water_juncs(Settings& settings,
	                      std::vector<std::vector<std::string>>& materials, 
	                      std::vector<std::vector<std::string>>& bound_cond,
						  std::vector<std::vector<std::string>>& init_cond);
    std::vector<JuncGeom>* get_geom_juncs() { return &juncs_geom; }
    std::vector<JuncWater>* get_water_juncs(){ return &juncs_water; }
};

#endif
