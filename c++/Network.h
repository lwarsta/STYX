#ifndef _NETWORK_H
#define _NETWORK_H
#include "Algorithms.h"
#include "Settings.h"
#include "GridBase.h"
#include "JuncGeom.h"
#include "JuncWater.h"
#include "LinkGeom.h"
#include "LinkWater.h"

class Network : GridBase {
private:
    void parseVTKData(
        const std::vector<std::vector<std::string>>& tokenized_data,
        std::vector<std::vector<double>>& points,
        std::vector<std::vector<int>>& cells,
        std::vector<int>& cell_types,
        std::vector<std::string>& field_data_int_names,
        std::vector<std::vector<std::vector<int>>>& field_data_int,
        std::vector<std::string>& field_data_double_names,
        std::vector<std::vector<std::vector<double>>>& field_data_double);
    void printVTKData(std::vector<std::vector<double>>& points,
        std::vector<std::vector<int>>& cells,
        std::vector<int>& cell_types,
        std::vector<std::string>& field_data_int_names,
        std::vector<std::vector<std::vector<int>>>& field_data_int,
        std::vector<std::string>& field_data_double_names,
        std::vector<std::vector<std::vector<double>>>& field_data_double);
    std::vector<JuncGeom> juncs_geom;
    std::vector<JuncWater> juncs_water;
	std::vector<LinkGeom> links_geom;
    std::vector<LinkWater> links_water;
public:
    void build_network(std::vector<std::vector<std::string>> &tokens_junc, 
                       std::vector<std::vector<std::string>> &tokens_link);
    void create_water_network_items();
    void init_water_network(Settings& settings,
	                        std::vector<std::vector<std::string>>& materials_net,
	                        std::vector<std::vector<std::string>>& bound_cond_net,
						    std::vector<std::vector<std::string>>& init_cond_net);
    std::vector<JuncGeom>* get_geom_juncs() { return &juncs_geom; }
    std::vector<JuncWater>* get_water_juncs(){ return &juncs_water; }
};

#endif
