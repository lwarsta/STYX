#ifndef _JUNCWATER_H
#define _JUNCWATER_H
#include "JuncGeom.h"
#include "CellBase.h"
#include "LinkWater.h"

class JuncWater : public CellBase
{
protected:
    JuncGeom * geom;
    std::vector<int> ids_link_end;
    std::vector<int> ids_link;
    std::vector<LinkWater*> links;
    std::vector<int> ids_junc;
    std::vector<JuncWater*> juncs;
    double water_depth;
    double water_depth_old;
    int type;
    double outfall_volume;
    std::vector<double> fluxes;
    std::vector<double> velocities;

public:
    JuncWater();
    void save_link(int ind_lnk_end, int ind_lnk, LinkWater* link_water, 
                   int ind_jnc, JuncWater* junc_water);
    void assign_geom(JuncGeom *geom_new){geom = geom_new;}
    JuncGeom * get_geom(){return geom;}
    std::vector<JuncWater*> get_juncs_neigh() { return juncs; }
    std::vector<LinkWater*> get_links() { return links; }
    std::vector<int> get_ids_lnk_end() { return ids_link_end; }
    void set_water_depth(double water_depth_new) { water_depth = water_depth_new; }
    double get_water_depth(){return water_depth;}
    double get_water_depth_old() { return water_depth_old; }
    void swap() { water_depth_old = water_depth; }
    void set_type(int type_new) { type = type_new; }
    int get_type() { return type; }
    double get_outfall_volume() { return outfall_volume; }
    void set_outfall_volume(double outfall_volume_new) { outfall_volume = outfall_volume_new;  }
    void assign_fluxes(std::vector<double> fluxes_new) { fluxes = fluxes_new;  }
    void assign_velocities(std::vector<double> velocities_new) { velocities = velocities_new; }
    std::vector<double> get_fluxes() { return fluxes; }
    std::vector<double> get_velocities() { return velocities; }
};

#endif
