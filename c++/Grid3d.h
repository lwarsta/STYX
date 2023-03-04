#ifndef _GRID3D_H
#define _GRID3D_H
#include "Settings.h"
#include "GridBase.h"
#include "CellGeom3d.h"
#include "CellWater3d.h"
#include "CellHeat3d.h"
#include "CellSolute3d.h"

class Grid3d : GridBase {
private:
    std::vector<CellGeom3d> cells_geom;
    std::vector<CellWater3d> cells_water;
    std::vector<CellHeat3d> cells_heat;
    std::vector<CellSolute3d> cells_solute;
public:
    void build_grid(std::vector < std::vector<std::string> > &tokens);
    void create_water_cells();
    void create_heat_cells();
    void create_solute_cells(size_t numOfSolutes);
    std::string parse_unstruct_vtk_mesh();
    void init_water_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond);
    void init_heat_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond);
    void init_solute_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond, std::vector < std::vector<std::string> >& soluteLib, std::vector <std::string>& species_names, std::vector<double> & conc_init);
    std::vector<std::string> find_solute_prop(std::string soluteName, std::vector<std::vector<std::string> > soluteLib); // MOVE TO COMMON PARENT
    std::vector<CellGeom3d> * get_geom_cells(){return &cells_geom;}
    std::vector<CellWater3d> * get_water_cells(){return &cells_water;}
    std::vector<CellHeat3d>* get_heat_cells() { return &cells_heat; }
    std::vector<CellSolute3d> * get_solute_cells(){return &cells_solute;}
};

#endif
