#ifndef _GRID2D_H
#define _GRID2D_H
#include "Settings.h"
#include "GridBase.h"
#include "CellGeom2d.h"
#include "CellWater2d.h"
#include "CellHeat2d.h"
#include "CellSolute2d.h"

class Grid2d : GridBase {
private:
    std::vector<CellGeom2d> cells_geom;
    std::vector<CellWater2d> cells_water;
    std::vector<CellHeat2d> cells_heat;
    std::vector<CellSolute2d> cells_solute;
public:
    void build_grid(std::vector < std::vector<std::string> > &tokens);
    void create_water_cells();
    void create_heat_cells();
    void create_solute_cells(size_t numOfSolutes);
    std::string parse_unstruct_vtk_mesh();
    void init_water_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond);
    void init_heat_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond);
    void init_solute_cells(Settings& settings, std::vector < std::vector<std::string> >& materials, std::vector < std::vector<std::string> >& bound_cond, std::vector < std::vector<std::string> >& init_cond, std::vector < std::vector<std::string> >& solute_lib, std::vector <std::string>& species_names);
    std::vector<std::string> find_solute_prop(std::string soluteName, std::vector<std::vector<std::string> > soluteLib); // MOVE TO COMMON PARENT
    std::vector<CellGeom2d>* get_geom_cells() { return &cells_geom; }
    std::vector<CellWater2d> * get_water_cells(){ return &cells_water; }
    std::vector<CellHeat2d>* get_heat_cells() { return &cells_heat; }
    std::vector<CellSolute2d> * get_solute_cells(){ return &cells_solute; }
};

#endif
