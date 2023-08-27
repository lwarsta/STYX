#include "Grid2d.h"

void Grid2d::build_grid(std::vector < std::vector<std::string> > &tokens)
{
    // Create a stream of tokens from loaded data.
    int itemCount = 0;
    bool pointsFound = false;
    int numOfPoints = 0;
    std::vector<double> coords;
    bool cellsFound = false;
    int numOfCells = 0;
    std::vector<int> vertIndices;
    bool cellTypesFound = false;
    int numOfCellTypes = 0;
    int cellInd = 0;
    int numOfFieldData = 0;
    bool fieldDataFound = false;
    int fieldCount = 0;
    int fieldItemCount = 0;
    std::string fieldName = "";
    int fieldDataCols = 0;
    std::string fieldDataType = "";
    std::vector<int> dataInt;
    std::vector<double> dataReal;
    bool processingField = false;
    for (size_t i = 0; i < tokens.size(); i++) {
        for (size_t j = 0; j < tokens.at(i).size(); j++) {
            if (tokens.at(i).at(j) == "POINTS" || pointsFound == true) {
                if (itemCount == 0) {
                    pointsFound = true;
                }
                else if (itemCount == 1) {
                    numOfPoints = atoi(tokens.at(i).at(j).c_str());
                }
                else if (itemCount > 2) {
                    coords.push_back(atof(tokens.at(i).at(j).c_str()));
                    if (coords.size() == 3) {
                        // warning C4267: 'argument': conversion from 'size_t' to 'int', possible loss of data
                        Vertex point((int)vertices.size(), coords.at(0), coords.at(1), coords.at(2));
                        coords.clear();
                        vertices.push_back(point);
                        if ((int)vertices.size() == numOfPoints) {
                            pointsFound = false;
                        }
                    }
                }
                if (pointsFound == true) {
                    itemCount++;
                }
                else {
                    itemCount = 0;
                }
            }
            else if (tokens.at(i).at(j) == "CELLS" || cellsFound == true) {
                if (itemCount == 0) {
                    cellsFound = true;
                }
                else if (itemCount == 1) {
                    numOfCells = atoi(tokens.at(i).at(j).c_str());
                }
                else if (itemCount > 2) {
                    vertIndices.push_back(atoi(tokens.at(i).at(j).c_str()));
                    if ((int)vertIndices.size() == (int)vertIndices.at(0) + 1) {
                        // Remove the number of points item in cell vertex data.
                        vertIndices.erase(vertIndices.begin());
                        std::vector<Vertex*> vertexPointers;
                        for (size_t l = 0; l < vertIndices.size(); l += 1) {
                            vertexPointers.push_back(&vertices.at(vertIndices.at(l)));
                        }
                        // Cell index and cell type are saved later in the process.
                        CellGeom2d cellGeom = CellGeom2d(-1, -1, vertIndices, vertexPointers);
                        cells_geom.push_back(cellGeom);
                        vertIndices.clear();
                        if ((int)cells_geom.size() == numOfCells) {
                            cellsFound = false;
                        }
                    }
                }
                if (cellsFound == true) {
                    itemCount++;
                }
                else {
                    itemCount = 0;
                }
            }
            else if (tokens.at(i).at(j) == "CELL_TYPES" || cellTypesFound == true) {
                if (itemCount == 0) {
                    cellTypesFound = true;
                }
                else if (itemCount == 1) {
                    numOfCellTypes = atoi(tokens.at(i).at(j).c_str());
                }
                else if (itemCount > 1) {
                    int cellType = atoi(tokens.at(i).at(j).c_str());
                    cells_geom.at(cellInd).setGeomType(cellType);
                    cellInd++;
                    if (cellInd == numOfCellTypes) {
                        cellTypesFound = false;
                    }
                }
                if (cellTypesFound == true) {
                    itemCount++;
                }
                else {
                    cellInd = 0;
                    itemCount = 0;
                }
            }
            else if (tokens.at(i).at(j) == "CELL_DATA") {
            }
            else if (tokens.at(i).at(j) == "FIELD" || fieldDataFound == true) {
                if (itemCount == 0) {
                    fieldDataFound = true;
                }
                else if (itemCount == 2) {
                    numOfFieldData = atoi(tokens.at(i).at(j).c_str());
                }
                else if (itemCount > 2) {
                    if (fieldItemCount == 0) {
                        fieldName = tokens.at(i).at(j);
                        processingField = true;
                    }
                    else if (fieldItemCount == 1) {
                        fieldDataCols = atoi(tokens.at(i).at(j).c_str());
                    }
                    else if (fieldItemCount == 3) {
                        fieldDataType = tokens.at(i).at(j);
                    }
                    else if (fieldItemCount > 3) {
                        // field column count should be here.
                        if (fieldDataType == "int") {
                            dataInt.push_back(atoi(tokens.at(i).at(j).c_str()));
                        }
                        else if (fieldDataType == "double") {
                            dataReal.push_back(atof(tokens.at(i).at(j).c_str()));
                        }
                        // Save field data into cells.
                        if (fieldName == "id" && (int)dataInt.size() == fieldDataCols) {
                            dataInt.clear();
                            cellInd++;
                        }
                        else if (fieldName == "grid_connection" && (int)dataInt.size() == fieldDataCols) {
                            cells_geom.at(cellInd).setGridConnection(dataInt.at(0));
                            dataInt.clear();
                            cellInd++;
                        }
                        else if (fieldName == "material" && (int)dataInt.size() == fieldDataCols) {
                            cells_geom.at(cellInd).setMaterial(dataInt.at(0));
                            dataInt.clear();
                            cellInd++;
                        }
                        else if (fieldName == "init_cond" && (int)dataInt.size() == fieldDataCols) {
                            cells_geom.at(cellInd).setInitCondInd(dataInt.at(0));
                            dataInt.clear();
                            cellInd++;
                        }
                        else if (fieldName == "bound_cond" && (int)dataInt.size() == fieldDataCols) {
                            cells_geom.at(cellInd).setBoundCondInd(dataInt.at(0));
                            dataInt.clear();
                            cellInd++;
                        }
                        else if ((int)dataInt.size() == fieldDataCols || (int)dataReal.size() == fieldDataCols) {
                            dataInt.clear();
                            dataReal.clear();
                            cellInd++;
                        }
                        // Move to the next field.
                        if (cellInd == numOfCells) { // should numOfCells be locally retrieved too?
                            fieldCount++;
                            fieldItemCount = 0;
                            cellInd = 0;
                            processingField = false;
                        }
                        if (fieldCount == numOfFieldData) {
                            fieldDataFound = false;
                        }
                    }
                    if (processingField == true) {
                        fieldItemCount++;
                    }
                }
                if (fieldDataFound == true) {
                    itemCount++;
                }
                else {
                    cellInd = 0; // is this needed here?
                    itemCount = 0;
                }
            }
        }
    }

    // Compute cell properties.
    for (size_t i = 0; i < cells_geom.size(); i++) {
        cells_geom.at(i).setId(i);
        cells_geom.at(i).computeCentrePoint();
        cells_geom.at(i).assignSideVertIndices();
        cells_geom.at(i).computeSideCentrePoints();
        cells_geom.at(i).computeSideLengths();
        cells_geom.at(i).computeAverageSlope();
        cells_geom.at(i).computeArea();
    }

    // Create a map of sorted cell side vertex indices.
    std::map< int, std::vector<std::vector<int> > > mapCellSideVertices;
    std::map< int, std::vector<int> > mapNeighIndices;
    std::map< std::vector<int>, std::vector<int> > mapSideCell;
    for (size_t i = 0; i < cells_geom.size(); i++) {
        std::vector<std::vector<int> > locSideVertIndices = cells_geom.at(i).getSideVertIndices();
        std::vector<int> vertIndices = cells_geom.at(i).getVertIndices();
        std::vector<std::vector<int> > sideVertIndices;
        for (size_t j = 0; j < locSideVertIndices.size(); j++) {
            std::vector<int> sideVerts;
            for (size_t k = 0; k < locSideVertIndices.at(j).size(); k++) {
                sideVerts.push_back(vertIndices.at(locSideVertIndices.at(j).at(k)));
            }
            std::sort(sideVerts.begin(), sideVerts.end());
            if (mapSideCell.find(sideVerts) == mapSideCell.end()) {
                std::vector<int> sideCell{ (int)i };
                mapSideCell[sideVerts] = sideCell;
            }
            else {
                mapSideCell[sideVerts].push_back(i);
            }
            sideVertIndices.push_back(sideVerts);
        }
        mapCellSideVertices[i] = sideVertIndices;
        std::vector<int> neighIndices(locSideVertIndices.size(), -1);
        mapNeighIndices[i] = neighIndices;
    }
    //
    for (size_t i = 0; i < cells_geom.size(); i++) {
        std::vector<std::vector<int> > sideVerts = mapCellSideVertices[i];
        for (size_t j = 0; j < sideVerts.size(); j++) {
            std::vector<int> sideCells = mapSideCell[sideVerts.at(j)];
            for (size_t k = 0; k < sideCells.size(); k++) {
                if (sideCells.at(k) != (int)i) {
                    mapNeighIndices[i].at(j) = sideCells.at(k);
                    //break;
                }
            }
        }
        // Save indices and pointers to cell.
        std::vector<CellGeom2d*> neighbourPointers;
        for (size_t j = 0; j < mapNeighIndices[i].size(); j++) {
            if (mapNeighIndices[i].at(j) != -1) {
                neighbourPointers.push_back(&cells_geom.at(mapNeighIndices[i].at(j)));
            }
            else {
                neighbourPointers.push_back(0);
            }
        }
        cells_geom.at(i).setNeighIndices(mapNeighIndices[i]);
        cells_geom.at(i).setNeighPointers(neighbourPointers);
        cells_geom.at(i).computeDistancesBetweenCells();
    }
}

void Grid2d::create_water_cells()
{
    cells_water.resize(cells_geom.size());
    for (size_t i = 0; i < cells_geom.size(); i++) {
        cells_water.at(i).setId(i);
        // Connect water cell to a geometry cell.
        cells_water.at(i).assignGeom(&cells_geom.at(i));
        // Connect water cell to water cells.
        std::vector<CellWater2d*> neighbourPointers;
        std::vector<int> neighIndices = cells_geom.at(i).getNeighIndices();
        for (size_t j = 0; j < neighIndices.size(); j++) {
            if (neighIndices.at(j) != -1) {
                neighbourPointers.push_back(&cells_water.at(neighIndices.at(j)));
            }
            else {
                neighbourPointers.push_back(0);
            }
        }
        cells_water.at(i).setNeighPointers(neighbourPointers);
    }
}

void Grid2d::create_heat_cells()
{
    cells_heat.resize(cells_geom.size());
    for (size_t i = 0; i < cells_geom.size(); i++)
    {
        cells_heat.at(i).setId(i);
        // Connect heat cell to a geometry cell.
        cells_heat.at(i).assignGeom(&cells_geom.at(i));
        // Connect heat cell to water cell.
        cells_heat.at(i).assignWater(&cells_water.at(i));
        // Connect heat cell to heat cells.
        std::vector<CellHeat2d*> neighbourPointers;
        std::vector<int> neighIndices = cells_geom.at(i).getNeighIndices();
        for (size_t j = 0; j < neighIndices.size(); j++) {
            if (neighIndices.at(j) != -1) {
                neighbourPointers.push_back(&cells_heat.at(neighIndices.at(j)));
            }
            else {
                neighbourPointers.push_back(0);
            }
        }
        cells_heat.at(i).setNeighPointers(neighbourPointers);
    }
}

void Grid2d::create_solute_cells(size_t numOfSolutes)
{
    cells_solute.resize(cells_geom.size() * numOfSolutes);
    
    for (int i = 0; i < numOfSolutes; i++)
    {
        for (size_t j = 0; j < cells_geom.size(); j++)
        {
            int index = i * cells_geom.size() + j;
            cells_solute.at(index).setId(index);
            // Connect solute cell to a geometry cell.
            cells_solute.at(index).assignGeom(&cells_geom.at(j));
            // Connect solute cell to water cell.
            cells_solute.at(index).assignWater(&cells_water.at(j));
            // Connect solute cell to solute cells.
            std::vector<CellSolute2d*> neighbourPointers;
            std::vector<int> neighIndices = cells_geom.at(j).getNeighIndices();
            for (size_t k = 0; k < neighIndices.size(); k++) {
                if (neighIndices.at(k) != -1) {
                    int indexNeigh = neighIndices.at(k) + i * cells_geom.size();
                    neighbourPointers.push_back(&cells_solute.at(indexNeigh));
                }
                else {
                    neighbourPointers.push_back(0);
                }
            }
            cells_solute.at(index).setNeighPointers(neighbourPointers);
        }
    }
}

void Grid2d::init_water_cells(
    Settings& settings,
    std::vector < std::vector<std::string> >& materials,
    std::vector < std::vector<std::string> >& bound_cond,
    std::vector < std::vector<std::string> >& init_cond)
{
    for (size_t i = 0; i < cells_water.size(); i++)
    {
        // Set material.
        CellGeom2d* geomCell = cells_water.at(i).getGeom();
        int mat = geomCell->getMaterial();
        cells_water.at(i).setMannN(atof(materials.at(mat + 1).at(0 + 1).c_str()));
        cells_water.at(i).setDeprStor(atof(materials.at(mat + 1).at(1 + 1).c_str()));
        // Set initial condition.
        int initInd = geomCell->getInitCondInd();
        cells_water.at(i).setWaterDepth(atof(init_cond.at(initInd + 1).at(0 + 1).c_str()));
        // Set boundary condition.
        int boundInd = geomCell->getBoundCondInd();
        cells_water.at(i).setOutletIndex(atoi(bound_cond.at(boundInd + 1).at(0 + 1).c_str()));
        cells_water.at(i).setDistToOutlet(atof(bound_cond.at(boundInd + 1).at(1 + 1).c_str()));
        cells_water.at(i).setSinkIndex(atoi(bound_cond.at(boundInd + 1).at(2 + 1).c_str()));
    }
}

void Grid2d::init_heat_cells(
    Settings& settings,
    std::vector < std::vector<std::string> > &materials,
    std::vector < std::vector<std::string> > &bound_cond,
    std::vector < std::vector<std::string> > &init_cond)
{
    for (size_t i = 0; i < cells_heat.size(); i++)
    {

    }
}

void Grid2d::init_solute_cells(
    Settings& settings,
    std::vector < std::vector<std::string> > &materials,
    std::vector < std::vector<std::string> > &bound_cond,
    std::vector < std::vector<std::string> > &init_cond,
    std::vector < std::vector<std::string> > &solute_lib,
    std::vector <std::string>& species_names)
{
    // Initialize 2d solute cells.
    size_t cellsPerSolute2d = 0;
    int num_of_species = settings.get_int("num_of_species");
    if (num_of_species > 0)
    {
        cellsPerSolute2d = cells_solute.size() / num_of_species;
    }
    for (size_t i = 0; i < num_of_species; i++)
    {
        for (size_t j = 0; j < cellsPerSolute2d; j++)
        {
            int index = i * cellsPerSolute2d + j;
            cells_solute.at(index).setId(index); // ALREADY DONE?
            std::vector<std::string> soluteProp = find_solute_prop(species_names.at(i), solute_lib);
            cells_solute.at(index).setDeposDryRate(atof(soluteProp.at(0 + 1).c_str()));
            cells_solute.at(index).setDeposWetRate(atof(soluteProp.at(1 + 1).c_str()));
        }
    }
}

std::vector<std::string> Grid2d::find_solute_prop(std::string soluteName, std::vector<std::vector<std::string> > soluteLib)
{
    // Bypass the the header line.
    for (size_t i = 1; i < soluteLib.size(); i++)
    {
        if (soluteName == soluteLib.at(i).at(0))
        {
            return soluteLib.at(i);
        }
    }
    // When properties are not found in the library, return generic parameters.
    // Remember to update this when new properties are added to the library.
    std::vector<std::string> soluteProp{ "generic", "0.0", "0.0", "-1", "0", "1", "0", "-1", "0", "-1" };
    return soluteProp;
}

std::string Grid2d::parse_unstruct_vtk_mesh()
{
    std::stringstream ss;
    ss << std::setprecision(4) << std::fixed;
    ss << "# vtk DataFile Version 2.0" << std::endl;
    ss << "STYX test mesh." << std::endl;
    ss << "ASCII" << std::endl;
    ss << "DATASET" << " " << "UNSTRUCTURED_GRID" << std::endl;
    ss << "POINTS" << " " << vertices.size() << " " << "double" << std::endl;
    for (size_t i = 0; i < vertices.size(); i++) {
        ss << vertices.at(i).x << " " << vertices.at(i).y << " " << vertices.at(i).z << std::endl;
    }
    ss << "CELLS" << " " << cells_geom.size() << " " << 4 * cells_geom.size() << std::endl; // 5
    for (size_t i = 0; i < cells_geom.size(); i++) {
        std::vector<int> vertIndices = cells_geom.at(i).getVertIndices();
        ss << vertIndices.size();
        for (size_t j = 0; j < vertIndices.size(); j++) {
            ss << " " << vertIndices.at(j);
        }
        ss << std::endl;
    }
    ss << "CELL_TYPES" << " " << cells_geom.size() << std::endl;
    for (size_t i = 0; i < cells_geom.size(); i++) {
        ss << cells_geom.at(i).getGeomType() << std::endl;
    }
    ss << "CELL_DATA" << " " << cells_geom.size() << std::endl;
    ss << "FIELD" << " " << "FieldData" << " " << "4" << std::endl;
    ss << "id" << " " << "1" << " " << cells_geom.size() << " " << "int" << std::endl;
    for (size_t i = 0; i < cells_geom.size(); i++) {
        ss << cells_geom.at(i).getId() << std::endl;
    }
    ss << "material" << " " << "1" << " " << cells_geom.size() << " " << "int" << std::endl;
    for (size_t i = 0; i < cells_geom.size(); i++) {
        ss << cells_geom.at(i).getMaterial() << std::endl;
    }
    ss << "water_depth" << " " << "1" << " " << cells_water.size() << " " << "double" << std::endl;
    for (size_t i = 0; i < cells_water.size(); i++) {
        ss << cells_water.at(i).getWaterDepth() << std::endl;
    }
    ss << "discharge" << " " << "1" << " " << cells_water.size() << " " << "double" << std::endl;
    for (size_t i = 0; i < cells_water.size(); i++) {
        ss << cells_water.at(i).getAvgDisch() << std::endl;
    }

    return ss.str();
}
