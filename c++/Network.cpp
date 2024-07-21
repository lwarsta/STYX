#include "Network.h"

void Network::build_network(std::vector<std::vector<std::string>>& tokens_jnc,
                            std::vector<std::vector<std::string>>& tokens_lnk) {    
    // Load junction data.
    std::vector<std::vector<double>> points_jnc;
    std::vector<std::vector<int>> cells_jnc;
    std::vector<int> cell_types_jnc;
    std::vector<std::string> field_data_int_names_jnc;
    std::vector<std::vector<std::vector<int>>> field_data_int_jnc;
    std::vector<std::string> field_data_double_names_jnc;
    std::vector<std::vector<std::vector<double>>> field_data_double_jnc;
    
    parseVTKData(tokens_jnc, points_jnc, cells_jnc, cell_types_jnc, field_data_int_names_jnc,
                 field_data_int_jnc, field_data_double_names_jnc, field_data_double_jnc);
    
    //printVTKData(points_jnc, cells_jnc, cell_types_jnc, field_data_int_names_jnc,
    //             field_data_int_jnc, field_data_double_names_jnc, 
    //             field_data_double_jnc);
    // Save points.
    for (const auto& point : points_jnc) {
        Vertex vert((int)vertices.size(), point.at(0), point.at(1), point.at(2));
        vertices.push_back(vert);
    }
    
    // Load link data.
    std::vector<std::vector<double>> points_lnk;
    std::vector<std::vector<int>> cells_lnk;
    std::vector<int> cell_types_lnk;
    std::vector<std::string> field_data_int_names_lnk;
    std::vector<std::vector<std::vector<int>>> field_data_int_lnk;
    std::vector<std::string> field_data_double_names_lnk;
    std::vector<std::vector<std::vector<double>>> field_data_double_lnk;
    parseVTKData(tokens_lnk, points_lnk, cells_lnk, cell_types_lnk, field_data_int_names_lnk,
        field_data_int_lnk, field_data_double_names_lnk, field_data_double_lnk);
    //printVTKData(points_lnk, cells_lnk, cell_types_lnk, field_data_int_names_lnk,
    //             field_data_int_lnk, field_data_double_names_lnk, 
    //             field_data_double_lnk);
    
    // Save points. All the points are saved into the same vector.
    // push_back seems to redo memory allocation?
    int cum_verts = (int)vertices.size();
    std::cout << "Number of vertices: " << vertices.size() << std::endl;

    for (const auto& point : points_lnk) {
        Vertex vert((int)vertices.size(), point.at(0), point.at(1), point.at(2));
        vertices.push_back(vert);
    }
    
    // Create geometry junction cells.
    for (size_t cell_ind = 0; cell_ind < cells_jnc.size(); cell_ind++) {
        // Create a vector of pointers to vertices.
        std::vector<Vertex*> vtx_ptrs;
        for (size_t coord_ind = 0; coord_ind < cells_jnc.at(cell_ind).size(); coord_ind += 1) {
            vtx_ptrs.push_back(&vertices.at(cells_jnc.at(cell_ind).at(coord_ind)));
        }
        // Create a geometry junction cell.
        JuncGeom cell_geom = JuncGeom(cell_ind, cell_types_jnc.at(cell_ind),
                                      cells_jnc.at(cell_ind), vtx_ptrs);
        juncs_geom.push_back(cell_geom);
    }
    
    // Save data to junction cells.
    for (size_t ser_ind = 0; ser_ind < field_data_int_names_jnc.size(); ser_ind++) {
        // Save grid connection index data.
        if (field_data_int_names_jnc.at(ser_ind) == "grid_connection") {
            for (size_t ind = 0; ind < juncs_geom.size(); ind++) {
                int data_item = field_data_int_jnc.at(ser_ind).at(ind).at(0);
                juncs_geom.at(ind).setGridConnection(data_item);
            }
        }
        // Save material index data.
        else if (field_data_int_names_jnc.at(ser_ind) == "material") {
            for (size_t ind = 0; ind < juncs_geom.size(); ind++) {
                int data_item = field_data_int_jnc.at(ser_ind).at(ind).at(0);
                juncs_geom.at(ind).setMaterial(data_item);
            }
        }
        // Save initial condition index data.
        else if (field_data_int_names_jnc.at(ser_ind) == "init_cond") {
            for (size_t ind = 0; ind < juncs_geom.size(); ind++) {
                int data_item = field_data_int_jnc.at(ser_ind).at(ind).at(0);
                juncs_geom.at(ind).setInitCondInd(data_item);
            }
        }
        // Save boundary condition index data.
        else if (field_data_int_names_jnc.at(ser_ind) == "bound_cond") {
            for (size_t ind = 0; ind < juncs_geom.size(); ind++) {
                int data_item = field_data_int_jnc.at(ser_ind).at(ind).at(0);
                juncs_geom.at(ind).setBoundCondInd(data_item);
            }
        }
    }
    
    // Create geometry link cells.
    for (size_t cell_ind = 0; cell_ind < cells_lnk.size(); cell_ind++) {
        // Create a vector of pointers to vertices.
        std::vector<Vertex*> vtx_ptrs;
        for (size_t coord_ind = 0; coord_ind < cells_lnk.at(cell_ind).size(); coord_ind += 1) {
            cells_lnk.at(cell_ind).at(coord_ind) += cum_verts;
            vtx_ptrs.push_back(&vertices.at(cells_lnk.at(cell_ind).at(coord_ind)));
        }
        // Create a geometry link cell.
        LinkGeom cell_geom = LinkGeom(cell_ind, cell_types_lnk.at(cell_ind),
                                      cells_lnk.at(cell_ind), vtx_ptrs);
        links_geom.push_back(cell_geom);
    }
    
    // Save data to link cells.
    for (size_t ser_ind = 0; ser_ind < field_data_int_names_lnk.size(); ser_ind++) {
        // Save grid connection index data.
        if (field_data_int_names_lnk.at(ser_ind) == "grid_connection") {
            for (size_t ind = 0; ind < links_geom.size(); ind++) {
                int data_item = field_data_int_jnc.at(ser_ind).at(ind).at(0);
                links_geom.at(ind).setGridConnection(data_item);
            }
        }
        // Save material index data.
        if (field_data_int_names_lnk.at(ser_ind) == "material") {
            for (size_t ind = 0; ind < links_geom.size(); ind++) {
                int data_item = field_data_int_lnk.at(ser_ind).at(ind).at(0);
                links_geom.at(ind).setMaterial(data_item);
            }
        }
        // Save initial condition index data.
        else if (field_data_int_names_lnk.at(ser_ind) == "init_cond") {
            for (size_t ind = 0; ind < links_geom.size(); ind++) {
                int data_item = field_data_int_lnk.at(ser_ind).at(ind).at(0);
                links_geom.at(ind).setInitCondInd(data_item);
            }
        }
        // Save boundary condition index data.
        else if (field_data_int_names_lnk.at(ser_ind) == "bound_cond") {
            for (size_t ind = 0; ind < links_geom.size(); ind++) {
                int data_item = field_data_int_lnk.at(ser_ind).at(ind).at(0);
                links_geom.at(ind).setBoundCondInd(data_item);
            }
        }
    }
    
    // Test print junctions.
    //for (size_t ind = 0; ind < juncs_geom.size(); ind++) {
    //    std::vector<Vertex*> vrts = juncs_geom.at(ind).getVertPointers();
    //    Vertex vtx = *vrts.at(0);
    //    std::cout << "Junction " << ind << " (" << vtx.x << ", " << vtx.y << ", " << vtx.z << ")" << std::endl;
    //}

    // Test print links.
    //for (size_t ind = 0; ind < links_geom.size(); ind++) {
    //    std::vector<Vertex*> vrts = links_geom.at(ind).getVertPointers();
    //    Vertex vtx = *vrts.at(0);
    //    std::cout << "Link " << ind << " (" << vtx.x << ", " << vtx.y << ", " << vtx.z << ")" << std::endl;
    //}
    
    // Connect junctions via links in the network.
    // Should some sort of threshold distance be used here?
    Algorithms algorithms;
    for (size_t ind_lnk = 0; ind_lnk < links_geom.size(); ind_lnk++) {
        std::vector<Vertex*> lnk_vrts = links_geom.at(ind_lnk).getVertPointers();
        // Compute closest junction to link end 0.
        int ind_jnc_0 = -1;
        double dist_0 = std::numeric_limits<double>::max();
        for (size_t ind_jnc = 0; ind_jnc < juncs_geom.size(); ind_jnc++) {
            std::vector<Vertex*> jnc_vrts = juncs_geom.at(ind_jnc).getVertPointers();
            Vertex vec = algorithms.create_vector(*lnk_vrts.at(0), *jnc_vrts.at(1));
            double vec_len = algorithms.compute_vector_length(vec);
            if (dist_0 > vec_len) {
                dist_0 = vec_len;
                ind_jnc_0 = ind_jnc;
            }
        }
        // Compute closest junction to link end 1.
        int ind_jnc_1 = -1;
        double dist_1 = std::numeric_limits<double>::max();
        for (size_t ind_jnc = 0; ind_jnc < juncs_geom.size(); ind_jnc++) {
            std::vector<Vertex*> jnc_vrts = juncs_geom.at(ind_jnc).getVertPointers();
            Vertex vec = algorithms.create_vector(*lnk_vrts.at(1), *jnc_vrts.at(1));
            double vec_len = algorithms.compute_vector_length(vec);
            if (dist_1 > vec_len) {
                dist_1 = vec_len;
                ind_jnc_1 = ind_jnc;
            }
        }
        // Save the result to the junctions.
        juncs_geom.at(ind_jnc_0).save_link(0, ind_lnk, &links_geom.at(ind_lnk), 
                      ind_jnc_1, &juncs_geom.at(ind_jnc_1));
        juncs_geom.at(ind_jnc_1).save_link(1, ind_lnk, &links_geom.at(ind_lnk),
                      ind_jnc_0, &juncs_geom.at(ind_jnc_0));
    }
    
}

void Network::parseVTKData(const std::vector<std::vector<std::string>>& tokens,
    std::vector<std::vector<double>>& points,
    std::vector<std::vector<int>>& cells,
    std::vector<int>& cell_types,
    std::vector<std::string>& field_data_int_names,
    std::vector<std::vector<std::vector<int>>>& field_data_int,
    std::vector<std::string>& field_data_double_names,
    std::vector<std::vector<std::vector<double>>>& field_data_double) {

    int itemCount = 0;
    int cellInd = 0;
    int numOfPoints = 0;
    int numOfCells = 0;
    int numOfVerts = 0;
    int numOfCellTypes = 0;
    int numOfFieldData = 0;
    bool pointsFound = false;
    bool cellsFound = false;
    bool cellTypesFound = false;
    bool fieldDataFound = false;
    bool processingField = false;
    int fieldItemCount = 0;
    std::string fieldName = "";
    int fieldDataCols = 0;
    std::string fieldDataType = "";
    int fieldCount = 0;
    std::vector<int> dataCellInt;
    std::vector<std::vector<int>> dataCellsInt;
    std::vector<double> dataCellReal;
    std::vector<std::vector<double>> dataCellsReal;
    
    for (size_t i = 0; i < tokens.size(); i++) {
        for (size_t j = 0; j < tokens.at(i).size(); j++) {
            // Process point data.
            if (tokens.at(i).at(j) == "POINTS" || pointsFound == true) {
                if (itemCount == 0) {
                    pointsFound = true;
                    points.clear();
                }
                else if (itemCount == 1) {
                    numOfPoints = atoi(tokens.at(i).at(j).c_str());
                    points.reserve(numOfPoints);

                    if (numOfPoints == 0) {
                        pointsFound = false;
                    }
                }
                else if (itemCount > 2) {
                    std::vector<double> point;
                    point.push_back(atof(tokens.at(i).at(j).c_str()));
                    j++;
                    point.push_back(atof(tokens.at(i).at(j).c_str()));
                    j++;
                    point.push_back(atof(tokens.at(i).at(j).c_str()));
                    points.push_back(point);

                    if (points.size() == static_cast<size_t>(numOfPoints)) {
                        pointsFound = false;
                    }
                }
                if (pointsFound == true) {
                    itemCount++;
                }
                else {
                    itemCount = 0;
                }
            }
            // Process cell vertex index data.
            else if (tokens.at(i).at(j) == "CELLS" || cellsFound == true) {
                if (itemCount == 0) {
                    cellsFound = true;
                    cells.clear();
                }
                else if (itemCount == 1) {
                    numOfCells = atoi(tokens.at(i).at(j).c_str());
                    cells.reserve(numOfCells);

                    if (numOfCells == 0) {
                        cellsFound = false;
                    }
                }
                else if (itemCount > 2) {
                    std::vector<int> vertIndices;
                    numOfVerts = atoi(tokens.at(i).at(j).c_str());
                    vertIndices.push_back(numOfVerts);
                    for (int k = 0; k < numOfVerts; k++) {
                        j++;
                        vertIndices.push_back(atoi(tokens.at(i).at(j).c_str()));
                    }
                    // Remove the number of points item in vertex index data.
                    vertIndices.erase(vertIndices.begin());
                    cells.push_back(vertIndices);

                    if (cells.size() == static_cast<size_t>(numOfCells)) {
                        cellsFound = false;
                    }
                }
                if (cellsFound == true) {
                    itemCount++;
                }
                else {
                    itemCount = 0;
                }
            }
            // Process cell type data.
            else if (tokens.at(i).at(j) == "CELL_TYPES" || cellTypesFound == true) {
                if (itemCount == 0) {
                    cellTypesFound = true;
                    cell_types.clear();
                }
                else if (itemCount == 1) {
                    numOfCellTypes = atoi(tokens.at(i).at(j).c_str());
                    cell_types.reserve(numOfCellTypes);

                    if (numOfCellTypes == 0) {
                        cellTypesFound = false;
                    }
                }
                else if (itemCount > 1) {
                    int cellType = atoi(tokens.at(i).at(j).c_str());
                    cell_types.push_back(cellType);

                    if (cell_types.size() == static_cast<size_t>(numOfCellTypes)) {
                        cellTypesFound = false;
                    }
                }
                if (cellTypesFound == true) {
                    itemCount++;
                }
                else {
                    itemCount = 0;
                }
            }
            // Process cell field data.
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
                        if (fieldDataType == "int") {
                            field_data_int_names.push_back(fieldName);
                        }
                        else if (fieldDataType == "double") {
                            field_data_double_names.push_back(fieldName);
                        }
                    }
                    else if (fieldItemCount > 3) {
                        // field column count should be here.
                        if (fieldDataType == "int") {
                            dataCellInt.push_back(atoi(tokens.at(i).at(j).c_str()));
                        }
                        else if (fieldDataType == "double") {
                            dataCellReal.push_back(atof(tokens.at(i).at(j).c_str()));
                        }

                        if (fieldDataType == "int" && (int)dataCellInt.size() == fieldDataCols)
                        {
                            dataCellsInt.push_back(dataCellInt);
                            dataCellInt.clear();
                            cellInd++;
                        }
                        else if (fieldDataType == "double" && (int)dataCellReal.size() == fieldDataCols)
                        {
                            dataCellsReal.push_back(dataCellReal);
                            dataCellReal.clear();
                            cellInd++;
                        }
                        else if ((int)dataCellInt.size() == fieldDataCols || (int)dataCellReal.size() == fieldDataCols) {
                            dataCellInt.clear();
                            dataCellReal.clear();
                            cellInd++;
                        }
                        // Move to the next field.
                        if (cellInd == numOfCells) { // should numOfCells be locally retrieved too?
                            if (fieldDataType == "int") {
                                field_data_int.push_back(dataCellsInt);
                            }
                            else if (fieldDataType == "double") {
                                field_data_double.push_back(dataCellsReal);
                            }
                            dataCellsInt.clear();
                            dataCellsReal.clear();
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
}

void Network::printVTKData(std::vector<std::vector<double>>& points,
    std::vector<std::vector<int>>& cells,
    std::vector<int>& cell_types,
    std::vector<std::string>& field_data_int_names,
    std::vector<std::vector<std::vector<int>>>& field_data_int,
    std::vector<std::string>& field_data_double_names,
    std::vector<std::vector<std::vector<double>>>& field_data_double) {

    std::cout << "Size of points vector: " << points.size() << "\n";
    std::cout << "Size of cells vector: " << cells.size() << "\n";
    std::cout << "Size of cell type vector: " << cell_types.size() << "\n";
    std::cout << "Size of int data name vector: " << field_data_int_names.size() << "\n";
    std::cout << "Size of int data vector: " << field_data_int.size() << "\n";
    std::cout << "Size of double data name vector: " << field_data_double_names.size() << "\n";
    std::cout << "Size of double data vector: " << field_data_double.size() << "\n";

    // Print points.
    std::cout << "Point coordinates" << std::endl;
    for (const auto& point : points) {
        for (double coordinate : point) {
            std::cout << coordinate << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print cells.
    std::cout << "Cell vertex ids" << std::endl;
    for (const auto& cell : cells) {
        for (double index : cell) {
            std::cout << index << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print cell types.
    std::cout << "Cell types" << std::endl;
    for (int cell_type : cell_types) {
            std::cout << cell_type << " ";
    }
    std::cout << std::endl;

    // Print int field data.
    size_t field_cnt = 0;
    for (const auto& cell_data_series : field_data_int) {
        std::cout << field_data_int_names.at(field_cnt).c_str() << std::endl;
        for (const auto& cell_data : cell_data_series) {
            for (int cell_data_item : cell_data) {
                std::cout << cell_data_item << " ";
            }
            std::cout << std::endl;
        }
        field_cnt++;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print double field data.
    field_cnt = 0;
    for (const auto& cell_data_series : field_data_double) {
        std::cout << field_data_double_names.at(field_cnt).c_str() << std::endl;
        for (const auto& cell_data : cell_data_series) {
            for (int cell_data_item : cell_data) {
                std::cout << cell_data_item << " ";
            }
            std::cout << std::endl;
        }
        field_cnt++;
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Network::create_water_network_items()
{
    // Create water links.
    links_water.resize(links_geom.size());

    for (size_t i = 0; i < links_water.size(); i++) {
        // Parametrise the water junction.
        links_water.at(i).setId(i);
        // Connect water junction to a geometry junction.
        links_water.at(i).assign_geom(&links_geom.at(i));
    }
    
    // Create water junctions.
    juncs_water.resize(juncs_geom.size());

    for (size_t i = 0; i < juncs_water.size(); i++) {
        // Parametrise the water junction.
        juncs_water.at(i).setId(i);
        // Connect water junction to a geometry junction.
        juncs_water.at(i).assign_geom(&juncs_geom.at(i));
        // Connect water junctions to water links and junctions.
        std::vector<int> ids_geom_link = juncs_geom.at(i).get_ids_lnk();
        std::vector<int> end_ids_geom_link = juncs_geom.at(i).get_ids_lnk_end();
        std::vector<int> ids_geom_junc = juncs_geom.at(i).get_ids_junc();
        for (size_t j = 0; j < ids_geom_link.size(); j++) {
            juncs_water.at(i).save_link(end_ids_geom_link.at(j), 
                                        ids_geom_link.at(j), 
                                        &links_water.at(ids_geom_link.at(j)),
                                        ids_geom_junc.at(j), 
                                        &juncs_water.at(ids_geom_junc.at(j))
                                        );
        }
    }
}

void Network::init_water_network(
    Settings& settings,
    std::vector<std::vector<std::string>>& materials_net_junc,
    std::vector<std::vector<std::string>>& bound_cond_net_junc,
    std::vector<std::vector<std::string>>& init_cond_net_junc,
    std::vector<std::vector<std::string>>& materials_new_link,
    std::vector<std::vector<std::string>>& bound_cond_net_link,
    std::vector<std::vector<std::string>>& init_cond_new_link)
{
    for (size_t i = 0; i < juncs_water.size(); i++)
    {
        // Set material.
        JuncGeom* geom_junc = juncs_water.at(i).get_geom();
        int mat = geom_junc->getMaterial();
        geom_junc->set_diameter(atof(materials_net_junc.at(mat + 1).at(0 + 1).c_str()));
        // Set initial condition.
        int initInd = geom_junc->getInitCondInd();
        // Set boundary condition.
        int bound_ind = geom_junc->getBoundCondInd();
        juncs_water.at(i).set_type(atoi(bound_cond_net_junc.at(bound_ind + 1).at(0 + 1).c_str()));
        // Compute geometric properties. IS THIS THE RIGHT PLACE TO DO THIS?
        geom_junc->comp_geom_properties();
        // TEMPORARILY test fill certain junctions with water to test the system.
        //if (i >= 0 && i <= 41) { // 34, 35, 36, 37, 38, 39
        //    juncs_water.at(i).set_water_depth(3.0);
        //    juncs_water.at(i).swap();
        //}
    }

    for (size_t i = 0; i < links_water.size(); i++)
    {
        // Set material.
        LinkGeom* geom_link = links_water.at(i).get_geom();
        int mat = geom_link->getMaterial();
        geom_link->set_diameter(atof(materials_new_link.at(mat + 1).at(0 + 1).c_str()));
        links_water.at(i).set_mann_n(atof(materials_new_link.at(mat + 1).at(1 + 1).c_str()));
        // Set initial condition.
        int initInd = geom_link->getInitCondInd();
        // Set boundary condition.
        int bound_ind = geom_link->getBoundCondInd();
        // Compute geometric properties. IS THIS THE RIGHT PLACE TO DO THESE?
        geom_link->comp_geom_properties();
        geom_link->compCenterPoint();
        // TEMPORARILY test fill certain link with water to test the system.
        //links_water.at(i).set_water_depth(geom_link->get_diameter());
        //links_water.at(i).swap();
    }
}