/*
for (int i = 0; i < cells2d.size(); i++) {
    std::vector<int> neighIndices = cells2d.at(i).getNeighIndices();
    std::cout << "Cell " << i << " neighbours:";
    for (int j = 0; j < neighIndices.size(); j++) {
        std::cout << " " << neighIndices.at(j);
    }
    std::cout << std::endl;
}
for (int i = 0; i < cells3d.size(); i++) {
    std::vector<int> neighIndices = cells3d.at(i).getNeighIndices();
    std::cout << "Cell " << i << " neighbours:";
    for (int j = 0; j < neighIndices.size(); j++) {
        std::cout << " " << neighIndices.at(j);
    }
    std::cout << std::endl;
}
*/
/*
// Set grid dimensions.
char nx = 2;
char ny = 2;
char nz = 2;
double lx = 10.0;
double ly = 10.0;
double lz = 10.0;
double dx = lx / nx;
double dy = ly / ny;
double dz = lz / nz;
// Create 2d vertices.
std::vector<Vertex> vertices2d;
int ind = 0;
for (int j = 0; j < ny + 1; j += 1) {
    for (int i = 0; i < nx + 1; i += 1) {
        double x = -0.5 * lx + i * dx;
        double y = -0.5 * ly + j * dy;
        double z = 0.0;
        vertices2d.push_back(Vertex(ind, x, y, z));
        ind += 1;
    }
}
// Create 2d cells and set pointers to vertices.
std::vector<CellGeom2d> cells2d;
ind = 0;
for (int j = 0; j < ny; j += 1) {
    for (int i = 0; i < nx; i += 1) {
        std::vector<int> vertIndices;
        vertIndices.push_back(i + j * (nx + 1) );
        vertIndices.push_back(i + 1 + j * (nx + 1));
        vertIndices.push_back(i + 1 + (j + 1) * (nx + 1) );
        vertIndices.push_back(i + (j + 1) * (nx + 1));
        std::vector<Vertex*> vertexPointers;
        for (int l = 0; l < vertIndices.size(); l += 1) {
            vertexPointers.push_back(&vertices2d.at(vertIndices.at(l)));
        }
        CellGeom2d cellGeom = CellGeom2d(ind, 9, vertIndices, vertexPointers);
        cells2d.push_back(cellGeom);
        ind += 1;
    }
}
// Set pointers to 2d cell neighbours.
ind = 0;
for (int j = 0; j < ny; j += 1) {
    for (int i = 0; i < nx; i += 1) {
        std::vector<int> neighIndices;
        if (j < ny - 1) {
            neighIndices.push_back(i + (j + 1) * nx);
        }
        else {
            neighIndices.push_back(-1);
        }
        if (i > 0) {
            neighIndices.push_back(i - 1 + j * nx);
        }
        else {
            neighIndices.push_back(-1);
        }
        if (j > 0) {
            neighIndices.push_back(i + (j - 1) * nx);
        }
        else {
            neighIndices.push_back(-1);
        }
        if (i < nx - 1) {
            neighIndices.push_back(i + 1 + j * nx);
        }
        else {
            neighIndices.push_back(-1);
        }
        std::vector<CellGeom2d*> neighbourPointers;
        for (int l = 0; l < neighIndices.size(); l += 1) {
            if (neighIndices.at(l) != -1) {
                neighbourPointers.push_back(&cells2d.at(neighIndices.at(l)));
            }
            else {
                neighbourPointers.push_back(0);
            }
        }
        cells2d.at(ind).setNeighIndices(neighIndices);
        cells2d.at(ind).setNeighPointers(neighbourPointers);
        ind += 1;
    }
}
// Create 3d vertices.
std::vector<Vertex> vertices3d;
ind = 0;
for (int k = 0; k < nz + 1; k += 1) {
    for (int j = 0; j < ny + 1; j += 1) {
        for (int i = 0; i < nx + 1; i += 1) {
            double x = -0.5 * lx + i * dx;
            double y = -0.5 * ly + j * dy;
            double z = k * -dz;
            vertices3d.push_back(Vertex(ind, x, y, z));
            ind += 1;
        }
    }
}
// Create 3d cells and set pointers to vertices.
std::vector<CellGeom3d> cells3d;
ind = 0;
for (int k = 0; k < nz; k += 1) {
    for (int j = 0; j < ny; j += 1) {
        for (int i = 0; i < nx; i += 1) {
            std::vector<int> vertIndices;
            vertIndices.push_back(i + j * (nx + 1) + (k + 1) * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + 1 + j * (nx + 1) + (k + 1) * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + 1 + (j + 1) * (nx + 1) + (k + 1) * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + (j + 1) * (nx + 1) + (k + 1) * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + j * (nx + 1) + k * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + 1 + j * (nx + 1) + k * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + 1 + (j + 1) * (nx + 1) + k * (nx + 1) * (ny + 1));
            vertIndices.push_back(i + (j + 1) * (nx + 1) + k * (nx + 1) * (ny + 1));
            std::vector<Vertex*> vertexPointers;
            for (int l = 0; l < vertIndices.size(); l += 1) {
                vertexPointers.push_back(&vertices3d.at(vertIndices.at(l)));
            }
            CellGeom3d cellGeom = CellGeom3d(ind, 12, vertIndices, vertexPointers);
            cells3d.push_back(cellGeom);
            ind += 1;
        }
    }
}
// Set pointers to 3d cell neighbours.
ind = 0;
for (int k = 0; k < nz; k += 1) {
    for (int j = 0; j < ny; j += 1) {
        for (int i = 0; i < nx; i += 1) {
            std::vector<int> neighIndices;
            if (k > 0) {
                neighIndices.push_back((k - 1) * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            if (j < ny - 1) {
                neighIndices.push_back(i + (j + 1) * nx + k * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            if (i > 0) {
                neighIndices.push_back(i - 1 + j * nx + k * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            if (j > 0) {
                neighIndices.push_back(i + (j - 1) * nx + k * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            if (i < nx - 1) {
                neighIndices.push_back(i + 1 + j * nx + k * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            if (k < nz - 1) {
                neighIndices.push_back((k + 1) * nx * ny);
            }
            else {
                neighIndices.push_back(-1);
            }
            std::vector<CellGeom3d*> neighbourPointers;
            for (int l = 0; l < neighIndices.size(); l += 1) {
                if (neighIndices.at(l) != -1) {
                    neighbourPointers.push_back(&cells3d.at(neighIndices.at(l)));
                }
                else {
                    neighbourPointers.push_back(0);
                }
            }
            cells3d.at(ind).setNeighIndices(neighIndices);
            cells3d.at(ind).setNeighPointers(neighbourPointers);
            ind += 1;
        }
    }
}
*/

/*
//std::vector<Vertex> vertices2d;
//std::vector<CellGeom2d> cells2d;
//bool isReadingCellData = false;
//bool isReadingFieldData = false;
//for (int i = 0; i < tokens2d.size(); i++) {
//    std::cout << "i: " << i;
//    for (int j = 0; j < tokens2d.at(i).size(); j++) {
//        std::cout << " " << tokens2d.at(i).at(j);
//    }
//    std::cout << std::endl;
//}
for (int i = 0; i < tokens2d.size(); i++) {
    if (tokens2d.at(i).at(0) == "POINTS") {
        int numOfPoints = atoi(tokens2d.at(i).at(1).c_str());
        int j;
        for (j = i + 1; j < i + 1 + numOfPoints; j++) {
            Vertex point;
            point.x = atof(tokens2d.at(j).at(0).c_str());
            point.y = atof(tokens2d.at(j).at(1).c_str());
            point.z = atof(tokens2d.at(j).at(2).c_str());
            vertices2d.push_back(point);
        }
        i = j;
    }
    if (tokens2d.at(i).at(0) == "CELLS") {
        int numOfCells = atoi(tokens2d.at(i).at(1).c_str());
        int j;
        for (j = i + 1; j < i + 1 + numOfCells; j++) {
            int pointsInCell = atoi(tokens2d.at(j).at(0).c_str());
            std::vector<int> vertIndices;
            for (int k = 1; k < pointsInCell + 1; k++) {
                vertIndices.push_back(atoi(tokens2d.at(j).at(k).c_str()));
            }
            std::vector<Vertex*> vertexPointers;
            for (int k = 0; k < vertIndices.size(); k += 1) {
                vertexPointers.push_back(&vertices2d.at(vertIndices.at(k)));
            }
            CellGeom2d cellGeom = CellGeom2d(-1, -1, vertIndices, vertexPointers);
            cells2d.push_back(cellGeom);
        }
        i = j;
    }
    if (tokens2d.at(i).at(0) == "CELL_TYPES") {
        int numOfCellTypes = atoi(tokens2d.at(i).at(1).c_str());
        int j;
        int ind = 0;
        for (j = i + 1; j < i + 1 + numOfCellTypes; j++) {
            char geomType = atoi(tokens2d.at(j).at(0).c_str());
            // Handle error.
            int result = cells2d.at(ind).setGeomType(geomType);
            ind++;
        }
        i = j;
    }
    if (tokens2d.at(i).at(0) == "CELL_DATA") {
        //std::cout << tokens2d.at(i).at(0) << std::endl;
        int numOfCellData = atoi(tokens2d.at(i).at(1).c_str()); // not used?
        isReadingCellData = true;
    }
    if (tokens2d.at(i).at(0) == "FIELD") {
        int numOfFields = atoi(tokens2d.at(i).at(2).c_str());
        std::cout << "Number of fields: " << numOfFields << std::endl;
        isReadingFieldData = true;
        int numOfCellData = 4;
        for (int j = 0; j < numOfFields; j++) {
            //std::cout << "index: " << i + 1 + j * numOfFields << std::endl;
            std::cout << "Number of cell data: " << numOfCellData << std::endl;
            int k;
            int ind = 0;
            for (k = i + 1 + 1 + j + j * numOfCellData; k < i + 1 + 1 + j + numOfCellData + j * numOfCellData; k++) {
                std::cout << "k: " << k << std::endl;
                ind++;
            }
        }
    }
    // Increment i here.
}
*/




/*
// Create a stream of tokens from loaded data.
itemCount = 0;
pointsFound = false;
numOfPoints = 0;
std::vector<Vertex> vertices3d;
coords.clear();
cellsFound = false;
numOfCells = 0;
vertIndices.clear();
std::vector<CellGeom3d> cells3d;
cellTypesFound = false;
numOfCellTypes = 0;
cellInd = 0;
numOfFieldData = 0;
fieldDataFound = false;
fieldCount = 0;
fieldItemCount = 0;
fieldName = "";
fieldDataCols = 0;
fieldDataType = "";
dataInt.clear();
dataReal.clear();
processingField = false;
for (int i = 0; i < tokens3d.size(); i++) {
    for (int j = 0; j < tokens3d.at(i).size(); j++) {
        if (tokens3d.at(i).at(j) == "POINTS" || pointsFound == true) {
            if (itemCount == 0) {
                std::cout << "POINTS found!" << std::endl;
                pointsFound = true;
            }
            else if (itemCount == 1) {
                numOfPoints = atoi(tokens3d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                coords.push_back(atof(tokens3d.at(i).at(j).c_str()));
                if (coords.size() == 3) {
                    Vertex point(-1, coords.at(0), coords.at(1), coords.at(2)); // index is missing
                    vertices3d.push_back(point);
                    coords.clear();
                    if (vertices3d.size() == numOfPoints) {
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
        else if (tokens3d.at(i).at(j) == "CELLS" || cellsFound == true) {
            if (itemCount == 0) {
                std::cout << "CELLS found" << std::endl;
                cellsFound = true;
            }
            else if (itemCount == 1) {
                numOfCells = atoi(tokens3d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                vertIndices.push_back(atoi(tokens3d.at(i).at(j).c_str()));
                if (vertIndices.size() == vertIndices.at(0) + 1) {
                    std::cout << "Saving vertex indices." << std::endl;
                    // Remove the number of points in cell item.
                    vertIndices.erase(vertIndices.begin());
                    std::vector<Vertex*> vertexPointers;
                    for (int l = 0; l < vertIndices.size(); l += 1) {
                        vertexPointers.push_back(&vertices3d.at(vertIndices.at(l)));
                    }
                    CellGeom3d cellGeom = CellGeom3d(-1, -1, vertIndices, vertexPointers); // index and cell are missing
                    cells3d.push_back(cellGeom);
                    vertIndices.clear();
                    if (cells3d.size() == numOfCells) {
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
        else if (tokens3d.at(i).at(j) == "CELL_TYPES" || cellTypesFound == true) {
            if (itemCount == 0) {
                std::cout << "CELL_TYPES found!" << std::endl;
                cellTypesFound = true;
            }
            else if (itemCount == 1) {
                numOfCellTypes = atoi(tokens3d.at(i).at(j).c_str());
            }
            else if (itemCount > 1) {
                int cellType = atoi(tokens3d.at(i).at(j).c_str());
                cells3d.at(cellInd).setGeomType(cellType);
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
        else if (tokens3d.at(i).at(j) == "CELL_DATA") {
            std::cout << "CELL_DATA found!" << std::endl;
        }
        else if (tokens3d.at(i).at(j) == "FIELD" || fieldDataFound == true) {
            if (itemCount == 0) {
                std::cout << "FIELD found!" << std::endl;
                fieldDataFound = true;
            }
            else if (itemCount == 2) {
                numOfFieldData = atoi(tokens3d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                if (fieldItemCount == 0) {
                    fieldName = tokens3d.at(i).at(j);
                    processingField = true;
                    std::cout << "Field name: " << fieldName << std::endl;
                }
                else if (fieldItemCount == 1) {
                    fieldDataCols = atoi(tokens3d.at(i).at(j).c_str());
                    std::cout << "Field data cols: " << fieldDataCols << std::endl;
                }
                else if (fieldItemCount == 3) {
                    fieldDataType = tokens3d.at(i).at(j);
                    std::cout << "Field data type: " << fieldDataType << std::endl;
                }
                else if (fieldItemCount > 3) {
                    std::cout << "Field item count: " << fieldItemCount << std::endl;
                    // field column count should be here.
                    if (fieldDataType == "int") {
                        dataInt.push_back(atoi(tokens3d.at(i).at(j).c_str()));
                    }
                    else if (fieldDataType == "double") {
                        dataReal.push_back(atof(tokens3d.at(i).at(j).c_str()));
                    }
                    // Save field data into cells.
                    if (fieldName == "id" && dataInt.size() == fieldDataCols) {
                        cells3d.at(cellInd).setId(dataInt.at(0));
                        dataInt.clear();
                        cellInd++;
                    }
                    else if (dataInt.size() == fieldDataCols || dataReal.size() == fieldDataCols) {
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
// Initialise cells.
for (int i = 0; i < cells3d.size(); i++) {
    std::cout << "Processing cell " << cells3d.at(i).getId() << std::endl;
    cells3d.at(i).computeCentrePoint();
    cells3d.at(i).assignFaceVertIndices();
    cells3d.at(i).computeFaceCentrePoints();
    cells3d.at(i).computeFaceAreas();
    cells3d.at(i).computeVolume();
}
// Create a map of sorted cell side vertex indicies.
std::map< int, std::vector<std::vector<int> > > mapCellFaceVertices;
for (int i = 0; i < cells3d.size(); i++) {
    std::vector<std::vector<int> > faceVertices = cells3d.at(i).getFaceVertIndices();
    for (int j = 0; j < faceVertices.size(); j++) {
        std::sort(faceVertices.at(j).begin(), faceVertices.at(j).end());
    }
    mapCellFaceVertices[i] = faceVertices;
}
// Find the cell neigbours and save the results to cells.
for (int i = 0; i < cells3d.size(); i++) {
    // Set current.
    std::vector<std::vector<int> > faceVertsSource = mapCellFaceVertices[i];
    std::vector<int> neighIndices;
    for (int j = 0; j < faceVertsSource.size(); j++) {
        bool neighFound = false;
        for (int k = 0; k < cells3d.size(); k++) {
            // Set target.
            std::vector<std::vector<int> > faceVertsTarget = mapCellFaceVertices[k];
            for (int l = 0; l < faceVertsTarget.size(); l++) {
                if (faceVertsSource.at(j).at(0) == faceVertsTarget.at(l).at(0) &&
                    faceVertsSource.at(j).at(1) == faceVertsTarget.at(l).at(1) &&
                    faceVertsSource.at(j).at(2) == faceVertsTarget.at(l).at(2)) {
                    neighIndices.push_back(k);
                    bool neighFound = true;
                    break;
                }
            }
            if (neighFound == true) {
                break;
            }
        }
        if (neighFound == false) {
            neighIndices.push_back(-1);
        }
    }
    // Save indices and pointers to cell.
    std::vector<CellGeom3d*> neighbourPointers;
    for (int j = 0; j < neighIndices.size(); j += 1) {
        if (neighIndices.at(j) != -1) {
            neighbourPointers.push_back(&cells3d.at(neighIndices.at(j)));
        }
        else {
            neighbourPointers.push_back(0);
        }
    }
    cells3d.at(i).setNeighIndices(neighIndices);
    cells3d.at(i).setNeighPointers(neighbourPointers);
}
// Compute intercell properties.
for (int i = 0; i < cells3d.size(); i++) {
    std::cout << "Processing cell " << cells3d.at(i).getId() << std::endl;
    cells3d.at(i).computeDistancesBetweenCells();
}
*/


/*
// Create a stream of tokens from loaded data.
int itemCount = 0;
bool pointsFound = false;
int numOfPoints = 0;
std::vector<Vertex> vertices2d;
std::vector<double> coords;
bool cellsFound = false;
int numOfCells = 0;
std::vector<int> vertIndices;
std::vector<CellGeom2d> cells2d;
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
for (int i = 0; i < tokens2d.size(); i++) {
    for (int j = 0; j < tokens2d.at(i).size(); j++) {
        if (tokens2d.at(i).at(j) == "POINTS" || pointsFound == true) {
            if (itemCount == 0) {
                std::cout << "POINTS found!" << std::endl;
                pointsFound = true;
            }
            else if (itemCount == 1) {
                numOfPoints = atoi(tokens2d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                coords.push_back( atof(tokens2d.at(i).at(j).c_str()) );
                if (coords.size() == 3) {
                    Vertex point(-1, coords.at(0), coords.at(1), coords.at(2)); // index is missing
                    vertices2d.push_back(point);
                    coords.clear();
                    if (vertices2d.size() == numOfPoints) {
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
        else if (tokens2d.at(i).at(j) == "CELLS" || cellsFound == true) {
            if (itemCount == 0) {
                std::cout << "CELLS found" << std::endl;
                cellsFound = true;
            }
            else if (itemCount == 1) {
                numOfCells = atoi(tokens2d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                vertIndices.push_back(atoi(tokens2d.at(i).at(j).c_str()));
                if (vertIndices.size() == vertIndices.at(0) + 1) {
                    std::cout << "Saving vertex indices." << std::endl;
                    // Remove the number of points in cell item.
                    vertIndices.erase(vertIndices.begin());
                    std::vector<Vertex*> vertexPointers;
                    for (int l = 0; l < vertIndices.size(); l += 1) {
                        vertexPointers.push_back(&vertices2d.at(vertIndices.at(l)));
                    }
                    CellGeom2d cellGeom = CellGeom2d(-1, -1, vertIndices, vertexPointers); // index and cell are missing
                    cells2d.push_back(cellGeom);
                    vertIndices.clear();
                    if (cells2d.size() == numOfCells) {
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
        else if (tokens2d.at(i).at(j) == "CELL_TYPES" || cellTypesFound == true) {
            if (itemCount == 0) {
                std::cout << "CELL_TYPES found!" << std::endl;
                cellTypesFound = true;
            }
            else if (itemCount == 1) {
                numOfCellTypes = atoi(tokens2d.at(i).at(j).c_str());
            }
            else if (itemCount > 1) {
                int cellType = atoi(tokens2d.at(i).at(j).c_str());
                cells2d.at(cellInd).setGeomType(cellType);
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
        else if (tokens2d.at(i).at(j) == "CELL_DATA") {
            std::cout << "CELL_DATA found!" << std::endl;
        }
        else if (tokens2d.at(i).at(j) == "FIELD" || fieldDataFound == true) {
            if (itemCount == 0) {
                std::cout << "FIELD found!" << std::endl;
                fieldDataFound = true;
            }
            else if (itemCount == 2) {
                numOfFieldData = atoi(tokens2d.at(i).at(j).c_str());
            }
            else if (itemCount > 2) {
                if (fieldItemCount == 0) {
                    fieldName = tokens2d.at(i).at(j);
                    processingField = true;
                    std::cout << "Field name: " << fieldName << std::endl;
                }
                else if (fieldItemCount == 1) {
                    fieldDataCols = atoi(tokens2d.at(i).at(j).c_str());
                    std::cout << "Field data cols: " << fieldDataCols << std::endl;
                }
                else if (fieldItemCount == 3) {
                    fieldDataType = tokens2d.at(i).at(j);
                    std::cout << "Field data type: " << fieldDataType << std::endl;
                }
                else if (fieldItemCount > 3) {
                    std::cout << "Field item count: " << fieldItemCount << std::endl;
                    // field column count should be here.
                    if (fieldDataType == "int") {
                        dataInt.push_back(atoi(tokens2d.at(i).at(j).c_str()));
                    }
                    else if (fieldDataType == "double") {
                        dataReal.push_back(atof(tokens2d.at(i).at(j).c_str()));
                    }
                    // Save field data into cells.
                    if (fieldName == "id" && dataInt.size() == fieldDataCols) {
                        cells2d.at(cellInd).setId(dataInt.at(0));
                        dataInt.clear();
                        cellInd++;
                    }
                    else if (dataInt.size() == fieldDataCols || dataReal.size() == fieldDataCols) {
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
// Initialise cells.
for (int i = 0; i < cells2d.size(); i++) {
    std::cout << "Processing cell " << cells2d.at(i).getId() << std::endl;
    cells2d.at(i).computeCentrePoint();
    cells2d.at(i).assignSideVertIndices();
    cells2d.at(i).computeSideCentrePoints();
    cells2d.at(i).computeSideLengths();
    cells2d.at(i).computeArea();
}
// Create a map of sorted cell side vertex indices.
std::map< int,std::vector<std::vector<int> > > mapCellSideVertices;
for (int i = 0; i < cells2d.size(); i++) {
    std::vector<std::vector<int> > sideVertices = cells2d.at(i).getSideVertIndices();
    for (int j = 0; j < sideVertices.size(); j++) {
        std::sort(sideVertices.at(j).begin(), sideVertices.at(j).end());
    }
    mapCellSideVertices[i] = sideVertices;
}
// Find the cell neigbours and save the results to cells.
for (int i = 0; i < cells2d.size(); i++) {
    // Set current.
    std::vector<std::vector<int> > sideVertsSource = mapCellSideVertices[i];
    std::vector<int> neighIndices;
    for (int j = 0; j < sideVertsSource.size(); j++) {
        bool neighFound = false;
        for (int k = 0; k < cells2d.size(); k++) {
            // Set target.
            std::vector<std::vector<int> > sideVertsTarget = mapCellSideVertices[k];
            for (int l = 0; l < sideVertsTarget.size(); l++) {
                if (sideVertsSource.at(j).at(0) == sideVertsTarget.at(l).at(0) &&
                    sideVertsSource.at(j).at(1) == sideVertsTarget.at(l).at(1)) {
                    neighIndices.push_back(k);
                    bool neighFound = true;
                    break;
                }
            }
            if (neighFound == true) {
                break;
            }
        }
        if (neighFound == false) {
            neighIndices.push_back(-1);
        }
    }
    // Save indices and pointers to cell.
    std::vector<CellGeom2d*> neighbourPointers;
    for (int j = 0; j < neighIndices.size(); j += 1) {
        if (neighIndices.at(j) != -1) {
            neighbourPointers.push_back(&cells2d.at(neighIndices.at(j)));
        }
        else {
            neighbourPointers.push_back(0);
        }
    }
    cells2d.at(i).setNeighIndices(neighIndices);
    cells2d.at(i).setNeighPointers(neighbourPointers);
}
// Compute intercell properties.
for (int i = 0; i < cells2d.size(); i++) {
    std::cout << "Processing cell " << cells2d.at(i).getId() << std::endl;
    cells2d.at(i).computeDistancesBetweenCells();
}
*/

/*
std::cout << "Vertex indices: " << std::endl;
for (int i = 0; i < cells2d.size(); i++) {
    std::cout << cells2d.at(i).getId() << ": ";
    std::vector<int> vertIndices = cells2d.at(i).getVertIndices();
    for (int j = 0; j < vertIndices.size(); j++) {
        std::cout << " " << vertIndices.at(j);
    }
    std::cout << std::endl;
}

std::cout << "Neighbour indices: " << std::endl;
for (int i = 0; i < cells2d.size(); i++) {
    std::cout << cells2d.at(i).getId() << ": ";
    std::vector<int> neighIndices = cells2d.at(i).getNeighIndices();
    for (int j = 0; j < neighIndices.size(); j++) {
        std::cout << " " << neighIndices.at(j);
    }
    std::cout << ", vector length: " << neighIndices.size();
    std::cout << std::endl;
}

std::cout << "Side vertex indices: " << std::endl;
for (int i = 0; i < cells2d.size(); i++) {
    std::vector<std::vector<int> > sideVertices = cells2d.at(i).getSideVertIndices();
    std::cout << "cell " << i;
    for (int j = 0; j < sideVertices.size(); j++) {
        //std::sort(sideVertices.at(j).begin(), sideVertices.at(j).end());
        std::cout << ", side " << j << ": ";
        for (int k = 0; k < sideVertices.at(j).size(); k++) {
            std::cout << " " << sideVertices.at(j).at(k);
        }
    }
    std::cout << std::endl;
}
*/

/*
std::cout << "Neighbour indices: " << std::endl;
for (int i = 0; i < cells3d.size(); i++) {
    std::cout << cells3d.at(i).getId() << ": ";
    std::vector<int> neighIndices = cells3d.at(i).getNeighIndices();
    for (int j = 0; j < neighIndices.size(); j++) {
        std::cout << " " << neighIndices.at(j);
    }
    std::cout << ", vector length: " << neighIndices.size();
    std::cout << std::endl;
}
*/
/*
// Create a map of sorted cell side vertex indicies.
std::map< int, std::vector<std::vector<int> > > mapCellFaceVertices;
for (int i = 0; i < cells.size(); i++) {
    std::vector<std::vector<int> > faceVertices = cells.at(i).getFaceVertIndices();
    for (int j = 0; j < faceVertices.size(); j++) {
        std::sort(faceVertices.at(j).begin(), faceVertices.at(j).end());
    }
    mapCellFaceVertices[i] = faceVertices;
}
*/

/*
// Create a map of sorted cell side vertex indices.
std::map< int, std::vector<std::vector<int> > > mapCellFaceVertices;
std::map< int, std::vector<int> > mapNeighIndices;
for (int i = 0; i < cells.size(); i++) {
    std::vector<std::vector<int> > locFaceVertIndices = cells.at(i).getFaceVertIndices();
    std::vector<int> vertIndices = cells.at(i).getVertIndices();
    std::vector<std::vector<int> > faceVertIndices;
    for (int j = 0; j < locFaceVertIndices.size(); j++) {
        std::vector<int> faceVerts;
        for (int k = 0; k < locFaceVertIndices.at(j).size(); k++) {
            faceVerts.push_back(vertIndices.at(locFaceVertIndices.at(j).at(k)));
        }
        std::sort(faceVerts.begin(), faceVerts.end());
        faceVertIndices.push_back(faceVerts);
    }
    mapCellFaceVertices[i] = faceVertIndices;
    std::vector<int> neighIndices(locFaceVertIndices.size(), -2);
    mapNeighIndices[i] = neighIndices;
}
// Find the cell neigbours and save the results to cells.
// This approach is extreamly slow. Optimization ideas below:
// - Initialise neighbour index to -2 for all sides/faces
// - Save the neighbour information to both tested cells.
// - Bypass fully processed cells (if no neighbour is found, the cell can still be assigned as processed).
// - Bypass processed sides/faces.
std::vector<bool> cellProcessed(cells.size(), false);
for (int i = 0; i < cells.size(); i++) {
    // Set current.
    std::vector<std::vector<int> > faceVertsCurrent = mapCellFaceVertices[i];
    for (int j = 0; j < faceVertsCurrent.size(); j++) {
        // Bypass a face that has already been processed.
        if (mapNeighIndices[i].at(j) != -2) {
            continue;
        }
        bool neighFound = false;
        for (int k = 0; k < cells.size(); k++) {
            // Bypass itself and fully processed cells.
            if (k == i || cellProcessed.at(k) == true) {
                continue;
            }
            // Set target.
            std::vector<std::vector<int> > faceVertsTarget = mapCellFaceVertices[k];
            for (int l = 0; l < faceVertsTarget.size(); l++) {
                // Face has already been connected.
                if (mapNeighIndices[k].at(l) != -2) {
                    continue;
                }
                if (faceVertsCurrent.at(j).at(0) == faceVertsTarget.at(l).at(0) &&
                    faceVertsCurrent.at(j).at(1) == faceVertsTarget.at(l).at(1) &&
                    faceVertsCurrent.at(j).at(2) == faceVertsTarget.at(l).at(2)) {
                    // Save the neigbour cell index to both tested cells.
                    mapNeighIndices[i].at(j) = k;
                    mapNeighIndices[k].at(l) = i;
                    neighFound = true;
                    break;
                }
            }
            if (neighFound == true) {
                break;
            }
        }
        if (neighFound == false) {
            mapNeighIndices[i].at(j) = -1;
        }
    }
    cellProcessed.at(i) = true;
    // Save indices and pointers to cell.
    std::vector<CellGeom3d*> neighbourPointers;
    for (int j = 0; j < mapNeighIndices[i].size(); j += 1) {
        if (mapNeighIndices[i].at(j) != -1) {
            neighbourPointers.push_back(&cells.at(mapNeighIndices[i].at(j)));
        }
        else {
            neighbourPointers.push_back(0);
        }
    }
    cells.at(i).setNeighIndices(mapNeighIndices[i]);
    cells.at(i).setNeighPointers(neighbourPointers);
    cells.at(i).computeDistancesBetweenCells();
}
*/
////////////////////////////////////////////////////////////////////////////////////

    /*
    // Create a map of sorted cell side vertex indices.
    std::map< int, std::vector<std::vector<int> > > mapCellSideVertices;
    std::map< int, std::vector<int> > mapNeighIndices;
    for (int i = 0; i < cells.size(); i++) {
        std::vector<std::vector<int> > locSideVertIndices = cells.at(i).getSideVertIndices();
        std::vector<int> vertIndices = cells.at(i).getVertIndices();
        std::vector<std::vector<int> > sideVertIndices;
        for (int j = 0; j < locSideVertIndices.size(); j++) {
            std::vector<int> sideVerts;
            for (int k = 0; k < locSideVertIndices.at(j).size(); k++) {
                sideVerts.push_back( vertIndices.at( locSideVertIndices.at(j).at(k) ) );
            }
            std::sort(sideVerts.begin(), sideVerts.end());
            sideVertIndices.push_back(sideVerts);
        }
        mapCellSideVertices[i] = sideVertIndices;
        std::vector<int> neighIndices(locSideVertIndices.size(), -2);
        mapNeighIndices[i] = neighIndices;
    }
    // Find the cell neigbours and save the results to cells.
    // This approach is extreamly slow. Optimization ideas below:
    // - Initialise neighbour index to -2 for all sides/faces
    // - Save the neighbour information to both tested cells.
    // - Bypass fully processed cells (if no neighbour is found, the cell can still be assigned as processed).
    // - Bypass processed sides/faces.
    std::vector<bool> cellProcessed(cells.size(), false);
    for (int i = 0; i < cells.size(); i++) {
        // Set current.
        std::vector<std::vector<int> > sideVertsCurrent = mapCellSideVertices[i];
        for (int j = 0; j < sideVertsCurrent.size(); j++) {
            // Bypass a side that has already been processed.
            if (mapNeighIndices[i].at(j) != -2) {
                continue;
            }
            bool neighFound = false;
            for (int k = 0; k < cells.size(); k++) {
                // Bypass itself and fully processed cells.
                if (k == i || cellProcessed.at(k) == true) {
                    continue;
                }
                // Set target.
                std::vector<std::vector<int> > sideVertsTarget = mapCellSideVertices[k];
                for (int l = 0; l < sideVertsTarget.size(); l++) {
                    // Side has already been connected.
                    if (mapNeighIndices[k].at(l) != -2) {
                        continue;
                    }
                    if (sideVertsCurrent.at(j).at(0) == sideVertsTarget.at(l).at(0) &&
                        sideVertsCurrent.at(j).at(1) == sideVertsTarget.at(l).at(1)) {
                        // Save the neigbour cell index to both tested cells.
                        mapNeighIndices[i].at(j) = k;
                        mapNeighIndices[k].at(l) = i;
                        neighFound = true;
                        break;
                    }
                }
                if (neighFound == true) {
                    break;
                }
            }
            if (neighFound == false) {
                mapNeighIndices[i].at(j) = -1;
            }
        }
        cellProcessed.at(i) = true;
        // Save indices and pointers to cell.
        std::vector<CellGeom2d*> neighbourPointers;
        for (int j = 0; j < mapNeighIndices[i].size(); j += 1) {
            if (mapNeighIndices[i].at(j) != -1) {
                neighbourPointers.push_back(&cells.at(mapNeighIndices[i].at(j)));
            }
            else {
                neighbourPointers.push_back(0);
            }
        }
        cells.at(i).setNeighIndices(mapNeighIndices[i]);
        cells.at(i).setNeighPointers(neighbourPointers);
        cells.at(i).computeDistancesBetweenCells();
    }
    */

    //////////////////////////////////////////////////////////////////////////////