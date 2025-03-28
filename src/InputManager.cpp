#include "InputManager.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

/**
 * @brief Loads a VTK file and parses the grid information including points, cells, cell types, and cell data.
 * 
 * @param filename Path to the VTK file.
 * @return Parsed VTKGrid structure containing points, cells, and associated metadata.
 * @throws std::runtime_error if file operations or parsing fails.
 */
VTKGrid InputManager::loadVTK(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open VTK file: " + filename);
    }

    VTKGrid grid;
    std::string line;

    std::cout << "Starting to parse VTK file: " << filename << std::endl;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "POINTS") {
            int numPoints;
            std::string dataType;
            iss >> numPoints >> dataType;

            std::cout << "Parsing POINTS, count: " << numPoints << std::endl;

            for (int i = 0; i < numPoints; ++i) {
                double x, y, z;
                if (!(ifs >> x >> y >> z)) {
                    throw std::runtime_error("Error reading points from: " + filename);
                }
                grid.points.push_back({x, y, z});
            }
            std::getline(ifs, line); // Clear newline
        }
        else if (token == "CELLS") {
            int numCells, totalInts;
            iss >> numCells >> totalInts;
            grid.cells.resize(numCells);

            std::cout << "Parsing CELLS, count: " << numCells << std::endl;

            for (int i = 0; i < numCells; ++i) {
                int numIndices;
                if (!(ifs >> numIndices)) {
                    throw std::runtime_error("Error reading cell vertices: " + filename);
                }
                std::vector<int> conn(numIndices);
                for (int j = 0; j < numIndices; ++j) {
                    if (!(ifs >> conn[j])) {
                        throw std::runtime_error("Error reading cell connectivity: " + filename);
                    }
                }
                grid.cells[i].connectivity = conn;
            }
            std::getline(ifs, line); // Clear newline
        }
        else if (token == "CELL_TYPES") {
            size_t numCellTypes;
            iss >> numCellTypes;

            if (numCellTypes != grid.cells.size()) {
                throw std::runtime_error("Mismatch CELL_TYPES and CELLS counts: " + filename);
            }

            std::cout << "Parsing CELL_TYPES, count: " << numCellTypes << std::endl;

            for (size_t i = 0; i < grid.cells.size(); ++i) {
                int ctype;
                if (!(ifs >> ctype)) {
                    throw std::runtime_error("Error reading cell type: " + filename);
                }
                grid.cells[i].cellType = ctype;
            }
            std::getline(ifs, line); // Clear newline
        }
        else if (token == "CELL_DATA") {
            size_t numCellData;
            iss >> numCellData;

            std::cout << "Parsing CELL_DATA, count: " << numCellData << std::endl;

            if (numCellData != grid.cells.size()) {
                throw std::runtime_error("Mismatch CELL_DATA and CELLS counts: " + filename);
            }

            // Read FIELD line
            std::string fieldLine;
            do {
                if (!std::getline(ifs, fieldLine)) {
                    throw std::runtime_error("EOF reached looking for FIELD line.");
                }
            } while (fieldLine.empty());

            std::istringstream fieldIss(fieldLine);
            std::string fieldToken, fieldLabel;
            int nFields;
            if (!(fieldIss >> fieldToken >> fieldLabel >> nFields) || fieldToken != "FIELD") {
                throw std::runtime_error("Invalid FIELD line format: " + fieldLine);
            }

            std::cout << "FIELD line found, fields count: " << nFields << std::endl;

            // Read each field
            for (int f = 0; f < nFields; ++f) {
                std::string fieldName, dataType;
                size_t numComponents, numTuples;

                // Read field header
                do {
                    if (!std::getline(ifs, line)) {
                        throw std::runtime_error("EOF while reading field header.");
                    }
                } while (line.empty());

                std::istringstream headerIss(line);
                if (!(headerIss >> fieldName >> numComponents >> numTuples >> dataType)) {
                    throw std::runtime_error("Invalid field header: " + line);
                }

                std::cout << "Reading field: " << fieldName 
                          << ", components: " << numComponents 
                          << ", tuples: " << numTuples << std::endl;

                std::vector<int> values(numTuples);
                for (size_t i = 0; i < numTuples; ++i) {
                    if (!(ifs >> values[i])) {
                        throw std::runtime_error("Error reading field values: " + fieldName);
                    }
                }

                // Assign field data to cells
                for (size_t i = 0; i < values.size(); ++i) {
                    if (fieldName == "id") grid.cells[i].id = values[i];
                    else if (fieldName == "grid_connection") grid.cells[i].gridConnection = values[i];
                    else if (fieldName == "material") grid.cells[i].material = values[i];
                    else if (fieldName == "init_cond") grid.cells[i].initCond = values[i];
                    else if (fieldName == "bound_cond") grid.cells[i].boundCond = values[i];
                }
            }
        }
    }

    // Print sample points and cells for quick inspection
    size_t pointsToPrint = std::min(grid.points.size(), size_t(3));
    for (size_t i = 0; i < pointsToPrint; ++i) {
        std::cout << "Point " << i << ": "
                  << grid.points[i].x << ", "
                  << grid.points[i].y << ", "
                  << grid.points[i].z << std::endl;
    }

    size_t cellsToPrint = std::min(grid.cells.size(), size_t(3));
    for (size_t i = 0; i < cellsToPrint; ++i) {
        std::cout << "Cell " << i << " - id: " << grid.cells[i].id
                  << ", type: " << grid.cells[i].cellType
                  << ", material: " << grid.cells[i].material
                  << ", init_cond: " << grid.cells[i].initCond
                  << ", bound_cond: " << grid.cells[i].boundCond
                  << ", connectivity:";
        size_t connToPrint = std::min(grid.cells[i].connectivity.size(), size_t(3));
        for (size_t j = 0; j < connToPrint; ++j)
            std::cout << " " << grid.cells[i].connectivity[j];
        std::cout << std::endl;
    }

    std::cout << "Finished parsing VTK file: " << filename << std::endl;
    return grid;
}

/**
 * @brief Loads a CSV file into a 2D vector of strings.
 * 
 * @param filename Path to the CSV file.
 * @return Data stored as a vector of vector of strings.
 * @throws std::runtime_error if the file cannot be opened.
 */
std::vector<std::vector<std::string>> InputManager::loadCSV(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open CSV file: " + filename);
    }

    std::vector<std::vector<std::string>> data;
    std::string line;
    while (std::getline(ifs, line)) {
        std::istringstream lineStream(line);
        std::string cell;
        std::vector<std::string> row;

        while (std::getline(lineStream, cell, ','))
            row.push_back(cell);
        
        data.push_back(row);
    }
    return data;
}

/**
 * @brief Loads an entire text file into a single string.
 * @param filename The path to the text file.
 * @return A string containing the file's contents.
 * @throws std::runtime_error if the file cannot be opened.
 */
std::string InputManager::loadTextFile(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open text file: " + filename);
    }
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    return buffer.str();
}