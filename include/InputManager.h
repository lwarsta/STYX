/**
 * @file InputManager.h
 * @brief Provides functionality to load VTK and CSV files.
 *
 * This file defines the InputManager class, which offers static methods for
 * loading VTK grid data and CSV files into appropriate data structures.
 */

#ifndef INPUT_MANAGER_H
#define INPUT_MANAGER_H

#include <string>
#include <vector>
#include "Common.h"

/**
 * @struct VTKCell
 * @brief Stores cell data parsed from a VTK file.
 *
 * The VTKCell structure holds the connectivity, type, id, grid connection,
 * material, initial condition, and boundary condition for a cell.
 */
struct VTKCell {
    std::vector<int> connectivity; ///< The indices of the cell’s vertices.
    int cellType;                  ///< VTK cell type (e.g., 3, 9, 12, etc.).
    int id;                        ///< Cell identifier (from CELL_DATA).
    int gridConnection;            ///< Connection data (to be updated).
    int material;                  ///< Material identifier.
    int initCond;                  ///< Initial condition.
    int boundCond;                 ///< Boundary condition.
};

/**
 * @struct VTKGrid
 * @brief Stores grid data parsed from a VTK file.
 *
 * The VTKGrid structure holds a vector of 3D points and a vector of VTK cells.
 */
struct VTKGrid {
    std::vector<Point3D> points; ///< List of points (vertex positions).
    std::vector<VTKCell> cells;  ///< List of cells.
};

/**
 * @class InputManager
 * @brief Loads VTK and CSV files.
 *
 * The InputManager class provides static functions to load VTK ASCII files
 * into a VTKGrid structure and CSV files into a vector of string rows.
 */
class InputManager {
public:
    /**
     * @brief Loads a VTK ASCII file.
     *
     * This function reads a VTK file and returns a VTKGrid structure containing
     * the points and cells. It throws a std::runtime_error if an error occurs.
     *
     * @param filename The name of the VTK file to load.
     * @return VTKGrid containing the parsed grid data.
     * @throws std::runtime_error if the file cannot be opened or parsed.
     */
    static VTKGrid loadVTK(const std::string &filename);

    /**
     * @brief Loads a CSV file.
     *
     * This function reads a CSV file and returns its contents as a vector of rows,
     * where each row is represented as a vector of strings.
     *
     * @param filename The name of the CSV file to load.
     * @return A vector of rows, each row being a vector of strings.
     * @throws std::runtime_error if the file cannot be opened.
     */
    static std::vector<std::vector<std::string>> loadCSV(const std::string &filename);
	
    /**
     * @brief Loads an entire text file into a single string.
     * @param filename The path to the text file.
     * @return A string containing the file's contents.
     * @throws std::runtime_error if the file cannot be opened.
     */
    static std::string loadTextFile(const std::string &filename);
};

#endif // INPUT_MANAGER_H
