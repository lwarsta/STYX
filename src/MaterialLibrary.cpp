/**
 * @file MaterialLibrary.cpp
 * @brief Implements functions for parsing material properties from CSV data.
 */

#include "MaterialLibrary.h"
#include <stdexcept>
#include <cstdlib>

/**
 * @brief Parses CSV data to extract 2D material properties.
 *
 * Each row (except the header row) should contain at least 7 columns.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of Material2D structures with parsed values.
 */
std::vector<Material2D> parseMaterial2D(const std::vector<std::vector<std::string>> &csvData) {
    std::vector<Material2D> materials;
    if (csvData.empty())
        return materials;
    
    // Skip header row (assumed row 0).
    for (size_t i = 1; i < csvData.size(); i++) {
        const auto &row = csvData[i];
        if (row.size() < 7)
            continue;
        Material2D m;
        m.name = row[0];
        m.mannings = std::stod(row[1]);
        m.depressionStorage = std::stod(row[2]);
        m.upperStorage = std::stod(row[3]);
        m.evaporationFraction = std::stod(row[4]);
        m.cropFactor = std::stod(row[5]);
        m.rootDepth = std::stod(row[6]);
        materials.push_back(m);
    }
    return materials;
}

/**
 * @brief Parses CSV data to extract 3D material properties.
 *
 * Each row (except the header row) should contain at least 13 columns.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of Material3D structures with parsed values.
 */
std::vector<Material3D> parseMaterial3D(const std::vector<std::vector<std::string>> &csvData) {
    std::vector<Material3D> materials;
    if (csvData.empty())
        return materials;
    
    // Skip header row (assumed row 0).
    for (size_t i = 1; i < csvData.size(); i++) {
        const auto &row = csvData[i];
        if (row.size() < 13)
            continue;
        Material3D m;
        m.name = row[0];
        m.saturatedHydraulicConductivity = std::stod(row[1]);
        m.compressibility = std::stod(row[2]);
        m.saturatedWaterContent = std::stod(row[3]);
        m.residualWaterContent = std::stod(row[4]);
        m.vanGenuchtenAlpha = std::stod(row[5]);
        m.vanGenuchtenN = std::stod(row[6]);
        m.dryWeight = std::stod(row[7]);
        m.heatCapacity = std::stod(row[8]);
        m.heatConductivity = std::stod(row[9]);
        m.heatConductivityMultiplier = std::stod(row[10]);
        m.dispersivityTrans = std::stod(row[11]);
        m.dispersivityLong = std::stod(row[12]);
        materials.push_back(m);
    }
    return materials;
}
