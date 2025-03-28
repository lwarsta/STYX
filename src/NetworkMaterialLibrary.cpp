/**
 * @file NetworkMaterialLibrary.cpp
 * @brief Implementation of CSV parsing functions for network material properties.
 */

#include "NetworkMaterialLibrary.h"
#include <stdexcept>
#include <cstdlib>

/**
 * @brief Parses CSV data to extract network junction material properties.
 *
 * Skips the header row (assumed row 0) and parses each subsequent row.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of NetworkJunctionMaterial structures.
 */
std::vector<NetworkJunctionMaterial> parseNetworkJunctionMaterials(
    const std::vector<std::vector<std::string>> &csvData)
{
    std::vector<NetworkJunctionMaterial> materials;
    if (csvData.empty())
        return materials;
    
    // Skip header row.
    for (size_t i = 1; i < csvData.size(); i++) {
        const auto &row = csvData[i];
        if (row.size() < 2)
            continue;
        NetworkJunctionMaterial m;
        m.name = row[0];
        m.diameter = std::stod(row[1]);
        materials.push_back(m);
    }
    return materials;
}

/**
 * @brief Parses CSV data to extract network link material properties.
 *
 * Skips the header row (assumed row 0) and parses each subsequent row.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of NetworkLinkMaterial structures.
 */
std::vector<NetworkLinkMaterial> parseNetworkLinkMaterials(
    const std::vector<std::vector<std::string>> &csvData)
{
    std::vector<NetworkLinkMaterial> materials;
    if (csvData.empty())
        return materials;
    
    // Skip header row.
    for (size_t i = 1; i < csvData.size(); i++) {
        const auto &row = csvData[i];
        if (row.size() < 3)
            continue;
        NetworkLinkMaterial m;
        m.name = row[0];
        m.diameter = std::stod(row[1]);
        m.roughness = std::stod(row[2]);
        materials.push_back(m);
    }
    return materials;
}
