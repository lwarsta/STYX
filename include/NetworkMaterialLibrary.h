/**
 * @file NetworkMaterialLibrary.h
 * @brief Material library definitions and CSV parsing functions for network components.
 *
 * This file defines the material property structures for network junctions (wells)
 * and network links (pipes) and declares functions to parse CSV data into these structures.
 */

#ifndef NETWORK_MATERIAL_LIBRARY_H
#define NETWORK_MATERIAL_LIBRARY_H

#include <string>
#include <vector>

/**
 * @brief Structure representing material properties for network junctions (wells).
 */
struct NetworkJunctionMaterial {
    std::string name;   ///< Material name.
    double diameter;    ///< Diameter in meters.
};

/**
 * @brief Structure representing material properties for network links (pipes).
 */
struct NetworkLinkMaterial {
    std::string name;   ///< Material name.
    double diameter;    ///< Diameter in meters.
    double roughness;   ///< Manning’s roughness coefficient.
};

/**
 * @brief Parses CSV data for network junction materials.
 *
 * This function assumes the first row of the CSV data is a header and skips it.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of NetworkJunctionMaterial structures.
 */
std::vector<NetworkJunctionMaterial> parseNetworkJunctionMaterials(
    const std::vector<std::vector<std::string>> &csvData);

/**
 * @brief Parses CSV data for network link materials.
 *
 * This function assumes the first row of the CSV data is a header and skips it.
 *
 * @param csvData A 2D vector of strings representing the CSV file content.
 * @return A vector of NetworkLinkMaterial structures.
 */
std::vector<NetworkLinkMaterial> parseNetworkLinkMaterials(
    const std::vector<std::vector<std::string>> &csvData);

#endif // NETWORK_MATERIAL_LIBRARY_H
