/**
 * @file MaterialLibrary.h
 * @brief Structures and functions for parsing material properties from CSV data.
 *
 * This file defines structures for 2D and 3D material properties and declares functions
 * to parse CSV data into these structures.
 */

#ifndef MATERIAL_LIBRARY_H
#define MATERIAL_LIBRARY_H

#include <string>
#include <vector>

/**
 * @brief Material properties for 2D cells.
 */
struct Material2D {
    std::string name;             ///< Material name.
    double mannings;              ///< Manning's roughness coefficient.
    double depressionStorage;     ///< Depression storage capacity.
    double upperStorage;          ///< Upper storage capacity.
    double evaporationFraction;   ///< Evaporation fraction.
    double cropFactor;            ///< Crop factor.
    double rootDepth;             ///< Root depth.
};

/**
 * @brief Material properties for 3D cells.
 */
struct Material3D {
    std::string name;                         ///< Material name.
    double saturatedHydraulicConductivity;    ///< Saturated hydraulic conductivity.
    double compressibility;                     ///< Compressibility.
    double saturatedWaterContent;             ///< Saturated water content.
    double residualWaterContent;              ///< Residual water content.
    double vanGenuchtenAlpha;                 ///< van Genuchten parameter alpha.
    double vanGenuchtenN;                     ///< van Genuchten parameter n (must be > 1).
    double dryWeight;                         ///< Dry weight.
    double heatCapacity;                      ///< Heat capacity.
    double heatConductivity;                  ///< Heat conductivity.
    double heatConductivityMultiplier;        ///< Heat conductivity multiplier.
    double dispersivityTrans;                 ///< Transverse dispersivity.
    double dispersivityLong;                  ///< Longitudinal dispersivity.
};

/**
 * @brief Parses CSV data for 2D material properties.
 *
 * The function assumes that the first row of csvData is a header row and skips it.
 *
 * @param csvData A 2D vector of strings representing the CSV data.
 * @return A vector of Material2D structures.
 */
std::vector<Material2D> parseMaterial2D(const std::vector<std::vector<std::string>> &csvData);

/**
 * @brief Parses CSV data for 3D material properties.
 *
 * The function assumes that the first row of csvData is a header row and skips it.
 *
 * @param csvData A 2D vector of strings representing the CSV data.
 * @return A vector of Material3D structures.
 */
std::vector<Material3D> parseMaterial3D(const std::vector<std::vector<std::string>> &csvData);

#endif // MATERIAL_LIBRARY_H
