/**
 * @file NetworkMeshFactory.h
 * @brief Factory for creating network meshes.
 *
 * This file defines the NetworkMeshFactory class, which provides static methods
 * to create and initialize a NetworkMesh object from VTK and CSV files.
 */

#ifndef NETWORK_MESH_FACTORY_H
#define NETWORK_MESH_FACTORY_H

#include "NetworkMesh.h"
#include "InputManager.h"
#include "PhysicsState.h"
#include "NetworkMaterialLibrary.h"
#include <memory>
#include <unordered_map>
#include <string>

/**
 * @brief A factory class to create a NetworkMesh.
 *
 * The NetworkMeshFactory loads junction and link data from VTK files and material,
 * initial condition, and boundary condition data from CSV files. It then constructs
 * a unified mesh with connectivity and assigns material properties and conditions.
 */
class NetworkMeshFactory {
public:
    /**
     * @brief Creates and initializes a NetworkMesh.
     *
     * Loads junction and link geometry from VTK files, maps the point indices to a
     * unified points vector, and creates cells (junctions and links). It then loads
     * material properties, initial conditions, and boundary conditions from CSV files,
     * and finally builds the mesh connectivity.
     *
     * @param juncVtkPath Path to the VTK file for network junctions.
     * @param linkVtkPath Path to the VTK file for network links.
     * @param materialsJuncPath Path to the CSV file for junction material properties.
     * @param materialsLinkPath Path to the CSV file for link material properties.
     * @param initCondJuncPath Path to the CSV file for junction initial conditions.
     * @param boundCondJuncPath Path to the CSV file for junction boundary conditions.
     * @param initCondLinkPath Path to the CSV file for link initial conditions.
     * @param boundCondLinkPath Path to the CSV file for link boundary conditions.
     * @return std::unique_ptr<NetworkMesh> A pointer to the fully initialized NetworkMesh.
     * @throws std::runtime_error if any file cannot be loaded.
     */
    static std::unique_ptr<NetworkMesh> createNetworkMesh(const std::string &juncVtkPath,
                                                            const std::string &linkVtkPath,
                                                            const std::string &materialsJuncPath,
                                                            const std::string &materialsLinkPath,
                                                            const std::string &initCondJuncPath,
                                                            const std::string &boundCondJuncPath,
                                                            const std::string &initCondLinkPath,
                                                            const std::string &boundCondLinkPath) {
        auto networkMesh = std::make_unique<NetworkMesh>();

        // Load junction data from VTK.
        VTKGrid vtkJunc = InputManager::loadVTK(juncVtkPath);

        // Create a mapping from original VTK indices to the unified points vector.
        std::unordered_map<int, int> pointIndexMap;

        for (long unsigned int i = 0; i < vtkJunc.points.size(); i++) {
            // Use the original index as the key.
            networkMesh->points.push_back(vtkJunc.points[i]);
            pointIndexMap[i] = i;
        }

        // Create junction cells.
        for (const auto &vtkCell : vtkJunc.cells) {
            auto cell = std::make_unique<NetworkJunctionCell>();

            // Set cell properties.
            cell->vertexIndices.push_back(pointIndexMap[vtkCell.connectivity[0]]);
            cell->vertexIndices.push_back(pointIndexMap[vtkCell.connectivity[1]]);
            cell->top = networkMesh->points[pointIndexMap[vtkCell.connectivity[0]]];
            cell->bottom = networkMesh->points[pointIndexMap[vtkCell.connectivity[1]]];
            cell->cellType = vtkCell.cellType;
            cell->id = vtkCell.id;
            cell->gridConnection = vtkCell.gridConnection;
            cell->material = vtkCell.material;
            cell->initCond = vtkCell.initCond;
            cell->boundCond = vtkCell.boundCond;

            // Initialize physics state.
            cell->waterState = std::make_shared<NetworkWaterState>();

            networkMesh->cells.push_back(std::move(cell));
        }

        // Load link data from VTK.
        int offset = vtkJunc.points.size();
        VTKGrid vtkLink = InputManager::loadVTK(linkVtkPath);

        for (long unsigned int i = 0; i < vtkLink.points.size(); i++) {
            networkMesh->points.push_back(vtkLink.points[i]);
            // Map the original VTK index to the new index.
            pointIndexMap[i] = i + offset;
        }

        // Create link cells.
        for (const auto &vtkCell : vtkLink.cells) {
            auto cell = std::make_unique<NetworkLinkCell>();

            // Set cell properties.
            int index1 = pointIndexMap[vtkCell.connectivity[0]];
            int index2 = pointIndexMap[vtkCell.connectivity[1]];
            cell->vertexIndices.push_back(index1);
            cell->vertexIndices.push_back(index2);
            cell->end1 = networkMesh->points[index1];
            cell->end2 = networkMesh->points[index2];
            cell->cellType = vtkCell.cellType;
            cell->id = vtkCell.id;
            cell->gridConnection = vtkCell.gridConnection;
            cell->material = vtkCell.material;
            cell->initCond = vtkCell.initCond;
            cell->boundCond = vtkCell.boundCond;

            // Initialize physics state.
            cell->waterState = std::make_shared<NetworkWaterState>();

            networkMesh->cells.push_back(std::move(cell));
        }

        // Build connectivity between cells.
        networkMesh->buildConnectivity();

        // Load junction material properties and assign them.
        auto materialsJuncCsv = InputManager::loadCSV(materialsJuncPath);
        auto materialsJunc = parseNetworkJunctionMaterials(materialsJuncCsv);
        for (auto &cellPtr : networkMesh->cells) {
            auto jnc = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
            if (!jnc)
                continue;
            int matIndex = jnc->material;
            if (matIndex >= 0 && matIndex < static_cast<int>(materialsJunc.size())) {
                jnc->diameter = materialsJunc[matIndex].diameter;
            }
        }

        // Load link material properties and assign them.
        auto materialsLinkCsv = InputManager::loadCSV(materialsLinkPath);
        auto materialsLink = parseNetworkLinkMaterials(materialsLinkCsv);
        for (auto &cellPtr : networkMesh->cells) {
            auto lnk = dynamic_cast<NetworkLinkCell*>(cellPtr.get());
            if (!lnk)
                continue;
            int matIndex = lnk->material;
            if (matIndex >= 0 && matIndex < static_cast<int>(materialsLink.size())) {
                lnk->diameter = materialsLink[matIndex].diameter;
                lnk->roughness = materialsLink[matIndex].roughness; // Ideally, roughness belongs in waterState.
            }
        }

        // Load junction initial condition values and assign them.
        auto initCondJunc = InputManager::loadCSV(initCondJuncPath);
        for (size_t i = 1; i < initCondJunc.size(); i++) {
            int juncId = std::stoi(initCondJunc[i][0]);
            double waterDepth = std::stod(initCondJunc[i][1]);
            for (auto &cellPtr : networkMesh->cells) {
                auto jnc = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
                if (!jnc || jnc->id != juncId || !jnc->waterState)
                    continue;
                auto netState = std::dynamic_pointer_cast<NetworkWaterState>(jnc->waterState);
                if (netState) {
                    netState->waterLevel = waterDepth;
                    netState->waterLevelOld = waterDepth;
                }
                break;
            }
        }

        // Load junction boundary condition values and assign them.
        auto boundCondJunc = InputManager::loadCSV(boundCondJuncPath);
        for (size_t i = 1; i < boundCondJunc.size(); i++) {
            int juncId = std::stoi(boundCondJunc[i][0]);
            int type = std::stoi(boundCondJunc[i][1]);
            for (auto &cellPtr : networkMesh->cells) {
                auto jnc = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
                if (!jnc || jnc->id != juncId)
                    continue;
                jnc->type = type;
                break;
            }
        }

        // Load link initial condition values and assign them.
        auto initCondLink = InputManager::loadCSV(initCondLinkPath);
        for (size_t i = 1; i < initCondLink.size(); i++) {
            int linkId = std::stoi(initCondLink[i][0]);
            for (auto &cellPtr : networkMesh->cells) {
                auto lnk = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
                if (!lnk || lnk->id != linkId || !lnk->waterState)
                    continue;
                // Currently, nothing is saved for link initial conditions.
                break;
            }
        }

        // Load link boundary condition values and assign them.
        auto boundCondLink = InputManager::loadCSV(boundCondLinkPath);
        for (size_t i = 1; i < boundCondLink.size(); i++) {
            int linkId = std::stoi(boundCondLink[i][0]);
            for (auto &cellPtr : networkMesh->cells) {
                auto lnk = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
                if (!lnk || lnk->id != linkId)
                    continue;
                // Currently, nothing is saved for link boundary conditions.
                break;
            }
        }

        return networkMesh;
    }
};

#endif // NETWORK_MESH_FACTORY_H
