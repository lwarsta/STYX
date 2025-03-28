/**
 * @file MeshFactory.h
 * @brief Factory for creating surface and subsurface meshes.
 *
 * This file defines the MeshFactory class, which provides static functions
 * to create SurfaceMesh and SubsurfaceMesh objects from VTK and CSV input files.
 */

#ifndef MESH_FACTORY_H
#define MESH_FACTORY_H

#include "SurfaceMesh.h"
#include "SubsurfaceMesh.h"
#include "InputManager.h"
#include "PhysicsState.h"
#include "MaterialLibrary.h"
#include <string>
#include <memory>

/**
 * @brief MeshFactory class.
 *
 * Provides static functions to create surface and subsurface meshes by loading VTK files
 * for geometry and CSV files for material properties, initial conditions, and boundary conditions.
 */
class MeshFactory {
public:
    /**
     * @brief Creates a SurfaceMesh from the provided VTK and CSV files.
     *
     * This function loads the surface grid from a VTK file, then reads material properties,
     * initial conditions, and boundary conditions from CSV files to initialize the mesh.
     *
     * @param vtkPath Path to the VTK file for the surface grid.
     * @param materialsPath Path to the CSV file containing 2D material properties.
     * @param initCondPath Path to the CSV file containing initial conditions for surface cells.
     * @param boundCondPath Path to the CSV file containing boundary conditions for surface cells.
     * @return std::unique_ptr<SurfaceMesh> The created and initialized surface mesh.
     */
    static std::unique_ptr<SurfaceMesh> createSurfaceMesh(const std::string &vtkPath,
                                                            const std::string &materialsPath,
                                                            const std::string &initCondPath,
                                                            const std::string &boundCondPath) {
        // Load the VTK grid.
        VTKGrid vtkGrid = InputManager::loadVTK(vtkPath);
        auto mesh = std::make_unique<SurfaceMesh>();
        mesh->points = vtkGrid.points;
		
        // Create a SurfaceCell for each VTK cell.
        for (const auto &vtkCell : vtkGrid.cells) {
            auto cell = std::make_unique<SurfaceCell>();
            
            // Set cell properties.
            cell->vertexIndices = vtkCell.connectivity;
            cell->cellType = vtkCell.cellType;
            cell->id = vtkCell.id;
            cell->gridConnection = vtkCell.gridConnection;
            cell->material = vtkCell.material;
            cell->initCond = vtkCell.initCond;
            cell->boundCond = vtkCell.boundCond;
			
            // Initialize the physics state.
            cell->waterState = std::make_shared<SurfaceWaterState>();
			
            mesh->cells.push_back(std::move(cell));
        }
		
        // Build connectivity between cells.
        mesh->buildConnectivity();
		
        // Load material properties for surface cells.
        auto material2DData = InputManager::loadCSV(materialsPath);
        auto materials2D = parseMaterial2D(material2DData);

        for (size_t i = 0; i < mesh->cells.size(); i++) {
            SurfaceCell* surfCell = dynamic_cast<SurfaceCell*>(mesh->cells[i].get());
            if (surfCell && surfCell->waterState) {
                int matIndex = surfCell->material;
                if (matIndex >= 0 && matIndex < static_cast<int>(materials2D.size())) {
                    const Material2D &mat = materials2D[matIndex];
                    auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(surfCell->waterState);
                    if (sws) {
                        sws->mannings = mat.mannings;
                        sws->depressionStorage = mat.depressionStorage;
                        sws->upperStorage = mat.upperStorage;
                        sws->evaporationFraction = mat.evaporationFraction;
                        sws->cropFactor = mat.cropFactor;
                        sws->rootDepth = mat.rootDepth;
                    }
                }
            }
        }
		
        // Load initial conditions for surface cells.
        auto initCond = InputManager::loadCSV(initCondPath);
        for (size_t i = 1; i < initCond.size(); i++) {
            int cellId = std::stoi(initCond[i][0]);
            double waterDepth = std::stod(initCond[i][1]);
            if (cellId >= 0 && cellId < static_cast<int>(mesh->cells.size())) {
                SurfaceCell* cell = dynamic_cast<SurfaceCell*>(mesh->cells[cellId].get());
                if (cell && cell->waterState) {
                    auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
                    sws->waterDepth = waterDepth;
                    sws->waterDepthOld = waterDepth;
                }
            }
        }
		
        // Load boundary conditions for surface cells.
        auto boundCond = InputManager::loadCSV(boundCondPath);
        for (size_t i = 1; i < boundCond.size(); i++) {
            int cellId = std::stoi(boundCond[i][0]);
            int outletId = std::stoi(boundCond[i][1]);
            double distanceToOutlet = std::stod(boundCond[i][2]);
            int sinkId = std::stoi(boundCond[i][3]);
            if (cellId >= 0 && cellId < static_cast<int>(mesh->cells.size())) {
                SurfaceCell* cell = dynamic_cast<SurfaceCell*>(mesh->cells[cellId].get());
                if (cell && cell->waterState) {
                    auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
                    sws->outletId = outletId;
                    sws->distanceToOutlet = distanceToOutlet;
                    sws->sinkId = sinkId;
                }
            }
        }
		
        // Initialize physics-related intercell variables (if applicable).
        // For surface mesh, this might be skipped if not needed.
        // (Note: This loop is here to match the original structure; adjust as needed.)
        for (auto &cellPtr : mesh->cells) {
            auto cell = dynamic_cast<VolumeCell*>(cellPtr.get());
            if (!cell) continue;
            auto waterState = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
            if (waterState)
                waterState->initializeIntercellVariables(cell->neighborIndices.size());
        }
		
        return mesh;
    }

    /**
     * @brief Creates a SubsurfaceMesh from the provided VTK and CSV files.
     *
     * This function loads the subsurface grid from a VTK file, then reads material properties,
     * initial conditions, and boundary conditions from CSV files to initialize the mesh.
     *
     * @param vtkPath Path to the VTK file for the subsurface grid.
     * @param materialsPath Path to the CSV file containing 3D material properties.
     * @param initCondPath Path to the CSV file containing initial conditions for subsurface cells.
     * @param boundCondPath Path to the CSV file containing boundary conditions for subsurface cells.
     * @return std::unique_ptr<SubsurfaceMesh> The created and initialized subsurface mesh.
     */
    static std::unique_ptr<SubsurfaceMesh> createSubsurfaceMesh(const std::string &vtkPath,
                                                                  const std::string &materialsPath,
                                                                  const std::string &initCondPath,
                                                                  const std::string &boundCondPath) {
        // Load the VTK grid.
        VTKGrid vtkGrid = InputManager::loadVTK(vtkPath);
        auto mesh = std::make_unique<SubsurfaceMesh>();
        mesh->points = vtkGrid.points;
		
        // Create a VolumeCell for each VTK cell.
        for (const auto &vtkCell : vtkGrid.cells) {
            auto cell = std::make_unique<VolumeCell>();
			
            // Set cell properties.
            cell->vertexIndices = vtkCell.connectivity;
            cell->cellType = vtkCell.cellType;
            cell->id = vtkCell.id;
            cell->gridConnection = vtkCell.gridConnection;
            cell->material = vtkCell.material;
            cell->initCond = vtkCell.initCond;
            cell->boundCond = vtkCell.boundCond;
			
            // Initialize the physics state.
            cell->waterState = std::make_shared<SubsurfaceWaterState>();
            
            mesh->cells.push_back(std::move(cell));
        }
		
        // Build connectivity between cells.
        mesh->buildConnectivity();
		
        // Load material properties for subsurface cells.
        auto material3DData = InputManager::loadCSV(materialsPath);
        auto materials3D = parseMaterial3D(material3DData);

        for (size_t i = 0; i < mesh->cells.size(); i++) {
            VolumeCell* volCell = dynamic_cast<VolumeCell*>(mesh->cells[i].get());
            if (volCell && volCell->waterState) {
                int matIndex = volCell->material;
                if (matIndex >= 0 && matIndex < static_cast<int>(materials3D.size())) {
                    const Material3D &mat = materials3D[matIndex];
                    auto subWater = std::dynamic_pointer_cast<SubsurfaceWaterState>(volCell->waterState);
                    if (subWater) {
                        subWater->saturatedHydraulicConductivity = mat.saturatedHydraulicConductivity;
                        subWater->compressibility = mat.compressibility;
                        subWater->saturatedWaterContent = mat.saturatedWaterContent;
                        subWater->residualWaterContent = mat.residualWaterContent;
                        subWater->vanGenuchtenAlpha = mat.vanGenuchtenAlpha;
                        subWater->vanGenuchtenN = mat.vanGenuchtenN;
                        subWater->dryWeight = mat.dryWeight;
                        subWater->heatCapacity = mat.heatCapacity;
                        subWater->heatConductivity = mat.heatConductivity;
                        subWater->heatConductivityMultiplier = mat.heatConductivityMultiplier;
                        subWater->dispersivityTrans = mat.dispersivityTrans;
                        subWater->dispersivityLong = mat.dispersivityLong;
						
                        // Compute van Genuchten m.
                        if (subWater->vanGenuchtenN > 0.0)
                            subWater->vanGenuchtenM = 1.0 - 1.0 / subWater->vanGenuchtenN;
                    }
                }
            }
        }
		
        // Load initial conditions for subsurface cells.
        auto initCond = InputManager::loadCSV(initCondPath);
        for (size_t i = 1; i < initCond.size(); i++) {
            int cellId = std::stoi(initCond[i][0]);
            double hydraulicHead = std::stod(initCond[i][1]);
            if (cellId >= 0 && cellId < static_cast<int>(mesh->cells.size())) {
                VolumeCell* cell = dynamic_cast<VolumeCell*>(mesh->cells[cellId].get());
                if (cell && cell->waterState) {
                    std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState)->hydraulicHead = hydraulicHead;
                }
            }
        }

        // Load boundary conditions for subsurface cells.
        auto boundCond = InputManager::loadCSV(boundCondPath);
        for (size_t i = 1; i < boundCond.size(); i++) {
            int cellId = std::stoi(boundCond[i][0]);
            if (cellId >= 0 && cellId < static_cast<int>(mesh->cells.size())) {
                VolumeCell* cell = dynamic_cast<VolumeCell*>(mesh->cells[cellId].get());
                if (cell && cell->waterState) {
                    // Boundary conditions can be processed here if needed.
                }
            }
        }
		
        return mesh;
    }
};

#endif // MESH_FACTORY_H
