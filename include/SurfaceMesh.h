/**
 * @file SurfaceMesh.h
 * @brief Defines the SurfaceMesh and SurfaceCell classes.
 *
 * These classes provide the geometry and connectivity for 2D surface meshes,
 * including methods to compute cell centers, areas, side lengths, and average slopes.
 */

#ifndef SURFACE_MESH_H
#define SURFACE_MESH_H

#include "Mesh.h"
#include "Timer.h"
#include "VectorHash.h"
#include "PhysicsState.h"
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <sstream>
#include <iomanip>

/**
 * @brief Represents a cell in a 2D surface mesh.
 *
 * A SurfaceCell is a specialization of the Cell base class that stores additional
 * information about the cell's edges (sides) and neighbor relationships.
 */
class SurfaceCell : public Cell {
public:
    std::vector<std::vector<int>> sideVertexIndices; ///< For each side (edge), stores the two vertex indices.

    /**
     * @brief Computes and assigns the side vertex indices.
     *
     * For triangular cells (vtk type 5) or quadrilateral cells (vtk type 9), the 
     * appropriate pairs of vertex indices are stored in sideVertexIndices. Also, 
     * the neighborIndices vector is initialized to -1 for each side.
     */
    void assignSideVertIndices() {
        sideVertexIndices.clear();
        // Triangle (vtk type 5)
        if (vertexIndices.size() == 3) {
            sideVertexIndices.push_back({ vertexIndices[2], vertexIndices[0] });
            sideVertexIndices.push_back({ vertexIndices[0], vertexIndices[1] });
            sideVertexIndices.push_back({ vertexIndices[1], vertexIndices[2] });
        }
        // Quadrilateral (vtk type 9)
        else if (vertexIndices.size() == 4) {
            sideVertexIndices.push_back({ vertexIndices[2], vertexIndices[3] });
            sideVertexIndices.push_back({ vertexIndices[3], vertexIndices[0] });
            sideVertexIndices.push_back({ vertexIndices[0], vertexIndices[1] });
            sideVertexIndices.push_back({ vertexIndices[1], vertexIndices[2] });
        }
        
        neighborIndices.assign(sideVertexIndices.size(), -1);
    }

    /**
     * @brief Computes the center (centroid) of the cell.
     * @param globalPoints The list of all points in the mesh.
     * @return The computed centroid as a Point3D.
     */
    Point3D computeCenter(const std::vector<Point3D>& globalPoints) const {
        Point3D center{0.0, 0.0, 0.0};
        for (int idx : vertexIndices) {
            center.x += globalPoints[idx].x;
            center.y += globalPoints[idx].y;
            center.z += globalPoints[idx].z;
        }
        double n = static_cast<double>(vertexIndices.size());
        center.x /= n;
        center.y /= n;
        center.z /= n;
        return center;
    }

    /**
     * @brief Computes the area of the cell as a planar polygon (ignoring z).
     * @param globalPoints The list of all points in the mesh.
     * @return The computed area.
     */
    double computeArea(const std::vector<Point3D>& globalPoints) const {
        double area = 0.0;
        size_t n = vertexIndices.size();
        for (size_t i = 0; i < n; i++) {
            const Point3D &p1 = globalPoints[vertexIndices[i]];
            const Point3D &p2 = globalPoints[vertexIndices[(i + 1) % n]];
            area += p1.x * p2.y - p2.x * p1.y;
        }
        return std::fabs(area) / 2.0;
    }

    /**
     * @brief Computes the lengths of each side (edge) of the cell.
     * @param globalPoints The list of all points in the mesh.
     * @return A vector containing the lengths of each side.
     */
    std::vector<double> computeSideLengths(const std::vector<Point3D>& globalPoints) const {
        std::vector<double> lengths;
        size_t n = vertexIndices.size();
        for (size_t i = 0; i < n; i++) {
            const Point3D &p1 = globalPoints[vertexIndices[i]];
            const Point3D &p2 = globalPoints[vertexIndices[(i + 1) % n]];
            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;
            double dz = p2.z - p1.z;
            lengths.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
        }
        return lengths;
    }

    /**
     * @brief Computes the average slope of the cell based on neighboring cells.
     *
     * The average slope is calculated from the slopes between the cell center and the centers of
     * all valid neighboring cells. If no valid neighbor exists, a minimum slope of 0.01 is returned.
     *
     * @param globalPoints The list of all points in the mesh.
     * @param cells The list of all cells in the mesh.
     * @return The average slope.
     */
    double computeAverageSlope(const std::vector<Point3D>& globalPoints,
                               const std::vector<std::unique_ptr<Cell>>& cells) const {
        double sumSlope = 0.0;
        int validCount = 0;
        
        // Loop through neighbor indices.
        for (size_t i = 0; i < neighborIndices.size(); i++) {
            int neighIdx = neighborIndices[i];
            if (neighIdx < 0 || static_cast<size_t>(neighIdx) >= cells.size())
                continue;
            const SurfaceCell* neighbor = dynamic_cast<const SurfaceCell*>(cells[neighIdx].get());
            if (!neighbor)
                continue;
            
            // Compute centers.
            Point3D myCenter = computeCenter(globalPoints);
            Point3D neighborCenter = neighbor->computeCenter(globalPoints);
            
            // Compute horizontal distance.
            double dx = neighborCenter.x - myCenter.x;
            double dy = neighborCenter.y - myCenter.y;
            double horizontalDist = std::sqrt(dx * dx + dy * dy);
            if (horizontalDist <= 0.0)
                continue;
            
            // Compute slope.
            double dz = std::fabs(neighborCenter.z - myCenter.z);
            double slope = dz / horizontalDist;
            sumSlope += slope;
            validCount++;
        }
        
        // Return minimum slope if no valid neighbor is found.
        if (validCount == 0)
            return 0.01;
        
        return sumSlope / validCount;
    }
};

/**
 * @brief Represents a 2D surface mesh.
 *
 * The SurfaceMesh class encapsulates a collection of SurfaceCell objects and the associated
 * geometric data. It provides methods to build connectivity and export the mesh in VTK format.
 */
class SurfaceMesh : public Mesh {
public:
    /**
     * @brief Builds the connectivity of the surface mesh.
     *
     * This method creates an unordered map of sorted edge indices to cell indices, and then
     * assigns neighbor indices for each cell based on shared edges.
     */
    void buildConnectivity() override {
        Timer timer("SurfaceMesh::buildConnectivity");
        std::unordered_map<std::vector<int>, std::vector<size_t>, VectorHash> edgeMap;
        size_t numCells = cells.size();

        // First pass: compute edges and populate the map.
        for (size_t i = 0; i < numCells; i++) {
            SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cells[i].get());
            if (!cell) continue;
            cell->assignSideVertIndices();
            cell->neighborIndices.clear();
            for (const auto &edge : cell->sideVertexIndices) {
                std::vector<int> sortedEdge = edge;
                std::sort(sortedEdge.begin(), sortedEdge.end());
                edgeMap[sortedEdge].push_back(i);
            }
        }

        // Second pass: assign neighbor indices for each cell.
        for (size_t i = 0; i < numCells; i++) {
            SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cells[i].get());
            if (!cell) continue;
            cell->neighborIndices.assign(cell->sideVertexIndices.size(), -1);
            for (size_t e = 0; e < cell->sideVertexIndices.size(); e++) {
                std::vector<int> sortedEdge = cell->sideVertexIndices[e];
                std::sort(sortedEdge.begin(), sortedEdge.end());
                auto it = edgeMap.find(sortedEdge);
                if (it != edgeMap.end()) {
                    const std::vector<size_t>& indices = it->second;
                    if (indices.size() == 1) {
                        cell->neighborIndices[e] = -1;
                    } else if (indices.size() == 2) {
                        cell->neighborIndices[e] = (indices[0] == i) ? static_cast<int>(indices[1]) : static_cast<int>(indices[0]);
                    } else {
                        for (size_t idx : indices) {
                            if (idx != i) {
                                cell->neighborIndices[e] = static_cast<int>(idx);
                                break;
                            }
                        }
                    }
                } else {
                    cell->neighborIndices[e] = -1;
                }
            }
        }
    }

    /**
     * @brief Exports the surface mesh in VTK format.
     *
     * The VTK output includes point coordinates, cell connectivity, cell types,
     * and cell data (such as material, upper storage water depth, and overland water depth).
     *
     * @return A string containing the VTK representation of the mesh.
     */
    std::string exportVTK() const {
        std::stringstream ss;
        ss << std::setprecision(4) << std::fixed;
        ss << "# vtk DataFile Version 2.0\n";
        ss << "STYX Surface Mesh.\n";
        ss << "ASCII\n";
        ss << "DATASET UNSTRUCTURED_GRID\n";
        ss << "POINTS " << points.size() << " double\n";
        for (const auto &pt : points) {
            ss << pt.x << " " << pt.y << " " << pt.z << "\n";
        }
        // Write cells.
        ss << "CELLS " << cells.size() << " ";
        int totalCellInts = 0;
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (cell)
                totalCellInts += cell->vertexIndices.size() + 1;
        }
        ss << totalCellInts << "\n";
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (!cell) continue;
            ss << cell->vertexIndices.size();
            for (int idx : cell->vertexIndices)
                ss << " " << idx;
            ss << "\n";
        }
        // Write cell types (VTK type 5 for triangle, 9 for quadrilateral).
        ss << "CELL_TYPES " << cells.size() << "\n";
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (!cell) continue;
            if (cell->vertexIndices.size() == 3)
                ss << "5\n";
            else if (cell->vertexIndices.size() == 4)
                ss << "9\n";
            else
                ss << "0\n"; // Unknown type.
        }
        // Write cell data.
        ss << "CELL_DATA " << cells.size() << "\n";
        ss << "FIELD FieldData 3\n";
        ss << "material 1 " << cells.size() << " int\n";
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (cell)
                ss << cell->material << "\n";
        }
        ss << "upper_storage_water_depth 1 " << cells.size() << " double\n";
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (!cell || !cell->waterState) {
                ss << "0.0\n";
                continue;
            }
            auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
            ss << (sws ? sws->upperStorageDepth : 0.0) << "\n";
        }
        ss << "overland_water_depth 1 " << cells.size() << " double\n";
        for (const auto &cellPtr : cells) {
            const SurfaceCell* cell = dynamic_cast<const SurfaceCell*>(cellPtr.get());
            if (!cell || !cell->waterState) {
                ss << "0.0\n";
                continue;
            }
            auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
            ss << (sws ? sws->waterDepth : 0.0) << "\n";
        }
        
        return ss.str();
    }
};

#endif // SURFACE_MESH_H
