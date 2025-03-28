/**
 * @file SubsurfaceMesh.h
 * @brief Defines the VolumeCell and SubsurfaceMesh classes for 3D grids.
 *
 * These classes provide methods for computing geometric properties of
 * volume cells and building connectivity in a 3D mesh.
 */

#ifndef SUBSURFACE_MESH_H
#define SUBSURFACE_MESH_H

#include "Mesh.h"
#include "Timer.h"
#include "VectorHash.h"
#include "PhysicsState.h"
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <iomanip>

/**
 * @brief Computes the Euclidean distance between two 3D points.
 * 
 * @param p1 The first point.
 * @param p2 The second point.
 * @return double The Euclidean distance.
 */
inline double euclideanDistance(const Point3D &p1, const Point3D &p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Represents a volume cell in a 3D subsurface mesh.
 *
 * A VolumeCell contains a list of faces (each a vector of vertex indices)
 * and methods to compute geometric properties such as cell center, face centers,
 * thickness, and distance to a neighboring cell.
 */
class VolumeCell : public Cell {
public:
    /// Each face is represented as a list of vertex indices.
    std::vector<std::vector<int>> faceVertexIndices;

    /**
     * @brief Computes the face indices based on the cell type.
     *
     * For hexahedra (8 vertices, VTK type 12) and wedges (6 vertices, VTK type 13),
     * the face indices are computed and stored in faceVertexIndices. The neighborIndices
     * vector is then resized (and initialized to -1) to match the number of faces.
     */
    void assignFaceVertIndices() {
        faceVertexIndices.clear();
        // Hexahedron (vtk type 12)
        if (vertexIndices.size() == 8) {
            faceVertexIndices.push_back({ vertexIndices[4], vertexIndices[5], vertexIndices[6], vertexIndices[7] }); // top
            faceVertexIndices.push_back({ vertexIndices[2], vertexIndices[3], vertexIndices[7], vertexIndices[6] }); // back
            faceVertexIndices.push_back({ vertexIndices[3], vertexIndices[0], vertexIndices[4], vertexIndices[7] }); // left
            faceVertexIndices.push_back({ vertexIndices[0], vertexIndices[1], vertexIndices[5], vertexIndices[4] }); // front
            faceVertexIndices.push_back({ vertexIndices[1], vertexIndices[2], vertexIndices[6], vertexIndices[5] }); // right
            faceVertexIndices.push_back({ vertexIndices[3], vertexIndices[2], vertexIndices[1], vertexIndices[0] }); // bottom
        }
        // Wedge (vtk type 13)
        else if (vertexIndices.size() == 6) {
            faceVertexIndices.push_back({ vertexIndices[3], vertexIndices[4], vertexIndices[5] }); // top
            faceVertexIndices.push_back({ vertexIndices[0], vertexIndices[3], vertexIndices[5], vertexIndices[2] }); // side 1
            faceVertexIndices.push_back({ vertexIndices[1], vertexIndices[4], vertexIndices[3], vertexIndices[0] }); // side 2
            faceVertexIndices.push_back({ vertexIndices[2], vertexIndices[5], vertexIndices[4], vertexIndices[0] }); // side 3
            faceVertexIndices.push_back({ vertexIndices[0], vertexIndices[2], vertexIndices[1] }); // bottom
        }

        neighborIndices.assign(faceVertexIndices.size(), -1);
    }

    /**
     * @brief Computes the center (centroid) of the cell.
     *
     * @param globalPoints The vector of all points in the mesh.
     * @return Point3D The centroid of the cell.
     */
    Point3D computeCenter(const std::vector<Point3D>& globalPoints) const {
        Point3D center{ 0.0, 0.0, 0.0 };
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
     * @brief Computes the center of a face given its vertex indices.
     *
     * @param faceIndices The indices of the vertices of the face.
     * @param globalPoints The vector of all points in the mesh.
     * @return Point3D The center of the face.
     */
    Point3D computeFaceCenter(const std::vector<int>& faceIndices, const std::vector<Point3D>& globalPoints) const {
        Point3D center{ 0.0, 0.0, 0.0 };
        for (int idx : faceIndices) {
            center.x += globalPoints[idx].x;
            center.y += globalPoints[idx].y;
            center.z += globalPoints[idx].z;
        }
        double n = static_cast<double>(faceIndices.size());
        center.x /= n;
        center.y /= n;
        center.z /= n;
        return center;
    }

    /**
     * @brief Computes the thickness of the cell.
     *
     * The thickness is defined as the absolute vertical distance between the centers
     * of the top and bottom faces. It is assumed that the first face is the top face and
     * the last face is the bottom face.
     *
     * @param globalPoints The vector of all points in the mesh.
     * @return double The cell thickness.
     */
    double computeThickness(const std::vector<Point3D>& globalPoints) const {
        if (faceVertexIndices.empty())
            return 0.0;
        // Assume the first face is the top and the last face is the bottom.
        Point3D topCenter = computeFaceCenter(faceVertexIndices.front(), globalPoints);
        Point3D bottomCenter = computeFaceCenter(faceVertexIndices.back(), globalPoints);
        return std::fabs(topCenter.z - bottomCenter.z);
    }

    /**
     * @brief Computes the distance to a neighbor cell.
     *
     * This function returns the distance between the centers of this cell and a neighbor,
     * taking into account common faces if available. If no common face is found, the direct
     * center-to-center distance is returned.
     *
     * @param neighborIdx The index in the neighborIndices vector.
     * @param globalPoints The vector of all points in the mesh.
     * @param cells The vector of all cells.
     * @return double The distance to the neighbor cell.
     */
    double getDistanceToNeighbor(size_t neighborIdx,
                                 const std::vector<Point3D>& globalPoints,
                                 const std::vector<std::unique_ptr<Cell>>& cells) const {
        if (neighborIdx >= neighborIndices.size())
            return std::numeric_limits<double>::max();
        size_t nIndex = neighborIndices[neighborIdx];
        if (nIndex >= cells.size())
            return std::numeric_limits<double>::max();
        const VolumeCell* neighborCell = dynamic_cast<const VolumeCell*>(cells[nIndex].get());
        if (!neighborCell)
            return std::numeric_limits<double>::max();
        
        // Compute centers of this cell and the neighbor.
        Point3D myCenter = computeCenter(globalPoints);
        Point3D neighborCenter = neighborCell->computeCenter(globalPoints);
        double minDistance = std::numeric_limits<double>::max();
        bool foundCommonFace = false;
        
        // Loop over each face of this cell.
        for (const auto &face : faceVertexIndices) {
            // Create a sorted copy of the face.
            std::vector<int> sortedFace = face;
            std::sort(sortedFace.begin(), sortedFace.end());
            // Check each face of the neighbor cell.
            for (const auto &nFace : neighborCell->faceVertexIndices) {
                std::vector<int> sortedNFace = nFace;
                std::sort(sortedNFace.begin(), sortedNFace.end());
                if (sortedFace == sortedNFace) {
                    foundCommonFace = true;
                    // Compute the center of the common face.
                    Point3D faceCenter = computeFaceCenter(face, globalPoints);
                    double d1 = euclideanDistance(myCenter, faceCenter);
                    double d2 = euclideanDistance(faceCenter, neighborCenter);
                    double total = d1 + d2;
                    if (total < minDistance)
                        minDistance = total;
                }
            }
        }
        if (!foundCommonFace) {
            // Fallback: return direct center-to-center distance.
            return euclideanDistance(myCenter, neighborCenter);
        }
        return minDistance;
    }

    /**
     * @brief Computes the area of a planar polygon defined by a set of 3D points.
     *
     * This function triangulates the polygon using its centroid as a reference and
     * sums the areas of the resulting triangles.
     *
     * @param pts A vector of 3D points representing the polygon vertices.
     * @return double The computed area of the polygon.
     */
    double computePolygonArea3D(const std::vector<Point3D>& pts) const {
        if (pts.size() < 3)
            return 0.0;
        // Compute the centroid of the polygon.
        Point3D centroid{0.0, 0.0, 0.0};
        for (const auto &p : pts) {
            centroid.x += p.x;
            centroid.y += p.y;
            centroid.z += p.z;
        }
        double n = static_cast<double>(pts.size());
        centroid.x /= n;
        centroid.y /= n;
        centroid.z /= n;

        double area = 0.0;
        // Triangulate the polygon: form triangles (centroid, pts[i], pts[(i+1)%n])
        for (size_t i = 0; i < pts.size(); i++) {
            const Point3D &p1 = pts[i];
            const Point3D &p2 = pts[(i + 1) % pts.size()];
            double ux = p1.x - centroid.x;
            double uy = p1.y - centroid.y;
            double uz = p1.z - centroid.z;
            double vx = p2.x - centroid.x;
            double vy = p2.y - centroid.y;
            double vz = p2.z - centroid.z;
            double crossX = uy * vz - uz * vy;
            double crossY = uz * vx - ux * vz;
            double crossZ = ux * vy - uy * vx;
            double triArea = 0.5 * std::sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
            area += triArea;
        }
        return area;
    }
    
    /**
     * @brief Computes the areas of all faces of the cell.
     *
     * For each face (represented by its vertex indices), the function extracts the corresponding
     * points from the global points vector and computes the face area.
     *
     * @param globalPoints The vector of all points in the mesh.
     * @return std::vector<double> A vector containing the area of each face.
     */
    std::vector<double> computeFaceAreas(const std::vector<Point3D>& globalPoints) const {
        std::vector<double> areas;
        for (const auto &faceIndices : faceVertexIndices) {
            std::vector<Point3D> facePts;
            for (int idx : faceIndices) {
                facePts.push_back(globalPoints[idx]);
            }
            areas.push_back(computePolygonArea3D(facePts));
        }
        return areas;
    }
    
    /**
     * @brief Computes the volume of the cell.
     *
     * For a convex polyhedron, the cell volume is computed by triangulating the cell into
     * tetrahedra formed by the cell's centroid and each triangle of a face (using the first
     * vertex as a common vertex) and summing the volumes of these tetrahedra.
     *
     * @param globalPoints The vector of all points in the mesh.
     * @return double The computed volume of the cell.
     */
    double computeVolume(const std::vector<Point3D>& globalPoints) const {
        Point3D cellCenter = computeCenter(globalPoints);
        double volume = 0.0;
        // For each face of the cell:
        for (const auto &faceIndices : faceVertexIndices) {
            std::vector<Point3D> facePts;
            for (int idx : faceIndices) {
                facePts.push_back(globalPoints[idx]);
            }
            // Skip faces with less than 3 vertices.
            if (facePts.size() < 3)
                continue;
            // Triangulate the face using the first vertex as a common vertex.
            for (size_t i = 1; i < facePts.size() - 1; i++) {
                // Compute the volume of the tetrahedron with vertices:
                // cellCenter, facePts[0], facePts[i], facePts[i+1]
                double Ax = facePts[0].x - cellCenter.x;
                double Ay = facePts[0].y - cellCenter.y;
                double Az = facePts[0].z - cellCenter.z;
                double Bx = facePts[i].x - cellCenter.x;
                double By = facePts[i].y - cellCenter.y;
                double Bz = facePts[i].z - cellCenter.z;
                double Cx = facePts[i+1].x - cellCenter.x;
                double Cy = facePts[i+1].y - cellCenter.y;
                double Cz = facePts[i+1].z - cellCenter.z;
                double crossX = By * Cz - Bz * Cy;
                double crossY = Bz * Cx - Bx * Cz;
                double crossZ = Bx * Cy - By * Cx;
                double dot = Ax * crossX + Ay * crossY + Az * crossZ;
                double tetraVolume = std::abs(dot) / 6.0;
                volume += tetraVolume;
            }
        }
        return volume;
    }
};

/**
 * @brief Represents a 3D subsurface mesh.
 *
 * The SubsurfaceMesh class manages a collection of VolumeCell objects and provides methods for
 * building connectivity and exporting the mesh in VTK format.
 */
class SubsurfaceMesh : public Mesh {
public:
    /**
     * @brief Builds connectivity for the subsurface mesh.
     *
     * This method computes the faces for each VolumeCell, uses an unordered_map to find cells sharing
     * the same face, and assigns neighbor indices accordingly.
     */
    void buildConnectivity() override {
        Timer timer("SubsurfaceMesh::buildConnectivity");
        std::unordered_map<std::vector<int>, std::vector<size_t>, VectorHash> faceMap;
        size_t numCells = cells.size();
        // First pass: for each cell, compute its faces and insert into the map.
        for (size_t i = 0; i < numCells; i++) {
            VolumeCell* cell = dynamic_cast<VolumeCell*>(cells[i].get());
            if (!cell) continue;
            cell->assignFaceVertIndices();
            // Insert each face into the map.
            for (const auto &face : cell->faceVertexIndices) {
                std::vector<int> sortedFace = face;
                std::sort(sortedFace.begin(), sortedFace.end());
                faceMap[sortedFace].push_back(i);
            }
        }
        // Second pass: For each cell, assign neighbor indices for each face.
        for (size_t i = 0; i < numCells; i++) {
            VolumeCell* cell = dynamic_cast<VolumeCell*>(cells[i].get());
            if (!cell) continue;
            // Resize neighborIndices to match the number of faces.
            cell->neighborIndices.assign(cell->faceVertexIndices.size(), -1);
            // For each face in the cell:
            for (size_t f = 0; f < cell->faceVertexIndices.size(); f++) {
                std::vector<int> sortedFace = cell->faceVertexIndices[f];
                std::sort(sortedFace.begin(), sortedFace.end());
                auto it = faceMap.find(sortedFace);
                if (it != faceMap.end()) {
                    const std::vector<size_t>& indices = it->second;
                    if (indices.size() == 1) {
                        // Only one cell has this face: no neighbor.
                        cell->neighborIndices[f] = -1;
                    } else if (indices.size() == 2) {
                        // Two cells share the face: assign the other cell.
                        if (indices[0] == i)
                            cell->neighborIndices[f] = static_cast<int>(indices[1]);
                        else if (indices[1] == i)
                            cell->neighborIndices[f] = static_cast<int>(indices[0]);
                    } else {
                        // More than two cells share the same face; choose the first neighbor that is not the current cell.
                        for (size_t idx : indices) {
                            if (idx != i) {
                                cell->neighborIndices[f] = static_cast<int>(idx);
                                break;
                            }
                        }
                    }
                } else {
                    cell->neighborIndices[f] = -1;
                }
            }
        }
    }
    
    /**
     * @brief Exports the subsurface mesh in VTK file format.
     *
     * The method generates a string in VTK format containing the points, cells, cell types, and cell data
     * (material and hydraulic head) of the subsurface mesh.
     *
     * @return std::string The VTK formatted mesh.
     */
    std::string exportVTK() const {
        std::stringstream ss;
        ss << std::setprecision(4) << std::fixed;
        ss << "# vtk DataFile Version 2.0\n";
        ss << "STYX Subsurface Mesh.\n";
        ss << "ASCII\n";
        ss << "DATASET UNSTRUCTURED_GRID\n";
        ss << "POINTS " << points.size() << " double\n";
        for (const auto &pt : points) {
            ss << pt.x << " " << pt.y << " " << pt.z << "\n";
        }
        // Determine cell type based on one cell's geometry type.
        int num_of_vert = 0;
        if (!cells.empty()) {
            const auto &cell0 = cells.front();
            if (cell0) {
                if (cell0->vertexIndices.size() == 8)
                    num_of_vert = 9;
                else if (cell0->vertexIndices.size() == 6)
                    num_of_vert = 7;
            }
        }
        ss << "CELLS " << cells.size() << " " << num_of_vert * cells.size() << "\n";
        for (const auto &cellPtr : cells) {
            std::vector<int> vertIndices = cellPtr->vertexIndices;
            ss << vertIndices.size();
            for (int idx : vertIndices)
                ss << " " << idx;
            ss << "\n";
        }
        ss << "CELL_TYPES " << cells.size() << "\n";
        for (const auto &cellPtr : cells) {
            if (cellPtr->vertexIndices.size() == 8)
                ss << "12\n";
            else if (cellPtr->vertexIndices.size() == 6)
                ss << "13\n";
            else
                ss << "0\n";
        }
        // Write cell data.
        ss << "CELL_DATA " << cells.size() << "\n";
        ss << "FIELD FieldData 2\n";
        ss << "material 1 " << cells.size() << " int\n";
        for (const auto &cellPtr : cells) {
            ss << cellPtr->material << "\n";
        }
        ss << "hydraulic_head 1 " << cells.size() << " double\n";
        for (const auto &cellPtr : cells) {
            auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cellPtr->waterState);
            ss << (state ? state->hydraulicHead : 0.0) << "\n";
        }
        return ss.str();
    }
};

#endif // SUBSURFACE_MESH_H
