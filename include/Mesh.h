/**
 * @file Mesh.h
 * @brief Contains the base classes for mesh representation.
 *
 * This file defines the Cell and Mesh classes used to represent a mesh. A mesh consists of
 * a collection of 3D points and cells, where each cell stores indices into the points vector,
 * neighbor cell indices, and pointers to physics state objects.
 */

#ifndef MESH_H
#define MESH_H

#include <vector>
#include <memory>
#include <iostream>
#include "Common.h"

/**
 * @brief Forward declaration for physics state objects.
 */
class PhysicsState;

/**
 * @brief Base class representing a cell in a mesh.
 *
 * A Cell stores vertex indices (into the Mesh::points vector), neighbor indices, and pointers to
 * associated physics state objects (for water, heat, and solute). It also stores several cell properties.
 */
class Cell {
public:
    /// Virtual destructor.
    virtual ~Cell() = default;

    /**
     * @brief Identify the cell.
     *
     * This virtual function prints a basic message indicating that this is a base cell.
     * Derived classes may override this function to provide more specific identification.
     */
    virtual void identify() const {
        std::cout << "Base Cell" << std::endl;
    }

    /// Indices into the Mesh::points array representing the cell's vertices.
    std::vector<int> vertexIndices;

    /// Neighboring cell indices.
    std::vector<int> neighborIndices;

    /// Pointer to the physics state for water in the cell.
    std::shared_ptr<PhysicsState> waterState;
    
    /// Pointer to the physics state for heat in the cell.
    std::shared_ptr<PhysicsState> heatState;
    
    /// Pointer to the physics state for solute in the cell.
    std::shared_ptr<PhysicsState> soluteState;

    /// Cell type identifier.
    int cellType = -1;
    
    /// Unique cell identifier.
    int id = -1;
    
    /// Grid connection identifier.
    int gridConnection = -1;
    
    /// Material identifier.
    int material = -1;
    
    /// Initial condition identifier.
    int initCond = -1;
    
    /// Boundary condition identifier.
    int boundCond = -1;
};

/**
 * @brief Base class representing a mesh.
 *
 * A Mesh is composed of a set of 3D points and a set of cells. Derived classes must implement
 * the buildConnectivity() function to establish connectivity (e.g., neighbor indices) between cells.
 */
class Mesh {
public:
    /// Virtual destructor.
    virtual ~Mesh() = default;

    /// A vector of 3D points representing the geometry of the mesh.
    std::vector<Point3D> points;

    /// A vector of cells that constitute the mesh.
    std::vector<std::unique_ptr<Cell>> cells;

    /**
     * @brief Build connectivity between cells.
     *
     * This pure virtual function must be implemented by derived classes to establish connectivity
     * information (such as neighbor indices) among the cells.
     */
    virtual void buildConnectivity() = 0;
};

#endif // MESH_H
