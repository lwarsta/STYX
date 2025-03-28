/**
 * @file NetworkMesh.h
 * @brief Defines the network mesh classes for representing junctions and links.
 *
 * This file contains the definitions for the NetworkJunctionCell, NetworkLinkCell, and
 * NetworkMesh classes. The NetworkMesh class extends the Mesh base class and provides
 * helper functions for computing geometric properties of network elements.
 */

#ifndef NETWORK_MESH_H
#define NETWORK_MESH_H

#include "Mesh.h"
#include "PhysicsState.h"
#include "Common.h"
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

/**
 * @brief Represents a network junction cell (e.g., a well) in the network mesh.
 */
class NetworkJunctionCell : public Cell {
public:
    /// The top point of the junction.
    Point3D top;
    /// The bottom point of the junction.
    Point3D bottom;
    /// Connected network link cells.
    std::vector<class NetworkLinkCell*> connectedLinks;

    /// Junction diameter in meters.
    double diameter = 0.0;
    /// Junction type (e.g., regular = 0, outfall = 1).
    int type = -1;

    /// Initial water depth for the junction.
    double initialWaterDepth = 0.0;

    /**
     * @brief Adds a network link cell to the connected links.
     *
     * @param link Pointer to the network link cell to add.
     */
    void addLink(class NetworkLinkCell* link) { 
        connectedLinks.push_back(link); 
    }
};

/**
 * @brief Represents a network link cell (e.g., a pipe) in the network mesh.
 */
class NetworkLinkCell : public Cell {
public:
    /// The first endpoint of the link.
    Point3D end1;
    /// The second endpoint of the link.
    Point3D end2;
    /// Pointer to the first connected junction cell.
    NetworkJunctionCell* junction1 = nullptr;
    /// Pointer to the second connected junction cell.
    NetworkJunctionCell* junction2 = nullptr;

    /// Link diameter in meters.
    double diameter = 0.0;
    /// Link roughness (Manning's n).
    double roughness = 0.0;
};

/**
 * @brief Represents the network mesh containing junction and link cells.
 *
 * The NetworkMesh class extends the Mesh base class and stores pointers to network
 * junction cells and link cells. It provides functions to build connectivity between
 * network elements and helper functions to compute geometric properties.
 */
class NetworkMesh : public Mesh {
public:
    /// Vector of pointers to network junction cells.
    std::vector<NetworkJunctionCell*> junctions;
    /// Vector of pointers to network link cells.
    std::vector<NetworkLinkCell*> links;

    /**
     * @brief Builds connectivity by assigning each link's endpoints to the nearest junction.
     *
     * The method separates junctions and links from the base Mesh::cells vector and then
     * matches each link to the closest junctions by comparing the link endpoints with the
     * junctions' top and bottom points.
     */
    void buildConnectivity() override {
        junctions.clear();
        links.clear();

        // Separate junctions and links from the base cells vector.
        for (auto &cellPtr : cells) {
            auto junc = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
            if (junc) {
                junctions.push_back(junc);
                continue;
            }
            auto link = dynamic_cast<NetworkLinkCell*>(cellPtr.get());
            if (link) {
                links.push_back(link);
            }
        }

        // For each link, find the closest junctions for each endpoint.
        for (auto &link : links) {
            double bestDist1 = std::numeric_limits<double>::max();
            double bestDist2 = std::numeric_limits<double>::max();
            NetworkJunctionCell* bestJunc1 = nullptr;
            NetworkJunctionCell* bestJunc2 = nullptr;

            // For each junction, compute the distance from both endpoints.
            for (auto &jnc : junctions) {
                // Compute distances from link->end1 to junction top and bottom.
                double d1_top = lateralDistance(link->end1, jnc->top);
                double d1_bottom = lateralDistance(link->end1, jnc->bottom);
                double d1 = std::min(d1_top, d1_bottom);

                // Compute distances from link->end2 to junction top and bottom.
                double d2_top = lateralDistance(link->end2, jnc->top);
                double d2_bottom = lateralDistance(link->end2, jnc->bottom);
                double d2 = std::min(d2_top, d2_bottom);

                if (d1 < bestDist1) {
                    bestDist1 = d1;
                    bestJunc1 = jnc;
                }
                if (d2 < bestDist2) {
                    bestDist2 = d2;
                    bestJunc2 = jnc;
                }
            }

            if (!bestJunc1 || !bestJunc2) {
                throw std::runtime_error("Could not connect a pipe to a junction (distance threshold exceeded?)");
            }

            // Distance threshold check to warn about large distances.
            const double DISTANCE_THRESHOLD = 1.0; // Magic number; adjust as needed.
            if (bestDist1 > DISTANCE_THRESHOLD || bestDist2 > DISTANCE_THRESHOLD) {
                std::cout << "Warning: Link " << link->id << " connected to junctions with large distances: "
                          << bestDist1 << ", " << bestDist2 << std::endl;
            }

            link->junction1 = bestJunc1;
            link->junction2 = bestJunc2;
            bestJunc1->addLink(link);
            bestJunc2->addLink(link);
        }
    }

    /**
     * @brief Computes the effective area of a junction.
     *
     * The effective area is approximated as the area of a circle with the given diameter.
     *
     * @param diameter The diameter of the junction in meters.
     * @return double The effective area.
     */
    inline double getJunctionArea(double diameter) const {
        double radius = 0.5 * diameter;
        return M_PI * radius * radius;
    }

    /**
     * @brief Computes the depth of a junction.
     *
     * The depth is calculated as the difference between the top and bottom z-coordinates.
     *
     * @param junction Pointer to a junction cell.
     * @return double The junction depth.
     */
    inline double getJunctionDepth(const NetworkJunctionCell* junction) const {
        return junction->top.z - junction->bottom.z;
    }

    /**
     * @brief Computes the slope of a network link.
     *
     * The slope is defined as the absolute vertical difference between the link's endpoints
     * divided by the horizontal distance.
     *
     * @param link Pointer to the network link cell.
     * @return double The slope of the link.
     */
    inline double getLinkSlope(const NetworkLinkCell* link) const {
        double dx = link->end2.x - link->end1.x;
        double dy = link->end2.y - link->end1.y;
        double horizontalDistance = std::sqrt(dx * dx + dy * dy);
        if (horizontalDistance == 0.0)
            return 0.0;
        return std::abs(link->end2.z - link->end1.z) / horizontalDistance;
    }

    /**
     * @brief Computes the lateral (horizontal) distance between the endpoints of a link.
     *
     * @param link Pointer to the network link cell.
     * @return double The lateral distance.
     */
    inline double getLinkLateralDistance(const NetworkLinkCell* link) const {
        return lateralDistance(link->end1, link->end2);
    }

    /**
     * @brief Computes the three-dimensional length of a link.
     *
     * @param link Pointer to the network link cell.
     * @return double The 3D length.
     */
    inline double getLinkLength(const NetworkLinkCell* link) const {
        double dx = link->end2.x - link->end1.x;
        double dy = link->end2.y - link->end1.y;
        double dz = link->end2.z - link->end1.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    /**
     * @brief Computes the cross-sectional area of a link.
     *
     * The area is approximated as the area of a circle with the link's diameter.
     *
     * @param link Pointer to the network link cell.
     * @return double The cross-sectional area.
     */
    inline double getLinkCrossectionArea(const NetworkLinkCell* link) const {
        double radius = link->diameter / 2.0;
        return M_PI * radius * radius;
    }

    /**
     * @brief Computes the centre point of a link.
     *
     * The centre is computed as the midpoint between the two endpoints.
     *
     * @param link Pointer to the network link cell.
     * @return Point3D The centre point.
     */
    inline Point3D getLinkCentrePoint(const NetworkLinkCell* link) const {
        Point3D cp;
        cp.x = (link->end1.x + link->end2.x) / 2.0;
        cp.y = (link->end1.y + link->end2.y) / 2.0;
        cp.z = (link->end1.z + link->end2.z) / 2.0;
        return cp;
    }

private:
    /**
     * @brief Computes the lateral Euclidean distance (ignoring z) between two points.
     *
     * @param p1 The first point.
     * @param p2 The second point.
     * @return double The lateral distance.
     */
    double lateralDistance(const Point3D& p1, const Point3D& p2) const {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        return std::sqrt(dx * dx + dy * dy);
    }
};

#endif // NETWORK_MESH_H
