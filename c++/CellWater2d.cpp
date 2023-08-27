#include "CellWater2d.h"

CellWater2d::CellWater2d()
{
    geom = 0;
    waterDepth = 0.0;
    waterDepthOld = 0.0;
    deprStor = 0.0;
    mannN = 0.0;
    outletIndex = -1;
    distToOutlet = 0.0;
    avgDisch = 0.0;
    sinkIndex = -1;
    watVolToSinks = 0.0;
}
/*
double CellWater2d::compFlowToManhole(double timeStep)
{
    double waterVolume = 0.0;

    if (geom != 0 && geom->getArea() > 0.0 && outletIndex != -1)
    {
        // The same radius is used for all manholes.
        double manholeRadius = 0.25; // could be 0.25 m
        // Compute slope.
        double slope = 0.0;
        Vertex cp = geom->getCentrePoint();

        for (size_t i = 0; i < neighPointers.size(); i++)
        {
            if (neighPointers.at(i) != 0)
            {
                CellGeom2d * neighGeom = neighPointers.at(i)->getGeom();
                Vertex neighCp = neighGeom->getCentrePoint();
                double elevation = cp.z - neighCp.z;
                double distance = sqrt((cp.x - neighCp.x) * 
                    (cp.x - neighCp.x) + (cp.y - neighCp.y) * (cp.y - neighCp.y));
                slope += (elevation / distance) * (elevation / distance);
            }
        }

        slope = sqrt(slope);

        // Compute water depth above depression storage.
        double waterDepthEff = waterDepth - deprStor;

        if (waterDepthEff < 0.0)
        {
            waterDepthEff = 0.0;
        }

        // Compute water volume entering the manhole during the time step.
        waterVolume = 2.0 * 3.14159265358979323846 * manholeRadius * 
            waterDepthEff * pow(waterDepthEff, 2.0/3.0) / mannN * sqrt(slope) * 
            timeStep;

        if (waterVolume > waterDepthEff * geom->getArea())
        {
            waterVolume = waterDepthEff * geom->getArea();
        }

        waterDepth -= waterVolume / geom->getArea();

        if (waterDepth < 0.0)
        {
            waterDepth = 0.0;
        }
    }
    return waterVolume;
}
*/
/*
void CellWater2d::addWatVolToStorm(double watVolStorm, double distToOutlet, 
    double simTime)
{
    // v = 1.273 q / D^2
    // q = 0.05 m^3 / s and D = 0.3 m -> 0.707 m / s
    //  It might be possible to estimate how much water is moving to the outlet to estimate discharge?
    StormWatItem stormWatItem;
    stormWatItem.watVol = watVolStorm;
    stormWatItem.dischTime = simTime + distToOutlet / 0.71;
    stormWatItem.distToOutlet = watVolStorm;
    //stormWatItems.push_back(stormWatItem);
    stormWatItems.insert(stormWatItems.begin(), stormWatItem);
}
*/
/*
void CellWater2d::removeWatVolFromStorm(double simTime)
{
    if (geom != 0 && geom->getArea() > 0.0)
    {
        while (!stormWatItems.empty())
        {
            StormWatItem stormWatItem = stormWatItems.back();

            if (stormWatItem.dischTime <= simTime)
            {
                waterDepth += stormWatItem.watVol / geom->getArea();
                stormWatItems.pop_back();
            }
            else {
                break;
            }
        }
    }
}
*/
void CellWater2d::removeWatBySinks(double timeStep)
{
    if (sinkIndex >= 0 && geom != 0 && geom->getArea() > 0.0)
    {
        /*
        double waterDepthEff = waterDepth - deprStor;
        if (waterDepthEff < 0.0)
        {
            waterDepthEff = 0.0;
        }
        double watDepthRemoved = timeStep * avgDisch / geom->getArea();
        if (watDepthRemoved > waterDepthEff)
        {
            watDepthRemoved = waterDepthEff;
        }
        waterDepth -= watDepthRemoved;
        if (waterDepth < 0.0)
        {
            waterDepth = 0.0;
        }
        watVolToSinks += watDepthRemoved * geom->getArea();
        */
        
        double waterDepthEff = waterDepth - deprStor;

        if (waterDepthEff < 0.0)
        {
            waterDepthEff = 0.0;
        }

        waterDepth -= waterDepthEff;

        if (waterDepth < 0.0)
        {
            waterDepth = 0.0;
        }

        watVolToSinks = waterDepthEff * geom->getArea();
    }
}