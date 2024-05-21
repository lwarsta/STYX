#include "CellWater2d.h"

CellWater2d::CellWater2d()
{
    geom = 0;
    waterDepth = 0.0;
    waterDepthOld = 0.0;
    up_stor_depth = 0.0;
    deprStor = 0.0;
    mannN = 0.0;
    up_stor_thresh = 0.0;
    evap_frac = 0.0;
    crop_factor = 0.0;
    root_depth = 0.0;
    outletIndex = -1;
    distToOutlet = 0.0;
    avgDisch = 0.0;
    sinkIndex = -1;
    watVolToSinks = 0.0;
}

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