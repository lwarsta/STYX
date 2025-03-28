#ifndef SURFACE_SNOW_MODEL_H
#define SURFACE_SNOW_MODEL_H

#include "IModel.h"
#include "SurfaceMesh.h"
#include "ForcingData.h"
#include "MeshCoupler.h"
#include <iostream>
#include <typeinfo>

class SurfaceSnowModel : public IModel {
public:
    /**
     * @brief Constructs a SurfaceSnowModel with the specified surface mesh.
     * @param mesh Pointer to a SurfaceMesh that represents the 2D surface domain.
     */
    explicit SurfaceSnowModel(SurfaceMesh* mesh)
        : surfaceMesh_(mesh),
          meshCoupler_(nullptr)
    {}

    /**
     * @brief Initializes the model with settings from a JSON object.
     * @param settings A JSON object containing the configuration settings.
     */
    virtual void initialize(const JsonValue &settings) override;

    /**
     * @brief Executes a single simulation step for the surface flow model.
     * @param timeStep The time step duration in seconds.
     * @param forcing A ForcingRecord containing the environmental forcing data.
     */
    virtual void runStep(double timeStep, const ForcingRecord &forcing) override;

    /**
     * @brief Sets the MeshCoupler to be used for heat exchange.
     * @param coupler Pointer to a MeshCoupler instance.
     */
    void setMeshCoupler(MeshCoupler* coupler) { meshCoupler_ = coupler; }

private:
    SurfaceMesh* surfaceMesh_;      ///< Pointer to the surface mesh.
    MeshCoupler* meshCoupler_;      ///< Pointer to the mesh coupler for heat exchange.

    // Common model parameters and variables.
    bool   enabled_;                ///< Flag indicating if the model is enabled.
    double weatherDataInterval_;    ///< Weather data interval (in seconds).
};

#endif // SURFACE_SNOW_MODEL_H
