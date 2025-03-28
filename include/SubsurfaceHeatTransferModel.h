#ifndef SUBSURFACE_HEAT_TRANSFER_MODEL_H
#define SUBSURFACE_HEAT_TRANSFER_MODEL_H

#include "IModel.h"
#include "SubsurfaceMesh.h"
#include "ForcingData.h"
#include "MeshCoupler.h"
#include <iostream>
#include <typeinfo>

class SubsurfaceHeatTransferModel : public IModel {
public:
    /**
     * @brief Constructs a SubsurfaceHeatTransferModel with the specified subsurface mesh.
     * @param mesh Pointer to a SubsurfaceMesh that represents the 3D surface domain.
     */
    explicit SubsurfaceHeatTransferModel(SubsurfaceMesh* mesh)
        : subsurfaceMesh_(mesh),
          meshCoupler_(nullptr)
    {}

    /**
     * @brief Initializes the model with settings from a JSON object.
     * @param settings A JSON object containing the configuration settings.
     */
    virtual void initialize(const JsonValue &settings) override;

    /**
     * @brief Executes a single simulation step for the subsurface heat transfer model.
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
    SubsurfaceMesh* subsurfaceMesh_;   ///< Pointer to the subsurface mesh.
    MeshCoupler* meshCoupler_;      ///< Pointer to the mesh coupler for heat exchange.

    // Model parameters and variables.
    bool   enabled_;                ///< Flag indicating if the model is enabled.
    double weatherDataInterval_;    ///< Weather data interval (in seconds).
};

#endif // SUBSURFACE_HEAT_TRANSFER_MODEL_H
