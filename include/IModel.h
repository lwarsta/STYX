#ifndef IMODEL_H
#define IMODEL_H

#include "json.h"
#include "ForcingData.h"

/**
 * @brief Abstract base class for simulation models.
 *
 * This interface defines the minimum functions that any simulation model must
 * implement in order to be used in the simulation framework.
 */
class IModel {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IModel() = default;

    /**
     * @brief Initialize the model with configuration settings.
     *
     * Derived classes should implement this function to initialize model-specific
     * parameters and any required pointers or references (e.g., to mesh data).
     *
     * @param settings A JSON object containing the configuration settings.
     */
    virtual void initialize(const JsonValue &settings) = 0;

    /**
     * @brief Execute a single simulation time step.
     *
     * This function performs the necessary computations for a single simulation
     * step using the provided time step duration and environmental forcing data.
     *
     * @param timeStep The simulation time step in seconds.
     * @param forcing The forcing record containing environmental data.
     */
    virtual void runStep(double timeStep, const ForcingRecord &forcing) = 0;
};

#endif // IMODEL_H
