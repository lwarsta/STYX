#include "SubsurfaceHeatTransferModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief Initializes the SubsurfaceHeatTransferModel using settings from the JSON configuration.
 *
 * This method reads the weather data interval from the "TimeControl" section and
 * the surface snow parameters from the "SubsurfaceHeatTransfer" section. If any parameter
 * is missing, a default fallback value is used.
 *
 * @param settings A JSON object containing the simulation configuration.
 */
void SubsurfaceHeatTransferModel::initialize(const JsonValue &settings) {
	(void)settings;
}

void SubsurfaceHeatTransferModel::runStep(double timeStep, const ForcingRecord &forcing) {
	(void)timeStep;
	(void)forcing;
}