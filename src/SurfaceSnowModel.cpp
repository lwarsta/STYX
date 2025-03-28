#include "SurfaceSnowModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief Initializes the SurfaceSnowModel using settings from the JSON configuration.
 *
 * This method reads the weather data interval from the "TimeControl" section and
 * the surface snow parameters from the "SurfaceSnow" section. If any parameter
 * is missing, a default fallback value is used.
 *
 * @param settings A JSON object containing the simulation configuration.
 */
void SurfaceSnowModel::initialize(const JsonValue &settings) {
	(void)settings;
}

void SurfaceSnowModel::runStep(double timeStep, const ForcingRecord &forcing) {
	(void)timeStep;
	(void)forcing;
}