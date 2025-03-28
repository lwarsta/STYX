/**
 * @file main.cpp
 * @brief Entry point for the simulation.
 *
 * This file loads simulation settings from a JSON file,
 * creates a Simulation instance, and runs the simulation.
 */

#include "json.h"
#include "InputManager.h"
#include "Simulation.h"
#include <iostream>
#include <exception>
#include <sstream>

/**
 * @brief Main function for the simulation.
 *
 * This function loads the simulation configuration, creates a Simulation object,
 * and runs the simulation. If any exception is thrown during initialization or
 * simulation, an error message is printed and the program exits with failure.
 *
 * @return int Exit status: 0 for success, non-zero for failure.
 */
int main() {
    try {
        // Load the settings file into a string using InputManager.
        std::string settingsStr = InputManager::loadTextFile("config/settings.json");
		
        // Parse the settings using the new JsonValue parser.
        JsonValue settings = JsonValue::parse(settingsStr);

        // Create and run the simulation.
        Simulation simulation(settings);
        simulation.run();
    } catch (const std::exception &ex) {
        std::cerr << "Fatal error during simulation: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
