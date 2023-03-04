#include <iostream>
#include <chrono>
#include "Framework.h"

void printBanner();

int main(int argc , char *argv[])
{
    // Start timing.
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    // Print banner and command line arguments for inspection.
    printBanner();
    std::cout << "Command line arguments:\n";
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Argument " << i << ": " << argv[i] << "\n";
    }

    // Initialise and run the simulation.
    Framework framework;
    int result;
    result = framework.initialize(argv[1]);
    if (result != 0)
    {
        std::cout << "-> Simulation initialization failed\n";
        return result;
    }
    result = framework.run();
    if (result != 0)
    {
        std::cout << "-> Simulation running failed\n";
        return result;
    }
    result = framework.finalize();
    if (result != 0)
    {
        std::cout << "-> Simulation finalization failed\n";
        return result;
    }

    // Print the execution time.
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << "Execution time: " << ms_int.count() << " ms\n\n";
    std::cout << "----------------------------------------------------------"
        << "\n";
    return 0;
}

void printBanner()
{
    std::cout << "---------------------------------------------------------\n";
    std::cout << " _____ _______   ____   __                      _      _ \n";
    std::cout << "/  ___|_   _\\ \\ / /\\ \\ / /                     | |    | |\n";
    std::cout << "\\ `--.  | |  \\ V /  \\ V /   _ __ ___   ___   __| | ___| |\n";
    std::cout << " `--. \\ | |   \\ /   /   \\  | '_ ` _ \\ / _ \\ / _` |/ _ \\ |\n";
    std::cout << "/\\__/ / | |   | |  / /^\\ \\ | | | | | | (_) | (_| |  __/ |\n";
    std::cout << "\\____/  \\_/   \\_/  \\/   \\/ |_| |_| |_|\\___/ \\__,_|\\___|_|\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "                   - STYX model -                        \n";
    std::cout << "         2d/3d hydrological modelling framework.         \n";
    std::cout << "                 Version 0.1 (2023).                     \n";
    std::cout << "            Developed under the MIT license.             \n";
    std::cout << "           Administrator: lassi.warsta@vtt.fi.           \n";
    std::cout << "---------------------------------------------------------\n";
}
