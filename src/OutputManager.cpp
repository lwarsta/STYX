#include "OutputManager.h"
#include "SurfaceMesh.h"
#include "SubsurfaceMesh.h"
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <fstream>

/*! \brief Constructs an OutputManager with the specified output file.
    \param filename The path to the output CSV file.
*/
OutputManager::OutputManager(const std::string &filename) {
    ofs_.open(filename);
    if (!ofs_.is_open()) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
}

/*! \brief Destructor that closes the output file stream if open.
*/
OutputManager::~OutputManager() {
    if (ofs_.is_open()) {
        ofs_.close();
    }
}

/*! \brief Writes the CSV header row.
*/
void OutputManager::writeHeader() {
    ofs_ << "Simulation Time (s),Excel Time (d),Precipitation Water Volume (m^3),"
         << "Potential PET Water Volume (m^3),Upper Storage Water Volume (m^3),"
         << "Evaporation Water Volume (m^3),Overland Water Volume (m^3),"
         << "Outfall Water Volume (m^3),Network Water Volume (m^3),"
         << "Infiltration Water Volume (m^3),Transpiration Water Volume (m^3),"
         << "Subsurface Water Volume (m^3)," << std::endl;
}

/*! \brief Writes a row of simulation results to the CSV.
*/
void OutputManager::writeRow(
    int simulationTime,
    double excelTime,
    double precipitationWaterVolume,
    double petWaterVolume,
    double upperStorageWaterVolume,
    double evaporationWaterVolume,
    double overlandWaterVolume, 
    double outfallWaterVolume,
    double networkWaterVolume,
    double infiltrationWaterVolume,
    double transpirationWaterVolume, 
    double subsurfaceWaterVolume
) {
    ofs_ << std::fixed << std::setprecision(4)
         << simulationTime << ","
         << excelTime << ","
         << precipitationWaterVolume << ","
         << petWaterVolume << ","
         << upperStorageWaterVolume << ","
         << evaporationWaterVolume << ","
         << overlandWaterVolume << ","
         << outfallWaterVolume << ","
         << networkWaterVolume << ","
         << infiltrationWaterVolume << ","
         << transpirationWaterVolume << ","
         << subsurfaceWaterVolume << ","
         << std::endl;
}

/*! \brief Closes the output file stream.
*/
void OutputManager::close() {
    if (ofs_.is_open()) {
        ofs_.close();
    }
}

/*! \brief Formats a step number into a three-digit string.
    \param stepNumber The simulation step number.
    \return A string representing the formatted step number.
*/
static std::string formatStepNumber(int stepNumber) {
    std::ostringstream ss;
    ss << std::setw(3) << std::setfill('0') << stepNumber;
    return ss.str();
}

/*! \brief Writes the network mesh as a VTK file.
    \param linkFolder Folder where link VTK files are saved.
    \param juncFolder Folder where junction VTK files are saved.
    \param stepNumber The current simulation step number.
    \param mesh The network mesh.
*/
void OutputManager::writeNetworkVTK(const std::string &linkFolder, const std::string &juncFolder, int stepNumber, const NetworkMesh &mesh) {
    // Currently a stub.
    (void)linkFolder;
    (void)juncFolder;
    (void)stepNumber;
    (void)mesh;
}

/*! \brief Writes the surface mesh as a VTK file.
    \param folder The folder to save the surface VTK file.
    \param stepNumber The current simulation step number.
    \param mesh The surface mesh.
*/
void OutputManager::writeSurfaceVTK(const std::string &folder, int stepNumber, const SurfaceMesh &mesh) {
    std::string vtkContent = mesh.exportVTK();
    std::string fileName = folder + formatStepNumber(stepNumber) + ".vtk";
    std::ofstream ofs(fileName);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file " + fileName + " for writing VTK output.");
    }
    ofs << vtkContent;
    ofs.close();
}

/*! \brief Writes the subsurface mesh as a VTK file.
    \param folder The folder to save the subsurface VTK file.
    \param stepNumber The current simulation step number.
    \param mesh The subsurface mesh.
*/
void OutputManager::writeSubsurfaceVTK(const std::string &folder, int stepNumber, const SubsurfaceMesh &mesh) {
    std::string vtkContent = mesh.exportVTK();
    std::string fileName = folder + formatStepNumber(stepNumber) + ".vtk";
    std::ofstream ofs(fileName);
    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file " + fileName + " for writing VTK output.");
    }
    ofs << vtkContent;
    ofs.close();
}
