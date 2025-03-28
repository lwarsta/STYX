#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

/*! 
 * \brief Forward declarations for mesh classes.
 */
class NetworkMesh;
class SurfaceMesh;
class SubsurfaceMesh;

/*! 
 * \brief Manages simulation output.
 *
 * The OutputManager class handles writing simulation results to CSV files as well as 
 * exporting VTK files for network, surface, and subsurface meshes.
 */
class OutputManager {
public:
    /*!
     * \brief Constructs an OutputManager with the given output file path.
     * \param filename The path to the output CSV file.
     */
    explicit OutputManager(const std::string &filename);

    /*!
     * \brief Destructor. Closes the output file if open.
     */
    ~OutputManager();

    /*!
     * \brief Writes the header row to the output CSV.
     */
    void writeHeader();

    /*!
     * \brief Writes one row of simulation results to the output CSV.
     * \param simulationTime Simulation time in seconds.
     * \param excelTime Excel serial date.
     * \param precipitationWaterVolume Precipitation water volume (m^3).
     * \param petWaterVolume Potential PET water volume (m^3).
     * \param upperStorageWaterVolume Upper storage water volume (m^3).
     * \param evaporationWaterVolume Evaporation water volume (m^3).
     * \param overlandWaterVolume Overland water volume (m^3).
     * \param outfallWaterVolume Outfall water volume (m^3).
     * \param networkWaterVolume Network water volume (m^3).
     * \param infiltrationWaterVolume Infiltration water volume (m^3).
     * \param transpirationWaterVolume Transpiration water volume (m^3).
     * \param subsurfaceWaterVolume Subsurface water volume (m^3).
     */
    void writeRow(
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
    );

    /*!
     * \brief Closes the output file.
     */
    void close();

    /*!
     * \brief Writes the network mesh in VTK format.
     * \param linkFolder Folder for link VTK files.
     * \param juncFolder Folder for junction VTK files.
     * \param stepNumber Simulation step number.
     * \param mesh The network mesh.
     */
    void writeNetworkVTK(const std::string &linkFolder, const std::string &juncFolder, int stepNumber, const NetworkMesh &mesh);

    /*!
     * \brief Writes the surface mesh in VTK format.
     * \param folder Folder for surface VTK files.
     * \param stepNumber Simulation step number.
     * \param mesh The surface mesh.
     */
    void writeSurfaceVTK(const std::string &folder, int stepNumber, const SurfaceMesh &mesh);

    /*!
     * \brief Writes the subsurface mesh in VTK format.
     * \param folder Folder for subsurface VTK files.
     * \param stepNumber Simulation step number.
     * \param mesh The subsurface mesh.
     */
    void writeSubsurfaceVTK(const std::string &folder, int stepNumber, const SubsurfaceMesh &mesh);

private:
    std::ofstream ofs_;  ///< Output file stream.
};

#endif // OUTPUT_MANAGER_H
