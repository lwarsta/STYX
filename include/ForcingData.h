#ifndef FORCING_DATA_H
#define FORCING_DATA_H

#include <vector>
#include <string>
#include <ctime>
#include <stdexcept>

/*!
 * \brief A simple record containing forcing data for a specific time offset.
 */
struct ForcingRecord {
    int offsetSeconds;      ///< Offset from simulation start time (in seconds)
    double precipitation;   ///< Precipitation in mm
    double pet;             ///< Potential evapotranspiration in mm
    double airTemperature;  ///< Air temperature in °C
};

/*!
 * \brief Manages time-based forcing data.
 *
 * This class loads and provides access to environmental forcing data 
 * (precipitation, PET, air temperature) from a CSV file.
 */
class ForcingData {
public:
    /*!
     * \brief Loads forcing data from a CSV.
     *
     * Each CSV row (except the header) should contain:
     *   - A datetime string,
     *   - Precipitation [mm],
     *   - PET [mm],
     *   - Air temperature [°C].
     *
     * The datetime is converted to an offset (in seconds) from the given start time.
     *
     * \param csvData The CSV data as a vector of string vectors.
     * \param startTime The simulation start time as a std::time_t.
     * \throw std::runtime_error if there is insufficient or invalid data.
     */
    void load(const std::vector<std::vector<std::string>> &csvData, std::time_t startTime);

    /*!
     * \brief Returns the forcing record for the given simulation time.
     *
     * The function finds the last record with an offset less than or equal to the given simulation time.
     * If simTimeSec exceeds the last record's offset, the last record is returned.
     *
     * \param simTimeSec Simulation time in seconds from start.
     * \return The corresponding ForcingRecord.
     */
    ForcingRecord getForcingAtTime(int simTimeSec) const;

private:
    std::vector<ForcingRecord> records_;  ///< List of forcing records.
};

#endif // FORCING_DATA_H
