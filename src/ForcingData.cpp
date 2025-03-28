#include "ForcingData.h"
#include "DateTimeUtils.h"  // Assumed to provide parseDateTime()
#include <algorithm>
#include <sstream>
#include <iomanip>

/*!
 * \brief Loads forcing data from a CSV.
 *
 * See the documentation in the header for details.
 */
void ForcingData::load(const std::vector<std::vector<std::string>> &csvData, std::time_t startTime) {
    records_.clear();
    if (csvData.size() < 2) {
        throw std::runtime_error("Forcing CSV is empty or has no data rows.");
    }

    // Skip header row (index 0)
    for (size_t i = 1; i < csvData.size(); i++) {
        const auto &row = csvData[i];
        if (row.size() < 4) {
            // Skip incomplete rows.
            continue;
        }
        // Parse the datetime in column 0.
        std::time_t rowTime = parseDateTime(row[0]);
        int offset = static_cast<int>(std::difftime(rowTime, startTime));

        // Parse precipitation, PET, and air temperature.
        double prec = std::stod(row[1]);
        double pet  = std::stod(row[2]);
        double temp = std::stod(row[3]);

        ForcingRecord rec;
        rec.offsetSeconds = offset;
        rec.precipitation = prec;
        rec.pet = pet;
        rec.airTemperature = temp;
        records_.push_back(rec);
    }
    if (records_.empty()) {
        throw std::runtime_error("No valid forcing rows after parsing CSV data.");
    }

    // Sort records chronologically if not already sorted.
    std::sort(records_.begin(), records_.end(),
              [](const ForcingRecord &a, const ForcingRecord &b) {
                  return a.offsetSeconds < b.offsetSeconds;
              });
}

/*!
 * \brief Returns the forcing record for the specified simulation time.
 *
 * See the documentation in the header for details.
 */
ForcingRecord ForcingData::getForcingAtTime(int simTimeSec) const {
    if (records_.empty()) {
        // Return a default record if empty.
        return {0, 0.0, 0.0, 0.0};
    }

    // Find the last record with offset <= simTimeSec.
    ForcingRecord result = records_.front();
    for (const auto &rec : records_) {
        if (rec.offsetSeconds <= simTimeSec) {
            result = rec;
        } else {
            break;
        }
    }
    return result;
}
