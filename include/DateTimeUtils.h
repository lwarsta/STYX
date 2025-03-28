/**
 * @file DateTimeUtils.h
 * @brief Utility functions for parsing and converting date/time values.
 *
 * This file contains inline functions for parsing a date/time string into a
 * std::time_t and for converting a std::time_t into an Excel serial date.
 */

#ifndef DATETIME_UTILS_H
#define DATETIME_UTILS_H

#include <string>
#include <ctime>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cmath>

/**
 * @brief Parses a date/time string into a std::time_t.
 *
 * The function expects the input string to be in the format "%Y-%m-%dT%H:%M:%S"
 * (e.g. "2023-02-14T15:30:00"). If parsing fails, a std::runtime_error is thrown.
 *
 * @param datetimeStr The date/time string to parse.
 * @return A std::time_t representing the parsed date and time.
 * @throws std::runtime_error if the date/time string cannot be parsed.
 */
inline std::time_t parseDateTime(const std::string &datetimeStr) {
    std::tm tm = {};
    std::istringstream ss(datetimeStr);
    ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S");
    if (ss.fail()) {
        throw std::runtime_error("Failed to parse datetime: " + datetimeStr);
    }
    tm.tm_isdst = -1;  // Let mktime decide DST
    return std::mktime(&tm);
}

/**
 * @brief Converts a std::time_t to an Excel serial date.
 *
 * The Excel serial date is defined as the number of days since 1900-01-00.
 * A time zone correction is applied by adding 0.125 days.
 *
 * @param t The time value as std::time_t.
 * @return The Excel serial date as a double.
 */
inline double timeToExcelDate(std::time_t t) {
    static const double SECONDS_PER_DAY = 86400.0;
    static const double EXCEL_EPOCH_OFFSET = 25569.0 + 0.125; // days between 1900-01-00 and 1970-01-01 + correction
    return EXCEL_EPOCH_OFFSET + (static_cast<double>(t) / SECONDS_PER_DAY);
}

#endif // DATETIME_UTILS_H
