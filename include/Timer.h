/**
 * @file Timer.h
 * @brief Simple RAII timer for measuring elapsed time.
 *
 * This header defines the Timer class, which starts timing upon construction and
 * prints the elapsed time (in milliseconds) upon destruction.
 */

#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>
#include <iostream>

/**
 * @brief RAII timer class that prints elapsed time on destruction.
 *
 * The Timer class records the start time when constructed and computes the elapsed time
 * when it is destroyed. This is useful for measuring the duration of code blocks.
 */
class Timer {
public:
    /**
     * @brief Constructs a Timer with a given name.
     * @param name The name of the timer (used for output).
     */
    explicit Timer(const std::string &name)
        : name_(name),
          start_(std::chrono::steady_clock::now())
    {}

    /**
     * @brief Destructor that prints the elapsed time.
     */
    ~Timer() {
        auto end = std::chrono::steady_clock::now();
        auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
        std::cout << "Timer [" << name_ << "] elapsed " << durationMs << " ms." << std::endl;
    }

private:
    std::string name_;                             ///< Name of the timer.
    std::chrono::steady_clock::time_point start_;  ///< Start time of the timer.
};

#endif // TIMER_H
