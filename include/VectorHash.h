/**
 * @file VectorHash.h
 * @brief Provides a hash function for std::vector<int> objects.
 *
 * This header defines the VectorHash structure, which computes a hash value
 * for a vector of integers. This can be used in unordered containers that require
 * a hash function for keys of type std::vector<int>.
 */

#ifndef VECTOR_HASH_H
#define VECTOR_HASH_H

#include <vector>
#include <functional>

/**
 * @brief Functor that computes a hash value for a std::vector<int>.
 */
struct VectorHash {
    /**
     * @brief Computes the hash value for a vector of integers.
     *
     * This operator iterates over the elements of the vector and combines their hash values
     * using bitwise operations.
     *
     * @param v The vector of integers to be hashed.
     * @return The computed hash value.
     */
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (int i : v) {
            seed ^= std::hash<int>()(i) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

#endif // VECTOR_HASH_H
