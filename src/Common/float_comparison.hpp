/**
 * @file float_comparison.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-08-31
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>

namespace uepm {

namespace utils {

/**
 * @brief Comparison of double numbers.
 *
 * @param A
 * @param B
 * @param maxUlps
 * @return true
 * @return false
 */
static inline bool doubles_are_equal(float A, float B, int maxUlps = 4) {
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    const int large_int = 0x80000000;
    assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    int aInt;
    std::memcpy(&aInt, &A, sizeof(int));
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0) {
        aInt = large_int - aInt;
    }
    // Make bInt lexicographically ordered as a twos-complement int
    int bInt;
    std::memcpy(&bInt, &B, sizeof(int));
    if (bInt < 0) {
        bInt = large_int - bInt;
    }
    int intDiff = abs(aInt - bInt);
    return (intDiff <= maxUlps);
}

}  //  namespace utils

} // namespace uepm