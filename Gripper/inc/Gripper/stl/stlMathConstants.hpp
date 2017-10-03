#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>    // _USE_MATH_DEFINES

// Standard C++ includes
#include <cmath>                        // M_E, M_PI, etc.


namespace math
{
    namespace constants
    {
        template <typename T> constexpr T e = static_cast<T>(M_E);
        template <typename T> constexpr T log2e = static_cast<T>(M_LOG2E);
        template <typename T> constexpr T log10e = static_cast<T>(M_LOG10E);
        template <typename T> constexpr T ln2 = static_cast<T>(M_LN2);
        template <typename T> constexpr T ln10 = static_cast<T>(M_LN10);
        template <typename T> constexpr T pi = static_cast<T>(M_PI);
        template <typename T> constexpr T pi_2 = static_cast<T>(M_PI_2);
        template <typename T> constexpr T pi_4 = static_cast<T>(M_PI_4);
        template <typename T> constexpr T pi_inv = static_cast<T>(M_1_PI);
        template <typename T> constexpr T pi_inv2 = static_cast<T>(M_2_PI);
        template <typename T> constexpr T sqrtpi_inv2 = static_cast<T>(M_2_SQRTPI);
        template <typename T> constexpr T sqrt2 = static_cast<T>(M_SQRT2);
        template <typename T> constexpr T sqrt2_inv = static_cast<T>(M_SQRT1_2);
    }
}
