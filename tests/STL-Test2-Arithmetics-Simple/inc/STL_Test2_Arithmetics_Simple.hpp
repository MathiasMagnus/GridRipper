#pragma once

// Gripper includes
#include <Gripper/stl/stlMultipoleTypes.hpp>

// Standard C++ includes
#include <cassert>          // assert
#include <cstdlib>          // EXIT_SUCCESS, EXIT_FAILURE
#include <type_traits>      // std::common_type_t
#include <cstdint>          // std::int32_t
#include <complex>          // std::complex
#include <iostream>         // std::cout
#include <chrono>           // std::chrono::high_resolution_clock
#include <cmath>            // std::abs


template < typename E1, typename E2 >
auto sum_diff(const E1& lhs, const E2& rhs)
{
    std::common_type_t<typename E1::value_type, typename E2::value_type> difference = 0;

    assert(lhs.extent() == rhs.extent());

    for (auto i = lhs.extent().initial(); lhs.extent().contains(i); ++i)
    {
        difference += std::abs(lhs.at(i) - rhs.at(i));
    }

    return difference;
}

template < typename E1, typename E2, typename T >
bool test(const E1& lhs, const E2& rhs, const T tolerance)
{
    if (lhs.extent() == rhs.extent())
        return sum_diff(lhs, rhs) > tolerance ? false : true;
    else
        return false;
}
