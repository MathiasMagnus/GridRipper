// Standard C++ includes
#include <cstdlib>          // EXIT_SUCCESS
#include <cstdint>          // std::int32_t
#include <iostream>         // std::cout
#include <chrono>           // std::chrono::high_resolution_clock
#include <cmath>            // std::exp

// Gripper includes
#include <Gripper/stl/stlMultipoleTypes.hpp>
#include <Gripper/PDE.hpp>

namespace std
{
    template <typename E>
    auto exp(const E& v)
    {
        return Multipole::stl::SpinWeightedSpherical::map(v, [](auto&& val) { return std::exp(val); });
    }
}