#ifndef TESTSTLARITHMETICS_HPP
#define TESTSTLARITHMETICS_HPP


// Gripper includes
#include <Gripper/stl/stlRuntime.hpp>

// Standard C++ includes
#include <algorithm>
#include <utility>


template<typename STLMultipoleType>
bool differ(const STLMultipoleType& a, const STLMultipoleType& b)
{
    if (a.size() != b.size()) return true;

    if (std::mismatch(a.cbegin(), a.cend(), b.cbegin(), [](const Multipole::stl::ValueType& c, const Multipole::stl::ValueType& d){return std::abs(c - d) > 1e-5; }) != std::make_pair(a.cend(), b.cend())) return false;;

    return false;
}


#endif // TESTSTLARITHMETICS_HPP