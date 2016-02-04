#ifndef TESTSTLARITHMETICS_HPP
#define TESTSTLARITHMETICS_HPP


// Gripper includes
#include <Gripper/stl/stlMultipoleTypes.hpp>
#include <Gripper/stl/stlRuntime.hpp>

// Standard C++ includes
#include <algorithm>
#include <utility>


template<typename STLMultipoleType_A, typename STLMultipoleType_B>
bool differ(const STLMultipoleType_A& a, const STLMultipoleType_B& b)
{
    ENTERING

    if (a.size() != b.size())
    {
        Gripper::stl::clog << Gripper::LogLevel::DEBUG << "a.size() = " << a.size() << "\tb.size() = " << b.size();

        LEAVING

        return true;
    }

    // STL-compatible version
    //if (std::mismatch(a.cbegin(), a.cend(), b.cbegin(), [](typename STLMultipoleType::reference_type c, typename STLMultipoleType::reference_type d){return std::abs(c - d) < 1e-5; }) != std::make_pair(a.cend(), b.cend())) return true;
    

    // Old-school version
    for (typename STLMultipoleType_A::size_type i = 0; i < a.size(); ++i)
        if (std::abs(a.at(i) - b.at(i)) > 1e-5)
        {
            Gripper::stl::clog << Gripper::LogLevel::DEBUG << "a.at(" << i << ") = " << a.at(i) << "\tb.at(" << i << ") = " << b.at(i);

            LEAVING

            return true;
        }

    LEAVING

    return false;
}


void logTest(bool differ, const char* test_name)
{
    ENTERING

    if (!differ)
        Gripper::stl::clog << Gripper::LogLevel::INFO << test_name << " PASSED";
    else
        Gripper::stl::clog << Gripper::LogLevel::CRIT << test_name << " FAILED";

    LEAVING
}


#endif // TESTSTLARITHMETICS_HPP