#pragma once

// Gripper includes
#include <Gripper/stl/stlMath.hpp>
#include <Gripper/stl/stlQuadrature.hpp>

// Structured exception handling
#include <Gripper/FloatingExceptions.hpp>

// Standard C++ includes
#include <cstdlib>          // EXIT_SUCCESS
#include <cstdint>          // std::int32_t
#include <iostream>         // std::cout
#include <fstream>          // std::ofstream
#include <chrono>           // std::chrono::high_resolution_clock
#include <cmath>            // std::exp
#include <complex>          // std::complex
#include <algorithm>        // std::transform
#include <vector>           // std::vector
#include <iterator>         // std::back_inserter


namespace stl
{
    template <typename ForwardIt,
              typename F,
              typename Pred>
    auto search_until(ForwardIt first, ForwardIt last, F f, Pred pred)
    {
        if (first == last) throw std::domain_error{ "stl::search_until cannot operate on empty ranges" };

        using result_type = decltype(f(std::declval<typename ForwardIt::value_type>()));

        result_type res;

        for (; first != last; ++first)
            if (pred(res = f(*first)))
                break;

        return res;
    }
}
