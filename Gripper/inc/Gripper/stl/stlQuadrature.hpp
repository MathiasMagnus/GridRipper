#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>
#include <Gripper/stl/stlMath.hpp>
#include <Gripper/stl/stlArithmeticProgression.hpp>

// Standard C++ includes
#include <cstdint>          // std::int_fast32_t
#include <numeric>          // std::accumulate

namespace math
{
    /// <summary>Numerical integrator from the Newton-Cotes family of formulas. Converges rapidly, when used on periodic functions over their periods.</summary>
    /// <param name="from">Lower bound of the integration interval.</param>
    /// <param name="to">Upper bound of the integration interval.</param>
    /// <param name="n">The number of intervals to approximate the result with. (Minimally 2)</param>
    /// <returns>The value of the numerical integral.</returns>
    /// <remarks>Uses <c>n + 1</c> function evaluations.</remarks>
    /// <cite>https://en.wikipedia.org/w/index.php?title=Trapezoidal_rule&oldid=771148256#Numerical_implementation</cite>
    ///
    template <typename Floating, typename F>
    auto periodic(const Floating from, const Floating to, const std::int_fast32_t n, const F f)
    {
        using int_type = std::int_fast32_t;
        using counter = stl::arithmetic_progression_iterator<int_type>;
        using result_type = decltype(f(std::declval<Floating>()));

        static constexpr Floating one = static_cast<Floating>(1);
        static constexpr Floating two = static_cast<Floating>(2);

        auto x_i_ab = [=](const int_type i) { return from + i * (to - from) / (n); };

        auto w_i0 = []() { return one; };
        auto w_in = []() { return one; };
        auto w_i = [](const int_type) { return two; };

        return (std::accumulate(++counter{ 0 },                                     // Left boundary condition treated in zero elem
                                counter{ n },                                       // Right boundary condition treated after std::accumulate
                                static_cast<result_type>(w_i0() * f(x_i_ab(0))),    // Zero elem is left boundary evaluation
                                [=](const result_type& sum,
                                    const int_type& i)
        {
            return sum + w_i(i) * f(x_i_ab(i));                                     // General case
        })
            + w_in() * f(x_i_ab(n)))                                                // Right boundary condition
            * (to - from) / (two * (n));                                            // Normalizing factor
    };


    template <typename Floating, typename F>
    auto chebysev(const Floating from, const Floating to, const std::int_fast32_t n, const F f)
    {
        if (n % 2) throw std::domain_error{ "math::chebysev only accepts even integration points" };

        using int_type = std::int_fast32_t;
        using counter = stl::arithmetic_progression_iterator<int_type>;
        using result_type = decltype(f(std::declval<Floating>()));
    
        static constexpr Floating one = static_cast<Floating>(1);
        static constexpr Floating two = static_cast<Floating>(2);
        static constexpr Floating one_per_two = one / two;
    
        auto x_i = [&](const int_type i) { return std::cos(pi<Floating> * i / n); };
        auto x_i_ab = [&](const int_type i) { return (x_i(i) + one) * one_per_two * (to - from) + from; };
    
        auto w_i0 = [&]() { return one; };
        auto w_in = [&]() { return one / (1 - n*n); };
        auto w_i  = [&](const int_type i) { return two / (1 - i*i); };
    
        auto w_ij0 = [&](const int_type i) { return two * std::cos(pi<Floating> * i * 0 / n) / n; };
        auto w_ijn = [&](const int_type i) { return two * std::cos(pi<Floating> * i * n / n) / n; };
        auto w_ij  = [&](const int_type i, const int_type j) { return std::cos(pi<Floating> * i * j / n) / n; };
    
        auto inner_sum = [=](const int_type i)
        {
            return std::accumulate(++counter{ 0 },              // Left boundary condition treated in zero elem
                                   counter{ n },                // Right boundary condition treated after std::accumulate
                                   w_ij0(i) * f(x_i_ab(0)),     // Zero elem is left boundary evaluation
                                   [=](const result_type& sum,
                                       const int_type& j)
            {
                return sum + w_ij(i, j) * f(x_i_ab(j));         // General case
            })
                + w_ij0(i) * f(x_i_ab(n));                      // Right boundary condition
        };
    
        auto outer_sum = [=](const auto inner)
        {
            return std::accumulate(++counter{ 0, 2 },           // Left boundary condition treated in zero elem
                                   counter{ n, 2 },             // Right boundary condition treated after std::accumulate
                                   w_i0() * inner(0),           // Zero elem is left boundary evaluation
                                   [=](const result_type& cum_sum,
                                       const int_type& i)
            {
                return cum_sum + w_i(i) * inner(i);             // General case
            })
                + w_in() * inner(n);                            // Right boundary condition
        };
    
        return outer_sum(inner_sum);
    };
}
