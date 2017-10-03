#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>
#include <Gripper/stl/stlMathConstants.hpp>
#include <Gripper/stl/stlYlms_dynamic.hpp>
#include <Gripper/stl/stlArithmeticProgression.hpp>

// Standard C++ includes
#include <cmath>            // std::tan, std::tgamma, std::beta (binom should be implemented via beta)
#include <numeric>          // std::accumulate
#include <type_traits>      // std::is_intergral
#include <exception>        // std::range_error
#include <string>           // std::to_string

namespace math
{
    template <typename Floating, typename Integral>
    auto fact(Integral n)
    {
        static_assert(std::is_integral<Integral>::value, "Value given to fact is not of integral type.");

        return std::tgamma(static_cast<Floating>(n + 1));
    }

    template <typename Integral>
    Integral binom(Integral n, Integral k)
    {
        static_assert(std::is_integral<Integral>::value, "Values given to binom are not of integral types.");

        // Shortcuts
        if (0 == k || n == k) {
            return 1;
        }
        if (k > n) {
            return 0;
        }
        if (k > (n - k)) {
            k = n - k;
        }
        if (1 == k) {
            return n;
        }

        Integral b = 1;

        for (Integral i = 1; i <= k; ++i)
        {
            b *= (n - (k - i));
            if (b < 0) throw std::range_error(std::string("Overflow in binom for values (") + std::to_string(n) + "," + std::to_string(k) + ")");
            b /= i;
        }

        return b;
    }

    template <typename Floating>
    auto cot(Floating x)
    {
        return std::tan(pi_2<Floating> - x);
    }

    template <typename Integral, typename F>
    auto sum(Integral from, Integral to, F f)
    {
        using result_type = decltype(f(std::declval<Integral>()));
        using counter_type = stl::arithmetic_progression_iterator<Integral>;

        return std::accumulate(counter_type{ from, ++to }, // Algorithms are non-inclusive, but math notation is, hence ++
                               counter_type{},
                               static_cast<result_type>(0),
                               [=](const result_type& cum_sum, const Integral& i) { return f(i) + cum_sum; });
    }
    /*
    template <typename Floating>
    auto Y_lms(Floating theta, Floating phi, int l, int m, int s)
    {
        using namespace std::complex_literals;

        // Helper functions
        //auto pown = [](Floating base, const int n)
        //{
        //    Floating result = 1.0;
        //
        //    for (int i = 0; i < n; ++i)
        //        result *= base;
        //
        //    return result;
        //};
        auto minus_one_pown = [](const int n) { return n % 2 ? static_cast<Floating>(-1) : static_cast<Floating>(1); };
        //auto cot_pown = [=](Floating radian, const int n) { return n < 0 ? pown(std::tan(radian), -n) : pown(cot(radian), n); };
        auto cot_pown = [=](Floating radian, const int n) { return n < 0 ? std::pow(std::tan(radian), -static_cast<Floating>(n)) : std::pow(cot(radian), static_cast<Floating>(n)); };
        auto trigon = [=](Floating radian, const int n)
        {
            return n < 0 ?
                std::pow(std::tan(radian), -static_cast<Floating>(n)) :
                std::pow(cot(radian), static_cast<Floating>(n));
        };

        // Equation
        auto sign = minus_one_pown(m);
        auto factor = std::sqrt((fact<Floating>(l + m) * fact<Floating>(l - m) * fact<Floating>(2 * l + 1)) /
                                (4 * pi<Floating> * fact<Floating>(l + s) * fact<Floating>(l - s)));
        //auto sinus = pown(std::sin(theta / 2), 2 * l);
        auto summation = sum(0, l - s, [=](const int r)
        {
            return
                binom(l - s, r) *
                binom(l + s, r + s - m) *
                minus_one_pown(l - r - s) *
                std::exp((phi * m) * std::complex<Floating>(0, 1)) *
                //cot_pown(theta / 2, 2 * r + s - m);
                cot_pown(theta / 2, 2 * r + s - m);
        });
        
        return
            sign *
            factor *
            sinus *
            summation;
    };

    template <typename Floating>
    auto Y_lm(Floating theta, Floating phi, int l, int m)
    {
        return Y_lms(theta, phi, l, m, 0);
    };
    */
}
