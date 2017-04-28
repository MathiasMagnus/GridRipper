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
		// TODO: Come up with something better for this test that is actually meaningful
		//
        //if (f(from) != f(to)) throw std::domain_error{ "Function provided to math::periodic is not periodic on the given interval." };

        using int_type = std::int_fast32_t;
        using counter_type = stl::arithmetic_progression_iterator<int_type>;
        using result_type = decltype(f(std::declval<Floating>()));

        static constexpr Floating one = static_cast<Floating>(1);
        static constexpr Floating two = static_cast<Floating>(2);

        auto x_i_ab = [=](const int_type i) { return from + i * (to - from) / (n); };

        auto w_i0 = []() { return one; };
        auto w_in = []() { return one; };
        auto w_i = [](const int_type) { return two; };

        return (std::accumulate(++counter_type{ 0, n },                             // Left boundary condition treated in zero elem
                                counter_type{},                                     // Right boundary condition treated after std::accumulate
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
	auto chebysev1(const Floating from, const Floating to, const std::int_fast32_t n, const F f)
	{
		//if (n % 2) throw std::domain_error{ "math::chebysev only accepts even integration points" };

		using int_type = std::int_fast32_t;
		using counter_type = stl::arithmetic_progression_iterator<int_type>;
		using result_type = decltype(f(std::declval<Floating>()));

		static constexpr Floating one = static_cast<Floating>(1);
		static constexpr Floating two = static_cast<Floating>(2);
		static constexpr Floating one_per_two = one / two;

		const Floating b_minus_a_per_two = (to - from) / two;
		const Floating b_plus_a_per_two = (to + from) / two;
		const Floating w_i = pi<Floating> / n;

		auto x_i = [=](const int_type i) { return std::cos((i - one_per_two) / (n) * pi<Floating>); };
		auto x_i_ab = [=](const int_type i) { return b_minus_a_per_two * x_i(i) + b_plus_a_per_two; };

		return (b_minus_a_per_two) * std::accumulate(counter_type{ 1, n },
													 counter_type{},
												     static_cast<result_type>(0),
													 [=](const result_type& sum,
														 const int_type& I)
		{
			return sum + w_i * std::sqrt(1 - x_i(I)*x_i(I)) * f(x_i_ab(I));
		});
	};


    template <typename Floating, typename F>
    auto chebysev2(const Floating from, const Floating to, const std::int_fast32_t n, const F f)
    {
        if (n % 2) throw std::domain_error{ "math::chebysev only accepts even integration points" };
	
        using int_type = std::int_fast32_t;
        using counter_type = stl::arithmetic_progression_iterator<int_type>;
        using result_type = decltype(f(std::declval<Floating>()));
    
        static constexpr Floating one = static_cast<Floating>(1);
        static constexpr Floating two = static_cast<Floating>(2);

        /// <notes>
        ///   The original integrator works for input function [-1:1]. We transform this into the input range.
        ///   <c>f_new</c> is indexed on the old range but return values from the input function over the input range.
        /// </notes>
        /// <cite>http://stackoverflow.com/a/929107</cite>
        ///
        Floating old_range = static_cast<Floating>(one - -one);
        Floating new_range = static_cast<Floating>(to - from);
        auto new_x = [=](const Floating& old_val) { return (((old_val - -one) * new_range) / old_range) + from; };
        auto f_new = [=](const Floating& old_x) { return f(new_x(old_x)); };

        auto x_i =   [=](const int_type i) { return std::cos(pi<Floating> * i / n); };
    
        auto w_i0 = [=]() { return one; };
        auto w_in = [=]() { return one / (1 - n*n); };
        auto w_i  = [=](const int_type i) { return two / (1 - i*i); };
    
        auto w_ij0 = [=](const int_type i) { return std::cos(pi<Floating> * i * 0 / n) / n; };
        auto w_ijn = [=](const int_type i) { return std::cos(pi<Floating> * i * n / n) / n; };
        auto w_ij  = [=](const int_type i, const int_type j) { return two * std::cos(pi<Floating> * i * j / n) / n; };
    
        auto inner_sum = [=](const int_type i)
        {
            return std::accumulate(++counter_type{ 0, n },      // Left boundary condition treated in zero elem
                                   counter_type{},              // Right boundary condition treated after std::accumulate
                                   w_ij0(i) * f_new(x_i(0)),    // Zero elem is left boundary evaluation
                                   [=](const result_type& sum,
                                       const int_type& j)
            {
                return sum + w_ij(i, j) * f_new(x_i(j));        // General case
            })
                + w_ijn(i) * f_new(x_i(n));                     // Right boundary condition
        };
    
        auto outer_sum = [=](const auto inner)
        {
            return std::accumulate(++counter_type{ 0, n, 2 },   // Left boundary condition treated in zero elem
                                   counter_type{},              // Right boundary condition treated after std::accumulate
                                   w_i0() * inner(0),           // Zero elem is left boundary evaluation
                                   [=](const result_type& cum_sum,
                                       const int_type& i)
            {
                return cum_sum + w_i(i) * inner(i);             // General case
            })
                + w_in() * inner(n);                            // Right boundary condition
        };
    
        return outer_sum(inner_sum) * (new_range / old_range);
    };


    template <typename Floating, typename F>
    auto chebysev3(const Floating from, const Floating to, const std::int_fast32_t N, const F f)
    {
        Floating old_range = static_cast<Floating>(1 - -1);
        Floating new_range = static_cast<Floating>(to - from);
        auto new_x = [=](const Floating& old_val) { return (((old_val - -1) * new_range) / old_range) + from; };
        auto f_new = [=](const Floating& old_x) { return f(new_x(old_x)); };

        //auto f = [](auto x) { return std::exp(x*0.5 + 0.5); };
        auto xi = [](int i, int n) { return std::cos(pi<Floating>*i / n); };
        auto wij = [](int i, int j, int n) { auto tmp = std::cos(pi<Floating>*i*j / n) / n; if (j == 0 || j == n) { return tmp; } else { return tmp*2.0; } };

        Floating sum = 0;
        //N even!
        for (int i = 0; i <= N; i += 2)
        {
            auto a = 0.0;
            for (int j = 0; j <= N; ++j)
            {
                a += wij(i, j, N) * f_new(xi(j, N));
            }

            double c = 1.0;
            if (i != 0) { c = 2.0 / (1 - i*i); }
            if (i == N) { c = 1.0 / (1 - i*i); }//ez hianyzott a multkor
            sum += c*a;
        }

        return sum * (new_range / old_range);
    }

    //template <typename Floating, typename F>
    //auto chebysev2(const Floating from, const Floating to, const std::int_fast32_t n, const F f)
    //{
    //    if (n % 2) throw std::domain_error{ "math::chebysev only accepts even integration points" };
    //
    //    using int_type = std::int_fast32_t;
    //    using counter_type = stl::arithmetic_progression_iterator<int_type>;
    //    using result_type = decltype(f(std::declval<Floating>()));
    //
    //    static constexpr Floating one = static_cast<Floating>(1);
    //    static constexpr Floating two = static_cast<Floating>(2);
    //    static constexpr Floating one_per_two = one / two;
    //
    //    auto x_i = [&](const int_type i) { return std::cos(pi<Floating> * i / n); };
    //    auto x_i_ab = [&](const int_type i) { return (x_i(i) + one) * one_per_two * (to - from) + from; };
    //
    //    auto w_i0 = [&]() { return one; };
    //    auto w_in = [&]() { return one / (1 - n*n); };
    //    auto w_i = [&](const int_type i) { return two / (1 - i*i); };
    //
    //    auto w_ij0 = [&](const int_type i) { return two * std::cos(pi<Floating> * i * 0 / n) / n; };
    //    auto w_ijn = [&](const int_type i) { return two * std::cos(pi<Floating> * i * n / n) / n; };
    //    auto w_ij = [&](const int_type i, const int_type j) { return std::cos(pi<Floating> * i * j / n) / n; };
    //
    //    auto inner_sum = [=](const int_type i)
    //    {
    //        return std::accumulate(++counter_type{ 0, n },      // Left boundary condition treated in zero elem
    //            counter_type{},              // Right boundary condition treated after std::accumulate
    //            w_ij0(i) * f(x_i_ab(0)),     // Zero elem is left boundary evaluation
    //            [=](const result_type& sum,
    //                const int_type& j)
    //        {
    //            return sum + w_ij(i, j) * f(x_i_ab(j));         // General case
    //        })
    //            + w_ij0(i) * f(x_i_ab(n));                      // Right boundary condition
    //    };
    //
    //    auto outer_sum = [=](const auto inner)
    //    {
    //        return std::accumulate(++counter_type{ 0, n, 2 },   // Left boundary condition treated in zero elem
    //            counter_type{},              // Right boundary condition treated after std::accumulate
    //            w_i0() * inner(0),           // Zero elem is left boundary evaluation
    //            [=](const result_type& cum_sum,
    //                const int_type& i)
    //        {
    //            return cum_sum + w_i(i) * inner(i);             // General case
    //        })
    //            + w_in() * inner(n);                            // Right boundary condition
    //    };
    //
    //    return outer_sum(inner_sum);
    //};
}
