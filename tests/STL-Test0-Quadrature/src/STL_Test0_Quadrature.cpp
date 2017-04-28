#include <STL_Test0_Quadrature.hpp>

int main()
{
	floating::scoped_exception_enabler fpe; // Floating point exceptions trigger C++ exceptions

    using integral = std::int_fast32_t;     // Type used to represent integral values
    using real = double;                    // Type used to represent real values
    using complex = std::complex<real>;     // Type used to represent complex values
    using counter = stl::arithmetic_progression_iterator<integral>; // Iterator used as series of numbers

    constexpr real convergence = 1e-6;		// Convergence threshold
	constexpr real precision = 1e-5;		// Precision requirement

    constexpr auto pi = math::pi<real>;
	constexpr auto e = math::e<real>;
	const auto const_one = [](const real&) { return (real)1; };
    const auto cos = [](const real& x) { return std::cos(x); };
	const auto sin_sq = [](const real& x) { auto a = std::sin(x); return a*a; };
	const auto exp = [](const real& x) { return std::exp(x); };

	auto convergence_test = [=](auto f, auto target, counter from, counter to)
	{
		counter it = adjacent_find(from, to, [=](const integral& n, const integral& m)
		{
			return std::abs(f(n) - f(m)) < convergence;
		});

		if (it == to)
			throw std::runtime_error{ "Function did not converge on the specified index range." };

		if (std::abs(f(*it) - target) > precision)
			throw std::runtime_error{ "Function did not converge to the specified target value." };

		return *it;
	};

    auto convergence_report = [](std::string function_name, integral count)
    {
        std::cout << function_name << " converged with " << count << " function evaluations.\n" << std::endl;
    };

	auto convergence_test_2d = [=](auto f, auto target, counter from1, counter to1, counter from2, counter to2)
	{
		counter it2 = to2;
		counter it1 = std::adjacent_find(from1, to1, [=, &it2](const integral& n1, const integral m1) mutable
		{
			it2 = std::find_if(from2, to2, [=](const integral& i2) { return std::abs(f(n1, i2) - f(m1, i2)) < convergence; });

			if (it2 != to2)
				return true;
			else
				return false;
		});

		if (it1 == to1)
			throw std::runtime_error{ "Function did not converge on the specified index ranges." };

		if (it2 == to2)
			throw std::logic_error{ "it2 is not expected to be end iterator if it1 was not." };

		if (std::abs(f(*it1, *it2) - target) > precision)
			throw std::runtime_error{ "Function did not converge to the specified target value." };

		return std::make_pair(*it1, *it2);
	};

    auto convergence_report_2d = [](std::string function_name, std::pair<integral, integral> counts)
    {
        std::cout << function_name << std::endl << "converged with M = " << counts.first << " and N = " << counts.second << " function evaluations.\n" << std::endl;
    };

    // Calulate pi by integrating a constant function
    auto calc_pi_trivial = [=](integral N) { return math::periodic<real>(0, pi, N, const_one); };

    convergence_report("math::periodic<real>(0, pi, N, const_one)",
                       convergence_test(calc_pi_trivial, pi, counter{ 2, 100 }, counter{}));

	// Calulate pi by integrating sin squared
	auto calc_pi = [=](integral N) { return math::periodic<real>(0, 2 * pi, N, sin_sq); };

    convergence_report("math::periodic<real>(0, 2 * pi, N, sin_sq)",
                       convergence_test(calc_pi, pi, counter{ 2, 100 }, counter{}));

	// Calculate e - 1 by integrating exp
	auto calc_e_minus_one = [=](integral N) { return math::chebysev2<real>(0, 1, N, exp); };

    convergence_report("math::chebysev<real>(0, 1, N, exp)",
                       convergence_test(calc_e_minus_one, e - 1, counter{ 4, 400, 2 }, counter{}));
    
    // Calculate 0 by integrating the cosine function from -pi to pi
    auto calc_zero = [=](integral N) { return math::chebysev2<real>(-pi, pi, N, cos); };

    convergence_report("math::chebysev<real>(-pi, pi, N, cos)",
                       convergence_test(calc_zero, 0, counter{ 4, 400, 2 }, counter{}));
    
    // Calculate 4pi by integrating the surface integral of the constant 1 function on a sphere of unit radius
    auto calc_4pi = [=](integral M, integral N)
    {
        return math::chebysev2<real>(0, pi, M, [=](const real theta)
        {
            return math::periodic<real>(0, 2 * pi, N, [=](const real phi)
            {
                return (real)1 * std::sin(theta);
            });
        });
    };

    convergence_report_2d(std::string{ "math::chebysev<real>(0, pi, M, [=](const real theta)\n" } +
                                           "\treturn math::periodic<real>(0, 2 * pi, N, [=](const real phi)\n" +
                                           "\t\treturn (real)1 * std::sin(theta);\n" +
                                           "\t});\n" +
                                       "})",
                          convergence_test_2d(calc_4pi, 4 * pi,
                                              counter{ 2, 400, 2 }, counter{},
                                              counter{ 2, 400 }, counter{}));
    
	return EXIT_SUCCESS;
}