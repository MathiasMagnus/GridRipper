#include <STL_Test0_Quadrature.hpp>

int main()
{
    using integral = std::int_fast32_t;     // Type used to represent integral values
    using real = double;                    // Type used to represent real values
    using complex = std::complex<real>;     // Type used to represent complex values
    using counter = stl::arithmetic_progression_iterator<integral>; // Iterator used as series of numbers

    constexpr real tolerance = 1e-6;

    constexpr auto pi = math::pi<real>;
    auto cos = [](const real& x) { return std::cos(x); };
    auto const_one = [](const real&) { return (real)1; };

    // Calulate pi by integrating a constant function
    auto calc_pi = [=](integral N)
    {
        return math::periodic<real>(0, pi, N, const_one);
    };

    std::cout <<
        "math::periodic<real>(0, pi, N, const_one) converged with " <<
        *std::adjacent_find(counter{ 2, 100 }, counter{}, [=](const integral& n, const integral& m) { return std::abs(calc_pi(n) - calc_pi(m)) < tolerance; }) <<
        " function evaluations" << std::endl;
    /*
    std::vector<real> almost_pi;

    std::transform(counter{3, 5},
                   counter{33, 5},
                   std::back_inserter(almost_pi),
                   [=](const integral& n)
    {
        return calc_pi(n);
    });

    auto is_any_too_far_from_pi = std::any_of(almost_pi.cbegin(),
                                              almost_pi.cend(),
                                              [=](const real& val)
    {
        return std::abs(val - pi) > 1e-6;
    });

    if (is_any_too_far_from_pi)
        std::exit(EXIT_FAILURE);
    */
    
    // Calculate 0 by integrating the cosine function from -pi to pi
    //auto calc_zero = [=](integral N) { return math::chebysev<real>(-pi, pi, N, cos); };
    auto calc_zero = [=](integral N) { return math::chebysev<real>(0, 2 * pi, N, cos); };

    std::vector<real> almost_zero;
    
    std::transform(counter{4, 54, 10},
                   counter{},
                   std::back_inserter(almost_zero),
                   [=](const integral& n)
    {
        return calc_zero(n);
    });
    
    auto is_any_too_far_from_zero = std::any_of(almost_zero.cbegin(),
                                                almost_zero.cend(),
                                                [](const real& val)
    {
        return std::abs(val) > 1e-6;
    });

    if (is_any_too_far_from_zero)
        std::exit(EXIT_FAILURE);
    
    // Calculate 4pi by integrating the surface integral of the constant 1 function on a sphere of unit radius
    auto calc_4pi = [](integral M, integral N)
    {
        return math::chebysev<real>(0, math::pi<real>, M, [=](const real theta)
        {
            return math::periodic<real>(0, 2 * math::pi<real>, N, [=](const real phi)
            {
                return (real)1 * std::sin(theta);
            });
        });
    };
    
    if (std::abs(calc_4pi(30, 30) - 4 * math::pi<real>) > 1e-6)
        std::exit(EXIT_FAILURE);
    
	return EXIT_SUCCESS;
}