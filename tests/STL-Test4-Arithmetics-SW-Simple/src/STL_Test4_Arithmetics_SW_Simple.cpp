#include <STL_Test4_Arithmetics_SW_Simple.hpp>

int main()
{
    // Simulation params
    using integral = std::int32_t;      // Type used to represent integral values
    using real = double;                // Type used to represent real values
    using complex = std::complex<real>; // Type used to represent complex values
    constexpr integral L_max = 7;       // Maximum multipole values for series expansion
    constexpr integral S_max = 3;       // Maximum spin values for series expansion
    const real neumann_length = 1e-4;   // Neumann-series expansion precision
    const real tolerance = 1e-3;        // Tolerance in sum of absoltue differences

                                        // Simulation type aliases
    using namespace Multipole::stl;
    using spherical_extent = SpinWeightedSpherical::Extent<L_max, S_max, integral>;
    using spherical_index = SpinWeightedSpherical::Index<L_max, S_max, integral>;
    using gaunt_matrix = SpinWeightedGaunt::Matrix<L_max, S_max, integral, real>;
    using real_coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, real>;
    using complex_coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, complex>;

    // Application type aliases
    using timer = std::chrono::high_resolution_clock;
    using time_point = std::chrono::high_resolution_clock::time_point;

    // Simulation variable declarations
    constexpr real one = 1, two = 2, three = 3, four = 4, eight = 8, one_per_two = one / two, three_per_two = three / two;
    spherical_extent lms_ext;
    real_coeff_vector result, serial_index, negative_serial_index, zeroes, a, b, c, d, e, f, up, down;
    spherical_index ai, bi, ci1, ci2, ci3, ci4, di, ei1, ei2, ei3, ei4, fi;
    gaunt_matrix gaunt;

    // Initialize states
    {
        lms_ext = spherical_extent({ 0, 0, 0 }, { L_max, L_max, S_max });
        serial_index = real_coeff_vector(lms_ext);
        negative_serial_index = real_coeff_vector(lms_ext);
        zeroes = real_coeff_vector(lms_ext);
        a = real_coeff_vector(lms_ext);
        b = real_coeff_vector(lms_ext);
        c = real_coeff_vector(lms_ext);
        d = real_coeff_vector(lms_ext);
        e = real_coeff_vector(lms_ext);
        f = real_coeff_vector(lms_ext);
        up = real_coeff_vector(lms_ext);
        down = real_coeff_vector(lms_ext);
        ai = { 3, -3, 1 };
        bi = { 3, 3, -1 };
        ci1 = { 0, 0, 0 };
        ci2 = { 2, 0, 0 };
        ci3 = { 4, 0, 0 };
        ci4 = { 6, 0, 0 };
        di = { 0, 0, 0 };
        ei1 = { 0, 0, 0 };
        ei2 = { 2, 0, 0 };
        ei3 = { 4, 0, 0 };
        ei4 = { 6, 0, 0 };
        fi = { 0, 0, 0 };

        // Unary/binary init
        integral count = 0;
        for (spherical_index i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            zeroes.at(i) = 0;
            serial_index.at(i) = count;
            negative_serial_index.at(i) = -count;

            ++count;
        }

        // Contraction init
        gaunt = gaunt_matrix(lms_ext);
        for (spherical_index i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            a.at(i) = (i == ai ? 1 : 0);
            b.at(i) = (i == bi ? 1 : 0);
            if (i == ci1)
                c.at(i) = 0.282095;
            else if (i == ci2)
                c.at(i) = -0.157696;
            else if (i == ci3)
                c.at(i) = 0.0128225;
            else if (i == ci4)
                c.at(i) = 0.0088908;
            else
                c.at(i) = 0;
        }

        for (spherical_index i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            d.at(i) = (i == di ? 1 : 0);
            if (i == ei1)
                e.at(i) = 0.5;
            else if (i == ei2)
                e.at(i) = 0.4;
            else if (i == ei3)
                e.at(i) = 0.2;
            else if (i == ei4)
                e.at(i) = 0.1;
            else
                e.at(i) = 0;
            f.at(i) = (i == fi ? 0.3352433995907 : 0);

        }

        // Edth init
        for (spherical_index i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            // Up
            if (i.s == std::min(i.l, up.s_max))
                up.at(i) = 0;
            else
            {
                auto factor_sq = (i.l - i.s) * (i.l + i.s + static_cast<integral>(1));

                if (factor_sq <= static_cast<integral>(0))
                    up.at(i) = 0;
                else
                {
                    auto index = i;

                    up.at(i) = std::sqrt(factor_sq) * serial_index.at(++index);
                }
            }

            // Down
            if (i.s == std::max(-i.l, -down.s_max))
                down.at(i) = 0;
            else
            {
                auto factor_sq = (i.l + i.s) * (i.l - i.s + static_cast<integral>(1));

                if (factor_sq <= static_cast<integral>(0))
                    down.at(i) = 0;
                else
                {
                    auto index = i;

                    down.at(i) = -std::sqrt(factor_sq) * serial_index.at(--index);
                }
            }
        }
    }

    // Unary operator tests
    if (!test(serial_index, +serial_index, tolerance)) return EXIT_FAILURE;
    if (!test(negative_serial_index, -serial_index, tolerance)) return EXIT_FAILURE;

    // Binary operator tests
    if (!test(zeroes, negative_serial_index + serial_index, tolerance)) return EXIT_FAILURE;
    if (!test(negative_serial_index, zeroes - serial_index, tolerance)) return EXIT_FAILURE;
    if (!test(negative_serial_index, -1 * serial_index, tolerance)) return EXIT_FAILURE;
    if (!test(negative_serial_index, serial_index / -1, tolerance)) return EXIT_FAILURE;

    // Contraction tests
    if (!test(c, SpinWeightedGaunt::contract(gaunt, a, b), tolerance)) return EXIT_FAILURE;
    if (!test(f, SpinWeightedGaunt::neumann(gaunt, d, e, neumann_length), tolerance)) return EXIT_FAILURE;

    // Spin tests
    if (!test(up, SpinWeightedSpherical::edth(serial_index), tolerance)) return EXIT_FAILURE;
    if (!test(down, SpinWeightedSpherical::edth_bar(serial_index), tolerance)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
