#include <STL_Test2_Arithmetics_Simple.hpp>

int main()
{
    // Simulation params
    using integral = std::int32_t;      // Type used to represent integral values
    using real = double;                // Type used to represent real values
    using complex = std::complex<real>; // Type used to represent complex values
    constexpr integral L_max = 7;       // Maximum multipole values for series expansion
    const real neumann_length = 1e-4;   // Neumann-series expansion precision
    const real tolerance = 1e-3;        // Tolerance in sum of absoltue differences

                                        // Simulation type aliases
    using namespace Multipole::stl;
    using spherical_extent = Spherical::Extent<L_max, integral>;
    using spherical_index = Spherical::Index<L_max, integral>;
    using gaunt_matrix = Gaunt::Matrix<L_max, integral, real>;
    using real_coeff_vector = Spherical::Vector<L_max, Parity::Even, integral, real>;
    using complex_coeff_vector = Spherical::Vector<L_max, Parity::Even, integral, complex>;

    // Application type aliases
    using timer = std::chrono::high_resolution_clock;
    using time_point = std::chrono::high_resolution_clock::time_point;

    // Simulation variable declarations
    constexpr real one = 1, two = 2, three = 3, four = 4, eight = 8, one_per_two = one / two, three_per_two = three / two;
    spherical_extent lm_ext;
    real_coeff_vector result, serial_index, negative_serial_index, zeroes, a, b, c, d, e, f, up, down;
    spherical_index ai, bi, ci1, ci2, ci3, ci4, di, ei1, ei2, ei3, ei4, fi;
    gaunt_matrix gaunt;

    // Initialize states
    {
        lm_ext = spherical_extent({ 0, 0 }, { L_max, L_max });
        serial_index = real_coeff_vector(lm_ext);
        negative_serial_index = real_coeff_vector(lm_ext);
        zeroes = real_coeff_vector(lm_ext);
        a = real_coeff_vector(lm_ext);
        b = real_coeff_vector(lm_ext);
        c = real_coeff_vector(lm_ext);
        d = real_coeff_vector(lm_ext);
        e = real_coeff_vector(lm_ext);
        f = real_coeff_vector(lm_ext);
        up = real_coeff_vector(lm_ext);
        down = real_coeff_vector(lm_ext);
        ai = { 3, -3 };
        bi = { 3, 3 };
        ci1 = { 0, 0 };
        ci2 = { 2, 0 };
        ci3 = { 4, 0 };
        ci4 = { 6, 0 };
        di = { 0, 0 };
        ei1 = { 0, 0 };
        ei2 = { 2, 0 };
        ei3 = { 4, 0 };
        ei4 = { 6, 0 };
        fi = { 0, 0 };

        // Unary/binary init
        integral count = 0;
        for (spherical_index i = lm_ext.initial(); lm_ext.contains(i); ++i)
        {
            zeroes.at(i) = 0;
            serial_index.at(i) = count;
            negative_serial_index.at(i) = -count;

            ++count;
        }

        // Contraction init
        gaunt = gaunt_matrix(lm_ext);
        for (spherical_index i = lm_ext.initial(); lm_ext.contains(i); ++i)
        {
            a.at(i) = (i == ai ? 1 : 0);
            b.at(i) = (i == bi ? 1 : 0);
            if (i == ci1)
                c.at(i) = -0.28209479177387797311;
            else if (i == ci2)
                c.at(i) = 0.21026104350167990065;
            else if (i == ci3)
                c.at(i) = -0.076934943211057260637;
            else if (i == ci4)
                c.at(i) = 0.011854396693264072568;
            else
                c.at(i) = 0;
        }

        for (spherical_index i = lm_ext.initial(); lm_ext.contains(i); ++i)
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
    if (!test(c, Gaunt::contract(gaunt, a, b), tolerance)) return EXIT_FAILURE;
    if (!test(f, Gaunt::neumann(gaunt, d, e, neumann_length), tolerance)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
