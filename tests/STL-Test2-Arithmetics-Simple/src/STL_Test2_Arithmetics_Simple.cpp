#include <STL_Test2_Arithmetics_Simple.hpp>

#include <fstream>

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
        gaunt.calculate();
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
    if (!test(up, SpinWeightedSpherical::eth(serial_index), tolerance)) return EXIT_FAILURE;
    if (!test(down, SpinWeightedSpherical::eth_bar(serial_index), tolerance)) return EXIT_FAILURE;

    /*
    Gripper::stl::initialize(argc, argv);

    auto L_max = Gripper::stl::getL();

    Gripper::stl::clog << Gripper::DEBUG << "L_max = " << L_max.i;

    const Multipole::stl::Radial::Extent radial(0, 128);
    const Multipole::stl::Spherical::Extent spherical(0, L_max);
    const Multipole::stl::Expansion::Extent combined(radial, spherical);

    using namespace Multipole::stl;
    using namespace Multipole::stl::Expansion;
    typedef Field<double>::index_type IndexType;

    auto start = std::chrono::high_resolution_clock::now();

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing self-contained operators";

        Radial::Vector<double> r_neg(radial);
        Radial::Vector<double> r_one(radial);
        Radial::Vector<double> r_two(radial);
        Spherical::Vector<double> s_neg(spherical);
        Spherical::Vector<double> s_one(spherical);
        Spherical::Vector<double> s_two(spherical);
        Field<double> f_neg(combined);
        Field<double> f_one(combined);
        Field<double> f_two(combined);

        for (auto& elem : r_neg) elem = -1.0;
        for (auto& elem : r_one) elem = 1.0;
        for (auto& elem : r_two) elem = 2.0;
        for (auto& elem : s_neg) elem = -1.0;
        for (auto& elem : s_one) elem = 1.0;
        for (auto& elem : s_two) elem = 2.0;
        for (auto& elem : f_neg) elem = -1.0;
        for (auto& elem : f_one) elem = 1.0;
        for (auto& elem : f_two) elem = 2.0;

        logTest(differ(r_one, +r_one), "operator+(Multipole::stl::Radial::Vector)");
        logTest(differ(r_neg, -r_one), "operator-(Multipole::stl::Radial::Vector)");

        logTest(differ(r_one + r_one, r_two), "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector)");
        logTest(differ(r_two - r_one, r_one), "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector)");
        logTest(differ(r_one * r_two, r_two), "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector)");
        logTest(differ(r_two / r_two, r_one), "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector)");

        logTest(differ(s_one, +s_one), "operator+(Multipole::stl::Spherical::Vector)");
        logTest(differ(s_neg, -s_one), "operator-(Multipole::stl::Spherical::Vector)");

        logTest(differ(s_one + s_one, s_two), "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector)");
        logTest(differ(s_two - s_one, s_one), "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector)");
        logTest(differ(s_one * s_two, s_two), "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector)");
        logTest(differ(s_two / s_two, s_one), "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector)");

        logTest(differ(f_one, +f_one), "operator+(Multipole::stl::Field)");
        logTest(differ(f_neg, -f_one), "operator-(Multipole::stl::Field)");

        logTest(differ(f_one + f_one, f_two), "operator+(Multipole::stl::Field, Multipole::stl::Field)");
        logTest(differ(f_two - f_one, f_one), "operator-(Multipole::stl::Field, Multipole::stl::Field)");
    }

    {
        Field<double> source1(combined);
        Field<double> source2(combined);
        Field<double> result(combined);

        for (Field<double>::radial_index_type r = 0; source1.extent().radial().contains(r); ++r)
            for (Field<double>::spherical_index_type i = 0; source1.extent().spherical().contains(i); ++i)
            {
                IndexType index(r, i);

                if (Spherical::IndexPair(i) == Spherical::IndexPair(1, -1))
                    source1.at(index) = 1;
                else 
                    source1.at(index) = 0;

                if (Spherical::IndexPair(i) == Spherical::IndexPair(1, +1))
                    source2.at(index) = 1;
                else
                    source2.at(index) = 0;

                if (Spherical::IndexPair(i) == Spherical::IndexPair(0, 0))
                    result.at(index) = -0.282095;
                else if 
                    (Spherical::IndexPair(i) == Spherical::IndexPair(2, 0))
                    result.at(index) = 0.126157;
                else
                    result.at(index) = 0;
            }

        logTest(differ(Field<double>(source1 * source2), result), "operator*(Multipole::stl::Field, Multipole::stl::Field)");
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Vector,Vector operators";

        Radial::Vector<double> Rsq(radial);
        Spherical::Vector<double> laplace_s2(spherical);
        Field<double> result(combined);
        
        for (Radial::Vector<double>::index_type r = 0; radial.contains(r); ++r)
            Rsq.at(r) = pown(r, 2);

        for (Spherical::Vector<double>::index_type i = 0; spherical.contains(i); ++i)
        {
            auto l = Spherical::IndexPair(i).l.i;
            laplace_s2.at(i) = l * (l + 1);
        }

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq * (l * (l + 1));
            }

        logTest(differ(Rsq * laplace_s2, result), "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq / (l * (l + 1));
            }

        logTest(differ(Rsq / laplace_s2, result), "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector)");

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) * rsq;
            }

        logTest(differ(laplace_s2 * Rsq, result), "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) / rsq;
            }

        logTest(differ(laplace_s2 / Rsq, result), "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector)");
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Radial,Field operators";

        using namespace Multipole::stl;

        Radial::Vector<double> Rsq(radial);
        Field<double> laplace_s2(combined);
        Field<double> result(combined);

        for (Radial::Vector<double>::index_type r = 0; radial.contains(r); ++r)
            Rsq.at(r) = pown(r, 2);

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                laplace_s2.at(IndexType(r, i)) = l * (l + 1);
            }

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq * (l * (l + 1));
            }

        logTest(differ(Rsq * laplace_s2, result), "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Field)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq / (l * (l + 1));
            }

        logTest(differ(Rsq / laplace_s2, result), "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Field)");

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) * rsq;
            }

        logTest(differ(laplace_s2 * Rsq, result), "operator*(Multipole::stl::Field, Multipole::stl::Radial::Vector)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) / rsq;
            }

        logTest(differ(laplace_s2 / Rsq, result), "operator/(Multipole::stl::Field, Multipole::stl::Radial::Vector)");
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Spherical,Field operators";

        using namespace Multipole::stl;

        Spherical::Vector<double> laplace_s2(spherical);
        Field<double> Rsq(combined);
        Field<double> result(combined);

        for (Gaunt::Index i = 0; spherical.contains(i); ++i)
        {
            auto l = Spherical::IndexPair(i).l.i;
            laplace_s2.at(i) = l * (l + 1);
        }

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
                Rsq.at(IndexType(r, i)) = pown(r, 2);

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) * rsq;
            }

        logTest(differ(laplace_s2 * Rsq, result), "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Field)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = (l * (l + 1)) / rsq;
            }

        logTest(differ(laplace_s2 / Rsq, result), "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Field)");

        // operator*

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq * (l * (l + 1));
            }

        logTest(differ(Rsq * laplace_s2, result), "operator*(Multipole::stl::Field, Multipole::stl::Spherical::Vector)");

        // operator/

        for (Radial::Index r = 0; radial.contains(r); ++r)
            for (Gaunt::Index i = 0; spherical.contains(i); ++i)
            {
                auto l = Spherical::IndexPair(i).l.i;
                auto rsq = pown(r, 2);
                result.at(IndexType(r, i)) = rsq / (l * (l + 1));
            }

        logTest(differ(Rsq / laplace_s2, result), "operator/(Multipole::stl::Field, Multipole::stl::Spherical::Vector)");
    }

    auto end = std::chrono::high_resolution_clock::now();

    Gripper::stl::clog << Gripper::TIMING << "Tests took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds.";
	*/
	return EXIT_SUCCESS;
}