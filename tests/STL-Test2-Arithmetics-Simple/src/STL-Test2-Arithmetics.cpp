#include <STL-Test2-Arithmetics.hpp>


int main(int argc, char* argv[])
{
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
	
	return EXIT_SUCCESS;
}