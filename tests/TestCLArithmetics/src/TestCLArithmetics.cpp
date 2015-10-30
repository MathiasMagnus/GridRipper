#include <TestSTLArithmetics.hpp>

const Multipole::stl::Radial::Index radial = 1024;
const Multipole::stl::Gaunt::Index combined = 1024;

int main(int argc, char* argv[])
{
    Gripper::stl::initialize(argc, argv);

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing self-contained operators";

        using namespace Multipole::stl;

        Radial::Vector r_neg(radial);
        Radial::Vector r_one(radial);
        Radial::Vector r_two(radial);
        Spherical::Vector s_neg(combined);
        Spherical::Vector s_one(combined);
        Spherical::Vector s_two(combined);
        Field f_neg(radial, combined);
        Field f_one(radial, combined);
        Field f_two(radial, combined);

        for (auto& elem : r_neg) elem = -1.0;
        for (auto& elem : r_one) elem = 1.0;
        for (auto& elem : r_two) elem = 2.0;
        for (auto& elem : s_neg) elem = -1.0;
        for (auto& elem : s_one) elem = 1.0;
        for (auto& elem : s_two) elem = 2.0;
        for (auto& elem : f_neg) elem = -1.0;
        for (auto& elem : f_one) elem = 1.0;
        for (auto& elem : f_two) elem = 2.0;

        if (!differ(r_one, +r_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Radial::Vector) FAILED";

        if (!differ(r_neg, -r_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Radial::Vector) FAILED";

        if (!differ(r_one + r_one, r_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) FAILED";

        if (!differ(r_two - r_one, r_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) FAILED";

        if (!differ(r_one * r_two, r_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) FAILED";

        if (!differ(r_two / r_two, r_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Radial::Vector) FAILED";

        if (!differ(s_one, +s_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(s_neg, -s_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(s_one + s_one, s_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(s_two - s_one, s_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(s_one * s_two, s_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(s_two / s_two, s_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Spherical::Vector) FAILED";

        if (!differ(f_one, +f_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Field) FAILED";

        if (!differ(f_neg, -f_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Field) FAILED";

        if (!differ(f_one + f_one, f_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Field, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Field, Multipole::stl::Field) FAILED";

        if (!differ(f_two - f_one, f_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Field, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Field, Multipole::stl::Field) FAILED";

        if (!differ(f_one * f_two, f_two))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Field, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Field, Multipole::stl::Field) FAILED";

        if (!differ(f_two / f_two, f_one))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Field, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Field, Multipole::stl::Field) FAILED";
    }

    {
        using namespace Multipole::stl;

        Field source(radial, combined);
        Field result(radial, combined);

        for (auto& elem : source) elem = 0.0;

        for(Gaunt::Index i = 0 ; i < source.sphericalSize() ; ++i)
            for (Radial::Index r = 0; r < source.radialSize(); ++r)
            {
                if (Spherical::IndexPair(i) == Spherical::IndexPair(1, -1)) source.at(r, i) = 1.0;
                if (Spherical::IndexPair(i) == Spherical::IndexPair(1, +1)) source.at(r, i) = 1.0;

                if (Spherical::IndexPair(i) == Spherical::IndexPair(0,  0)) result.at(r, i) = 1.0;
            }

        if (!differ(pown(source, 2), result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "Multipole::stl::pown(Multipole::stl::Field, int) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "Multipole::stl::pown(Multipole::stl::Field, int) FAILED";
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Vector,Vector operators";

        using namespace Multipole::stl;

        Radial::Vector inv_Rsq(radial.r);
        Spherical::Vector laplace_s2(combined.i);
        Field result(radial.r, combined.i);
        
        for (Radial::Index r = 0; r < radial; ++r) inv_Rsq.at(r) = pown(r, -2);
        for (Gaunt::Index i = 0; i < combined; ++i)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            laplace_s2.at(i) = l * (l + 1);
        }

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq + (l * (l + 1));
        }

        if (!differ(inv_Rsq + laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq - (l * (l + 1));
        }

        if (!differ(inv_Rsq - laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq * (l * (l + 1));
        }

        if (!differ(inv_Rsq * laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq / (l * (l + 1));
        }

        if (!differ(inv_Rsq / laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Spherical::Vector) FAILED";

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) + inv_rsq;
        }

        if (!differ(laplace_s2 + inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) - inv_rsq;
        }

        if (!differ(laplace_s2 - inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) * inv_rsq;
        }

        if (!differ(laplace_s2 * inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) / inv_rsq;
        }

        if (!differ(laplace_s2 / inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Radial::Vector) FAILED";
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Radial,Field operators";

        using namespace Multipole::stl;

        Radial::Vector inv_Rsq(radial.r);
        Field result(radial.r, combined.i);
        Field laplace_s2(radial.r, combined.i);

        for (Radial::Index r = 0; r < radial; ++r) inv_Rsq.at(r) = pown(r, -2);
        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            laplace_s2.at(r, i) = l * (l + 1);
        }

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq + (l * (l + 1));
        }

        if (!differ(inv_Rsq + laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Radial::Vector, Multipole::stl::Field) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq - (l * (l + 1));
        }

        if (!differ(inv_Rsq - laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Radial::Vector, Multipole::stl::Field) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq * (l * (l + 1));
        }

        if (!differ(inv_Rsq * laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Radial::Vector, Multipole::stl::Field) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq / (l * (l + 1));
        }

        if (!differ(inv_Rsq / laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Radial::Vector, Multipole::stl::Field) FAILED";

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) + inv_rsq;
        }

        if (!differ(laplace_s2 + inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Field, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Field, Multipole::stl::Radial::Vector) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) - inv_rsq;
        }

        if (!differ(laplace_s2 - inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Field, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Field, Multipole::stl::Radial::Vector) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) * inv_rsq;
        }

        if (!differ(laplace_s2 * inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Field, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Field, Multipole::stl::Radial::Vector) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) / inv_rsq;
        }

        if (!differ(laplace_s2 / inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Field, Multipole::stl::Radial::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Field, Multipole::stl::Radial::Vector) FAILED";
    }

    {
        Gripper::stl::clog << Gripper::LogLevel::INFO << "Testing mixed Spherical,Field operators";

        using namespace Multipole::stl;

        Spherical::Vector laplace_s2(combined);
        Field result(radial, combined);
        Field inv_Rsq(radial, combined);

        for (Gaunt::Index i = 0; i < combined; ++i)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            laplace_s2.at(i) = l * (l + 1);
        }
        for (Gaunt::Index i = 0; i < combined; ++i)
            for (Radial::Index r = 0; r < radial; ++r)
                inv_Rsq.at(r, i) = pown(r, -2);

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) + inv_rsq;
        }

        if (!differ(laplace_s2 + inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Spherical::Vector, Multipole::stl::Field) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) - inv_rsq;
        }

        if (!differ(laplace_s2 - inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Spherical::Vector, Multipole::stl::Field) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) * inv_rsq;
        }

        if (!differ(laplace_s2 * inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Spherical::Vector, Multipole::stl::Field) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = (l * (l + 1)) / inv_rsq;
        }

        if (!differ(laplace_s2 / inv_Rsq, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Field) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Spherical::Vector, Multipole::stl::Field) FAILED";

        // operator+

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq + (l * (l + 1));
        }

        if (!differ(inv_Rsq + laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator+(Multipole::stl::Field, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator+(Multipole::stl::Field, Multipole::stl::Spherical::Vector) FAILED";

        // operator-

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq - (l * (l + 1));
        }

        if (!differ(inv_Rsq - laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator-(Multipole::stl::Field, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator-(Multipole::stl::Field, Multipole::stl::Spherical::Vector) FAILED";

        // operator*

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq * (l * (l + 1));
        }

        if (!differ(inv_Rsq * laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator*(Multipole::stl::Field, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator*(Multipole::stl::Field, Multipole::stl::Spherical::Vector) FAILED";

        // operator/

        for (Gaunt::Index i = 0; i < combined; ++i) for (Radial::Index r = 0; r < radial; ++r)
        {
            ValueType l = Spherical::IndexPair(i).l.i;
            ValueType inv_rsq = pown(r, -2);
            result.at(r, i) = inv_rsq / (l * (l + 1));
        }

        if (!differ(inv_Rsq / laplace_s2, result))
            Gripper::stl::clog << Gripper::LogLevel::INFO << "operator/(Multipole::stl::Field, Multipole::stl::Spherical::Vector) PASSED";
        else
            Gripper::stl::clog << Gripper::LogLevel::CRIT << "operator/(Multipole::stl::Field, Multipole::stl::Spherical::Vector) FAILED";
    }
	
	return EXIT_SUCCESS;
}