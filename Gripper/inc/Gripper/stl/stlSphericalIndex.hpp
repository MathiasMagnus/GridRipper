#ifndef STLSPHERICALINDEX_HPP
#define STLSPHERICALINDEX_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlGauntIndex.hpp>        // Spherical::IndexPair can construct itself from Gaunt::Index

#include <Gripper/MultipoleTypes.hpp>           // Index is just an alias of the basic template type

// Standard C++ includes
#include <cmath>


namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            class EXPORT Index
            {
            public:

                // Common typedefs

                typedef int32_t  value_type;

                // Common interface

                Index();
                Index(const Index& in);
                Index(Index&& src);
                ~Index();

                Index& operator=(const Index& rhs);

                // Lattice interface

                Index(const value_type& in);

                explicit operator value_type&();

                Index& operator+=(const Index& rhs);
                Index& operator-=(const Index& rhs);
                Index& operator*=(const Index& rhs);
                Index& operator/=(const Index& rhs);

                Index operator++(int rhs);
                Index& operator++();
                Index operator--(int rhs);
                Index& operator--();

                value_type i;
            };

            class EXPORT IndexPair
            {
            public:

                // Common typedefs

                typedef Index  value_type;

                // Common interface

                IndexPair();
                IndexPair(const IndexPair& in);
                IndexPair(IndexPair&& src);
                ~IndexPair();

                // Lattice Interface

                IndexPair(const value_type& l_in, const value_type& m_in);
                IndexPair(const Gaunt::Index& i_in);

                value_type l;
                value_type m;
            };

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename ArithmeticType = std::int8_t>
            using Index = Multipole::SpinWeightedSpherical::Index<ArithmeticType>;

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


///////////////////////////////////////////
// Spherical::Index non-member operators //
///////////////////////////////////////////

// Unary 
EXPORT Multipole::stl::Spherical::Index operator+(const Multipole::stl::Spherical::Index& rhs);
EXPORT Multipole::stl::Spherical::Index operator-(const Multipole::stl::Spherical::Index& rhs);

// Binary
EXPORT Multipole::stl::Spherical::Index operator+(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT Multipole::stl::Spherical::Index operator-(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT Multipole::stl::Spherical::Index operator*(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT Multipole::stl::Spherical::Index operator/(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT Multipole::stl::Spherical::Index operator%(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);

EXPORT bool operator< (const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT bool operator> (const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT bool operator<=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT bool operator>=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT bool operator==(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);
EXPORT bool operator!=(const Multipole::stl::Spherical::Index& lhs, const Multipole::stl::Spherical::Index& rhs);

///////////////////////////////////////////////
// Spherical::IndexPair non-member operators //
///////////////////////////////////////////////

// Binary
EXPORT bool operator==(const Multipole::stl::Spherical::IndexPair& lhs, const Multipole::stl::Spherical::IndexPair& rhs);
EXPORT bool operator!=(const Multipole::stl::Spherical::IndexPair& lhs, const Multipole::stl::Spherical::IndexPair& rhs);

#endif // STLSPHERICALINDEX_HPP