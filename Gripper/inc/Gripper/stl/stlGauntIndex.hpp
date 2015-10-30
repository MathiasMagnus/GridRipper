#ifndef STLGAUNTINDEX_HPP
#define STLGAUNTINDEX_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlSphericalIndex.hpp>    // Gaunt::Index can construct itself from Spherical::IndexPair
#include <Gripper/MultipoleTypes.hpp>           // Instantiate Gaunt::stl::index from generic types


namespace Multipole
{
    namespace stl
    {
        namespace Gaunt
        {
            class EXPORT Index
            {
            public:

                // Lattice typedefs

                typedef uint32_t  value_type;

                // Common interface

                Index();
                Index(const Index& in);
                Index(Index&& src);
                ~Index();

                Index& operator=(const Index& rhs);

                // Lattice interface

                Index(const value_type& i_in);
                Index(const Multipole::stl::Spherical::IndexPair& lm_in);
                Index(const Multipole::Spherical::IndexPair& lm_in);

                explicit operator const value_type&();

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

        } // namespace Gaunt

    } // namespace stl

} // namespace Multipole


///////////////////////////////////////
// Gaunt::Index non-member operators //
///////////////////////////////////////

// Binary
EXPORT Multipole::stl::Gaunt::Index operator+(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT Multipole::stl::Gaunt::Index operator-(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT Multipole::stl::Gaunt::Index operator*(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT Multipole::stl::Gaunt::Index operator/(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT Multipole::stl::Gaunt::Index operator%(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);

EXPORT bool operator< (const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT bool operator> (const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT bool operator<=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT bool operator>=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT bool operator==(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);
EXPORT bool operator!=(const Multipole::stl::Gaunt::Index& lhs, const Multipole::stl::Gaunt::Index& rhs);

#endif // STLGAUNTINDEX_HPP