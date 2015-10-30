#ifndef STLRADIALINDEX_HPP
#define STLRADIALINDEX_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types

// Standard C++ includes
#include <cmath>


namespace Multipole
{
    namespace stl
    {
        namespace Radial
        {
            class EXPORT Index
            {
            public:

                // Common typedef

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

                value_type r;
            };

        } // namespace Radial

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////
// Radial::Index non-member operators //
////////////////////////////////////////

// Unary 
EXPORT Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& rhs);
namespace Multipole { namespace stl { namespace Radial { EXPORT double pown(const Index& base, int power); } } }

// Binary
EXPORT Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator*(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator/(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator%(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);

EXPORT bool operator< (const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator> (const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator<=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator>=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator==(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator!=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);

#endif // STLRADIALINDEX_HPP