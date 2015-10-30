#ifndef STLGAUNTCOEFFICIENT_HPP
#define STLGAUNTCOEFFICIENT_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlGauntIndex.hpp>
#include <Gripper/MultipoleTypes.hpp>           // Up-conversion from generic Multipole::Gaunt::Coefficient


namespace Multipole
{
    namespace stl
    {
        namespace Gaunt
        {
            class EXPORT Coefficient
            {
            public:

                // Lattice typedefs

                typedef double  value_type;
                typedef Index   index_type;

                // Common interface

                Coefficient();
                Coefficient(const Coefficient& in);
                Coefficient(Coefficient&& src);
                ~Coefficient();

                Coefficient& operator=(const Coefficient& rhs)/* = default*/;

                // Lattice interface

                Coefficient(const Multipole::Gaunt::Coefficient& in);
                Coefficient(const index_type& i1_in, const index_type& i2_in, const index_type& i3_in, const value_type& value_in);
                Coefficient(const Spherical::IndexPair& l1m1_in, const Spherical::IndexPair& l2m2_in, const Spherical::IndexPair& l3m3_in, const value_type& value_in);

                index_type i1;
                index_type i2;
                index_type i3;
                value_type value;
            };

        } // namespace Gaunt

    } // namespace stl

} // namespace Multipole


/////////////////////////////////////////////
// Gaunt::Coefficient non-member operators //
/////////////////////////////////////////////

EXPORT bool operator< (const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);
EXPORT bool operator> (const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);
EXPORT bool operator<=(const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);
EXPORT bool operator>=(const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);
EXPORT bool operator==(const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);
EXPORT bool operator!=(const Multipole::stl::Gaunt::Coefficient& lhs, const Multipole::stl::Gaunt::Coefficient& rhs);

#endif // STLGAUNTCOEFFICIENT_HPP