#ifndef STLEXPANSIONEXTENT_HPP
#define STLEXPANSIONEXTENT_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlExpansionIndex.hpp>
#include <Gripper/stl/stlRadialExtent.hpp>
#include <Gripper/stl/stlSphericalExtent.hpp>


namespace Multipole
{
    namespace stl
    {
        namespace Expansion
        {
            class EXPORT Extent
            {
            public:

                // Lattice typedefs

                typedef Radial::Extent                      radial_extent_type;
                typedef Spherical::Extent                   spherical_extent_type;
                typedef radial_extent_type::index_type      radial_index_type;
                typedef spherical_extent_type::index_type   spherical_index_type;
                typedef Index                               index_type;

                // Common interface

                Extent();
                Extent(const Extent& in);
                Extent(Extent&& src);
                ~Extent();

                Extent& operator=(const Extent& rhs);

                // Lattice interface

                Extent(const radial_extent_type& radial, const spherical_extent_type& spherical);

                const radial_extent_type& radial() const;
                const spherical_extent_type& spherical() const;

                bool contains(const index_type& index) const;

            private:

                radial_extent_type m_radial;
                spherical_extent_type m_spherical;
            };

        } // namespace Expansion

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Expansion::Extent non-member operators //
////////////////////////////////////////////

// Binary
EXPORT bool operator==(const Multipole::stl::Expansion::Extent& lhs, const Multipole::stl::Expansion::Extent& rhs);
EXPORT bool operator!=(const Multipole::stl::Expansion::Extent& lhs, const Multipole::stl::Expansion::Extent& rhs);

#endif // STLEXPANSIONEXTENT_HPP