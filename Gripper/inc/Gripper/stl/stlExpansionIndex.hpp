#ifndef STLEXPANSIONINDEX_HPP
#define STLEXPANSIONINDEX_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlRadialIndex.hpp>
#include <Gripper/stl/stlSphericalIndex.hpp>


namespace Multipole
{
    namespace stl
    {
        namespace Expansion
        {
            class EXPORT Index
            {
            public:

                // Lattice typedefs

                typedef Radial::Index   radial_index_type;
                typedef Gaunt::Index    spherical_index_type;

                // Common interface

                Index();
                Index(const Index& in);
                Index(Index&& src);
                ~Index();

                // Lattice interface

                Index(const radial_index_type& radial, const spherical_index_type& spherical);

                const radial_index_type& radial() const;
                const spherical_index_type& spherical() const;

            private:

                radial_index_type m_radial;
                spherical_index_type m_spherical;
            };

        } // Expansion

    } // namespace stl

} // namespace Multipole


/////////////////////////////////////////
// Radial::Extent non-member operators //
/////////////////////////////////////////

// Binary
EXPORT bool operator==(const Multipole::stl::Expansion::Index& lhs, const Multipole::stl::Expansion::Index& rhs);
EXPORT bool operator!=(const Multipole::stl::Expansion::Index& lhs, const Multipole::stl::Expansion::Index& rhs);

#endif // STLEXPANSIONINDEX_HPP