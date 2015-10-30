#ifndef STLRADIALEXTENT_HPP
#define STLRADIALEXTENT_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlRadialIndex.hpp>
#include <Gripper/stl/stlLogger.hpp>            // Trace logging


namespace Multipole
{
    namespace stl
    {
        namespace Radial
        {
            class EXPORT Extent
            {
            public:

                // Lattice typedefs

                typedef Index  index_type;

                // Common interface

                Extent();
                Extent(const Extent& in);
                Extent(Extent&& src);
                ~Extent();

                Extent& operator=(const Extent& rhs);

                // Lattice interface

                Extent(const index_type& initial, const index_type& final);

                const index_type& initial() const;
                const index_type& final() const;

                bool contains(const index_type& index) const;

            private:

                index_type m_initial;
                index_type m_final;
            };

        } // namespace Radial

    } // namespace stl

} // namespace Multipole


/////////////////////////////////////////
// Radial::Extent non-member operators //
/////////////////////////////////////////

// Binary
EXPORT bool operator==(const Multipole::stl::Radial::Extent& lhs, const Multipole::stl::Radial::Extent& rhs);
EXPORT bool operator!=(const Multipole::stl::Radial::Extent& lhs, const Multipole::stl::Radial::Extent& rhs);

#endif // STLRADIALEXTENT_HPP