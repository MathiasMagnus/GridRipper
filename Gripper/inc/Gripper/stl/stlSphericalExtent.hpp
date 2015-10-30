#ifndef STLSPHERICALEXTENT_HPP
#define STLSPHERICALEXTENT_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlGauntIndex.hpp>        // Spherical::Extent uses Gaunt::Index as representation
#include <Gripper/stl/stlLogger.hpp>            // Trace logging


namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            class EXPORT Extent
            {
            public:

                // Lattice typedefs

                typedef Gaunt::Index  index_type;

                // Common interface

                Extent();
                Extent(const Extent& in);
                Extent(Extent&& src);
                ~Extent();

                Extent& operator=(const Extent& rhs);

                // Lattice interface

                Extent(const index_type& initial, const index_type& final);
                Extent(const Index& initial, const Index& final);

                const index_type& initial() const;
                const index_type& final() const;

                bool contains(const index_type& index) const;

            private:

                index_type m_initial;
                index_type m_final;
            };

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename ArithmeticType>
            struct Extent
            {
                // Lattice typedefs

                using index_type = Index<ArithmeticType>;

                // Common interface

                Extent() = default;
                Extent(const Extent& in) = default;
                Extent(Extent&& src) = default;
                ~Extent() = default;

                Extent& operator=(const Extent& rhs) { this->initial = rhs.initial; this->final = rhs.final; return *this; }

                // Lattice interface

                Extent(const index_type& initial_in, const index_type& final_in) : initial(initial_in), final(final_in) {}

                bool contains(const index_type& index) const { return (index >= m_initial) && (index < m_final) ? true : false; }

                index_type initial;
                index_type final;
            };

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Spherical::Extent non-member operators //
////////////////////////////////////////////

// Binary
EXPORT bool operator==(const Multipole::stl::Spherical::Extent& lhs, const Multipole::stl::Spherical::Extent& rhs);
EXPORT bool operator!=(const Multipole::stl::Spherical::Extent& lhs, const Multipole::stl::Spherical::Extent& rhs);

////////////////////////////////////////////////////////
// SpinWeightedSpherical::Extent non-member operators //
////////////////////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::stl::SpinWeightedSpherical::Extent<AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<AT>& rhs) { return (lhs.final == rhs.final) && (lhs.initial == rhs.initial); }
template <typename AT> bool operator!=(const Multipole::stl::SpinWeightedSpherical::Extent<AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<AT>& rhs) { return (lhs.final != rhs.final) || (lhs.initial != rhs.initial); }

#endif // STLSPHERICALEXTENT_HPP