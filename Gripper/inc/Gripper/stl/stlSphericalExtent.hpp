#pragma once

// Gripper includes
#include <Gripper/stl/stlSphericalIndex.hpp>    // Multipole::stl::Spherical::Index, Multipole::stl::SpinWeightedSpherical::Index

// Standard C++ includes
#include <initializer_list>                     // std::initializer_list


namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            template <typename ArithemticType = Index<>::value_type>
            class Extent
            {
            public:

                // Lattice typedefs

                typedef Index<ArithemticType> index_type;

                // Common interface

                Extent() = default;
                Extent(const Extent&) = default;
                Extent(Extent&&) = default;
                ~Extent() = default;

                Extent& operator=(const Extent&) = default;
                Extent& operator=(Extent&&) = default;

                // Lattice interface

                Extent(index_type&& initial, index_type&& final) : m_initial(initial), m_final(final) {}
                Extent(std::initializer_list<index_type> init) : m_initial(*(init.begin())), m_final(*(init.begin() + 1)) { static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::Spherical::Extent must be 2"); }

                const index_type& initial() const { return m_initial; }
                const index_type& final() const { return m_final; }

                bool contains(const index_type& index) const { return (index >= m_initial) && (index <= m_final) ? true : false; }

            private:

                index_type m_initial;
                index_type m_final;
            };

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename ArithemticType = Index<>::value_type>
            class Extent
            {
            public:

                // Lattice typedefs

                typedef Index<ArithemticType> index_type;

                // Common interface

                Extent() = default;
                Extent(const Extent&) = default;
                Extent(Extent&&) = default;
                ~Extent() = default;

                Extent& operator=(const Extent&) = default;
                Extent& operator=(Extent&&) = default;

                // Lattice interface

                Extent(index_type&& initial, index_type&& final) : m_initial(initial), m_final(final) {}
                Extent(std::initializer_list<index_type> init) : m_initial(*(init.begin())), m_final(*(init.begin() + 1)) { static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::SpinWeightedSpherical::Extent must be 2"); }

                const index_type& initial() const { return m_initial; }
                const index_type& final() const { return m_final; }

                bool contains(const index_type& index) const { return (index >= m_initial) && (index <= m_final) ? true : false; }

            private:

                index_type m_initial;
                index_type m_final;
            };

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Spherical::Extent non-member operators //
////////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::stl::Spherical::Extent<AT>& lhs, const Multipole::stl::Spherical::Extent<AT>& rhs) { return (lhs.final == rhs.final) && (lhs.initial == rhs.initial); }
template <typename AT> bool operator!=(const Multipole::stl::Spherical::Extent<AT>& lhs, const Multipole::stl::Spherical::Extent<AT>& rhs) { return (lhs.final != rhs.final) || (lhs.initial != rhs.initial); }

////////////////////////////////////////////////////////
// SpinWeightedSpherical::Extent non-member operators //
////////////////////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::stl::SpinWeightedSpherical::Extent<AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<AT>& rhs) { return (lhs.final == rhs.final) && (lhs.initial == rhs.initial); }
template <typename AT> bool operator!=(const Multipole::stl::SpinWeightedSpherical::Extent<AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<AT>& rhs) { return (lhs.final != rhs.final) || (lhs.initial != rhs.initial); }