#pragma once

// Gripper includes
#include <Gripper/stl/stlSphericalIndex.hpp>    // Multipole::stl::Spherical::Index, Multipole::stl::SWS::Index

// Standard C++ includes
#include <initializer_list>                     // std::initializer_list


namespace Multipole
{
    namespace stl
    {
        namespace Radial
        {
            template <typename ArithemticType = typename std::int32_t>
            class Extent
            {
            public:

                // Lattice typedefs

                typedef ArithemticType index_type;

                // Common interface

                Extent() = default;
                Extent(const Extent&) = default;
                Extent(Extent&&) = default;
                ~Extent() = default;

                Extent& operator=(const Extent&) = default;
                Extent& operator=(Extent&&) = default;

                // Lattice interface

                Extent(const index_type& initial, const index_type& final) : m_initial(initial), m_final(final) {}
                Extent(std::initializer_list<index_type> init) : m_initial(*(init.begin())), m_final(*(init.begin() + 1)) { static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::Radial::Extent must be 2"); }

                const index_type& initial() const { return m_initial; }
                const index_type& final() const { return m_final; }

                bool contains(const index_type& index) const { return (index >= m_initial) && (index <= m_final) ? true : false; }

            private:

                index_type m_initial;
                index_type m_final;
            };

            template <typename AT>
            std::ostream& operator<<(std::ostream& os, const Extent<AT>& extent)
            {
                //////////////////////////////////////////////////////////////////////////////////////
                //                                                                                  //
                //                              !!!!!! WARNING !!!!!!!                              //
                //                                                                                  //
                //                              STL NON-SENSE DETECTED                              //
                //                                                                                  //
                //////////////////////////////////////////////////////////////////////////////////////
                //
                // The STL does not provide a mechanism to query the open-mode of a stream.
                // Therefor there is no way to distinguish between formatted and binary streams.
                //
                // ASSUMPTION: here we make the assumption that std::ostream& is some derivate of
                //             a formatted console entity.

                os << "[ " << static_cast<int>(extent.initial()) << ", " << static_cast<int>(extent.final()) << " ]";

                return os;
            }

        } // namespace Radial

    } // namespace stl

} // namespace Multipole

/////////////////////////////////////////
// Radial::Extent non-member operators //
/////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::stl::Radial::Extent<AT>& lhs, const Multipole::stl::Radial::Extent<AT>& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
template <typename AT> bool operator!=(const Multipole::stl::Radial::Extent<AT>& lhs, const Multipole::stl::Radial::Extent<AT>& rhs) { return (lhs.final() != rhs.final()) || (lhs.initial() != rhs.initial()); }