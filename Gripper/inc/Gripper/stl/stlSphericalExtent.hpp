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
            template <std::size_t L_Max, typename ArithemticType = typename Index<L_Max>::value_type>
            class Extent
            {
            public:

                // Lattice typedefs

                typedef Index<L_Max, ArithemticType> index_type;
                typedef typename index_type::value_type index_internal_type;

                // Common interface

                Extent() = default;
                Extent(const Extent&) = default;
                Extent(Extent&&) = default;
                ~Extent() = default;

                Extent& operator=(const Extent&) = default;
                Extent& operator=(Extent&&) = default;

                static index_internal_type l_max = static_cast<index_internal_type>(L_Max);

                // Lattice interface

                Extent(const index_type& initial, const index_type& final) : m_initial(initial), m_final(final) {}
                Extent(std::initializer_list<index_type> init) : m_initial(*(init.begin())), m_final(*(init.begin() + 1)) { static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::Spherical::Extent must be 2"); }

                const index_type& initial() const { return m_initial; }
                const index_type& final() const { return m_final; }

                bool contains(const index_type& index) const { return (index >= m_initial) && (index <= m_final) ? true : false; }

            private:

                index_type m_initial;
                index_type m_final;
            };

            template <std::size_t L_Max, typename AT>
            std::ostream& operator<<(std::ostream& os, const Extent<L_Max, AT>& extent)
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

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <std::size_t L_Max, std::size_t S_Max, typename ArithemticType = typename Index<L_Max, S_Max>::value_type>
            class Extent
            {
            public:

                // Lattice typedefs

                typedef Index<L_Max, S_Max, ArithemticType> index_type;
                typedef typename index_type::value_type index_internal_type;

                // Common interface

                Extent() = default;
                Extent(const Extent&) = default;
                Extent(Extent&&) = default;
                ~Extent() = default;

                Extent& operator=(const Extent&) = default;
                Extent& operator=(Extent&&) = default;

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const index_internal_type s_max = static_cast<index_internal_type>(S_Max);

                // Lattice interface

                Extent(const index_type& initial, const index_type& final) : m_initial(initial), m_final(final) {}
                Extent(std::initializer_list<index_type> init) : m_initial(*(init.begin())), m_final(*(init.begin() + 1)) { static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::SpinWeightedSpherical::Extent must be 2"); }

                const index_type& initial() const { return m_initial; }
                const index_type& final() const { return m_final; }

                bool contains(const index_type& index) const { return (index >= m_initial) && (index <= m_final) ? true : false; }

            private:

                index_type m_initial;
                index_type m_final;
            };

            template <std::size_t L_Max, std::size_t S_Max, typename AT>
            std::ostream& operator<<(std::ostream& os, const Extent<L_Max, S_Max, AT>& extent)
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

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Spherical::Extent non-member operators //
////////////////////////////////////////////

// Binary
template <std::size_t L_Max, typename AT> bool operator==(const Multipole::stl::Spherical::Extent<L_Max, AT>& lhs, const Multipole::stl::Spherical::Extent<L_Max, AT>& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
template <std::size_t L_Max, typename AT> bool operator!=(const Multipole::stl::Spherical::Extent<L_Max, AT>& lhs, const Multipole::stl::Spherical::Extent<L_Max, AT>& rhs) { return (lhs.final() != rhs.final()) || (lhs.initial() != rhs.initial()); }

////////////////////////////////////////////////////////
// SpinWeightedSpherical::Extent non-member operators //
////////////////////////////////////////////////////////

// Binary
template <std::size_t L_Max, std::size_t S_Max, typename AT> bool operator==(const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
template <std::size_t L_Max, std::size_t S_Max, typename AT> bool operator!=(const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& rhs) { return (lhs.final() != rhs.final()) || (lhs.initial() != rhs.initial()); }