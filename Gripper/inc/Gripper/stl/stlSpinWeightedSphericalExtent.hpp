#pragma once

// Gripper includes
#include <Gripper/stl/stlSpinWeightedSphericalIndex.hpp>    // Multipole::stl::SpinWeightedSpherical::Index

// Standard C++ includes
#include <cassert>              // assert
#include <cstddef>              // std::size_t
#include <initializer_list>     // std::initializer_list
#include <ostream>              // std::ostream


namespace Multipole
{
    namespace stl
    {
        namespace SpinWeightedSpherical
        {
            /// <summary>Class template representing the extent of spherical expansions.</summary>
            ///
            template <std::size_t L_Max, std::size_t S_Max, typename ArithemticType = typename Index<L_Max, S_Max>::value_type>
            class Extent
            {
            public:

                // Lattice type aliases

                typedef Index<L_Max, S_Max, ArithemticType> index_type;
                typedef typename index_type::value_type index_internal_type;

                // Lattice static members

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const index_internal_type s_max = static_cast<index_internal_type>(S_Max);

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                Extent() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                Extent(const Extent&) = default;

                /// <summary>Default move constructor.</summary>
                ///
                Extent(Extent&&) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~Extent() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                Extent& operator=(const Extent&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                Extent& operator=(Extent&&) = default;

                // Lattice interface

                /// <summary>Constructor that validates contents.</summary>
                ///
                Extent(const index_type& initial, const index_type& final) : _initial(initial), _final(final) { assert(_initial <= _final); }

                /// <summary>Initializer list constructor that validates contents.</summary>
                ///
                Extent(std::initializer_list<index_type> init) : _initial(*(init.begin())), _final(*(init.begin() + 1))
                {
                    //static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::SpinWeightedSpherical::Extent must be 2");
                
                    assert(_initial <= _final);
                }

                /// <summary>Returns the origin of the extent.</summary>
                ///
                const index_type& initial() const { return _initial; }

                /// <summary>Returns the end of the extent.</summary>
                ///
                const index_type& final() const { return _final; }

                /// <summary>Tests whether an index is inside the extent.</summary>
                ///
                bool contains(const index_type& index) const { return (index >= _initial) && (index <= _final) ? true : false; }

            private:

                index_type _initial;
                index_type _final;
            };


            /// <summary>STL stream operator overload for formatted console output.</summary>
            ///
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

                os << "[ " << extent.initial() << ", " << extent.final() << " ]";

                return os;
            }

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////////////////
// SpinWeightedSpherical::Extent non-member operators //
////////////////////////////////////////////////////////

// Binary
template <std::size_t L_Max, std::size_t S_Max, typename AT> bool operator==(const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& rhs) { return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial()); }
template <std::size_t L_Max, std::size_t S_Max, typename AT> bool operator!=(const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& lhs, const Multipole::stl::SpinWeightedSpherical::Extent<L_Max, S_Max, AT>& rhs) { return (lhs.final() != rhs.final()) || (lhs.initial() != rhs.initial()); }
