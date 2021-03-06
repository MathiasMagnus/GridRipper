#pragma once

// Standard C++ includes
#include <cassert>              // assert
#include <cstdlib>              // std::size_t
#include <cstdint>              // std::int8_t
#include <initializer_list>     // std::initializer_list
#include <cmath>                // std::round, std::sqrt
#include <ostream>              // std::ostream


namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            /// <summary>Traits class consisting of type aliases describing spherical expansion indeces.</summary>
            ///
            template <typename VT>
            struct IndexTraits
            {
                // Common type aliases

                using value_type = VT;

                // STL type aliases

                using size_type = std::size_t;
            };


            /// <summary>Class template representing spherical expansion indeces.</summary>
            ///
            template <std::size_t L_Max, typename ArithmeticType = std::int8_t>
            struct Index : IndexTraits<ArithmeticType>
            {
                // Common type aliases

                using typename IndexTraits<ArithmeticType>::value_type;

                // STL type aliases

                using typename IndexTraits<ArithmeticType>::size_type;

                // Lattice static members

                static const value_type l_max = static_cast<value_type>(L_Max);

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                Index() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                Index(const Index&) = default;

                /// <summary>Default move constructor.</summary>
                ///
                Index(Index&&) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~Index() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                Index& operator=(const Index& in) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                Index& operator=(Index&& in) = default;

                /// <summary>Templated constructor that validates contents.</summary>
                ///
                template <typename TT>
                Index(const TT l_in, const TT m_in) : l(l_in), m(m_in) { assert(l <= l_max && std::abs(m) <= l_max); }

                /// <summary>Initializer list constructor that validates contents.</summary>
                ///
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 2, "Initializer-list of Index must have 2 elements: {l, m}");
                    l = *init.begin();
                    m = *(init.begin() + 1);

                    assert(l <= l_max && std::abs(m) <= l_max);
                }

                // Member operators

                /// <summary>Strictly weak ordering comparator.</summary>
                ///
                bool operator<(const Index& rhs) const
                {
                    if (l < rhs.l) return true;
                    if (rhs.l < l) return false;

                    if (m < rhs.m) return true;
                    return false;
                }

                /// <summary>Strictly weak ordering comparator.</summary>
                ///
                bool operator>(const Index& rhs) const
                {
                    if (l > rhs.l) return true;
                    if (rhs.l > l) return false;

                    if (m > rhs.m) return true;
                    return false;
                }

                /// <summary>Equality test operator.</summary>
                ///
                bool operator==(const Index& rhs) const
                {
                    return (l == rhs.l) && (m == rhs.m);
                }

                /// <summary>Negated equality test operator.</summary>
                ///
                bool operator!=(const Index& rhs) const
                {
                    return (l != rhs.l) || (m != rhs.m);
                }

                /// <summary>Less than or equal to comparator.</summary>
                ///
                bool operator<=(const Index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                /// <summary>Freater than or equal to comparator.</summary>
                ///
                bool operator>=(const Index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
                }

                /// <summary>Prefix increment operator.</summary>
                ///
                Index& operator++()
                {
                    bool rewind_m_and_step_l = m == l;

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_m_and_step_l ? -l : m + 1;

                    return *this;
                }

                /// <summary>Prefix decrement operator.</summary>
                ///
                Index& operator--()
                {
                    bool rewind_m_and_step_l = -m == l;

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_m_and_step_l ? l : m - 1;

                    return *this;
                }

                /// <summary>Postfix increment operator.</summary>
                ///
                Index operator++(int)
                {
                    Index result = *this;

                    bool rewind_m_and_step_l = m == l;

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_m_and_step_l ? -l : m + 1;

                    return result;
                }

                /// <summary>Postfix decrement operator.</summary>
                ///
                Index operator--(int)
                {
                    Index result = *this;

                    bool rewind_m_and_step_l = -m == l;

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_m_and_step_l ? l : m - 1;

                    return result;
                }

                value_type l;
                value_type m;
            };

            /// <summary>STL stream operator overload for formatted console output.</summary>
            ///
            template <std::size_t L, typename AT>
            std::ostream& operator<<(std::ostream& os, const Index<L, AT>& index)
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

                os << "{ " << static_cast<int>(index.l) << ", " << static_cast<int>(index.m) << " }";

                return os;
            }

        } // namespace Spherical

    } // namespace stl

} // namespace Multipole
