#pragma once

// Standard C++ includes
#include <cstdlib>          // std::size_t
#include <cstdint>          // std::int8_t
#include <cassert>          // assert
#include <initializer_list> // std::initializer_list
#include <cmath>            // std::round, std::sqrt


namespace Multipole
{
    namespace stl
    {
        namespace SpinWeightedSpherical
        {
            /// <summary>Traits class consisting of type aliases describing spin-weighted spherical expansion indeces.</summary>
            ///
            template <typename VT>
            struct IndexTraits
            {
                // Common type aliases

                using value_type = VT;

                // STL type aliases

                using size_type = std::size_t;
            };


            /// <summary>Class template representing spin-weighted spherical expansion indeces.</summary>
            ///
            template <std::size_t L_Max, std::size_t S_Max, typename ArithmeticType = std::int8_t>
            struct Index : IndexTraits<ArithmeticType>
            {
                // Common type aliases

                using typename IndexTraits<ArithmeticType>::value_type;

                // STL type aliases

                using typename IndexTraits<ArithmeticType>::size_type;

                // Lattice static members

                static const value_type l_max = static_cast<value_type>(L_Max);
                static const value_type s_max = static_cast<value_type>(S_Max);

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
                Index& operator=(const Index&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                Index& operator=(Index&&) = default;

                /// <summary>Templated constructor that validates contents.</summary>
                ///
                template <typename TT>
                Index(const TT l_in, const TT m_in, const TT s_in) : l(l_in), m(m_in), s(s_in) { assert(l <= l_max && std::abs(m) <= l_max && std::abs(s) <= std::min(l, s_max)); }

                /// <summary>Initializer list constructor that validates contents.</summary>
                ///
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 3, "Initializer-list of Multipole::stl::SpinWeightedSpherical::Index must have 3 elements: {l, m, s}");
                    assert(init.size() == 3);

                    l = *init.begin();
                    m = *(init.begin() + 1);
                    s = *(init.begin() + 2);

                    assert(l <= l_max && std::abs(m) <= l_max && std::abs(s) <= std::min(l, s_max));
                }

                // Member operators

                /// <summary>Strictly weak ordering comparator.</summary>
                ///
                bool operator<(const Index& rhs) const
                {
                    if (l < rhs.l) return true;
                    if (rhs.l < l) return false;

                    if (m < rhs.m) return true;
                    if (rhs.m < m) return false;

                    if (s < rhs.s) return true;
                    return false;
                }

                /// <summary>Strictly weak ordering comparator.</summary>
                ///
                bool operator>(const Index& rhs) const
                {
                    if (l > rhs.l) return true;
                    if (rhs.l > l) return false;

                    if (m > rhs.m) return true;
                    if (rhs.m > m) return false;

                    if (s > rhs.s) return true;
                    return false;
                }

                /// <summary>Equality test operator.</summary>
                ///
                bool operator==(const Index& rhs) const
                {
                    return (l == rhs.l) && (m == rhs.m) && (s == rhs.s);
                }

                /// <summary>Negated equality test operator.</summary>
                ///
                bool operator!=(const Index& rhs) const
                {
                    return (l != rhs.l) || (m != rhs.m) || (s != rhs.s);
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
                    bool rewind_s_and_step_m = s == std::min(s_max, l);
                    bool rewind_m_and_step_l = rewind_s_and_step_m && (m == l);

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_s_and_step_m ? (rewind_m_and_step_l ? -l : m + 1) : m;
                    s = rewind_s_and_step_m ? -std::min(s_max, l) : s + 1;

                    return *this;
                }

                /// <summary>Prefix decrement operator.</summary>
                ///
                Index& operator--()
                {
                    bool rewind_s_and_step_m = s == std::max(-l, -s_max);
                    bool rewind_m_and_step_l = rewind_s_and_step_m && (-m == l);

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_s_and_step_m ? (rewind_m_and_step_l ? -l : m - 1) : m;
                    s = rewind_s_and_step_m ? std::min(s_max, l) : s - 1;

                    return *this;
                }

                /// <summary>Postfix increment operator.</summary>
                ///
                Index operator++(int)
                {
                    Index result = *this;

                    this->operator++();

                    return result;
                }

                /// <summary>Postfix decrement operator.</summary>
                ///
                Index operator--(int)
                {
                    Index result = *this;

                    this->operator--();

                    return result;
                }

                value_type l;
                value_type m;
                value_type s;
            };

            /// <summary>STL stream operator overload for formatted console output.</summary>
            ///
            template <std::size_t L, std::size_t S, typename AT>
            std::ostream& operator<<(std::ostream& os, const Index<L, S, AT>& index)
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

                os << "{ " << static_cast<int>(index.l) << ", " << static_cast<int>(index.m) << ", " << static_cast<int>(index.s) << " }";

                return os;
            }

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole
