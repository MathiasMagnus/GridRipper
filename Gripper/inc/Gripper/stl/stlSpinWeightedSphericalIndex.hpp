#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>

// Standard C++ includes
#include <cassert>              // assert
#include <cstdlib>              // std::size_t
#include <cstdint>              // std::int8_t
#include <initializer_list>     // std::initializer_list
#include <cmath>                // std::abs
#include <ostream>              // std::ostream


namespace math
{
    namespace sws
    {
        /// <summary>Traits class consisting of type aliases describing spin-weighted spherical expansion indicies.</summary>
        ///
        template <typename VT>
        struct index_traits
        {
            // Common type aliases

            using value_type = VT;
        };


        /// <summary>Class template representing spin-weighted spherical expansion indicies.</summary>
        ///
        template <typename ArithmeticType = std::int8_t>
        struct index : index_traits<ArithmeticType>
        {
            // Common type aliases

            using typename index_traits<ArithmeticType>::value_type;

            // Constructors / Destructors / Assignment operators

            /// <summary>Default constructor.</summary>
            /// <remarks>Default constructed objects are in an invalid state.</remarks>
            ///
            index() = default;

            /// <summary>Default copy constructor.</summary>
            ///
            index(const index&) = default;

            /// <summary>Default move constructor.</summary>
            ///
            index(index&&) = default;

            /// <summary>Default destructor.</summary>
            ///
            ~index() = default;

            /// <summary>Default copy assignment operator.</summary>
            ///
            index& operator=(const index&) = default;

            /// <summary>Default move assignment operator.</summary>
            ///
            index& operator=(index&&) = default;

            /// <summary>Templated constructor that validates contents.</summary>
            ///
            template <typename TT>
            index(const TT l_in, const TT m_in, const TT s_in) : l(l_in), m(m_in), s(s_in)
            {
                // Invalid indices will be instantiated
                //
                //assert(s < 0 ? -s <= l : s <= l); // SpinWeightedSphericalindex must have |S| <= L_Max.
                //assert(l <= l_max && std::abs(m) <= l_max);
            }

            /// <summary>Initializer list constructor that validates contents.</summary>
            ///
            index(std::initializer_list<value_type> init)
            {
                //static_assert(init.size() == 3, "Initializer-list of Multipole::stl::SWS::index must have 3 elements: {l, m, s}");
                assert(init.size() == 3);

                l = *(init.begin() + 0);
                m = *(init.begin() + 1);
                s = *(init.begin() + 2);
            }

            // Member operators
            /*
                /// <summary>Strictly weak ordering comparator.</summary>
                ///
                bool operator<(const index& rhs) const
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
                bool operator>(const index& rhs) const
                {
                    if (l > rhs.l) return true;
                    if (rhs.l > l) return false;

                    if (m > rhs.m) return true;
                    if (rhs.m > m) return false;

                    if (s > rhs.s) return true;
                    return false;
                }
                */
            /// <summary>Equality test operator.</summary>
            ///
            bool operator==(const index& rhs) const
            {
                return (l == rhs.l) && (m == rhs.m) && (s == rhs.s);
            }

            /// <summary>Negated equality test operator.</summary>
            ///
            bool operator!=(const index& rhs) const
            {
                return (l != rhs.l) || (m != rhs.m) || (s != rhs.s);
            }
            /*
                /// <summary>Less than or equal to comparator.</summary>
                ///
                bool operator<=(const index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                /// <summary>Greater than or equal to comparator.</summary>
                ///
                bool operator>=(const index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
                }
                
                /// <summary>Prefix increment operator.</summary>
                ///
                index& operator++()
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
                index& operator--()
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
                index operator++(int)
                {
                    index result = *this;

                    this->operator++();

                    return result;
                }

                /// <summary>Postfix decrement operator.</summary>
                ///
                index operator--(int)
                {
                    index result = *this;

                    this->operator--();

                    return result;
                }
                */
            value_type l;
            value_type m;
            value_type s;
        };



        /// <summary>STL stream operator overload for formatted console output.</summary>
        ///
        template <typename AT>
        std::ostream& operator<<(std::ostream& os, const index<AT>& index)
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

            os <<
                "{ " << static_cast<int>(index.l) <<
                ", " << static_cast<int>(index.m) <<
                ", " << static_cast<int>(index.s) <<
                " }";

            return os;
        }

    } // namespace sws

} // namespace math
