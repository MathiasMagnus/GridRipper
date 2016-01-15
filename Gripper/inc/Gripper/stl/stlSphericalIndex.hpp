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
        namespace Spherical
        {
            template <typename VT>
            class IndexTraits
            {
            public:

                // Common typedefs

                typedef VT          value_type;

                // STL typedefs

                typedef std::size_t size_type;
            };

            template <std::size_t L_Max, typename ArithmeticType = std::int8_t>
            struct Index : IndexTraits<ArithmeticType>
            {
                // Common typedefs

                using typename IndexTraits<ArithmeticType>::value_type;
                using typename IndexTraits<ArithmeticType>::size_type;

                Index() = default;
                Index(const Index&) = default;
                Index(Index&&) = default;
                ~Index() = default;

                Index& operator=(const Index& in) = default;
                Index& operator=(Index&& in) = default;

                static value_type l_max = static_cast<value_type>(L_Max);

                template <typename TT>
                Index(const TT l_in, const TT m_in) : l(l_in), m(m_in) { assert(l <= l_max && std::abs(m) <= l_max); }
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 2, "Initializer-list of Index must have 2 elements: {l, m}");
                    l = *init.begin();
                    m = *(init.begin() + 1);

                    assert(l <= l_max && std::abs(m) <= l_max);
                }

                bool operator<(const Index& rhs) const
                {
                    if (l < rhs.l) return true;
                    if (rhs.l < l) return false;

                    if (m < rhs.m) return true;
                    return false;
                }

                bool operator>(const Index& rhs) const
                {
                    if (l > rhs.l) return true;
                    if (rhs.l > l) return false;

                    if (m > rhs.m) return true;
                    return false;
                }

                bool operator==(const Index& rhs) const
                {
                    return (l == rhs.l) && (m == rhs.m);
                }

                bool operator!=(const Index& rhs) const
                {
                    return (l != rhs.l) || (m != rhs.m);
                }

                bool operator<=(const Index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                bool operator>=(const Index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
                }

                Index& operator++()
                {
                    bool rewind_m_and_step_l = m == l;

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_m_and_step_l ? -l : m + 1;

                    return *this;
                }

                Index& operator--()
                {
                    bool rewind_m_and_step_l = -m == l;

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_m_and_step_l ? l : m - 1;

                    return *this;
                }

                Index operator++(int)
                {
                    Index result = *this;

                    bool rewind_m_and_step_l = m == l;

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_m_and_step_l ? -l : m + 1;

                    return result;
                }

                Index operator--(int)
                {
                    Index result = *this;

                    bool rewind_m_and_step_l = -m == l;

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_m_and_step_l ? l : m - 1;

                    return result;
                }

                static Index convert(const std::size_t& i)
                {
                    value_type l = static_cast<value_type>(std::round((-1. + std::sqrt(1. + 4 * i)) / 2.));
                    value_type m = i - l*(l + 1);

                    return Index{ l, m };
                }

                static size_type convert(const Index& i)
                {
                    return i.l * (i.l + 1) + i.m;
                }

                value_type l;
                value_type m;
            };

            template <std::size_t L, typename  AT>
            Index<L, AT> next(const Index<L, AT>& index)
            {
                bool rewind_m_and_step_l = index.m == index.l;

                return Index<AT>
                {
                    rewind_m_and_step_l ? index.l + 1 : index.l,
                        rewind_m_and_step_l ? -(index.l + 1) : index.m + 1
                };
            }

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename VT>
            class IndexTraits
            {
            public:

                // Common typedefs

                typedef VT          value_type;

                // STL typedefs

                typedef std::size_t size_type;
            };

            template <std::size_t L_Max, std::size_t S_Max, typename ArithmeticType = std::int8_t>
            struct Index : IndexTraits<ArithmeticType>
            {
                // Common typedefs

                using typename IndexTraits<ArithmeticType>::value_type;
                using typename IndexTraits<ArithmeticType>::size_type;

                Index() = default;
                Index(const Index&) = default;
                Index(Index&&) = default;
                ~Index() = default;

                Index& operator=(const Index&) = default;
                Index& operator=(Index&&) = default;

                static const value_type l_max = static_cast<value_type>(L_Max);
                static const value_type s_max = static_cast<value_type>(S_Max);

                template <typename TT>
                Index(const TT l_in, const TT m_in, const TT s_in) : l(l_in), m(m_in), s(s_in) { assert(l <= l_max && std::abs(m) <= l_max && std::abs(s) <= std::min(l, s_max)); }
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 3, "Initializer-list of Multipole::stl::SpinWeightedSpherical::Index must have 3 elements: {l, m, s}");
                    assert(init.size() == 3);

                    l = *init.begin();
                    m = *(init.begin() + 1);
                    s = *(init.begin() + 2);

                    assert(l <= l_max && std::abs(m) <= l_max && std::abs(s) <= std::min(l, s_max));
                }

                bool operator<(const Index& rhs) const
                {
                    if (l < rhs.l) return true;
                    if (rhs.l < l) return false;

                    if (m < rhs.m) return true;
                    if (rhs.m < m) return false;

                    if (s < rhs.s) return true;
                    return false;
                }

                bool operator>(const Index& rhs) const
                {
                    if (l > rhs.l) return true;
                    if (rhs.l > l) return false;

                    if (m > rhs.m) return true;
                    if (rhs.m > m) return false;

                    if (s > rhs.s) return true;
                    return false;
                }

                bool operator==(const Index& rhs) const
                {
                    return (l == rhs.l) && (m == rhs.m) && (s == rhs.s);
                }

                bool operator!=(const Index& rhs) const
                {
                    return (l != rhs.l) || (m != rhs.m) || (s != rhs.s);
                }

                bool operator<=(const Index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                bool operator>=(const Index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
                }

                Index& operator++()
                {
                    bool rewind_s_and_step_m = s == std::min(s_max, l);
                    bool rewind_m_and_step_l = rewind_s_and_step_m && (m == l);

                    l = rewind_m_and_step_l ? l + 1 : l;
                    m = rewind_s_and_step_m ? (rewind_m_and_step_l ? -l : m + 1) : m;
                    s = rewind_s_and_step_m ? -std::min(s_max, l) : s + 1;

                    return *this;
                }

                Index& operator--()
                {
                    bool rewind_s_and_step_m = s == std::max(-l, -s_max);
                    bool rewind_m_and_step_l = rewind_s_and_step_m && (-m == l);

                    l = rewind_m_and_step_l ? l - 1 : l;
                    m = rewind_s_and_step_m ? (rewind_m_and_step_l ? -l : m - 1) : m;
                    s = rewind_s_and_step_m ? std::min(s_max, l) : s - 1;

                    return *this;
                }

                Index operator++(int)
                {
                    Index result = *this;

                    this->operator++();

                    return result;
                }

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