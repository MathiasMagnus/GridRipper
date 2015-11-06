#pragma once

// Standard C++ includes
#include <cstdint>          // std::int8_t
#include <initializer_list> // std::initializer_list
#include <cmath>            // std::round, std::sqrt


namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            template <typename VT>
            class Traits
            {
            public:

                // Common typedefs

                typedef VT          value_type;

                // STL typedefs

                typedef std::size_t size_type;
            };

            template <typename ArithmeticType = std::int8_t>
            struct Index : Traits<ArithmeticType>
            {
                Index() = default;
                Index(const Index&) = default;
                Index(Index&&) = default;
                ~Index() = default;

                Index& operator=(const Index& in) = default;
                Index& operator=(Index&& in) = default;

                template <typename TT>
                Index(const TT l_in, const TT m_in) : l(l_in), m(m_in) {}
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 2, "Initializer-list of Index must have 2 elements: {l, m}");
                    l = *init.begin();
                    m = *(init.begin() + 1);
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

                bool operator<=(const Index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                bool operator>=(const Index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
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

            template <typename AT>
            Index<AT> next(const Index<AT>& index)
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
            class Traits
            {
            public:

                // Common typedefs

                typedef VT          value_type;

                // STL typedefs

                typedef std::size_t size_type;
            };

            template <typename ArithmeticType = std::int8_t>
            struct Index : Traits<ArithmeticType>
            {
                typedef ArithmeticType value_type;

                Index() = default;
                Index(const Index&) = default;
                Index(Index&&) = default;
                ~Index() = default;

                Index& operator=(const Index& in) = default;

                template <typename TT>
                Index(const TT l_in, const TT m_in, const TT s_in) : l(l_in), m(m_in), s(s_in) {}
                Index(std::initializer_list<value_type> init)
                {
                    //static_assert(init.size() == 3, "Initializer-list of Multipole::stl::SpinWeightedSpherical::Index must have 3 elements: {l, m, s}");
                    l = *init.begin();
                    m = *(init.begin() + 1);
                    s = *(init.begin() + 2);
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

                bool operator<=(const Index& rhs) const
                {
                    return (*this < rhs) || (*this == rhs);
                }

                bool operator>=(const Index& rhs) const
                {
                    return (*this > rhs) || (*this == rhs);
                }

                static Index convert(const std::size_t& i, const value_type& max_s)
                {
                    value_type i_corr = static_cast<value_type>(i / (2 * max_s + 1)); // truncation on division is intensional (omit std::floor and double up-cast)
                    value_type l = static_cast<value_type>(std::round((-1. + std::sqrt(1. + 4 * i_corr)) / 2.));
                    value_type m = i_corr - l*(l + 1);
                    value_type s = (i % (2 * max_s + 1)) - max_s;

                    return Index{ l, m, s };
                }

                static size_type convert(const Index& i, const value_type& max_s)
                {
                    return (i.l * (i.l + 1) + i.m) * (2 * max_s + 1) + i.s + max_s;
                }

                value_type l;
                value_type m;
                value_type s;
            };

            template <typename AT>
            Index<AT> next(const Index<AT>& index, const AT& max_l, const AT& max_s)
            {
                bool rewind_s_and_step_m = index.s == max_s;
                bool rewind_m_and_step_l = rewind_s_and_step_m && (index.m == index.l);

                return Index<AT>
                {
                    rewind_m_and_step_l ? index.l + 1 : index.l,
                        rewind_s_and_step_m ? (rewind_m_and_step_l ? -(index.l + 1) : index.m + 1) : index.m,
                        rewind_s_and_step_m ? -max_s : index.s + 1
                };
            }

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole