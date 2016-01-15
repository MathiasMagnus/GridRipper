#pragma once

// Standard C++ includes
#include <cstdlib>          // std::size_t
#include <cstdint>          // std::int8_t


namespace Multipole
{
    namespace stl
    {
        namespace Radial
        {/*
            template <typename VT>
            class IndexTraits
            {
            public:

                // Common typedefs

                typedef VT          value_type;

                // STL typedefs

                typedef std::size_t size_type;
            };

            template <typename ArithmeticType = std::int8_t>
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

                template <typename TT>
                Index(const TT r_in) : r(r_in) {}

                bool operator<(const Index& rhs) const { return r < rhs.r ? true : false; }
                bool operator>(const Index& rhs) const { return r > rhs.r ? true : false; }

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

                value_type r;
            };

            class EXPORT Index
            {
            public:

                // Common typedef

                typedef int32_t  value_type;

                // Common interface

                Index();
                Index(const Index& in);
                Index(Index&& src);
                ~Index();

                Index& operator=(const Index& rhs);

                // Lattice interface

                Index(const value_type& in);

                explicit operator value_type&();

                Index& operator+=(const Index& rhs);
                Index& operator-=(const Index& rhs);
                Index& operator*=(const Index& rhs);
                Index& operator/=(const Index& rhs);

                Index operator++(int rhs);
                Index& operator++();
                Index operator--(int rhs);
                Index& operator--();

                value_type r;
            };
            */
        } // namespace Radial

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////
// Radial::Index non-member operators //
////////////////////////////////////////
/*
// Unary 
EXPORT Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& rhs);
namespace Multipole { namespace stl { namespace Radial { EXPORT double pown(const Index& base, int power); } } }

// Binary
EXPORT Multipole::stl::Radial::Index operator+(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator-(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator*(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator/(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT Multipole::stl::Radial::Index operator%(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);

EXPORT bool operator< (const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator> (const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator<=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator>=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator==(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
EXPORT bool operator!=(const Multipole::stl::Radial::Index& lhs, const Multipole::stl::Radial::Index& rhs);
*/