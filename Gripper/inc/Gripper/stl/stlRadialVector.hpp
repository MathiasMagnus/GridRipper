#pragma once

// Gripper includes
#include <Gripper/stl/stlParity.hpp>
#include <Gripper/stl/stlRadialIndex.hpp>
#include <Gripper/stl/stlRadialExtent.hpp>

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <utility>                  // Needed for std::pair, std::declval

namespace Multipole
{
    namespace stl
    {
        namespace Radial
        {
            template <typename DT, typename IT, typename VT>
            class VectorTraits
            {
            public:

                // Common typedefs

                typedef VT value_type;
                typedef DT distance_type;

                // STL typedefs

                typedef std::vector<value_type>                         container_type;
                typedef typename container_type::size_type              size_type;
                typedef typename container_type::iterator               iterator_type;
                typedef typename container_type::const_iterator         const_iterator_type;
                typedef typename container_type::reverse_iterator       reverse_iterator_type;
                typedef typename container_type::const_reverse_iterator const_reverse_iterator_type;

                // Lattice typedefs

                typedef Extent<IT>                              extent_type;
                typedef typename extent_type::index_type        index_type;
            };


            template <typename ET, typename DT, typename IT, typename VT>
            class Expression : public VectorTraits<DT, IT, VT>
            {
            public:

                // Common typedefs

                using typename VectorTraits<DT, IT, VT>::value_type;
                using typename VectorTraits<DT, IT, VT>::distance_type;

                // STL typedefs

                using typename VectorTraits<DT, IT, VT>::container_type;
                using typename VectorTraits<DT, IT, VT>::size_type;
                using typename VectorTraits<DT, IT, VT>::iterator_type;
                using typename VectorTraits<DT, IT, VT>::const_iterator_type;
                using typename VectorTraits<DT, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<DT, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs

                using typename VectorTraits<DT, IT, VT>::extent_type;
                using typename VectorTraits<DT, IT, VT>::index_type;
                
                // STL interface

                size_type   size()                  const { return static_cast<ET const&>(*this).size(); }
                value_type  operator[](size_type i) const { return static_cast<ET const&>(*this)[i]; }
                value_type  at(size_type i)         const { return static_cast<ET const&>(*this).at(i); }

                // Lattice interface

                distance_type separation()          const { return static_cast<ET const&>(*this).separation(); }
                extent_type extent()                const { return static_cast<ET const&>(*this).extent(); }
                value_type  at(const index_type& i) const { return static_cast<ET const&>(*this).at(i); }

                // Expression interface

                operator ET&()             { return static_cast<ET&>(*this); }
                operator ET const&() const { return static_cast<const ET&>(*this); }
            };


            template <typename DT, typename IT, typename VT>
            class Vector : public Expression<Vector<DT, IT, VT>, DT, IT, VT>
            {
            public:

                // Common typedefs

                using typename VectorTraits<DT, IT, VT>::value_type;
                using typename VectorTraits<DT, IT, VT>::distance_type;

                // STL typedefs

                using typename VectorTraits<DT, IT, VT>::container_type;
                using typename VectorTraits<DT, IT, VT>::size_type;
                using typename VectorTraits<DT, IT, VT>::iterator_type;
                using typename VectorTraits<DT, IT, VT>::const_iterator_type;
                using typename VectorTraits<DT, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<DT, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs

                using typename VectorTraits<DT, IT, VT>::extent_type;
                using typename VectorTraits<DT, IT, VT>::index_type;

                // Common interface

                Vector() = default;
                Vector(const Vector&) = default;
                Vector(Vector&&) = default;
                ~Vector() = default;

                Vector& operator=(Vector&) = default;
                Vector& operator=(Vector&&) = default;

                // STL interface

                iterator_type begin() { return m_data.begin(); }
                const_iterator_type begin() const { return m_data.begin(); }
                const_iterator_type cbegin() const { return m_data.cbegin(); }

                iterator_type end() { return m_data.end(); }
                const_iterator_type end() const { return m_data.end(); }
                const_iterator_type cend() const { return m_data.cend(); }

                reverse_iterator_type rbegin() { return m_data.rbegin(); }
                const_reverse_iterator_type rbegin() const { return m_data.rbegin(); }
                const_reverse_iterator_type crbegin() const { return m_data.crbegin(); }

                reverse_iterator_type rend() { return m_data.rend(); }
                const_reverse_iterator_type rend() const { return m_data.rend(); }
                const_reverse_iterator_type crend() const { return m_data.crend(); }

                value_type& at(size_type pos) { return m_data.at(pos); }                        // Returns specified element of the Vector with bounds checking
                const value_type& at(size_type pos) const { return m_data.at(pos); }            // Returns specified const element of the Vector with bounds checking
                value_type& operator[](size_type pos) { return m_data[pos]; }                   // Returns specified element of the Vector
                const value_type& operator[](size_type pos) const { return m_data[pos]; }       // Returns specified const element of the Vector

                value_type* data() { return m_data.data(); }                                    // Returns pointer to the underlying memory
                const value_type* data() const { return m_data.data(); }                        // Returns pointer to the underlying memory

                size_type size() const { return m_data.size(); }                                // Returns the number of elements inside the Vector
                void clear() { m_data.clear(); }                                                // Clear the contents of the Vector

                // Lattice interface

                Vector(const distance_type& separation, const extent_type& ext) : m_separation(separation), m_extent(ext), m_data(static_cast<size_type>(ext.final() - ext.initial() + 1)) {}

                //
                // Note: the const and non-const interface differs on purpose. While it would be possible to return a proxy object of a mirrored
                //       value, it is prohibited intentionally. It is an algorithmic mistake if one tries to modified values with no real storage.
                //

                value_type& at(const index_type& pos)
                {
                    if (!m_extent.contains(pos)) assert("Radial::Vector: index out of range");

                    return m_data.at(convert(pos));
                }

                const value_type& at(const index_type& pos) const
                {
                    if (pos >= 0)
                    {
                        if (!m_extent.contains(pos)) assert("Radial::Vector: index out of range");

                        return m_data.at(convert(pos));
                    }
                    else
                    {
                        if (!m_extent.contains(-pos)) assert("Radial::Vector: index out of range");

                        const auto& elem = m_data.at(convert(-pos));

                        return static_cast<value_type>(value_type::parity * value_type::l_parity(elem.extent()) * elem);
                    }
                }

                distance_type separation() const { return m_separation; }
                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, typename IndexType, typename ValueType>
                Vector(Expression<VecExpr, distance_type, IndexType, ValueType> const& expr)
                {
                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        //std::cout << "Resizing from " << m_extent << " to " << v.extent() << std::endl;
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    // Unconditional copy of separation
                    m_separation = v.separation();

                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i)
                        this->at(i) = v.at(i);
                }

            private:

                size_type convert(const index_type& i) const
                {
                    return static_cast<size_type>(i - m_extent.initial());
                }

                distance_type m_separation;
                extent_type m_extent;
                container_type m_data;
            };


            template <typename DT, typename ET, typename VT, typename F>
            struct Func : public Expression<Func<DT, ET, VT, F>, DT, typename ET::index_type, VT>
            {
            public:

                // Common typedefs

                using typename VectorTraits<DT, typename ET::index_type, VT>::value_type;
                using typename VectorTraits<DT, typename ET::index_type, VT>::distance_type;

                // STL typedefs

                using typename VectorTraits<DT, typename ET::index_type, VT>::size_type;

                // Lattice typedefs

                using typename VectorTraits<DT, typename ET::index_type, VT>::extent_type;
                using typename VectorTraits<DT, typename ET::index_type, VT>::index_type;

                Func(const distance_type& sep, const extent_type& ext, const F& f) : m_sep(sep), m_ext(ext), m_f(f) {}

                // STL interface

                size_type  size()                  const { return static_cast<size_type>(extent().final() - extent().initial() + 1); }
                value_type operator[](size_type i) const { return m_f[static_cast<size_type>(i - extent().initial())]; }
                value_type at(size_type i)         const { return m_f.at(static_cast<size_type>(i - extent().initial())); }

                // Lattice interface

                distance_type separation()         const { return m_sep; }
                extent_type extent()               const { return m_ext; }
                value_type at(index_type i)        const { return m_f(i); }

            private:

                distance_type m_sep;
                extent_type m_ext;
                F m_f;
            };


            template <typename E, typename F>
            class Map : public Expression<Map<E, F>, typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>::value_type;
                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>::distance_type;

                // STL typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>::extent_type;
                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename std::result_of<F(typename E::value_type)>::type>::index_type;

                Map(Expression<E, typename E::distance_type, typename E::index_type, typename E::value_type> const& u, const F& f) : _u(u), _f(f) {}

                // STL interface

                size_type  size()                  const { return _u.size(); }
                value_type operator[](size_type i) const { return _f(u[i]); }
                value_type at(size_type i)         const { return _f(_u.at(i)); }

                // Lattice interface

                distance_type separation()         const { return _u.separation(); }
                extent_type extent()               const { return _u.extent(); }
                value_type at(index_type i)        const { return _f(_u.at(i)); }

            private:

                const E& _u;
                F _f;
            };


            template <typename E1, typename E2, typename F>
            class Zip : public Expression<Zip<E1, E2, F>,
                                          decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()),
                                          typename E1::index_type,
                                          typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()), typename E1::index_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::value_type;
                using typename VectorTraits<decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()), typename E1::index_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::distance_type;

                // STL typedefs

                using typename VectorTraits<decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()), typename E1::index_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()), typename E1::index_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::extent_type;
                using typename VectorTraits<decltype(std::declval<typename E1::distance_type>() * std::declval<typename E2::distance_type>()), typename E1::index_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::index_type;

                Zip(Expression<E1, typename E1::distance_type, typename E1::index_type, typename E1::value_type> const& u, Expression<E2, typename E2::distance_type, typename E2::index_type, typename E2::value_type> const& v, const F& f) : _u(u), _v(v), _f(f)
                {
                    assert(u.separation() == v.separation());
                    assert(u.extent() == v.extent());
                }

                // STL interface

                size_type  size()                  const { return _u.size(); }
                value_type operator[](size_type i) const { return _f(u[i], _v[i]); }
                value_type at(size_type i)         const { return _f(_u.at(i), _v.at(i)); }

                // Lattice interface

                distance_type separation()         const { return _u.separation(); }
                extent_type extent()               const { return _u.extent(); }
                value_type at(index_type i)        const { return _f(_u.at(i), _v.at(i)); }

            private:

                const E1& _u;
                const E2& _v;
                F _f;
            };


            template <typename E>
            class Derivate : public Expression<Derivate<E>, typename E::distance_type, typename E::index_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename E::value_type>::value_type;
                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename E::value_type>::distance_type;

                // STL typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename E::value_type>::extent_type;
                using typename VectorTraits<typename E::distance_type, typename E::index_type, typename E::value_type>::index_type;

                Derivate(Expression<E, typename E::distance_type, typename E::index_type, typename E::value_type> const& u) : m_u(u) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { static_assert(false, "This function is not yet implemented"); return m_u(u[i]); }
                value_type at(size_type i)         const { static_assert(false, "This function is not yet implemented"); return m_u.at(i); }

                // Lattice interface

                distance_type separation()         const { return m_u.separation(); }
                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const
                {
                    if (i < extent().final() - 2)
                        return static_cast<value_type>((8 * (m_u.at(i+1) - m_u.at(i-1)) + m_u.at(i-2) - m_u.at(i+2)) / (12 * m_u.separation()));
                    else
                    {
                        if (i == extent().final() - 1)
                            return static_cast<value_type>((m_u.at(i+1) - m_u.at(i-1)) / (2 * m_u.separation()));
                        else
                            return static_cast<value_type>((m_u.at(i) - m_u.at(i-1)) / m_u.separation());
                    }
                }

            private:

                const E& m_u;
            };

            template <typename DT, typename ET, typename F>
            auto func(const DT& sep, const ET& ext, const F& f) { return Func<DT, ET, typename std::result_of<F(typename ET::index_type)>::type, F>(sep, ext, f); }

            template <typename E, typename F>
            auto map(const Expression<E, typename E::distance_type, typename E::index_type, typename E::value_type>& u, const F& f) { return Map<E, F>(u, f); }

            template <typename E1, typename E2, typename F>
            auto zip(const Expression<E1, typename E1::distance_type, typename E1::index_type, typename E1::value_type>& u, const Expression<E2, typename E2::distance_type, typename E2::index_type, typename E2::value_type>& v, const F& f) { return Zip<E1, E2, F>(u, v, f); }

            template <typename E>
            auto derivate(const Expression<E, typename E::distance_type, typename E::index_type, typename E::value_type>& u) { return Derivate<E>(u); }

        } // namespace Radial

    } // namespace stl

} // namespace Multipole


/////////////////////////////////////////
// Radial::Vector non-member operators //
/////////////////////////////////////////

// Unary 
template <typename E, typename D, typename I, typename T> auto operator+(Multipole::stl::Radial::Expression<E, D, I, T> const& v) { return Multipole::stl::Radial::map(v, [](auto&& val) { return static_cast<T>(1) * val; }); }
template <typename E, typename D, typename I, typename T> auto operator-(Multipole::stl::Radial::Expression<E, D, I, T> const& v) { return Multipole::stl::Radial::map(v, [](auto&& val) { return static_cast<T>(-1) * val; }); }

// Binary
template <typename E1, typename E2, typename D, typename I1, typename I2, typename T1, typename T2> auto operator+(Multipole::stl::Radial::Expression<E1, D, I1, T1> const& u, Multipole::stl::Radial::Expression<E2, D, I2, T2> const& v) { return Multipole::stl::Radial::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, typename D, typename I1, typename I2, typename T1, typename T2> auto operator-(Multipole::stl::Radial::Expression<E1, D, I1, T1> const& u, Multipole::stl::Radial::Expression<E2, D, I2, T2> const& v) { return Multipole::stl::Radial::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs - rhs; }); }
template <typename E1, typename E2, typename D, typename I1, typename I2, typename T1, typename T2> auto operator*(Multipole::stl::Radial::Expression<E1, D, I1, T1> const& u, Multipole::stl::Radial::Expression<E2, D, I2, T2> const& v) { return Multipole::stl::Radial::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs * rhs; }); }

template <typename E, typename D, typename I, typename T, typename Scalar> auto operator*(const Scalar& alpha, Multipole::stl::Radial::Expression<E, D, I, T> const& v) { return Multipole::stl::Radial::map(v, [=](auto&& val) { return alpha * val; }); }
template <typename E, typename D, typename I, typename T, typename Scalar> auto operator*(Multipole::stl::Radial::Expression<E, D, I, T> const& v, const Scalar& alpha) { return Multipole::stl::Radial::map(v, [=](auto&& val) { return alpha * val; }); }