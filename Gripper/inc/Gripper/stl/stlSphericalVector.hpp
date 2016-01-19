#pragma once

// Gripper includes
#include <Gripper/stl/stlParity.hpp>
#include <Gripper/stl/stlSphericalIndex.hpp>
#include <Gripper/stl/stlSphericalExtent.hpp>

// GSL includes
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf.h>

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <future>                   // Needed for guided parallel processing
#include <array>                    // Needed for guided parallel processing
#include <utility>                  // Needed for std::pair, std::declval
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            template <std::size_t L_Max, Parity P, typename IT, typename VT>
            class VectorTraits
            {
            public:

                // Common typedefs

                typedef VT value_type;

                // STL typedefs

                typedef std::vector<value_type>                         container_type;
                typedef typename container_type::size_type              size_type;
                typedef typename container_type::iterator               iterator_type;
                typedef typename container_type::const_iterator         const_iterator_type;
                typedef typename container_type::reverse_iterator       reverse_iterator_type;
                typedef typename container_type::const_reverse_iterator const_reverse_iterator_type;

                // Lattice typedefs

                typedef Extent<L_Max, IT>                       extent_type;
                typedef typename extent_type::index_type        index_type;
                typedef typename index_type::value_type         index_internal_type;

                // Lattice static members

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const Parity parity = P;
            };


            template <typename ET, std::size_t L_Max, Parity P, typename IT, typename VT>
            class Expression : public VectorTraits<L_Max, P, IT, VT>
            {
            public:

                // Common typedefs
                
                using expression_type = ET;
                using typename VectorTraits<L_Max, P, IT, VT>::value_type;

                // STL typedefs

                using typename VectorTraits<L_Max, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, P, IT, VT>::size_type;
                using typename VectorTraits<L_Max, P, IT, VT>::iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::const_iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs
 
                using typename VectorTraits<L_Max, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, P, IT, VT>::index_internal_type;
                
                // STL interface

                size_type   size()                  const { return static_cast<ET const&>(*this).size(); }
                value_type  operator[](size_type i) const { return static_cast<ET const&>(*this)[i]; }
                value_type  at(size_type i)         const { return static_cast<ET const&>(*this).at(i); }

                // Lattice interface

                extent_type extent()                const { return static_cast<ET const&>(*this).extent(); }
                value_type  at(const index_type& i) const { return static_cast<ET const&>(*this).at(i); }

                // Expression interface

                operator ET&()             { return static_cast<ET&>(*this); }
                operator ET const&() const { return static_cast<const ET&>(*this); }
            };

            template <typename E> auto l_parity(const typename E::extent_type&); // FIXME: this forward declaration should not really exist

            template <std::size_t L_Max, Parity P, typename IT, typename VT>
            class Vector : public Expression<Vector<L_Max, P, IT, VT>, L_Max, P, IT, VT>
            {
            public:

                // Common typedefs

                using typename VectorTraits<L_Max, P, IT, VT>::value_type;

                // STL typedefs

                using typename VectorTraits<L_Max, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, P, IT, VT>::size_type;
                using typename VectorTraits<L_Max, P, IT, VT>::iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::const_iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<L_Max, P, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs

                using typename VectorTraits<L_Max, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, P, IT, VT>::index_internal_type;

                template <typename E>
                static auto l_parity(const typename E::extent_type& ext) { return Multipole::stl::Spherical::l_parity(ext); }

                // Common interface

                Vector() = default;
                Vector(const Vector& in) = default;
                Vector(Vector&& in) = default;
                ~Vector() = default;

                Vector& operator=(const Vector&) = default;
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

                Vector(const extent_type& ext) : m_data(static_cast<size_type>(ext.final() - ext.initial())), m_extent(ext) {}

                value_type& at(const index_type& pos) { return m_data.at(index_type::convert(pos)); }
                const value_type& at(const index_type& pos) const { return m_data.at(index_type::convert(pos)); }

                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, std::size_t L_Max, typename IndexType, typename ValueType>
                Vector(Expression<VecExpr, L_Max, P, IndexType, ValueType> const& expr)
                {
                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i)
                        this->at(i) = v.at(i);
                }

            private:

                container_type m_data;
                extent_type m_extent;
            };


            template <typename E>
            struct Id : public Expression<Id<E>, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Id(Expression<E, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u) : m_u(u) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_u[i]; }
                value_type at(size_type i)         const { return m_u.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_u.at(i); }

            private:

                const E& m_u;
            };


            template <typename E, typename F>
            struct Func : public Expression<Func<E, F>, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Func(const extent_type& ext, const F& f) : m_ext(ext), m_f(f) {}

                // STL interface

                size_type  size()                  const { static_assert(false, "This function is not yet implemented"); return static_cast<size_type>(0); }
                value_type operator[](size_type i) const { static_assert(false, "This function is not yet implemented"); return m_f[i]; }
                value_type at(size_type i)         const { static_assert(false, "This function is not yet implemented"); return m_f.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_ext; }
                value_type at(index_type i)        const { return m_f(i); }

            private:

                extent_type m_ext;
                F m_f;
            };


            template <typename E, typename F>
            class Map : public Expression<Map<E, F>, E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::extent_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::index_type;
                using typename VectorTraits<E::l_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::index_internal_type;

                Map(Expression<E, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u, const F& f) : m_u(u), m_f(f) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_f(m_u[i]); }
                value_type at(size_type i)         const { return m_f(m_u.at(i)); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_f(m_u.at(i)); }

            private:

                const E& m_u;
                F m_f;
            };


            template <typename E1, typename E2, typename F>
            class Zip : public Expression<Zip<E1, E2, F>, E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::value_type;

                // STL typedefs

                using typename VectorTraits<E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::extent_type;
                using typename VectorTraits<E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::index_type;
                using typename VectorTraits<E1::l_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::index_internal_type;

                Zip(Expression<E1, E1::l_max, E1::parity, typename E1::index_internal_type, typename E1::value_type> const& u, Expression<E2, E2::l_max, E2::parity, typename E2::index_internal_type, typename E2::value_type> const& v, const F& f) : _u(u), _v(v), _f(f)
                {
                    static_assert(E1::l_max == E2::l_max, "Spherical::Zip(l_max value mismatch)");
                    static_assert(E1::parity == E2::parity, "Spherical::Zip(parity value mismatch)");

                    assert(u.extent() == v.extent());
                }

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_f(m_u[i], m_v[i]); }
                value_type at(size_type i)         const { return m_f(m_u.at(i), m_v.at(i)); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_f(m_u.at(i), m_v.at(i)); }

            private:

                const E1& m_u;
                const E2& m_v;
                F m_f;
            };


            template <typename E, Parity P>
            struct Par : public Expression<Par<E, P>, E::l_max, P, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, P, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, P, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, P, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, P, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, P, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Par(Expression<E, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u) : m_u(u) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_u[i]; }
                value_type at(size_type i)         const { return m_u.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_u.at(i); }

            private:

                const E& m_u;
            };

            template <typename E, typename F>
            auto map(const Expression<E, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F& f) { return Map<E, F>(u, f); }

            template <typename E1, typename E2, typename F>
            auto zip(const Expression<E1, E1::l_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u, const Expression<E2, E2::l_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v, const F& f) { return Zip<E1, E2, F>(u, v, f); }

            template <Parity P, typename E>
            auto parity(const Expression<E, E::l_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Par<E, P>(u); }

            namespace impl
            {
                template <typename E, typename F>
                auto l_parity(const typename E::extent_type& ext, const F& f) { return Func<E, F>(ext, f); }
            }

            template <typename E>
            auto l_parity(const typename E::extent_type& ext) {
                return impl::l_parity(ext, [](const typename E::index_type& i) { return i.l % 2 ? -1 : 1; });
            }

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename E> auto l_parity(const typename E::extent_type&); // FIXME: this forward declaration should not really exist

            

            template <std::size_t L_Max, std::size_t S_Max, Parity P, typename IT, typename VT>
            class VectorTraits
            {
            public:

                // Common typedefs

                typedef VT value_type;

                // STL typedefs

                typedef std::vector<value_type>                         container_type;
                typedef typename container_type::size_type              size_type;
                typedef typename container_type::iterator               iterator_type;
                typedef typename container_type::const_iterator         const_iterator_type;
                typedef typename container_type::reverse_iterator       reverse_iterator_type;
                typedef typename container_type::const_reverse_iterator const_reverse_iterator_type;

                // Lattice typedefs

                typedef Extent<L_Max, S_Max, IT>                        extent_type;
                typedef typename extent_type::index_type                index_type;
                typedef typename index_type::value_type                 index_internal_type;

                // Lattice static members

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const index_internal_type s_max = static_cast<index_internal_type>(S_Max);
                static const Parity parity = P;
            };


            template <typename ET, std::size_t L_Max, std::size_t S_Max, Parity P, typename IT, typename VT>
            class Expression : public VectorTraits<L_Max, S_Max, P, IT, VT>
            {
            public:

                // Common typedefs

                using expression_type = ET;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::value_type;
 
                // STL typedefs

                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::size_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::const_iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs

                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::index_internal_type;

                // STL interface

                size_type   size()                  const { return static_cast<ET const&>(*this).size(); }
                value_type  operator[](size_type i) const { return static_cast<ET const&>(*this)[i]; }
                value_type  at(size_type i)         const { return static_cast<ET const&>(*this).at(i); }

                // Lattice interface

                extent_type extent()                const { return static_cast<ET const&>(*this).extent(); }
                value_type  at(const index_type& i) const { return static_cast<ET const&>(*this).at(i); }

                // Expression interface

                operator ET&() { return static_cast<ET&>(*this); }
                operator ET const&() const { return static_cast<const ET&>(*this); }
            };


            template <typename E>
            struct Id : public Expression<Id<E>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Id(const Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) : m_u(u) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_u[i]; }
                value_type at(size_type i)         const { return m_u.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_u.at(i); }

            private:

                const E& m_u;
            };


            template <typename E, typename F>
            struct Func : public Expression<Func<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Func(const extent_type& ext, const F& f) : m_ext(ext), m_f(f) {}

                // STL interface

                size_type  size()                  const { static_assert(false, "This function is not yet implemented"); return static_cast<size_type>(0); }
                value_type operator[](size_type i) const { static_assert(false, "This function is not yet implemented"); return m_f[i]; }
                value_type at(size_type i)         const { static_assert(false, "This function is not yet implemented"); return m_f.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_ext; }
                value_type at(index_type i)        const { return m_f(i); }

            private:

                extent_type m_ext;
                F m_f;
            };


            template <typename E, typename F>
            class Map : public Expression<Map<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::extent_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::index_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>::index_internal_type;

                Map(Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u, const F& f) : _u(u), _f(f) {}

                // STL interface

                size_type  size()                  const { return _u.size(); }
                value_type operator[](size_type i) const { return _f(u[i]); }
                value_type at(size_type i)         const { return _f(_u.at(i)); }

                // Lattice interface

                extent_type extent()               const { return _u.extent(); }
                value_type at(index_type i)        const { return _f(_u.at(i)); }

            private:

                const E& _u;
                F _f;
            };


            template <typename E1, typename E2, typename F>
            class Zip : public Expression<Zip<E1, E2, F>, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::value_type;

                // STL typedefs

                using typename VectorTraits<E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::extent_type;
                using typename VectorTraits<E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::index_type;
                using typename VectorTraits<E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>::index_internal_type;

                Zip(const Expression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u, const Expression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v, const F& f) : _u(u), _v(v), _f(f)
                {
                    static_assert(E1::l_max == E2::l_max, "Spherical::Zip(l_max value mismatch)");
                    static_assert(E1::s_max == E2::s_max, "Spherical::Zip(s_max value mismatch)");
                    static_assert(E1::parity == E2::parity, "Spherical::Zip(parity value mismatch)");

                    assert(u.extent() == v.extent());
                }

                // STL interface

                size_type  size()                  const { return _u.size(); }
                value_type operator[](size_type i) const { return _f(u[i], _v[i]); }
                value_type at(size_type i)         const { return _f(_u.at(i), _v.at(i)); }

                // Lattice interface

                extent_type extent()               const { return _u.extent(); }
                value_type at(index_type i)        const { return _f(_u.at(i), _v.at(i)); }

            private:

                const E1& _u;
                const E2& _v;
                F _f;
            };


            template <typename E, Parity P>
            struct Par : public Expression<Par<E, P>, E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::s_max, P, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Par(Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u) : m_u(u) {}

                // STL interface

                size_type  size()                  const { return m_u.size(); }
                value_type operator[](size_type i) const { return m_u[i]; }
                value_type at(size_type i)         const { return m_u.at(i); }

                // Lattice interface

                extent_type extent()               const { return m_u.extent(); }
                value_type at(index_type i)        const { return m_u.at(i); }

            private:

                const E& m_u;
            };

            
            namespace Spin
            {
                struct Up {};
                struct Down {};
            }

            namespace impl
            {
                template <typename S> struct SpinStepper;

                template <> struct SpinStepper<Spin::Up>
                {
                    template <typename E>
                    static auto at(const E& e, const typename E::index_type& i)
                    {
                        typename E::value_type result = 0;

                        if (i.s == std::min(i.l, E::s_max)) return result;

                        auto factor_sq = static_cast<typename E::index_internal_type>((i.l - i.s) * (i.l + i.s + static_cast<typename E::index_internal_type>(1)));

                        if (factor_sq <= static_cast<typename E::index_internal_type>(0)) return result;

                        auto index = i;

                        result = static_cast<decltype(result)>(std::sqrt(factor_sq) * e.at(++index));

                        return result;
                    };
                };

                template <> struct SpinStepper<Spin::Down>
                {
                    template <typename E>
                    static auto at(const E& e, const typename E::index_type& i)
                    {
                        typename E::value_type result = 0;

                        if (i.s == std::max(-i.l, -E::s_max)) return result;

                        auto factor_sq = static_cast<typename E::index_internal_type>((i.l + i.s) * (i.l - i.s + static_cast<typename E::index_internal_type>(1)));

                        if (factor_sq <= static_cast<typename E::index_internal_type>(0)) return result;

                        auto index = i;

                        result = static_cast<decltype(result)>(-std::sqrt(factor_sq) * e.at(--index));

                        return result;
                    };
                };
            }

            template <typename E, typename S>
            class Eth : public Expression<Eth<E, S>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Common typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::value_type;

                // STL typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::size_type;

                // Lattice typedefs

                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::extent_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_type;
                using typename VectorTraits<E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>::index_internal_type;

                Eth(Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type> const& u) : _u(u) {}

                // STL interface

                size_type  size()                  const { return _u.size(); }
                //value_type operator[](size_type i) const { return impl::SpinStepper<S>::at(_u, index_type::convert(i)); }
                //vaue_type at(size_type i)         const { return impl::SpinStepper<S>::at(_u, index_type::convert(i)); }

                // Lattice interface

                extent_type extent()               const { return _u.extent(); }
                value_type at(index_type i)        const { return impl::SpinStepper<S>::at(_u, i); }

            private:

                const E& _u;
            };
            
            template <typename E, typename F>
            auto map(const Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F& f) { return Map<E, F>(u, f); }

            template <typename E1, typename E2, typename F>
            auto zip(const Expression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u, const Expression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v, const F& f) { return Zip<E1, E2, F>(u, v, f); }

            template <Parity P, typename E>
            auto parity(const Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Par<E, P>(u); }

            namespace impl
            {
                template <typename E, typename F>
                auto l_parity(const typename E::extent_type& ext, const F& f) { return Func<E, F>(ext, f); }
            }

            template <typename E>
            auto l_parity(const typename E::extent_type& ext) { return impl::l_parity<E>(ext, [](const typename E::index_type& i) { return i.l % 2 ? -1 : 1; }); }
            
            template <typename E>
            auto eth(const Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Eth<E, Spin::Down>(u); }

            template <typename E>
            auto eth_bar(const Expression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Eth<E, Spin::Up>(u); }


            template <std::size_t L_Max, std::size_t S_Max, Parity P, typename IT, typename VT>
            class Vector : public Expression<Vector<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>
            {
            public:

                // Common typedefs
;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::value_type;

                // STL typedefs

                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::size_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::const_iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::reverse_iterator_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::const_reverse_iterator_type;

                // Lattice typedefs

                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, S_Max, P, IT, VT>::index_internal_type;

                static auto l_parity(const extent_type& ext) { return Multipole::stl::SpinWeightedSpherical::l_parity<Expression<Vector<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>>(ext); }

                // Common interface

                Vector() = default;
                Vector(const Vector& in) = default;
                Vector(Vector&& in) = default;
                ~Vector() = default;

                Vector& operator=(const Vector&) = default;
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

                Vector(const extent_type ext) : m_extent(ext) // TODO: when a proper indexing function has been devised, replace loop with initializer
                {
                    std::size_t counter = 0;

                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i)
                        ++counter;

                    m_data.resize(counter);
                }

                value_type& at(const index_type& pos)
                {
                    if (!m_extent.contains(pos)) assert("SpinWeightedSpherical::Vector: index out of range");

                    return m_data.at(convert(pos));
                }
                const value_type& at(const index_type& pos) const
                {
                    if (!m_extent.contains(pos)) assert("SpinWeightedSpherical::Vector: index out of range");

                    return m_data.at(convert(pos));
                }

                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, typename IndexType, typename ValueType>
                Vector(Expression<VecExpr, L_Max, S_Max, P, IndexType, ValueType> const& expr) : Vector()
                {
                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    // Serial evaluator
                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i)
                        this->at(i) = v.at(i);
                    /*
                    // Parallel evaluator
                    std::size_t count = 0;
                    std::array<std::future<void>, std::thread::hardware_concurrency()> futures;
                    std::array<extent_type, std::thread::hardware_concurrency()> extents;

                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i) ++count;

                    for (auto I = 0; I < extents.size(); ++I)
                        for (index_type i = m_extent.initial() + i * (count / extents.size()); ++i)
                            */
                }

                index_type convert(const size_type& i) const
                {
                    static_assert(false, "This function currently does not work");

                    IT i_corr = static_cast<IT>(i / (2 * m_extent.final().s + 1)); // truncation on division is intentional (omit std::floor and double up-cast)
                    IT l = static_cast<IT>(std::round((-1. + std::sqrt(1. + 4 * i_corr)) / 2.));
                    IT m = i_corr - l*(l + 1);
                    IT s = (i % (2 * m_extent().final.s + 1)) - m_extent.final().s;

                    return index_type{ l, m, s };
                }

            private:

                size_type convert(const index_type& i) const
                {
                    // FIXME: I am a terribly slow indexing funcion work-around!!!

                    std::size_t counter = 0;

                    for (index_type I = m_extent.initial(); I != i; ++I)
                        ++counter;

                    return counter;
                }

                extent_type m_extent;
                container_type m_data;
            };
            
        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Spherical::Vector non-member operators //
////////////////////////////////////////////

// Unary 
template <typename E, std::size_t L, Multipole::stl::Parity P, typename I, typename T> auto operator+(Multipole::stl::Spherical::Expression<E, L, P, I, T> const& v) { return Multipole::stl::Spherical::map(v, [](auto&& val) { return static_cast<T>(1) * val; }); }
template <typename E, std::size_t L, Multipole::stl::Parity P, typename I, typename T> auto operator-(Multipole::stl::Spherical::Expression<E, L, P, I, T> const& v) { return Multipole::stl::Spherical::map(v, [](auto&& val) { return static_cast<T>(-1) * val; }); }

// Binary
template <typename E1, typename E2, std::size_t L, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(Multipole::stl::Spherical::Expression<E1, L, P, I1, T1> const& u, Multipole::stl::Spherical::Expression<E2, L, P, I2, T2> const& v) { return Multipole::stl::Spherical::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, std::size_t L, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(Multipole::stl::Spherical::Expression<E1, L, P, I1, T1> const& u, Multipole::stl::Spherical::Expression<E2, L, P, I2, T2> const& v) { return Multipole::stl::Spherical::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs - rhs; }); }

template <typename E, std::size_t L, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar& alpha, Multipole::stl::Spherical::Expression<E, L, P, I, T> const& v) { return Multipole::stl::Spherical::map(v, [=](auto&& val) { return alpha * val; }); }
template <typename E, std::size_t L, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(Multipole::stl::Spherical::Expression<E, L, P, I, T> const& v, const Scalar& alpha) { return Multipole::stl::Spherical::map(v, [=](auto&& val) { return alpha * val; }); }

////////////////////////////////////////////////////////
// SpinWeightedSpherical::Vector non-member operators //
////////////////////////////////////////////////////////

// Unary 
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator+(Multipole::stl::SpinWeightedSpherical::Expression<E, L, S, P, I, T> const& v) { return Multipole::stl::SpinWeightedSpherical::map(v, [](auto&& val) { return static_cast<T>(1) * val; }); }
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator-(Multipole::stl::SpinWeightedSpherical::Expression<E, L, S, P, I, T> const& v) { return Multipole::stl::SpinWeightedSpherical::map(v, [](auto&& val) { return static_cast<T>(-1) * val; }); }

// Binary
template <typename E1, typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const Multipole::stl::SpinWeightedSpherical::Expression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SpinWeightedSpherical::Expression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SpinWeightedSpherical::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const Multipole::stl::SpinWeightedSpherical::Expression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SpinWeightedSpherical::Expression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SpinWeightedSpherical::zip(u, v, [](auto&& lhs, auto&& rhs) { return lhs - rhs; }); }

template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar& alpha, const Multipole::stl::SpinWeightedSpherical::Expression<E, L, S, P, I, T>& v) { return Multipole::stl::SpinWeightedSpherical::map(v, [=](auto&& val) { return alpha * val; }); }
//template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Multipole::stl::SpinWeightedSpherical::Expression<E, L, S, P, I, T>& v, const Scalar& alpha) { return Multipole::stl::SpinWeightedSpherical::map(v, [=](auto&& val) { return alpha * val; }); }