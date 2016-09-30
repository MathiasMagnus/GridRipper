#pragma once

// Gripper includes
#include <Gripper/stl/stlParity.hpp>                        // Multipole::stl::Parity
#include <Gripper/stl/stlSpinWeightedSphericalIndex.hpp>    // Multipole::stl::SWS::Index
#include <Gripper/stl/stlSpinWeightedSphericalExtent.hpp>   // Multipole::stl::SWS::Extent

// Standard C++ includes
#include <cstddef>                  // std::size_t
#include <functional>               // std::reference_wrapper, std::ref, std::cref
#include <vector>                   // std::vector
#include <type_traits>              // std::result_of
#include <complex>                  // std::complex


namespace Multipole
{
    namespace stl
    {
        namespace SWS
        {
            // Forward decalaration of l_parity helper
            template <typename E> auto l_parity(const typename E::extent_type&); // FIXME: this forward declaration should not really exist

            /// <summary>Traits class consisting of type aliases describing spin-weighted spherical expansion coefficient vectors.</summary>
            ///
            template <std::size_t L_Max, std::size_t S, Parity P, typename IT, typename VT>
            struct VectorTraits
            {
                // Common type aliases

                using value_type = VT;

                // STL type aliases

                using container_type = std::vector<value_type>;
                using size_type = typename container_type::size_type;

                // Lattice type aliases

                using extent_type = Extent<L_Max, S, IT>;
                using index_type = typename extent_type::index_type;
                using index_internal_type = typename index_type::value_type;

                // Lattice static members

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const index_internal_type s = static_cast<index_internal_type>(S);
                static const Parity parity = P;
            };


            /// <summary>Read-only Expression Template base class of spin-weighted spherical expansion coefficient vectors to be implemented statically.</summary>
            ///
            template <typename ET, std::size_t L_Max, std::size_t S, Parity P, typename IT, typename VT>
            struct ConstExpression : public VectorTraits<L_Max, S, P, IT, VT>
            {
                // ConstExpression type aliases

                using expression_type = ET;

                // Common type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::value_type;

                // STL type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::size_type;

                // Lattice type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using VectorTraits<L_Max, S, P, IT, VT>::l_max;
                using VectorTraits<L_Max, S, P, IT, VT>::s;
                using VectorTraits<L_Max, S, P, IT, VT>::parity;

                // ConstExpression interface

                operator const expression_type&()   const { return static_cast<const expression_type&>(*this); }

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type   size()                  const { return static_cast<const expression_type&>(*this).size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type  operator[](size_type i) const { return static_cast<const expression_type&>(*this)[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type  at(size_type i)         const { return static_cast<const expression_type&>(*this).at(i); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()                const { return static_cast<const expression_type&>(*this).extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type  at(const index_type& i) const { return static_cast<const expression_type&>(*this).at(i); }
            };


            /// <summary>Read-write Expression Template base class of spherical expansion coefficient vectors to be implemented statically.</summary>
            ///
            template <typename ET, std::size_t L_Max, std::size_t S, Parity P, typename IT, typename VT>
            struct Expression : public ConstExpression<ET, L_Max, S, P, IT, VT>
            {
                // Expression type aliases

                using expression_type = ET;

                // Common type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::value_type;

                // STL type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::size_type;

                // Lattice type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using VectorTraits<L_Max, S, P, IT, VT>::l_max;
                using VectorTraits<L_Max, S, P, IT, VT>::s;
                using VectorTraits<L_Max, S, P, IT, VT>::parity;

                // Expression interface

                operator expression_type&() { return static_cast<expression_type&>(*this); }

                // STL interface

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type& operator[](size_type i) { return static_cast<expression_type&>(*this)[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type& at(size_type i) { return static_cast<expression_type&>(*this).at(i); }

                // Lattice interface

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type& at(const index_type& i) { return static_cast<expression_type&>(*this).at(i); }
            };


            /// <summary>Class storing the coefficients of a series expansion over spin-weighted spherical harmonics.</summary>
            ///
            /// <remarks>The Vector class template is not an Expression Template. There are seperate View classes for that.</remarks>
            /// <remarks>Also note that the Vector class slightly overallocates, due to the complexity of inverting the indexing function.</remarks>
            ///
            template <std::size_t L_Max, std::size_t S, Parity P, typename IT, typename VT>
            class Vector : public VectorTraits<L_Max, S, P, IT, VT>
            {
            public:

                // Common type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::value_type;

                // STL type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::container_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::size_type;

                // Lattice type aliases

                using typename VectorTraits<L_Max, S, P, IT, VT>::extent_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_type;
                using typename VectorTraits<L_Max, S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using VectorTraits<L_Max, S, P, IT, VT>::l_max;
                using VectorTraits<L_Max, S, P, IT, VT>::s;
                using VectorTraits<L_Max, S, P, IT, VT>::parity;

                static auto l_parity(const extent_type& ext) { return Multipole::stl::SWS::l_parity<Expression<Vector<L_Max, S, P, IT, VT>, L_Max, S, P, IT, VT>>(ext); }

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                Vector() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                Vector(const Vector& in) = default;

                /// <summary>Default move constructor.</summary>
                ///
                Vector(Vector&& in) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~Vector() = default;
                //~Vector() { std::cout << "Vector DTOR" << std::endl; };

                /// <summary>Default copy assignment operator.</summary>
                ///
                Vector& operator=(const Vector&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                Vector& operator=(Vector&&) = default;

                ///<summary>Constructs a series expansion vector of <c>ext</c> size.</summary>
                ///
                Vector(const extent_type& ext) : m_extent(ext), m_data() { m_data.resize(distance(ext.initial(), ext.final())); }

                ///<summary>Constructs a series expansion vector from the expression <c>expr</c>.</summary>
                ///
                template <typename ConstVecExpr, std::size_t L_Max, std::size_t S, typename IndexType, typename ValueType>
                Vector(const ConstExpression<ConstVecExpr, L_Max, S, P, IndexType, ValueType>& expr)
                {
                    serial_evaluator(expr);
                }

                ///<summary>Constructs a series expansion vector from the expression <c>expr</c>.</summary>
                ///
                template <typename ConstVecExpr, std::size_t L_Max, std::size_t S, typename IndexType, typename ValueType>
                Vector& operator=(const ConstExpression<ConstVecExpr, L_Max, S, P, IndexType, ValueType>& expr)
                {
                    serial_evaluator(expr);

                    return *this;
                }

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type size() const { return m_data.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                const value_type& operator[](size_type i) const { return m_data[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                const value_type&  at(size_type i)         const { return m_data.at(i); }

                // STL interface

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type& operator[](size_type i) { return m_data[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type& at(size_type i) { return m_data.at(i); }

                // Extra STL functions

                /// <summary>Returns pointer to the underlying memory.</summary>
                ///
                value_type* data() { return m_data.data(); }

                /// <summary>Returns pointer to the underlying memory.</summary>
                ///
                const value_type* data() const { return m_data.data(); }

                /// <summary>Clears the contents of the Vector.</summary>
                ///
                void clear() { m_data.clear(); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent() const { return m_extent; }

                /// <summary>Returns the coefficient corresponding to index <c>pos</c>.</summary>
                ///
                const value_type& at(const index_type& pos) const
                {
#ifndef NDEBUG
                    if (!m_extent.contains(pos)) assert("SWS::Vector: index out of range");
#endif
                    return m_data.at(convert(pos));
                }

                // Lattice interface

                /// <summary>Returns the coefficient corresponding to index <c>pos</c>.</summary>
                ///
                value_type& at(const index_type& pos)
                {
#ifndef NDEBUG
                    if (!m_extent.contains(pos)) assert("SWS::Vector: index out of range");
#endif
                    return m_data.at(convert(pos));
                }

            private:

                /// <summary>Converts from size_type to index_type, when loopong over STL indices but need Multipole index values.</summary>
                ///
                index_type convert(const size_type& i) const
                {
                    index_type result{ 0, 0, 0 };

                    result.l = std::lround((std::sqrt(1 + 4 * i) - 1) / 2);
                    result.m = i - result.l * (result.l + 1);
                    result.s = s;

                    return index_type;
                }

                /// <summary>Converts from index_type to size_type for indexing into own container.</summary>
                ///
                size_type convert(const index_type& pos) const
                {
                    return pos.l * (pos.l + static_cast<index_internal_type>(1u)) + pos.m;
                }

                /// <summary>Converts from index_type to size_type for indexing into own container.</summary>
                ///
                /// <note>The distance of two indicies only has meaning in the context of a vector or matrix, but not the extent itself.</note>
                ///
                size_type distance(const index_type& from, const index_type to) const
                {
                    return convert(to) - convert(from);
                }

                /// <summary>Initializes internal states and evaluates elements of the expression <c>expr</c> in a serial manner.</summary>
                ///
                template <typename ConstVecExpr, std::size_t L_Max, std::size_t S_Max, typename IndexType, typename ValueType>
                void serial_evaluator(const ConstExpression<ConstVecExpr, L_Max, S_Max, P, IndexType, ValueType>& expr)
                {
                    // Extract type from encapsulating expression
                    const ConstVecExpr& v = expr;

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    for (index_type i = m_extent.initial(); m_extent.contains(i); ++i)
                        this->at(i) = v.at(i);
                }

                extent_type m_extent;
                container_type m_data;
            };


            /// <summary>Expression Template providing read-only access with reference semantics to the elements of a series expansion with storage.</summary>
            ///
            template <std::size_t L_Max, std::size_t S_Max, Parity P, typename IT, typename VT>
            class ConstView : public ConstExpression<ConstView<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<ConstView<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::container_type;
                using typename expression_type::size_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

            private:

                using vector_type = Vector<l_max, s_max, parity, index_internal_type, value_type>;
                using reference_type = std::reference_wrapper<const vector_type>;

            public:

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
                ///
                ConstView() = delete;

                /// <summary>Default copy constructor.</summary>
                ///
                ConstView(const ConstView&) = default;

                /// <summary>Default move constructor.</summary>
                /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
                ///
                ConstView(ConstView&&) = default;//delete;

                /// <summary>Default destructor.</summary>
                ///
                ~ConstView() = default;

                /// <summary>Default copy assign operator.</summary>
                ///
                ConstView& operator=(const ConstView&) = default;

                /// <summary>Default move assign operator.</summary>
                /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
                ///
                ConstView& operator=(ConstView&&) = delete;

                /// <summary>Constructs a <c>View</c> from a <c>Vector</c>.</summary>
                ///
                ConstView(const vector_type& v) : _v(std::cref(v)) {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                   const { return _v.get().size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                const value_type& operator[](size_type i) const { return _v.get()[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                const value_type& at(size_type i)         const { return _v.get().at(i); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent() const { return _v.get().extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                const value_type& at(const index_type& i) const { return _v.get().at(i); }

            private:

                const reference_type _v;
            };


            /// <summary>Expression Template providing read-write access with reference semantics to the elements of a series expansion with storage.</summary>
            ///
            template <std::size_t L_Max, std::size_t S_Max, Parity P, typename IT, typename VT>
            class View : public Expression<View<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = Expression<View<L_Max, S_Max, P, IT, VT>, L_Max, S_Max, P, IT, VT>;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::container_type;
                using typename expression_type::size_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

            private:

                using vector_type = Vector<l_max, s_max, parity, index_internal_type, value_type>;
                using reference_type = std::reference_wrapper<vector_type>;

            public:

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
                ///
                View() = delete;

                /// <summary>Default copy constructor.</summary>
                ///
                View(const View&) = default;

                /// <summary>Default move constructor.</summary>
                /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
                ///
                View(View&&) = default;//delete;

                /// <summary>Default destructor.</summary>
                ///
                ~View() = default;

                /// <summary>Default copy assign operator.</summary>
                ///
                View& operator=(const View&) = default;

                /// <summary>Default move assign operator.</summary>
                /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
                ///
                View& operator=(View&&) = default;//delete;

                /// <summary>Constructs a <c>View</c> from a <c>Vector</c>.</summary>
                ///
                explicit View(vector_type& v) : _v(std::ref(v)) {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                  const { return _v.get().size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                const value_type& operator[](size_type i) const { return _v.get()[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                const value_type& at(size_type i)         const { return _v.get().at(i); }

                // STL interface

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type& operator[](size_type i) { return _v.get()[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type& at(size_type i) { return _v.get().at(i); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent() const { return _v.get().extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                const value_type& at(const index_type& i) const { return _v.get().at(i); }

                // Lattice interface

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type& at(const index_type& i) { return _v.get().at(i) }

            private:

                reference_type _v;
            };


            /// <summary>Expression Template performing Identity transformation on a series expansion.</summary>
            ///
            template <typename E>
            class Id : public ConstExpression<Id<E>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<Id<E>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>;
                using outer_expression_type = E;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::size_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs an <c>Id</c> from a <c>ConstExpression</c>.</summary>
                ///
                Id(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u)
                    : _u(static_cast<const E&>(u))
                {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                  const { return _u.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type operator[](size_type i) const { return _u[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type at(size_type i)         const { return _u.at(i); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _u.extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const { return _u.at(i); }

            private:

                const E _u;
            };


            /// <summary>Expression Template that evaluates a function object to yield series expansion coefficients.</summary>
            ///
            template <typename E, typename F>
            struct Func : public ConstExpression<Func<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<Func<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>;

                // Common type alises

                using typename expression_type::value_type;

                // STL type alises

                using typename expression_type::size_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Func</c> from an extent <c>ext</c> and a function object <c>f</c>.</summary>
                ///
                Func(const extent_type ext, const F f) : m_ext(ext), m_f(f) {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                size_type  size()                  const { static_assert(false, "This function is not yet implemented"); return static_cast<size_type>(0); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                value_type operator[](size_type i) const { static_assert(false, "This function is not yet implemented"); return _f[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                value_type at(size_type i)         const { static_assert(false, "This function is not yet implemented"); return _f.at(i); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _ext; }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const { return _f(i); }

            private:

                extent_type _ext;
                F _f;
            };


            /// <summary>Expression Template that applies a unary function object to every element of a series expansion.</summary>
            ///
            template <typename E, typename F>
            class Map : public ConstExpression<Map<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<Map<E, F>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>;
                using outer_expression_type = E;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::size_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Map</c> from an expression <c>u</c> and a function object <c>f</c>.</summary>
                ///
                Map(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F f)
                    : _u(static_cast<const E&>(u))
                    , _f(f)
                {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                  const { return _u.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type operator[](size_type i) const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u[i]); }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type at(size_type i)         const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u.at(i)); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _u.extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const { return _f(_u.at(i)); }

            private:

                const E _u;
                const F _f;
            };


            /// <summary>Expression Template that applies a binary function object to every element of two series expansions to yield the result.</summary>
            ///
            template <typename E1, typename E2, typename F>
            class Zip : public ConstExpression<Zip<E1, E2, F>, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<Zip<E1, E2, F>, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>;
                using outer_expression_type1 = E1;
                using outer_expression_type2 = E2;

                // Common typedefs

                using typename expression_type::value_type;

                // STL typedefs

                using typename expression_type::size_type;

                // Lattice typedefs

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Zip</c> from two possibly varying expressions <c>u</c> and <c>v</c> and a function object <c>f</c>.</summary>
                ///
                Zip(const ConstExpression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u,
                    const ConstExpression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v,
                    const F f)
                    : _u(static_cast<const E1&>(u))
                    , _v(static_cast<const E2&>(v))
                    , _f(f)
                {
                    //std::cout << "sw::Zip _u = " << _u.extent() << " _v = " << _v.extent() << std::endl;

                    static_assert(E1::l_max == E2::l_max, "SWS::Zip(l_max value mismatch)");
                    static_assert(E1::s_max == E2::s_max, "SWS::Zip(s_max value mismatch)");
                    static_assert(E1::parity == E2::parity, "SWS::Zip(parity value mismatch)");

                    assert(_u.extent() == _v.extent());
                }

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                  const { return _u.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type operator[](size_type i) const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u[i], _v[i]); }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type at(size_type i)         const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u.at(i), _v.at(i)); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _u.extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const { return _f(_u.at(i), _v.at(i)); }

            private:

                const E1 _u;
                const E2 _v;
                F _f;
            };
            
            /// <summary>Namespace for tags of the spin stepper operator.</summary>
            ///
            namespace Spin
            {
                struct Up {};
                struct Down {};
            }

            namespace impl
            {
                template <typename VT> struct realify;
                template <> struct realify<float> { using type = float; };
                template <> struct realify<double> { using type = double; };
                template <> struct realify<std::complex<float>> { using type = float; };
                template <> struct realify<std::complex<double>> { using type = double; };

                /// <summary>Unimplemented Exression Template of the spin stepping edth operator.</summary>
                ///
                template <typename S> struct SpinStepper;

                /// <summary>Specialization of the spin stepping helper function for the Expression Template implementing the spin-up operator.</summary>
                /// <note>J. N. Goldberg, A. J. Macfarlane, E. T. Newman, F. Rohrlich, and E. C. G. Sudarshan Spin‐s Spherical Harmonics and ð, eq. (2.7a)</note>
                ///
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
                        
                        result = static_cast<decltype(result)>(static_cast<typename realify<decltype(result)>::type>(std::sqrt(factor_sq)) * e.at(++index));

                        return result;
                    }
                };


                /// <summary>Specialization of the spin stepping helper function for the Expression Template implementing the spin-down operator.</summary>
                /// <note>J. N. Goldberg, A. J. Macfarlane, E. T. Newman, F. Rohrlich, and E. C. G. Sudarshan Spin‐s Spherical Harmonics and ð, eq. (2.7b)</note>
                ///
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

                        result = static_cast<decltype(result)>(-static_cast<typename realify<decltype(result)>::type>(std::sqrt(factor_sq)) * e.at(--index));

                        return result;
                    }
                };

            } // namespace impl


            /// <summary>Expression Template that applies applies the spin stepping operator (edth) to a series expansion.</summary>
            ///
            template <typename E, typename S>
            class Edth : public ConstExpression<Edth<E, S>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Expression type aliases

                using expression_type = ConstExpression<Edth<E, S>, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>;
                using outer_expression_type = E;

                // Common typedefs

                using typename expression_type::value_type;

                // STL typedefs

                using typename expression_type::size_type;

                // Lattice typedefs

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::l_max;
                using expression_type::s_max;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Zip</c> from two possibly varying expressions <c>u</c> and <c>v</c> and a function object <c>f</c>.</summary>
                ///
                Edth(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u)
                    : _u(static_cast<const E&>(u))
                {}

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                ///
                size_type  size()                  const { return _u.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type operator[](size_type i) const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u[i], _v[i]); }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                /// <remarks>This function cannot be implemented with a single function object.</remarks>
                ///
                value_type at(size_type i)         const { static_assert(false, "This function cannot be implemented with a single function object"); return _f(_u.at(i), _v.at(i)); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _u.extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const { return impl::SpinStepper<S>::at(_u, i); }

            private:

                const E _u;
            };
            
            template <typename E, typename F>
            auto map(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F f) { return Map<E, F>(u, f); }

            template <typename E1, typename E2, typename F>
            auto zip(const ConstExpression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u,
                     const ConstExpression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v,
                     const F f)
            {
                return Zip<E1, E2, F>(u, v, f);
            }

            //template <Parity P, typename E>
            //auto parity(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Par<E, P>(u); }

            namespace impl
            {
                template <typename E, typename F>
                auto l_parity(const typename E::extent_type& ext, const F& f) { return Func<E, F>(ext, f); }
            }

            template <typename E>
            auto l_parity(const typename E::extent_type& ext) { return impl::l_parity<E>(ext, [](const typename E::index_type& i) { return i.l % 2 ? -1 : 1; }); }
            
        } // namespace SWS

    } // namespace stl

} // namespace Multipole

////////////////////////////////////////////////////////
// SWS::Vector non-member operators //
////////////////////////////////////////////////////////

// Unary 
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator+(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v) { return Multipole::stl::SWS::map(v, [](const auto& val) { return static_cast<T>(1) * val; }); }
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator-(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v) { return Multipole::stl::SWS::map(v, [](const auto& val) { return static_cast<T>(-1) * val; }); }
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto conjugate(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v) { return Multipole::stl::SWS::map(v, [](const auto& val) { return std::conj(val); }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator+(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [](const auto& val) { return static_cast<T>(1) * val; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator-(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [](const auto& val) { return static_cast<T>(-1) * val; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto conjugate(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [](const auto& val) { return std::conj(val); }); }

// Binary
template <typename E1, typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const Multipole::stl::SWS::ConstExpression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SWS::ConstExpression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const Multipole::stl::SWS::ConstExpression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SWS::ConstExpression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const Multipole::stl::SWS::Vector<L, S, P, I1, T1>& u, const Multipole::stl::SWS::ConstExpression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(Multipole::stl::SWS::ConstView<L, S, P, I1, T1>(u), v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E2, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const Multipole::stl::SWS::Vector<L, S, P, I1, T1>& u, const Multipole::stl::SWS::ConstExpression<E2, L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(Multipole::stl::SWS::ConstView<L, S, P, I1, T1>(u), v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E1, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const Multipole::stl::SWS::ConstExpression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SWS::Vector<L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(u, Multipole::stl::SWS::ConstView<L, S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const Multipole::stl::SWS::ConstExpression<E1, L, S, P, I1, T1>& u, const Multipole::stl::SWS::Vector<L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(u, Multipole::stl::SWS::ConstView<L, S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const Multipole::stl::SWS::Vector<L, S, P, I1, T1>& u, const Multipole::stl::SWS::Vector<L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(Multipole::stl::SWS::ConstView<L, S, P, I1, T1>(u), Multipole::stl::SWS::ConstView<L, S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const Multipole::stl::SWS::Vector<L, S, P, I1, T1>& u, const Multipole::stl::SWS::Vector<L, S, P, I2, T2>& v) { return Multipole::stl::SWS::zip(Multipole::stl::SWS::ConstView<L, S, P, I1, T1>(u), Multipole::stl::SWS::ConstView<L, S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }

template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar alpha, const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v) { return Multipole::stl::SWS::map(v, [=](const auto& val) { return alpha * val; }); }
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v, const Scalar alpha) { return Multipole::stl::SWS::map(v, [=](const auto& val) { return alpha * val; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar alpha, const Multipole::stl::SWS::Vector<L, S, P, I, T>& v) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [=](const auto& val) { return alpha * val; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator*(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v, const Scalar alpha) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [=](const auto& val) { return alpha * val; }); }

template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator/(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v, const Scalar alpha) { return Multipole::stl::SWS::map(v, [=](const auto& val) { return val / alpha; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T, typename Scalar> auto operator/(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v, const Scalar alpha) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [=](const auto& val) { return val / alpha; }); }

// FIXME: should not be restricted to float, when using Scalar, SWS::ConstExpressions will also match
template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator-(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v, const float alpha) { return Multipole::stl::SWS::map(v, [=](const auto& val) { return val - alpha; }); }
template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto operator-(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v, const float alpha) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [=](const auto& val) { return val - alpha; }); }

namespace Multipole
{
    namespace stl
    {
        namespace SWS
        {
            template <typename E, std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto real_cast(const Multipole::stl::SWS::ConstExpression<E, L, S, P, I, T>& v) { return Multipole::stl::SWS::map(v, [](auto&& val) { return val.real(); }); }
            template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto real_cast(const Multipole::stl::SWS::Vector<L, S, P, I, T>& v) { return Multipole::stl::SWS::map(Multipole::stl::SWS::ConstView<L, S, P, I, T>(v), [](auto&& val) { return val.real(); }); }

            template <typename E> auto edth(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Edth<E, Spin::Up>(u); }
            template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto edth(const Vector<L, S, P, I, T>& u) { return edth(ConstView<L, S, P, I, T>(u)); }

            template <typename E> auto edth_bar(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Edth<E, Spin::Down>(u); }
            template <std::size_t L, std::size_t S, Multipole::stl::Parity P, typename I, typename T> auto edth_bar(const Vector<L, S, P, I, T>& u) { return edth_bar(ConstView<L, S, P, I, T>(u)); }

        } // namespace SWS

    } // namespace stl

} // namespace Multipole
