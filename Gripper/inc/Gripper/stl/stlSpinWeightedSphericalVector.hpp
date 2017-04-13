#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>
#include <Gripper/stl/stlParity.hpp>                        // math::sws::parity
#include <Gripper/stl/stlSpinWeightedSphericalIndex.hpp>    // math::sws::index
#include <Gripper/stl/stlSpinWeightedSphericalExtent.hpp>   // math::sws::extent
//#include <Gripper/stl/stlIndexableIterator.hpp>             // indexable_iterator

// Custom C++ includes
#include <Gripper/stl/stlArithmeticProgression.hpp>           // stl::arithmetic_progression_iterator

// Standard C++ includes
#include <cstddef>                  // std::size_t
#include <functional>               // std::reference_wrapper, std::ref, std::cref
#include <vector>                   // std::vector
#include <type_traits>              // std::result_of
#include <complex>                  // std::complex


namespace math
{
    namespace sws
    {
            // Forward decalaration of l_parity helper
            template <typename E> auto l_parity(const typename E::extent_type&); // FIXME: this forward declaration should not really exist

            /// <summary>Traits class consisting of type aliases describing spin-weighted spherical expansion coefficient vectors.</summary>
            ///
            template <std::size_t S, parity P, typename IT, typename VT>
            struct vector_traits
            {
                // Common type aliases

                using value_type = VT;

                // STL type aliases

                using container_type = std::vector<value_type>;
                using size_type = typename container_type::size_type;

                // ConstIterator aliases

                using const_iterator_type = typename container_type::const_iterator;
                using const_reverse_iterator_type = typename container_type::const_reverse_iterator;

                // Iterator aliases

                using iterator_type = typename container_type::iterator;
                using reverse_iterator_type = typename container_type::reverse_iterator;

                // Lattice type aliases

                using extent_type = extent<IT>;
                using index_type = typename extent_type::index_type;
                using index_internal_type = typename index_type::value_type;

                // Lattice static members

                static const index_internal_type s = static_cast<index_internal_type>(S);
                static const parity parity = P;
            };


            /// <summary>Read-only Expression Template base class of spin-weighted spherical expansion coefficient vectors to be implemented statically.</summary>
            ///
            template <typename ET, std::size_t S, parity P, typename IT, typename VT>
            struct const_expression : public vector_traits<S, P, IT, VT>
            {
                // ConstExpression type aliases

                using expression_type = ET;

                // Common type aliases

                using typename vector_traits<S, P, IT, VT>::value_type;

                // STL type aliases

                using typename vector_traits<S, P, IT, VT>::container_type;
                using typename vector_traits<S, P, IT, VT>::size_type;

                // ConstIterator aliases

                using const_iterator_type = typename vector_traits<S, P, IT, VT>::const_iterator_type;
                using const_reverse_iterator_type = typename vector_traits<S, P, IT, VT>::const_reverse_iterator_type;

                // Lattice type aliases

                using typename vector_traits<S, P, IT, VT>::extent_type;
                using typename vector_traits<S, P, IT, VT>::index_type;
                using typename vector_traits<S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using vector_traits<S, P, IT, VT>::s;
                using vector_traits<S, P, IT, VT>::parity;

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return static_cast<const expression_type&>(*this).cbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return static_cast<const expression_type&>(*this).cend(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return static_cast<const expression_type&>(*this).crbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return static_cast<const expression_type&>(*this).crend(); }

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
            template <typename ET, std::size_t S, parity P, typename IT, typename VT>
            struct expression : public const_expression<ET, S, P, IT, VT>
            {
                // Expression type aliases

                using expression_type = ET;

                // Common type aliases

                using typename vector_traits<S, P, IT, VT>::value_type;

                // STL type aliases

                using typename vector_traits<S, P, IT, VT>::container_type;
                using typename vector_traits<S, P, IT, VT>::size_type;

                // Iterator aliases

                using iterator_type = typename expression_type::iterator_type;
                using reverse_iterator_type = typename expression_type::reverse_iterator_type;

                // Lattice type aliases

                using typename vector_traits<S, P, IT, VT>::extent_type;
                using typename vector_traits<S, P, IT, VT>::index_type;
                using typename vector_traits<S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using vector_traits<S, P, IT, VT>::s;
                using vector_traits<S, P, IT, VT>::parity;

                // Expression interface

                operator expression_type&() { return static_cast<expression_type&>(*this); }

                // STL interface

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                ///
                value_type& operator[](size_type i) { return static_cast<expression_type&>(*this)[i]; }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                ///
                value_type& at(size_type i) { return static_cast<expression_type&>(*this).at(i); }

                // Iterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                iterator_type begin()             { return static_cast<const expression_type&>(*this).begin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                iterator_type end()               { return static_cast<const expression_type&>(*this).end(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                reverse_iterator_type rbegin()    { return static_cast<const expression_type&>(*this).rbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                reverse_iterator_type rend()      { return static_cast<const expression_type&>(*this).rend(); }

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
            template <std::size_t S, parity P, typename IT, typename VT>
            class vector : public vector_traits<S, P, IT, VT>
            {
            public:

                // Common type aliases

                using typename vector_traits<S, P, IT, VT>::value_type;

                // STL type aliases

                using typename vector_traits<S, P, IT, VT>::container_type;
                using typename vector_traits<S, P, IT, VT>::size_type;

                // ConstIterator aliases

                using typename vector_traits<S, P, IT, VT>::const_iterator_type;
                using typename vector_traits<S, P, IT, VT>::const_reverse_iterator_type;

                // Iterator aliases

                using typename vector_traits<S, P, IT, VT>::iterator_type;
                using typename vector_traits<S, P, IT, VT>::reverse_iterator_type;

                // Lattice type aliases

                using typename vector_traits<S, P, IT, VT>::extent_type;
                using typename vector_traits<S, P, IT, VT>::index_type;
                using typename vector_traits<S, P, IT, VT>::index_internal_type;

                // Lattice static members

                using vector_traits<S, P, IT, VT>::s;
                using vector_traits<S, P, IT, VT>::parity;

                static auto l_parity(const extent_type& ext) { return math::sws::l_parity<expression<vector<S, P, IT, VT>, S, P, IT, VT>>(ext); }

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                vector() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                vector(const vector& in) = default;

                /// <summary>Default move constructor.</summary>
                ///
                vector(vector&& in) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~vector() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                vector& operator=(const vector&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                vector& operator=(vector&&) = default;

                ///<summary>Constructs a series expansion vector of <c>ext</c> size.</summary>
                ///
                //vector(const extent_type& ext) : m_extent(ext), m_data() { m_data.resize(distance(ext.initial(), ext.final())); }

                ///<summary>Constructs a series expansion vector of <c>ext</c> size.</summary>
                ///
                vector(index_internal_type L) : m_extent({ 0, 0, s }, { L, L, s }), m_data() { m_data.resize(distance(m_extent.initial(), m_extent.final()) + 1); }

                ///<summary>Constructs a series expansion vector from the expression <c>expr</c>.</summary>
                ///
                template <typename ConstVecExpr, std::size_t SS, typename IndexType, typename ValueType>
                vector(const const_expression<ConstVecExpr, SS, P, IndexType, ValueType>& expr)
                {
                    static_assert(s == expr.s, "Spin weights of variables differ on the two sides of assignment.");

                    serial_evaluator(expr);
                }

                ///<summary>Constructs a series expansion vector from the expression <c>expr</c>.</summary>
                ///
                template <typename ConstVecExpr, std::size_t SS, typename IndexType, typename ValueType>
                vector& operator=(const const_expression<ConstVecExpr, SS, P, IndexType, ValueType>& expr)
                {
                    static_assert(s == expr.s, "Spin weights of variables differ on the two sides of assignment.");

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return m_data.cbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return m_data.cend(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return m_data.crbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return m_data.crend(); }

                // Iterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                iterator_type begin() { return m_data.begin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                iterator_type end() { return m_data.end(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                reverse_iterator_type rbegin() { return m_data.rbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                reverse_iterator_type rend() { return m_data.rend(); }

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
                template <typename ConstVecExpr, std::size_t SS, typename IndexType, typename ValueType>
                void serial_evaluator(const const_expression<ConstVecExpr, SS, P, IndexType, ValueType>& expr)
                {
                    // Extract type from encapsulating expression
                    const ConstVecExpr& v = expr;

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    // After extents have been set to match, we can index with container_type::(const_)iterator
                    std::copy(v.cbegin(), v.cend(), this->begin());

                    // After resize, it's safe to utilize container_type::size_type for indexing.
                    //std::for_each(stl::arithmetic_progression_iterator<size_type>(convert(m_extent.initial())),
                    //              stl::arithmetic_progression_iterator<size_type>(convert(m_extent.final())),
                    //              [this, &v](const size_type& i) mutable
                    //{
                    //    this->at(i) = v.at(i);
                    //});
                }

                extent_type m_extent;
                container_type m_data;
            };


            /// <summary>Expression Template providing read-only access with reference semantics to the elements of a series expansion with storage.</summary>
            ///
            template <std::size_t S, parity P, typename IT, typename VT>
            class const_view : public const_expression<const_view<S, P, IT, VT>, S, P, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = const_expression<const_view<S, P, IT, VT>, S, P, IT, VT>;
                using trait_type = vector_traits<S, P, IT, VT>;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::container_type;
                using typename expression_type::size_type;

                // ConstIterator aliases

                using const_iterator_type = typename vector_traits<S, P, IT, VT>::const_iterator_type;
                using const_reverse_iterator_type = typename vector_traits<S, P, IT, VT>::const_reverse_iterator_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::s;
                using expression_type::parity;

            private:

                using vector_type = vector<s, parity, index_internal_type, value_type>;
                using reference_type = std::reference_wrapper<const vector_type>;

            public:

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
                ///
                const_view() = delete;

                /// <summary>Default copy constructor.</summary>
                ///
                const_view(const const_view&) = default;

                /// <summary>Default move constructor.</summary>
                /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
                ///
                const_view(const_view&&) = default;//delete;

                /// <summary>Default destructor.</summary>
                ///
                ~const_view() = default;

                /// <summary>Default copy assign operator.</summary>
                ///
                const_view& operator=(const const_view&) = default;

                /// <summary>Default move assign operator.</summary>
                /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
                ///
                const_view& operator=(const_view&&) = delete;

                /// <summary>Constructs a <c>View</c> from a <c>Vector</c>.</summary>
                ///
                const_view(const vector_type& v) : _v(std::cref(v)) {}

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return _v.get().cbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return _v.get().cend(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return _v.get().crbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return _v.get().crend(); }

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
            template <std::size_t S, parity P, typename IT, typename VT>
            class view : public expression<view<S, P, IT, VT>, S, P, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = expression<view<S, P, IT, VT>, S, P, IT, VT>;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::container_type;
                using typename expression_type::size_type;

                // ConstIterator aliases

                using typename expression_type::const_iterator_type;
                using typename expression_type::const_reverse_iterator_type;

                // Iterator aliases

                using typename expression_type::iterator_type;
                using typename expression_type::reverse_iterator_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::s;
                using expression_type::parity;

            private:

                using vector_type = vector<s, parity, index_internal_type, value_type>;
                using reference_type = std::reference_wrapper<vector_type>;

            public:

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
                ///
                view() = delete;

                /// <summary>Default copy constructor.</summary>
                ///
                view(const view&) = default;

                /// <summary>Default move constructor.</summary>
                /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
                ///
                view(view&&) = default;//delete;

                /// <summary>Default destructor.</summary>
                ///
                ~view() = default;

                /// <summary>Default copy assign operator.</summary>
                ///
                view& operator=(const view&) = default;

                /// <summary>Default move assign operator.</summary>
                /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
                ///
                view& operator=(view&&) = default;

                /// <summary>Constructs a <c>view</c> from a <c>Vector</c>.</summary>
                ///
                explicit view(vector_type& v) : _v(std::ref(v)) {}

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return _v.get().cbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return _v.get().cend(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return _v.get().crbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return _v.get().crend(); }

                // Iterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                iterator_type begin() { return _v.get().begin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                iterator_type end() { return _v.get().end(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                iterator_type rbegin() { return _v.get().rbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                iterator_type rend() { return _v.get().rend(); }

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
            class id : public const_expression<id<E>, E::s, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // Expression type aliases

                using expression_type = const_expression<id<E>, E::s, E::parity, typename E::index_internal_type, typename E::value_type>;
                using outer_expression_type = E;

                // Common type aliases

                using typename expression_type::value_type;

                // STL type aliases

                using typename expression_type::size_type;

                // ConstIterator aliases

                using typename expression_type::const_iterator_type;
                using typename expression_type::const_reverse_iterator_type;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::s;
                using expression_type::parity;

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs an <c>id</c> from a <c>const_epxression</c>.</summary>
                ///
                id(const const_expression<E, E::s, E::parity, typename E::index_internal_type, typename E::value_type>& u)
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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return _u.cbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return _u.cend(); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return _u.crbegin(); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return _u.crend(); }

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
            class func : public const_expression<func<E, F>, E::s, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // ConstIterator forward-declarations

                class bidir_const_iter;
                class reverse_bidir_const_iter;

                // Expression type aliases

                using expression_type = const_expression<func<E, F>, E::s, E::parity, typename E::index_internal_type, typename E::value_type>;

                // Common type alises

                using typename expression_type::value_type;

                // STL type alises

                using typename expression_type::size_type;

                // ConstIterator aliases

                using const_iterator_type = bidir_const_iter;
                using const_reverse_iterator_type = reverse_bidir_const_iter;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::s;
                using expression_type::parity;

                // ConstIterator member classes

            protected:

                using iter_base = std::iterator<std::bidirectional_iterator_tag,
                                                typename E::value_type,
                                                std::ptrdiff_t,
                                                const typename E::value_type*,
                                                const typename E::value_type&>;

            public:

                /// <summary>Iterator class that invokes a function object when dereferencing.</summary>
                ///
                class bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using index_iter = typename extent_type::const_iterator_type;
                    using index_type = typename index_iter::value_type;
                    using func_type = F;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
                    bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
                    bidir_const_iter(const bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
                    bidir_const_iter(bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(const bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
                    bidir_const_iter(const index_type& index, const F& f) : _it(index), _f(f) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return _f(*_it); }

                    /// <summary>Prefix increment operator.</summary>
                    ///
                    bidir_const_iter& operator++()
                    {
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const bidir_const_iter& rhs)
                    {
                        return _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const bidir_const_iter& rhs)
                    {
                        return _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const { return &_f(*_it); }

                    /// <summary>Postfix increment operator.</summary>
                    ///
                    bidir_const_iter operator++(int)
                    {
                        bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
                    bidir_const_iter& operator--()
                    {
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
                    bidir_const_iter operator--(int)
                    {
                        bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    index_iter _it;
                    func_type _f;
                };

                /// <summary>Iterator class that invokes a function object when dereferencing.</summary>
                ///
                class reverse_bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using index_iter = typename extent_type::reverse_const_iterator_type;
                    using index_type = typename index_iter::value_type;
                    using func_type = F;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
                    reverse_bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
                    reverse_bidir_const_iter(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
                    reverse_bidir_const_iter(reverse_bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~reverse_bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator=(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator=(reverse_bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
                    reverse_bidir_const_iter(const index_type& index, const F& f) : _it(index), _f(f) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return _f(*_it); }

                    /// <summary>Prefix increment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator++()
                    {
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const reverse_bidir_const_iter& rhs)
                    {
                        return _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const reverse_bidir_const_iter& rhs)
                    {
                        return _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const { return &_f(*_it); }

                    /// <summary>Postfix increment operator.</summary>
                    ///
                    reverse_bidir_const_iter operator++(int)
                    {
                        reverse_bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator--()
                    {
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
                    reverse_bidir_const_iter operator--(int)
                    {
                        reverse_bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    index_iter _it;
                    func_type _f;
                };


                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>func</c> from an extent <c>ext</c> and a function object <c>f</c>.</summary>
                ///
                func(const extent_type ext, const F f) : m_ext(ext), m_f(f) {}

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return const_iterator_type(_ext.cbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return const_iterator_type(_ext.cend(), _f); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return const_reverse_iterator_type(_ext.crbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return const_reverse_iterator_type(_ext.crend(), _f); }

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
            class map : public const_expression<map<E, F>, E::s, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>
            {
            public:

				// ConstIterator forward-declarations

				class bidir_const_iter;
				class reverse_bidir_const_iter;

				// Expression type aliases

				using expression_type = const_expression<map<E, F>, E::s, E::parity, typename E::index_internal_type, typename std::result_of<F(typename E::value_type)>::type>;
				using outer_expression_type = E;

				// Common type aliases

				using typename expression_type::value_type;

				// STL type aliases

				using typename expression_type::size_type;

				// ConstIterator aliases

				using const_iterator_type = bidir_const_iter;
				using const_reverse_iterator_type = reverse_bidir_const_iter;

				// Lattice type aliases

				using typename expression_type::extent_type;
				using typename expression_type::index_type;
				using typename expression_type::index_internal_type;

				// Lattice static members

				using expression_type::s;
				using expression_type::parity;

                // ConstIterator member classes

            protected:

                using iter_base = std::iterator<std::bidirectional_iterator_tag,
                                               typename std::result_of<F(typename E::value_type)>::type,
                                               std::ptrdiff_t,
                                               const typename std::result_of<F(typename E::value_type)>::type*,
                                               const typename std::result_of<F(typename E::value_type)>::type&>;

            public:

                /// <summary>Iterator class that invokes a function object when dereferencing.</summary>
                ///
                class bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using iter_type = typename expression_type::const_iterator_type;
                    using func_type = F;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
                    bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
                    bidir_const_iter(const bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
                    bidir_const_iter(bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(const bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
                    bidir_const_iter(const iter_type& it, const F& f) : _it(it), _f(f) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return _f(*_it); }

                    /// <summary>Prefix increment operator.</summary>
                    ///
                    bidir_const_iter& operator++()
                    {
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const bidir_const_iter& rhs)
                    {
                        return _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const bidir_const_iter& rhs)
                    {
                        return _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const { return &_f(*_it); }

                    /// <summary>Postfix increment operator.</summary>
                    ///
                    bidir_const_iter operator++(int)
                    {
                        bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
                    bidir_const_iter& operator--()
                    {
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
                    bidir_const_iter operator--(int)
                    {
                        bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    iter_type _it;
                    const func_type _f;
                };

				/// <summary>Iterator class that invokes a function object when dereferencing.</summary>
				///
                class reverse_bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

					using iter_type = typename expression_type::const_reverse_iterator_type;
                    using func_type = F;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
					reverse_bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
					reverse_bidir_const_iter(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
					reverse_bidir_const_iter(reverse_bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~reverse_bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
					reverse_bidir_const_iter& operator=(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
					reverse_bidir_const_iter& operator=(reverse_bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
					reverse_bidir_const_iter(const iter_type& it, const F& f) : _it(it), _f(f) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return _f(*_it); }

                    /// <summary>Prefix increment operator.</summary>
                    ///
					reverse_bidir_const_iter& operator++()
                    {
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const bidir_const_iter& rhs)
                    {
                        return _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const bidir_const_iter& rhs)
                    {
                        return _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const { return &_f(*_it); }

                    /// <summary>Postfix increment operator.</summary>
                    ///
					reverse_bidir_const_iter operator++(int)
                    {
                        bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
					reverse_bidir_const_iter& operator--()
                    {
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
					reverse_bidir_const_iter operator--(int)
                    {
                        bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    iter_type _it;
                    const func_type _f;
                };

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Map</c> from an expression <c>u</c> and a function object <c>f</c>.</summary>
                ///
                map(const const_expression<E, E::s, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F f)
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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return const_iterator_type(_u.cbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return const_iterator_type(_u.cend(), _f); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return const_reverse_iterator_type(_u.crbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return const_reverse_iterator_type(_u.crend(), _f); }

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
            class zip : public const_expression<zip<E1, E2, F>, E1::s, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>
            {
            public:

				// ConstIterator forward-declarations

				class bidir_const_iter;
				class reverse_bidir_const_iter;

				// Expression type aliases

				using expression_type = const_expression<zip<E1, E2, F>, E1::s, E1::parity, typename E1::index_internal_type, typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type>;
				using outer_expression_type1 = E1;
				using outer_expression_type2 = E2;

				// Common typedefs

				using typename expression_type::value_type;

				// STL typedefs

				using typename expression_type::size_type;

				// ConstIterator aliases

				using const_iterator_type = bidir_const_iter;
				using const_reverse_iterator_type = reverse_bidir_const_iter;

				// Lattice typedefs

				using typename expression_type::extent_type;
				using typename expression_type::index_type;
				using typename expression_type::index_internal_type;

				// Lattice static members

				using expression_type::s;
				using expression_type::parity;

				// ConstIterator member classes

            protected:

                using iter_base = std::iterator<std::bidirectional_iterator_tag,
                                                typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type,
                                                std::ptrdiff_t,
                                                const typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type*,
                                                const typename std::result_of<F(typename E1::value_type, typename E2::value_type)>::type&>;

            public:

				/// <summary>Iterator class that invokes a function object when dereferencing.</summary>
				///
				class bidir_const_iter : public iter_base
				{
				public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

					using iter_type1 = typename E1::expression_type::const_iterator_type;
					using iter_type2 = typename E2::expression_type::const_iterator_type;
					using func_type = F;

                public:

					// Constructors / Destructors / Assignment operators

					/// <summary>Default constructor.</summary>
					/// <remarks>Default constructed objects are in an invalid state.</remarks>
					///
					bidir_const_iter() = default;

					/// <summary>Default copy constructor.</summary>
					///
					bidir_const_iter(const bidir_const_iter&) = default;

					/// <summary>Default move constructor.</summary>
					///
					bidir_const_iter(bidir_const_iter&&) = default;

					/// <summary>Default destructor.</summary>
					///
					~bidir_const_iter() = default;

					/// <summary>Default copy assignment operator.</summary>
					///
					bidir_const_iter& operator=(const bidir_const_iter&) = default;

					/// <summary>Default move assignment operator.</summary>
					///
					bidir_const_iter& operator=(bidir_const_iter&&) = default;

					/// <summary>Index-functor pair constructor.</summary>
					///
					bidir_const_iter(const iter_type1& it1, const iter_type2& it2, const F& f) : _it1(it1), _it2(it2), _f(f) {}

					// Iterator concept

					/// <summary>Dereference operator.</summary>
					///
					const value_type operator*() const { return _f(*_it1, *_it2); }

					/// <summary>Prefix increment operator.</summary>
					///
					bidir_const_iter& operator++()
					{
						++_it1;
                        ++_it2;

						return *this;
					}

					// InputIterator concept

					/// <summary>Equality operator.</summary>
					/// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
					///
					inline bool operator==(const bidir_const_iter& rhs)
					{
						return (_it1 == rhs._it1) &&
                               (_it2 == rhs._it2);
					}

					/// <summary>Unequality operator.</summary>
					/// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
					///
					inline bool operator!=(const bidir_const_iter& rhs)
					{
						return (_it1 != rhs._it1) ||
                               (_it2 != rhs._it2);
					}

					/// <summary>Arrow operator.</summary>
					///
					const pointer operator->() const { return &_f(*_it1, *_it2); }

					/// <summary>Postfix increment operator.</summary>
					///
					bidir_const_iter operator++(int)
					{
						bidir_const_iter tmp = *this;

						++*this;

						return tmp;
					}

					// BidirectionalIterator concept

					/// <summary>Prefix decrement operator.</summary>
					///
					bidir_const_iter& operator--()
					{
						--_it1;
                        --_it2;

						return *this;
					}

					/// <summary>Postfix decrement operator.</summary>
					///
					bidir_const_iter operator--(int)
					{
						bidir_const_iter tmp = *this;

						--*this;

						return tmp;
					}

				private:

					iter_type1 _it1;
					iter_type2 _it2;
					const func_type _f;
				};

				/// <summary>Iterator class that invokes a function object when dereferencing.</summary>
				///
				class reverse_bidir_const_iter : public iter_base
				{
				public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using iter_type1 = typename E1::expression_type::const_reverse_iterator_type;
                    using iter_type2 = typename E2::expression_type::const_reverse_iterator_type;
					using func_type = F;

                public:

					// Constructors / Destructors / Assignment operators

					/// <summary>Default constructor.</summary>
					/// <remarks>Default constructed objects are in an invalid state.</remarks>
					///
					reverse_bidir_const_iter() = default;

					/// <summary>Default copy constructor.</summary>
					///
					reverse_bidir_const_iter(const reverse_bidir_const_iter&) = default;

					/// <summary>Default move constructor.</summary>
					///
					reverse_bidir_const_iter(reverse_bidir_const_iter&&) = default;

					/// <summary>Default destructor.</summary>
					///
					~reverse_bidir_const_iter() = default;

					/// <summary>Default copy assignment operator.</summary>
					///
					reverse_bidir_const_iter& operator=(const reverse_bidir_const_iter&) = default;

					/// <summary>Default move assignment operator.</summary>
					///
					reverse_bidir_const_iter& operator=(reverse_bidir_const_iter&&) = default;

					/// <summary>Index-functor pair constructor.</summary>
					///
					reverse_bidir_const_iter(const iter_type1& it1, const iter_type2& it2, const F& f) : _it1(it1), _it2(it2), _f(f) {}

					// Iterator concept

					/// <summary>Dereference operator.</summary>
					///
					const reference operator*() const { return _f(*_it1, *_it2); }

					/// <summary>Prefix increment operator.</summary>
					///
					reverse_bidir_const_iter& operator++()
					{
						++_it1;
                        ++_it2;

						return *this;
					}

					// InputIterator concept

					/// <summary>Equality operator.</summary>
					/// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
					///
					inline bool operator==(const bidir_const_iter& rhs)
					{
						return (_it1 == rhs._it1) &&
                               (_it2 == rhs._it2);
					}

					/// <summary>Unequality operator.</summary>
					/// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
					///
					inline bool operator!=(const bidir_const_iter& rhs)
					{
						return (_it1 != rhs._it1) ||
                               (_it2 != rhs._it2);
					}

					/// <summary>Arrow operator.</summary>
					///
					const pointer operator->() const { return &_f(*_it1, *_it2); }

					/// <summary>Postfix increment operator.</summary>
					///
					reverse_bidir_const_iter operator++(int)
					{
						bidir_const_iter tmp = *this;

						++*this;

						return tmp;
					}

					// BidirectionalIterator concept

					/// <summary>Prefix decrement operator.</summary>
					///
					reverse_bidir_const_iter& operator--()
					{
						--_it1;
                        --_it2;

						return *this;
					}

					/// <summary>Postfix decrement operator.</summary>
					///
					reverse_bidir_const_iter operator--(int)
					{
						bidir_const_iter tmp = *this;

						--*this;

						return tmp;
					}

				private:

					iter_type1 _it1;
                    iter_type2 _it2;
					const func_type _f;
				};

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Zip</c> from two possibly varying expressions <c>u</c> and <c>v</c> and a function object <c>f</c>.</summary>
                ///
                zip(const const_expression<E1, E1::s, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u,
                    const const_expression<E2, E2::s, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v,
                    const F f)
                    : _u(static_cast<const E1&>(u))
                    , _v(static_cast<const E2&>(v))
                    , _f(f)
                {
                    //std::cout << "sw::Zip _u = " << _u.extent() << " _v = " << _v.extent() << std::endl;

                    static_assert(E1::s == E2::s, "sws::zip(s value mismatch)");
                    static_assert(E1::parity == E2::parity, "sws::zip(parity value mismatch)");

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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return const_iterator_type(_u.cbegin(), _v.cbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return const_iterator_type(_u.cend(), _v.cend(), _f); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return const_reverse_iterator_type(_u.crbegin(), _v.crbegin(), _f); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return const_reverse_iterator_type(_u.crend(), _v.crend(), _f); }

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
            namespace spin
            {
                struct up {};
                struct down {};
            }

            namespace impl
            {
                /// <summary>Unimplemented helper to workaround unconditional widening of results due to double return types of common math functions with Integral params.</summary>
                ///
                template <typename VT> struct realify;

                template <> struct realify<float> { using type = float; };
                template <> struct realify<double> { using type = double; };
                template <> struct realify<std::complex<float>> { using type = float; };
                template <> struct realify<std::complex<double>> { using type = double; };

                /// <summary>Unimplemented helper of the spin stepping edth operator Exression Template for function partial specialization.</summary>
                ///
                template <typename E, typename S> struct spin_stepper;

                /// <summary>Specialization of the spin stepping helper function for the Expression Template implementing the spin-up operator.</summary>
                /// <note>J. N. Goldberg, A. J. Macfarlane, E. T. Newman, F. Rohrlich, and E. C. G. Sudarshan Spin‐s Spherical Harmonics and ð, eq. (2.7a)</note>
                ///
                template <typename E> struct spin_stepper<E, spin::up>
                {
                    static typename E::index_internal_type s = E::s + 1;

                    static auto factor(const typename E::index_type& i)
                    {
                        using ii_type = typename E::index_internal_type;

                        ii_type result = (i.l - i.s) * (i.l + i.s + static_cast<ii_type>(1));

                        if (result < static_cast<ii_type>(0))
                            return static_cast<ii_type>(0);
                        else
                            return result;
                    }
                    /*
                    static auto at(const E& e, const typename E::index_type& i)
                    {
                        using result_type = typename E::value_type;
                        result_type result = 0;

                        if (i.s == std::min(i.l, E::s_max)) return result;

                        auto factor_sq = static_cast<typename E::index_internal_type>((i.l - i.s) * (i.l + i.s + static_cast<typename E::index_internal_type>(1)));

                        if (factor_sq <= static_cast<typename E::index_internal_type>(0)) return result;

                        auto index = i;
                        
                        result = static_cast<result_type>(static_cast<typename realify<result_type>::type>(std::sqrt(factor_sq)) * e.at(++index));

                        return result;
                    }*/
                };


                /// <summary>Specialization of the spin stepping helper function for the Expression Template implementing the spin-down operator.</summary>
                /// <note>J. N. Goldberg, A. J. Macfarlane, E. T. Newman, F. Rohrlich, and E. C. G. Sudarshan Spin‐s Spherical Harmonics and ð, eq. (2.7b)</note>
                ///
                template <typename E> struct spin_stepper<E, spin::down>
                {
                    static typename E::index_internal_type s = E::s - 1;

                    static auto factor(const typename E::index_type& i)
                    {
                        using ii_type = typename E::index_internal_type;

                        ii_type result = (i.l + i.s) * (i.l - i.s + static_cast<ii_type>(1));

                        if (result < static_cast<ii_type>(0))
                            return static_cast<ii_type>(0);
                        else
                            return result;
                    }
                    /*
                    static auto at(const E& e, const typename E::index_type& i)
                    {
                        using result_type = typename E::value_type;
                        result_type result = 0;

                        if (i.s == std::max(-i.l, -E::s_max)) return result;

                        auto factor_sq = static_cast<typename E::index_internal_type>((i.l + i.s) * (i.l - i.s + static_cast<typename E::index_internal_type>(1)));

                        if (factor_sq <= static_cast<typename E::index_internal_type>(0)) return result;

                        auto index = i;

                        result = static_cast<result_type>(-static_cast<typename realify<result_type>::type>(std::sqrt(factor_sq)) * e.at(--index));

                        return result;
                    }*/
                };

            } // namespace impl


            /// <summary>Expression Template that applies applies the spin stepping operator (edth) to a series expansion.</summary>
            ///
            template <typename E, typename S>
            class edth : public const_expression<edth<E, S>, impl::spin_stepper<E, S>::s, E::parity, typename E::index_internal_type, typename E::value_type>
            {
            public:

                // ConstIterator forward-declarations

                class bidir_const_iter;
                class reverse_bidir_const_iter;

                // Expression type aliases 

                using expression_type = const_expression<edth<E, S>, impl::spin_stepper<E, S>::s, E::parity, typename E::index_internal_type, typename E::value_type>;
                using outer_expression_type = E;

                // Common typedefs

                using typename expression_type::value_type;

                // STL typedefs

                using typename expression_type::size_type;

                // ConstIterator aliases

                using const_iterator_type = bidir_const_iter;
                using const_reverse_iterator_type = reverse_bidir_const_iter;

                // Lattice typedefs

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // Lattice static members

                using expression_type::s;
                using expression_type::parity;

                // ConstIterator member classes

            protected:

                using iter_base = std::iterator<std::bidirectional_iterator_tag,
                                                value_type,
                                                std::ptrdiff_t,
                                                value_type*,
                                                value_type&>;

            public:

                /// <summary>Iterator class that invokes a function object when dereferencing.</summary>  
                ///
                class bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using index_iter_type = typename extent_type::const_iterator_type;
                    using iter_type = typename expression_type::const_iterator_type;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
                    bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
                    bidir_const_iter(const bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
                    bidir_const_iter(bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(const bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
                    bidir_const_iter& operator=(bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
                    bidir_const_iter(const index_iter_type& idx, const iter_type& it) : _idx(idx), _it(it) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return std::sqrt(impl::spin_stepper<E, S>::factor(static_cast<typename impl::realify<value_type>::type>(*idx))) * (*it); }

                    /// <summary>Prefix increment operator.</summary>
                    ///
                    bidir_const_iter& operator++()
                    {
                        ++_idx;
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const bidir_const_iter& rhs)
                    {
                        return _idx == rhs._idx && _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const bidir_const_iter& rhs)
                    {
                        return _idx != rhs._idx || _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const
                    {
                        static_assert(false, "This function prototype is invalid to instantiate. Returning pointer to temporary has no meaning.");

                        return nullptr;
                    }

                    /// <summary>Postfix increment operator.</summary>
                    ///
                    bidir_const_iter operator++(int)
                    {
                        bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
                    bidir_const_iter& operator--()
                    {
                        --_idx;
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
                    bidir_const_iter operator--(int)
                    {
                        bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    index_iter_type _idx;
                    iter_type _it;
                };

                /// <summary>Iterator class that invokes a function object when dereferencing.</summary>
                ///
                class reverse_bidir_const_iter : public iter_base
                {
                public:

                    // Iterator aliases

                    using typename iter_base::iterator_category;
                    using typename iter_base::value_type;
                    using typename iter_base::difference_type;
                    using typename iter_base::pointer;
                    using typename iter_base::reference;

                private:

                    using index_iter_type = typename extent_type::const_reverse_iterator_type;
                    using iter_type = typename expression_type::const_reverse_iterator_type;

                public:

                    // Constructors / Destructors / Assignment operators

                    /// <summary>Default constructor.</summary>
                    /// <remarks>Default constructed objects are in an invalid state.</remarks>
                    ///
                    reverse_bidir_const_iter() = default;

                    /// <summary>Default copy constructor.</summary>
                    ///
                    reverse_bidir_const_iter(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move constructor.</summary>
                    ///
                    reverse_bidir_const_iter(reverse_bidir_const_iter&&) = default;

                    /// <summary>Default destructor.</summary>
                    ///
                    ~reverse_bidir_const_iter() = default;

                    /// <summary>Default copy assignment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator=(const reverse_bidir_const_iter&) = default;

                    /// <summary>Default move assignment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator=(reverse_bidir_const_iter&&) = default;

                    /// <summary>Index-functor pair constructor.</summary>
                    ///
                    reverse_bidir_const_iter(const index_iter_type& idx, const iter_type& it) : _idx(idx), _it(it) {}

                    // Iterator concept

                    /// <summary>Dereference operator.</summary>
                    ///
                    const value_type operator*() const { return std::sqrt(impl::spin_stepper<E, S>::factor(static_cast<typename impl::realify<value_type>::type>(*idx))) * (*it);
                    }

                    /// <summary>Prefix increment operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator++()
                    {
                        ++_idx;
                        ++_it;

                        return *this;
                    }

                    // InputIterator concept

                    /// <summary>Equality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator==(const bidir_const_iter& rhs)
                    {
                        return _idx == rhs._idx && _it == rhs._it;
                    }

                    /// <summary>Unequality operator.</summary>
                    /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                    ///
                    inline bool operator!=(const bidir_const_iter& rhs)
                    {
                        return _idx != rhs._idx || _it != rhs._it;
                    }

                    /// <summary>Arrow operator.</summary>
                    ///
                    const pointer operator->() const
                    {
                        static_assert(false, "This function prototype is invalid to instantiate. Returning pointer to temporary has no meaning.");

                        return nullptr;
                    }

                    /// <summary>Postfix increment operator.</summary>
                    ///
                    reverse_bidir_const_iter operator++(int)
                    {
                        bidir_const_iter tmp = *this;

                        ++*this;

                        return tmp;
                    }

                    // BidirectionalIterator concept

                    /// <summary>Prefix decrement operator.</summary>
                    ///
                    reverse_bidir_const_iter& operator--()
                    {
                        --_idx;
                        --_it;

                        return *this;
                    }

                    /// <summary>Postfix decrement operator.</summary>
                    ///
                    reverse_bidir_const_iter operator--(int)
                    {
                        bidir_const_iter tmp = *this;

                        --*this;

                        return tmp;
                    }

                private:

                    index_iter_type _idx;
                    iter_type _it;
                };

                // Constructors / Destructors / Assignment operators

                /// <summary>Constructs a <c>Zip</c> from two possibly varying expressions <c>u</c> and <c>v</c> and a function object <c>f</c>.</summary>
                ///
                edth(const expression_type& u)
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

                // ConstIterator interface

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_iterator_type cbegin()             const { return const_iterator_type(_u.extent().cbegin(), _u.cbegin()); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_iterator_type cend()               const { return const_iterator_type(_u.extent().cend(), _u.cend()); }

                /// <summary>Returns a constant iterator to the beginning of the array.</summary>
                ///
                const_reverse_iterator_type crbegin()    const { return const_reverse_iterator_type(_u.extent().crbegin(), _u.crbegin()); }

                /// <summary>Returns a constant iterator to the one-past-the-end of the array.</summary>
                ///
                const_reverse_iterator_type crend()      const { return const_reverse_iterator_type(_u.extent().crend(), _u.crend()); }

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
            auto make_map(const const_expression<E, E::s, E::parity, typename E::index_internal_type, typename E::value_type>& u, const F f)
            {
                return map<E, F>(u, f);
            }

            template <typename E1, typename E2, typename F>
            auto make_zip(const const_expression<E1, E1::s, E1::parity, typename E1::index_internal_type, typename E1::value_type>& u,
                          const const_expression<E2, E2::s, E2::parity, typename E2::index_internal_type, typename E2::value_type>& v,
                          const F f)
            {
                return zip<E1, E2, F>(u, v, f);
            }

            //template <Parity P, typename E>
            //auto parity(const ConstExpression<E, E::l_max, E::s_max, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return Par<E, P>(u); }

            namespace impl
            {
                template <typename E, typename F>
                auto l_parity(const typename E::extent_type& ext, const F& f) { return func<E, F>(ext, f); }
            }

            template <typename E>
            auto l_parity(const typename E::extent_type& ext) { return impl::l_parity<E>(ext, [](const typename E::index_type& i) { return i.l % 2 ? -1 : 1; }); }

    } // namespace sws

} // namespace math

////////////////////////////////////////////////////////
// math::sws::vector non-member operators //
////////////////////////////////////////////////////////

// Unary 
template <typename E, std::size_t S, math::sws::parity P, typename I, typename T> auto operator+(const math::sws::const_expression<E, S, P, I, T>& v) { return math::sws::make_map(v, [](const auto& val) { return static_cast<T>(1) * val; }); }
template <typename E, std::size_t S, math::sws::parity P, typename I, typename T> auto operator-(const math::sws::const_expression<E, S, P, I, T>& v) { return math::sws::make_map(v, [](const auto& val) { return static_cast<T>(-1) * val; }); }
template <typename E, std::size_t S, math::sws::parity P, typename I, typename T> auto conjugate(const math::sws::const_expression<E, S, P, I, T>& v) { return math::sws::make_map(v, [](const auto& val) { return std::conj(val); }); }
template <std::size_t S, math::sws::parity P, typename I, typename T> auto operator+(const math::sws::vector<S, P, I, T>& v) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [](const auto& val) { return static_cast<T>(1) * val; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T> auto operator-(const math::sws::vector<S, P, I, T>& v) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [](const auto& val) { return static_cast<T>(-1) * val; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T> auto conjugate(const math::sws::vector<S, P, I, T>& v) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [](const auto& val) { return std::conj(val); }); }

// Binary
template <typename E1, typename E2, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const math::sws::const_expression<E1, S, P, I1, T1>& u, const math::sws::const_expression<E2, S, P, I2, T2>& v) { return math::sws::make_zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const math::sws::const_expression<E1, S, P, I1, T1>& u, const math::sws::const_expression<E2, S, P, I2, T2>& v) { return math::sws::make_zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E2, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const math::sws::vector<S, P, I1, T1>& u, const math::sws::const_expression<E2, S, P, I2, T2>& v) { return math::sws::make_zip(math::sws::const_view<S, P, I1, T1>(u), v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E2, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const math::sws::vector<S, P, I1, T1>& u, const math::sws::const_expression<E2, S, P, I2, T2>& v) { return math::sws::make_zip(math::sws::const_view<S, P, I1, T1>(u), v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E1, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const math::sws::const_expression<E1, S, P, I1, T1>& u, const math::sws::vector<S, P, I2, T2>& v) { return math::sws::make_zip(u, math::sws::const_view<S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const math::sws::const_expression<E1, S, P, I1, T1>& u, const math::sws::vector<S, P, I2, T2>& v) { return math::sws::make_zip(u, math::sws::const_view<S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator+(const math::sws::vector<S, P, I1, T1>& u, const math::sws::vector<S, P, I2, T2>& v) { return math::sws::make_zip(math::sws::const_view<S, P, I1, T1>(u), math::sws::const_view<S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <std::size_t S, math::sws::parity P, typename I1, typename I2, typename T1, typename T2> auto operator-(const math::sws::vector<S, P, I1, T1>& u, const math::sws::vector<S, P, I2, T2>& v) { return math::sws::make_zip(math::sws::const_view<S, P, I1, T1>(u), math::sws::const_view<S, P, I2, T2>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }

template <typename E, std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar alpha, const math::sws::const_expression<E, S, P, I, T>& v) { return math::sws::make_map(v, [=](const auto& val) { return alpha * val; }); }
template <typename E, std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator*(const math::sws::const_expression<E, S, P, I, T>& v, const Scalar alpha) { return math::sws::make_map(v, [=](const auto& val) { return alpha * val; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator*(const Scalar alpha, const math::sws::vector<S, P, I, T>& v) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [=](const auto& val) { return alpha * val; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator*(const math::sws::vector<S, P, I, T>& v, const Scalar alpha) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [=](const auto& val) { return alpha * val; }); }

template <typename E, std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator/(const math::sws::const_expression<E, S, P, I, T>& v, const Scalar alpha) { return math::sws::make_map(v, [=](const auto& val) { return val / alpha; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T, typename Scalar> auto operator/(const math::sws::vector<S, P, I, T>& v, const Scalar alpha) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [=](const auto& val) { return val / alpha; }); }

// FIXME: should not be restricted to float, when using Scalar, SWS::ConstExpressions will also match
template <typename E, std::size_t S, math::sws::parity P, typename I, typename T> auto operator-(const math::sws::const_expression<E, S, P, I, T>& v, const float alpha) { return math::sws::make_map(v, [=](const auto& val) { return val - alpha; }); }
template <std::size_t S, math::sws::parity P, typename I, typename T> auto operator-(const math::sws::vector<S, P, I, T>& v, const float alpha) { return math::sws::make_map(math::sws::const_view<S, P, I, T>(v), [=](const auto& val) { return val - alpha; }); }

namespace math
{
    namespace sws
    {
        template <typename E, std::size_t S, parity P, typename I, typename T> auto real_cast(const const_expression<E, S, P, I, T>& v) { return make_map(v, [](auto&& val) { return val.real(); }); }
        template <std::size_t S, parity P, typename I, typename T> auto real_cast(const vector<S, P, I, T>& v) { returnmake_map(const_view<S, P, I, T>(v), [](auto&& val) { return val.real(); }); }

        template <typename E> auto eth(const const_expression<E, E::s, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return edth<E, spin::up>(u); }
        template <std::size_t S, parity P, typename I, typename T> auto eth(const vector<S, P, I, T>& u) { return edth(const_view<S, P, I, T>(u)); }

        template <typename E> auto eth_bar(const const_expression<E, E::s, E::parity, typename E::index_internal_type, typename E::value_type>& u) { return edth<E, spin::down>(u); }
        template <std::size_t S, parity P, typename I, typename T> auto edth_bar(const vector<S, P, I, T>& u) { return edth_bar(const_view<S, P, I, T>(u)); }

    } // namespace sws

} // namespace math
