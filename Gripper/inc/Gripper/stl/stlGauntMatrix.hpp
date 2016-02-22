#pragma once

// GSL includes
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf.h>

// Gripper includes
#include <Gripper/stl/stlMath.hpp>              // Multipole::stl::pi
#include <Gripper/stl/stlParity.hpp>            // Multipole::stl::Parity
#include <Gripper/stl/stlSphericalIndex.hpp>    // Multipole::stl::Spherical::Index
#include <Gripper/stl/stlSphericalExtent.hpp>   // Multipole::stl::Spherical::Extent
#include <Gripper/stl/stlSphericalVector.hpp>   // Multipole::stl::Spherical::Vector
#include <Gripper/stl/stlRadialVector.hpp>      // Multipole::stl::Radial::Vector

// Standard C++ includes
#include <cassert>                  // assert
#include <cstddef>                  // std::size_t
#include <cstdint>                  // std::int32_t
#include <functional>               // std::reference_wrapper, std::ref, std::cref
#include <map>                      // std::map
#include <vector>                   // std::vector
#include <array>                    // std::array
#include <ostream>                  // std::ostream
#include <numeric>                  // std::accumulate
#include <algorithm>                // std::find_if, std::copy
#include <iterator>                 // std::back_inserter
#include <cmath>                    // std::fabs, std::sqrt
#include <future>                   // std::future, std::async


namespace Multipole
{
    namespace stl
    {      
        namespace Gaunt
        {
            /// <summary>Traits class consisting of type aliases describing the cross integrals of spherical harmonics (aka. Gaunt-matrix of Gaunt-coefficients).</summary>
            ///
            template <std::size_t L_Max, typename IT, typename VT>
            struct MatrixTraits
            {
                // Lattice type aliases

                using extent_type = Spherical::Extent<L_Max, IT>;
                using index_type = typename extent_type::index_type;
                using index_internal_type = typename index_type::value_type;

                // STL type aliases

                using mapped_type = VT;
                using key_type = std::array<index_type, 3>;
                using value_type = std::pair<const key_type, mapped_type>;
                using container_type = std::vector<value_type>;
                using size_type = typename container_type::size_type;

                // ConstIterator aliases

                using const_iterator_type = typename container_type::const_iterator;
                using const_reverse_iterator_type = typename container_type::const_reverse_iterator;

                // Iterator aliases

                using iterator_type = typename container_type::iterator;
                using reverse_iterator_type = typename container_type::reverse_iterator;

                // Matrix type aliases

                using marker_type = std::pair<const_iterator_type, const_iterator_type>;
                using accelerator_structure_type = std::map<const index_type, marker_type>;

                // Lattice static members

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
            };


            /// <summary>Read-only Expression Template base class of spherical expansion coefficient vectors to be implemented statically.</summary>
            /// <remarks>Gaunt matrices only have ConstExpression interfaces, as their contents are not meant to change after construction.</remarks>
            ///
            template <typename ET, std::size_t L_Max, typename IT, typename VT>
            struct ConstExpression : public MatrixTraits<L_Max, IT, VT>
            {
                // ConstExpression type aliases

                using expression_type = ET;

                // Lattice type aliases

                using typename MatrixTraits<L_Max, IT, VT>::extent_type;
                using typename MatrixTraits<L_Max, IT, VT>::index_type;
                using typename MatrixTraits<L_Max, IT, VT>::index_internal_type;

                // STL type aliases

                using typename MatrixTraits<L_Max, IT, VT>::mapped_type;
                using typename MatrixTraits<L_Max, IT, VT>::key_type;
                using typename MatrixTraits<L_Max, IT, VT>::value_type;
                using typename MatrixTraits<L_Max, IT, VT>::container_type;
                using typename MatrixTraits<L_Max, IT, VT>::size_type;

                // ConstIterator aliases

                using typename MatrixTraits<L_Max, IT, VT>::const_iterator_type;
                using typename MatrixTraits<L_Max, IT, VT>::const_reverse_iterator_type;

                // Matrix type aliases

                using typename MatrixTraits<L_Max, IT, VT>::marker_type;
                using typename MatrixTraits<L_Max, IT, VT>::accelerator_structure_type;

                // Lattice static members

                using MatrixTraits<L_Max, IT, VT>::l_max;

                // ConstExpression interface

                operator const expression_type&()   const { return static_cast<const expression_type&>(*this); }

                // ConstAssocSTL interface

                /// <summary>Returns the number of coefficients inside the matrix.</summary>
                ///
                size_type size() const { return static_cast<const expression_type&>(*this).size(); }

                /// <summary>Returns the element corresponding to key <c>pos</c> without bounds checking.</summary>
                ///
                const mapped_type& operator[](const key_type& pos) const { return static_cast<const expression_type&>(*this)[pos]; }

                /// <summary>Returns the element corresponding to key <c>pos</c> with bounds checking.</summary>
                ///
                const mapped_type& at(const key_type& pos) const { return static_cast<const expression_type&>(*this).at(pos); }

                // ConstIterator interface

                /// <summary>Returns an immutable iterator to the first coefficient in the matrix.</summary>
                ///
                const_iterator_type cbegin() const { return static_cast<const expression_type&>(*this).cbegin(); }

                /// <summary>Returns an immutable iterator to the coefficient after the last one in the matrix.</summary>
                ///
                const_iterator_type cend() const { return static_cast<const expression_type&>(*this).cend(); }

                /// <summary>Returns an immutable reverse iterator to the first coefficient in the matrix.</summary>
                ///
                const_iterator_type crbegin() const { return static_cast<const expression_type&>(*this).crbegin(); }

                /// <summary>Returns an immutable reverse iterator to the coefficient after the last one in the matrix.</summary>
                ///
                const_iterator_type crend() const { return static_cast<const expression_type&>(*this).crend(); }

                // GauntMatrix interface

                /// <summary>Returns a pair of iterators to all coefficients with first index <c>index</c>.</summary>
                ///
                const marker_type& get_marker(const index_type& index) const { return static_cast<const expression_type&>(*this).get_marker(index); }

                /// <summary>Returns a pair of indecies representing the span of the Gaunt-matrix.</summary>
                ///
                const extent_type& extent() const { return static_cast<const expression_type&>(*this).extent(); }
            };


            /// <summary>Expression Template that applies applies the spin stepping operator (edth) to a series expansion.</summary>
            ///
            template <std::size_t L_Max, typename IT, typename VT>
            class Matrix : public MatrixTraits<L_Max, IT, VT>
            {
            public:

                // Lattice type aliases

                using typename MatrixTraits<L_Max, IT, VT>::extent_type;
                using typename MatrixTraits<L_Max, IT, VT>::index_type;
                using typename MatrixTraits<L_Max, IT, VT>::index_internal_type;

                // STL type aliases

                using typename MatrixTraits<L_Max, IT, VT>::mapped_type;
                using typename MatrixTraits<L_Max, IT, VT>::key_type;
                using typename MatrixTraits<L_Max, IT, VT>::value_type;
                using typename MatrixTraits<L_Max, IT, VT>::container_type;
                using typename MatrixTraits<L_Max, IT, VT>::size_type;

                // Iterator aliases

                using typename MatrixTraits<L_Max, IT, VT>::iterator_type;
                using typename MatrixTraits<L_Max, IT, VT>::const_iterator_type;
                using typename MatrixTraits<L_Max, IT, VT>::reverse_iterator_type;
                using typename MatrixTraits<L_Max, IT, VT>::const_reverse_iterator_type;

                // Matrix type aliases

                using typename MatrixTraits<L_Max, IT, VT>::marker_type;
                using typename MatrixTraits<L_Max, IT, VT>::accelerator_structure_type;

                // Lattice static members

                using MatrixTraits<L_Max, IT, VT>::l_max;

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                Matrix() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                Matrix(const Matrix& in) = default;

                /// <summary>Default move constructor.</summary>
                ///
                Matrix(Matrix&& in) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~Matrix() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                Matrix& operator=(const Matrix&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                Matrix& operator=(Matrix&&) = default;

                ///<summary>Constructs a series expansion vector of <c>ext</c> size.</summary>
                ///
                Matrix(const extent_type& ext) : _extent(ext), _data(), _accel() { compute(); }

                // Iterator interface

                iterator_type begin() { return _data.begin(); }
                const_iterator_type begin() const { return _data.begin(); }
                const_iterator_type cbegin() const { return _data.cbegin(); }

                iterator_type end() { return _data.end(); }
                const_iterator_type end() const { return _data.end(); }
                const_iterator_type cend() const { return _data.cend(); }

                reverse_iterator_type rbegin() { return _data.rbegin(); }
                const_reverse_iterator_type rbegin() const { return _data.rbegin(); }
                const_reverse_iterator_type crbegin() const { return _data.crbegin(); }

                reverse_iterator_type rend() { return _data.rend(); }
                const_reverse_iterator_type rend() const { return _data.rend(); }
                const_reverse_iterator_type crend() const { return _data.crend(); }

                // ConstAssocSTL interface

                /// <summary>Returns the number of coefficients inside the matrix.</summary>
                ///
                size_type size() const { return _data.size(); }

                /// <summary>Returns the element corresponding to key <c>pos</c> without bounds checking.</summary>
                ///
                const mapped_type& operator[](const key_type& pos) const { return _data[pos]; }

                /// <summary>Returns the element corresponding to key <c>pos</c> with bounds checking.</summary>
                ///
                const mapped_type& at(const key_type& pos) const { return _data.at(pos); }

                // GauntMatrix interface
                
                /// <summary>Returns a pair of iterators to all coefficients with first index <c>index</c>.</summary>
                ///
                const marker_type& get_marker(const index_type& index) const { return _accel.at(index); }

                /// <summary>Returns a pair of indecies representing the span of the Gaunt-matrix.</summary>
                ///
                const extent_type& extent() const { return _extent; }

            private:

                extent_type _extent;                    // Extent of l1m1 indices
                container_type _data;                   // Container to store values
                accelerator_structure_type _accel;      // Container to store accelerator markers

                /// <summary>Computes matrix coefficients with first index in the extent of the matrix.</summary>
                ///
                void compute()
                {
                    std::vector<std::future<container_type>> jobs;

                    // Launch threads
                    for (index_type i = extent().initial() ; extent().contains(i); ++i)
                        jobs.push_back(std::async(std::launch::async, &Multipole::stl::Gaunt::Matrix<l_max, index_internal_type, mapped_type>::computeLM, this, i.l, i.m));

                    // Wait for jobs to complete
                    for (auto& job : jobs) job.wait();

                    // Merge results
                    for (auto& job : jobs)
                    {
                        auto temp = job.get();
                        std::copy(temp.cbegin(), temp.cend(), std::back_inserter(_data));
                    }

                    // Build accelerating structure
                    for (index_type i = extent().initial(); extent().contains(i); ++i)
                    {
                        auto first = std::find_if(_data.cbegin(), _data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                        auto last = std::find_if_not(first, _data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                        _accel.insert(typename accelerator_structure_type::value_type(i, marker_type(first, last)));
                    }
                }

                /// <summary>Computes matrix coefficients with first index <c>L1</c> and <c>M1</c>.</summary>
                ///
                container_type computeLM(const index_internal_type L1, const index_internal_type M1)
                {
                    container_type result;

                    // Compute all coefficients with designated L1,M1 value
                    for (index_internal_type L2 = 0; L2 <= l_max; ++L2)
                        for (index_internal_type M2 = -L2; M2 <= L2; ++M2)
                            for (index_internal_type L3 = 0; L3 <= l_max; ++L3)
                                for (index_internal_type M3 = -L3; M3 <= L3; ++M3)
                                {
                                    // Shortcuts
                                    if ((L1>(L2 + L3)) || (L2>(L1 + L3)) || (L3>(L1 + L2))) continue;
                                    if ((L1 + L2 + L3) % 2) continue;
                                    if ((M1 + M2 + M3) != 0) continue;

                                    // NOTE: * 2 is because of GSL calling convention
                                    mapped_type temp1 = static_cast<mapped_type>(gsl_sf_coupling_3j(L1 * 2, L2 * 2, L3 * 2, M1 * 2, M2 * 2, M3 * 2));
                                    mapped_type temp2 = static_cast<mapped_type>(gsl_sf_coupling_3j(L1 * 2, L2 * 2, L3 * 2, 0, 0, 0));

                                    mapped_type res = static_cast<mapped_type>(std::sqrt(static_cast<mapped_type>((2 * L1 + 1) * (2 * L2 + 1) * (2 * L3 + 1))) * temp1 * temp2 / (2 * std::sqrt(pi<mapped_type>)));

                                    if (std::fabs(res) > 1e-6)
                                        result.push_back(value_type(key_type{ {index_type{ L1, M1 },
                                                                             index_type{ L2, M2 },
                                                                             index_type{ L3, M3 } } },
                                                                    res));
                                }

                    return result;
                }
            };


            /// <summary>Expression template providing read-only access with reference semantics to the elements of a Gaunt-matrix with storage.</summary>
            ///
            template <std::size_t L_Max, typename IT, typename VT>
            class ConstView : public ConstExpression<ConstView<L_Max, IT, VT>, L_Max, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = typename ConstExpression<ConstView<L_Max, IT, VT>, L_Max, IT, VT>;

                // Lattice type aliases

                using typename expression_type::extent_type;
                using typename expression_type::index_type;
                using typename expression_type::index_internal_type;

                // STL type aliases

                using typename expression_type::mapped_type;
                using typename expression_type::key_type;
                using typename expression_type::value_type;
                using typename expression_type::container_type;
                using typename expression_type::size_type;

                // ConstIterator aliases

                using typename expression_type::const_iterator_type;
                using typename expression_type::const_reverse_iterator_type;

                // Matrix type aliases

                using typename expression_type::marker_type;
                using typename expression_type::accelerator_structure_type;

                // Lattice static members

                using expression_type::l_max;

            private:

                using matrix_type = Matrix<l_max, index_internal_type, mapped_type>;
                using reference_type = std::reference_wrapper<const matrix_type>;

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
                ConstView(ConstView&&) = delete;

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
                ConstView(const matrix_type& m) : _m(std::cref(m)) {}

                // ConstAssocSTL interface

                /// <summary>Returns the number of coefficients inside the matrix.</summary>
                ///
                size_type size() const { return _m.size(); }

                /// <summary>Returns the element corresponding to key <c>pos</c> without bounds checking.</summary>
                ///
                const mapped_type& operator[](const key_type& pos) const { return _m[pos]; }

                /// <summary>Returns the element corresponding to key <c>pos</c> with bounds checking.</summary>
                ///
                const mapped_type& at(const key_type& pos) const { return _m.at(pos); }

                // ConstIterator interface

                /// <summary>Returns an immutable iterator to the first coefficient in the matrix.</summary>
                ///
                const_iterator_type cbegin() const { return _m.cbegin(); }

                /// <summary>Returns an immutable iterator to the coefficient after the last one in the matrix.</summary>
                ///
                const_iterator_type cend() const { return _m.cend(); }

                /// <summary>Returns an immutable reverse iterator to the first coefficient in the matrix.</summary>
                ///
                const_iterator_type crbegin() const { return _m.crbegin(); }

                /// <summary>Returns an immutable reverse iterator to the coefficient after the last one in the matrix.</summary>
                ///
                const_iterator_type crend() const { return _m.crend(); }

                // GauntMatrix interface

                /// <summary>Returns a pair of iterators to all coefficients with first index <c>index</c>.</summary>
                ///
                const marker_type& get_marker(const index_type& index) const { return _m.get_marker(index); }

                /// <summary>Returns a pair of indecies representing the span of the Gaunt-matrix.</summary>
                ///
                const extent_type& extent() const { return _m.extent(); }

            private:

                const reference_type _m;
            };


            /// <summary>Expression Template performing a contraction of two spherical series expansion vectors with a Gaunt-matrix.</summary>
            ///
            template <typename E1,
                typename E2,
                typename G3,
                Parity P = static_cast<Parity>(E1::parity * E2::parity),
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
            class Contraction : public Spherical::ConstExpression<Contraction<E1, E2, G3>, E1::l_max, P, IT, VT>
            {
            public:

                // Expression type aliases

                using expression_type = typename Spherical::ConstExpression<Contraction<E1, E2, G3>, E1::l_max, P, IT, VT>;

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
                using expression_type::parity;

                /// <summary>Constructs a <c>Contraction</c> from an extent <c>ext</c> and a function object <c>f</c>.</summary>
                ///
                Contraction(const Spherical::ConstExpression<E1, E1::l_max, E1::parity, typename E1::index_type::value_type, typename E1::value_type>& u,
                            const Spherical::ConstExpression<E2, E2::l_max, E2::parity, typename E2::index_type::value_type, typename E2::value_type>& v,
                            const G3& g) :
                    _u(u), _v(v), _g(g)
                {
                    // Assert that u and v have identical extents AND that they are maximal
                    assert(u.extent() == v.extent() && (v.extent() == extent_type({ 0, 0 }, { G3::l_max, G3::l_max })));
                };

                // ConstSTL interface

                /// <summary>Returns the number of coefficients inside the vector.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                size_type  size()                  const { return _u.size(); }

                /// <summary>Returns the <c>i</c>th element without bounds checking.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                value_type operator[](size_type i) const { static_assert(false, "This function requires size_type >> index_type conversion to be implemented"); return static_cast<value_type>(0); }

                /// <summary>Returns the <c>i</c>th element with bounds checking.</summary>
                /// <remarks>This function is not yet implemented.</remarks>
                ///
                value_type at(size_type i)         const { static_assert(false, "This function requires size_type >> index_type conversion to be implemented"); return static_cast<value_type>(0); }

                // ConstLattice interface

                /// <summary>Returns a pair of indecies representing the span of the series expansion.</summary>
                ///
                extent_type extent()               const { return _g.extent(); }

                /// <summary>Returns the coefficient corresponding to index <c>i</c>.</summary>
                ///
                value_type at(index_type i)        const
                {
                    return std::accumulate(_g.get_marker(i).first,
                                           _g.get_marker(i).second,
                                           static_cast<value_type>(0),
                                           [&](const value_type& accum, const typename G3::value_type& elem)
                    {
                        return accum + _u.at(elem.first.at(1)) * _v.at(elem.first.at(2)) * elem.second;
                    });
                }

            private:

                const E1 _u;
                const E2 _v;
                const G3 _g;
            };


            /// <summary>STL stream operator overload for formatted console output.</summary>
            ///
            template <std::size_t L_Max, typename IT, typename VT>
            std::ostream& operator<<(std::ostream& os, const Matrix<L_Max, IT, VT>& mat)
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

                for (auto& elem : mat)
                {
                    for (auto& index : elem.first)
                        os << "{ " << index.l << " }{ " << index.m << " }\t";
                    os << "= " << elem.second << "\n";
                }

                os.flush();

                return os;
            }

        } // namespace Gaunt

    } // namespace stl

} // namespace Multipole

namespace Multipole
{
    namespace stl
    {
        namespace Gaunt
        {
            /// <summary>Contraction helper function with template type deduction.</summary>
            ///
            template <typename E1,
                typename E2,
                typename G3>
                auto contract(const G3& gaunt,
                              const Spherical::ConstExpression<E1, E1::l_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& lhs,
                              const Spherical::ConstExpression<E2, E2::l_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& rhs)
            {
                return Contraction<E1, E2, G3>(lhs, rhs, gaunt);
            };

            /// <summary>Contraction helper function with template type deduction.</summary>
            ///
            template <std::size_t L1, Parity P1, typename I1, typename T1,
                typename E2,
                typename G3>
                auto contract(const G3& gaunt,
                              const Spherical::Vector<L1, P1, I1, T1>& lhs,
                              const Spherical::ConstExpression<E2, E2::l_max, E2::parity, typename E2::index_internal_type, typename E2::value_type>& rhs)
            {
                using view_type = Spherical::ConstView<L1, P1, I1, T1>;

                return Contraction<view_type, E2, G3>(view_type(lhs), rhs, gaunt);
            };

            /// <summary>Contraction helper function with template type deduction.</summary>
            ///
            template <typename E1,
                std::size_t L2, Multipole::stl::Parity P2, typename I2, typename T2,
                typename G3>
                auto contract(const G3& gaunt,
                              const Spherical::ConstExpression<E1, E1::l_max, E1::parity, typename E1::index_internal_type, typename E1::value_type>& lhs,
                              const Spherical::Vector<L2, P2, I2, T2>& rhs)
            {
                using view_type = Spherical::ConstView<L2, P2, I2, T2>;

                return Contraction<E1, view_type, G3>(lhs, view_type(rhs), gaunt);
            };

            /// <summary>Contraction helper function with template type deduction.</summary>
            ///
            template <std::size_t L1, Parity P1, typename I1, typename T1,
                std::size_t L2, Parity P2, typename I2, typename T2,
                typename G3>
                auto contract(const G3& gaunt,
                              const Spherical::Vector<L1, P1, I1, T1>& lhs,
                              const Spherical::Vector<L2, P2, I2, T2>& rhs)
            {
                using view_type1 = Spherical::ConstView<L1, P1, I1, T1>;
                using view_type2 = Spherical::ConstView<L2, P2, I2, T2>;

                return Contraction<view_type1, view_type2, G3>(view_type1(lhs), view_type2(rhs), gaunt);
            };

            /// <summary>Map helper function with template type deduction.</summary>
            ///
            template <typename RE1,
                typename RE2,
                typename G3,
                typename D,
                typename I,
                typename VT1,
                typename VT2>
                auto contract(const G3& gaunt,
                              const Multipole::stl::Radial::Expression<RE1, D, I, VT1>& lhs,
                              const Multipole::stl::Radial::Expression<RE2, D, I, VT2>& rhs)
            {
                return Radial::zip(lhs, rhs, [&](auto&& l, auto&& r) { return Contraction<typename VT1::expression_type, typename VT2::expression_type, G3>(l, r, gaunt); });
            };


            namespace impl
            {
                template <typename E1,
                          typename E2,
                          typename G3,
                          typename PT>
                const auto neumann(const G3& gaunt,
                                   const E1& lhs,
                                   const E2& rhs,
                                   const PT percentile)
                {
                    Spherical::Vector<E2::l_max, E2::parity, typename E2::index_type::value_type, typename E2::value_type> running_x = rhs, temp;
                    auto f_00_Y_00_inv = static_cast<typename E2::value_type>(1) / rhs.at(typename E2::index_type{ 0, 0, 0 });
                    typename E2::value_type x_prime, sum, prev_sum;
                    PT deviation;

                    // 0th order term in summation
                    x_prime = static_cast<typename E2::value_type>(1);
                    sum = x_prime;

                    // 1st order term in summation
                    x_prime = 0;

                    for (auto i = ++typename E2::index_type(running_x.extent().initial()); running_x.extent().contains(i); ++i)
                    {
                        x_prime += running_x.at(i);
                    }
                    x_prime *= f_00_Y_00_inv;

                    prev_sum = sum;
                    sum -= x_prime;
                    deviation = std::abs(prev_sum - sum) / std::max(prev_sum, sum);

                    // 2nd order and onward
                    for (std::int32_t I = 2; deviation > percentile; ++I)
                    {
                        x_prime = 0;
                        temp = contract(gaunt, running_x, rhs);

                        for (auto i = ++typename E2::index_type(running_x.extent().initial()); running_x.extent().contains(i); ++i)
                        {
                            x_prime += temp.at(i);
                        }
                        x_prime *= f_00_Y_00_inv;

                        prev_sum = sum;
                        sum += static_cast<typename E2::value_type>((I % 2) ? -1.f : 1.f) * x_prime;
                        deviation = std::abs(prev_sum - sum) / std::max(prev_sum, sum);
                        running_x = temp;
                    }

                    return lhs * (sum * f_00_Y_00_inv);
                }

            } // namespace impl
            
            /// <summary>Contraction helper function with template type deduction performing Neumann-series division.</summary>
            ///
            template <typename E1,
                      typename E2,
                      typename G3,
                      typename PT>
            auto neumann(const G3& gaunt,
                         const Spherical::ConstExpression<E1, E1::l_max, E1::parity, typename E1::internal_index_type, typename E1::value_type>& lhs,
                         const Spherical::ConstExpression<E2, E2::l_max, E2::parity, typename E2::internal_index_type, typename E2::value_type>& rhs,
                         const PT percentile)
            {
                return impl::neumann(gaunt, lhs, rhs, percentile);
            };

            /// <summary>Contraction helper function with template type deduction performing Neumann-series division.</summary>
            ///
            template <std::size_t L1, Multipole::stl::Parity P1, typename I1, typename T1,
                typename E2,
                typename G3,
                typename PT>
                auto neumann(const G3& gaunt,
                             const Spherical::Vector<L1, P1, I1, T1>& lhs,
                             const Spherical::ConstExpression<E2, E2::l_max, E2::parity, typename E2::internal_index_type, typename E2::value_type>& rhs,
                             const PT percentile)
            {
                return impl::neumann(gaunt,
                                     Multipole::stl::Spherical::ConstView<L1, P1, I1, T1>(lhs),
                                     rhs,
                                     percentile);
            };

            /// <summary>Contraction helper function with template type deduction performing Neumann-series division.</summary>
            ///
            template <typename E1,
                std::size_t L2, Multipole::stl::Parity P2, typename I2, typename T2,
                typename G3,
                typename PT>
                auto neumann(const G3& gaunt,
                             const Spherical::ConstExpression<E1, E1::l_max, E1::parity, typename E1::internal_index_type, typename E1::value_type>& lhs,
                             const Spherical::Vector<L2, P2, I2, T2>& rhs,
                             const PT percentile)
            {
                return impl::neumann(gaunt,
                                     lhs,
                                     Multipole::stl::Spherical::ConstView<L2, P2, I2, T2>(rhs),
                                     percentile);
            };

            /// <summary>Contraction helper function with template type deduction performing Neumann-series division.</summary>
            ///
            template <std::size_t L1, Multipole::stl::Parity P1, typename I1, typename T1,
                std::size_t L2, Multipole::stl::Parity P2, typename I2, typename T2,
                typename G3,
                typename PT>
                auto neumann(const G3& gaunt,
                             const Spherical::Vector<L1, P1, I1, T1>& lhs,
                             const Spherical::Vector<L2, P2, I2, T2>& rhs,
                             const PT percentile)
            {
                return impl::neumann(gaunt,
                                     Multipole::stl::Spherical::ConstView<L1, P1, I1, T1>(lhs),
                                     Multipole::stl::Spherical::ConstView<L2, P2, I2, T2>(rhs),
                                     percentile);
            };
            
        } // namespace Gaunt

    } // namespace stl

} // namespace Multipole
