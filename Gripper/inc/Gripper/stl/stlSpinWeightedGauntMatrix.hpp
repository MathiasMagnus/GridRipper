#pragma once

// GSL includes
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf.h>

// Gripper includes
#include <Gripper/stl/stlSphericalIndex.hpp>
#include <Gripper/stl/stlSphericalExtent.hpp>
#include <Gripper/stl/stlSphericalVector.hpp>
#include <Gripper/stl/stlRadialVector.hpp>

// Standard C++ includes
#include <map>                      // Needed for internal storage
#include <vector>                   // Needed for internal storage
#include <array>                    // Needed for internal storage
#include <ostream>                  // Needed for output
#include <sstream>                  // Needed for output
#include <numeric>                  // std::accumulate
#include <algorithm>                // std::find_if
#include <cassert>                  // assert


namespace Multipole
{
    namespace stl
    {      
        namespace SpinWeightedGaunt
        {
            template <std::size_t L_Max, std::size_t S_Max, typename IndexInternalType, typename ValueType>
            class Matrix;

            template <std::size_t L_Max, std::size_t S_Max, typename IndexInternalType, typename ValueType>
            std::ostream& operator<<(std::ostream& os, const Matrix<L_Max, S_Max, IndexInternalType, ValueType>& mat)
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
                        os << "{ " << index.l << " }{ " << index.m << " }{ " << index.s << " }\t";
                    os << "= " << elem.second << "\n";
                }

                os.flush();

                return os;
            }

            template <std::size_t L_Max, std::size_t S_Max, typename IndexInternalType = typename SpinWeightedSpherical::Index<L_Max, S_Max>::value_type, typename ValueType = double>
            class Matrix
            {
            public:

                typedef IndexInternalType                                               index_internal_type;
                typedef SpinWeightedSpherical::Index<L_Max, S_Max, index_internal_type> index_type;
                typedef ValueType                                                       mapped_type;
                typedef std::array<index_type, 3>                                       key_type;
                typedef std::pair<const key_type, mapped_type>                          value_type;
                typedef std::vector<value_type>                                         container_type;
                typedef typename container_type::size_type                              size_type;
                typedef typename container_type::iterator                               iterator_type;
                typedef typename container_type::const_iterator                         const_iterator_type;
                typedef typename container_type::reverse_iterator                       reverse_iterator_type;
                typedef typename container_type::const_reverse_iterator                 const_reverse_iterator_type;
                typedef std::pair<const_iterator_type, const_iterator_type>             marker_type;
                typedef std::map<const index_type, marker_type>                         accelerator_structure_type;

                Matrix() = default;
                Matrix(const Matrix&) = default;
                Matrix(Matrix&&) = default;
                ~Matrix() = default;

                Matrix& operator=(Matrix&) = default;
                Matrix& operator=(Matrix&&) = default;

                static const index_internal_type l_max = static_cast<index_internal_type>(L_Max);
                static const index_internal_type s_max = static_cast<index_internal_type>(S_Max);

                void calculate() { compute(); }

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

                mapped_type& at(const key_type& pos) { return m_data.at(pos); }                  // Returns specified element of the Matrix with bounds checking
                const mapped_type& at(const key_type& pos) const { return m_data.at(pos); }      // Returns specified const element of the Matrix with bounds checking
                mapped_type& operator[](const key_type& pos) { return m_data[pos]; }             // Returns specified element of the Matrix
                const mapped_type& operator[](const key_type& pos) const { return m_data[pos]; } // Returns specified const element of the Matrix

                size_type size() const { return m_data.size(); }    // Returns the number of elements inside the Matrix
                const marker_type& get_marker(const index_type& index) const { return m_accel.at(index); }

            private:

                container_type m_data;                              // Container to store values
                accelerator_structure_type m_accel;                 // Container to store accelerator markers

                void compute()                                      // Compute entire matrix
                {
                    std::vector<std::future<container_type>> jobs;

                    for (index_internal_type l1 = 0; l1 <= l_max; ++l1)
                        for (index_internal_type m1 = -l_max; m1 <= l_max; ++m1)
                            jobs.push_back(std::async(std::launch::async, &Multipole::stl::SpinWeightedGaunt::Matrix<l_max, s_max, IndexInternalType, ValueType>::computeLM, this, l1, m1));
                    for (auto& job : jobs)
                        job.wait();

                    for (auto& job : jobs)
                    {
                        auto temp = job.get();
                        std::copy(temp.cbegin(), temp.cend(), std::back_inserter(m_data));
                    }

                    for (auto i = index_type{ 0, 0, 0 }; i <= index_type{ l_max, l_max, s_max }; ++i)
                    {
                        auto first = std::find_if(m_data.cbegin(), m_data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                        auto last = std::find_if_not(first, m_data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                        m_accel.insert(typename accelerator_structure_type::value_type(i, marker_type(first, last)));
                    }
                }

                container_type computeLM(const index_internal_type L1, const index_internal_type M1) // Compute coefficients with given L1 value
                {
                    container_type result;

                    // Compute all coefficients with designated L1,M1 value

                    for (index_internal_type S1 = -std::min(s_max, L1); S1 <= std::min(s_max, L1); ++S1)
                        for (index_internal_type L2 = 0; L2 <= l_max; ++L2)
                            for (index_internal_type M2 = -L2; M2 <= L2; ++M2)
                                for (index_internal_type S2 = -std::min(s_max, L2); S2 <= std::min(s_max, L2); ++S2)
                                    for (index_internal_type L3 = 0; L3 <= l_max; ++L3)
                                        for (index_internal_type M3 = -L3; M3 <= L3; ++M3)
                                            for (index_internal_type S3 = -std::min(s_max, L3); S3 <= std::min(s_max, L3); ++S3)
                                            {
                                                // Shortcuts
                                                if ((L1>(L2 + L3)) || (L2>(L1 + L3)) || (L3>(L1 + L2))) continue;
                                                if ((L1 + L2 + L3) % 2) continue;
                                                if ((M1 + M2 + M3) != 0) continue;

                                                // NOTE: * 2 is because of GSL calling convention
                                                mapped_type temp1 = static_cast<mapped_type>(gsl_sf_coupling_3j(L1 * 2, L2 * 2, L3 * 2, M1 * 2, M2 * 2, M3 * 2));
                                                mapped_type temp2 = static_cast<mapped_type>(gsl_sf_coupling_3j(L1 * 2, L2 * 2, L3 * 2, -S1 * 2, -S2 * 2, -S3 * 2));

                                                mapped_type res = static_cast<mapped_type>(std::sqrt((2 * L1 + 1) * (2 * L2 + 1) * (2 * L3 + 1) / (4.0 * 3.1415926535897931159979634)) * temp1 * temp2);

                                                if (std::fabs(res) > 1e-6)
                                                    result.push_back(value_type(key_type{ {index_type{ L1, M1, S1 },
                                                                                         index_type{ L2, M2, S2 },
                                                                                         index_type{ L3, M3, S3 } } },
                                                                                res));
                                            }

                    return result;
                }
            };
            
            template <typename E1,
                typename E2,
                typename G3,
                Parity P = static_cast<Parity>(E1::parity * E2::parity),
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
            class Contraction : public SpinWeightedSpherical::Expression<Contraction<E1, E2, G3>, E1::l_max, E1::s_max, P, IT, VT>
            {
            public:

                // Common typedefs

                using value_type = VT;

                // STL typedefs

                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::size_type;

                // Lattice typedefs

                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::extent_type;
                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::index_type;
                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::index_internal_type;

                Contraction() = delete;
                Contraction(const Contraction&) = default;
                Contraction(Contraction&&) = default;
                ~Contraction() = default;

                Contraction& operator=(const Contraction&) = default;
                Contraction& operator=(Contraction&&) = default;

                Contraction(const SpinWeightedSpherical::Expression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_type::value_type, typename E1::value_type>& u,
                            const SpinWeightedSpherical::Expression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_type::value_type, typename E2::value_type>& v,
                            const G3& g) :
                    _u(u), _v(v), _g(g)
                {
                    assert(u.extent() == v.extent() && (v.extent() == extent_type({ 0, 0, 0 }, { G3::l_max, G3::l_max, G3::s_max })));
                };

                // STL interface

                size_type   size()                  const { return _v.size(); }
                value_type  operator[](size_type i) const { return this->at(index_type::convert(i)); }
                value_type  at(size_type i)         const { return this->at(index_type::convert(i)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(index_type i)        const
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

                const E1& _u;
                const E2& _v;
                const G3& _g;
            };

            /*
            template <typename E1,
                typename E2,
                typename G3,
                typename PT,
                Parity P = static_cast<Parity>(E1::parity * E2::parity),
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>()),
                typename TT = SpinWeightedSpherical::Vector<E2::l_max, E2::s_max, P, typename E2::index_type::value_type, typename E2::value_type >>
            class Neumann : public SpinWeightedSpherical::Expression<Neumann<E1, E2, G3, PT>, E1::l_max, E1::s_max, P, IT, VT>
            {
            public:

                // Common typedefs

                using value_type = VT;

                // STL typedefs

                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::size_type;

                // Lattice typedefs

                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::extent_type;
                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::index_type;
                using typename SpinWeightedSpherical::VectorTraits<E1::l_max, E1::s_max, P, IT, VT>::index_internal_type;

                Neumann(const SpinWeightedSpherical::Expression<E1, E1::l_max, E1::s_max, E1::parity, typename E1::index_type::value_type, typename E1::value_type>& u,
                        const SpinWeightedSpherical::Expression<E2, E2::l_max, E2::s_max, E2::parity, typename E2::index_type::value_type, typename E2::value_type>& v,
                        const G3& g,
                        const PT& percentile) :
                    _u(u), _v(v), _g(g), _sum()
                {
                    assert(u.extent() == v.extent() && (v.extent() == extent_type({ 0, 0, 0 }, { G3::l_max, G3::l_max, G3::s_max })));

                    PT actual = static_cast<PT>(1.0);

                    std::vector<TT> series;

                    for (std::size_t count = 1; actual > percentile; ++count)
                    {
                        series.push_back(Contraction<E2, E2, G3>(series.at(count), v, g));

                        PT max_dev = static_cast<PT>(1.0);
                        
                        for (typename TT::index_type i = series.at(count).extent().initial(); series.at(count).extent().contains(i); ++i)
                        {
                            auto act_dev = std::abs(series.at(count - 1).at(i) - series.at(count).at(i) / series.at(count - 1).at(i));
                            
                            max_dev = std::max(max_dev, act_dev);
                        }

                        actual = max_dev; // The largest relative pair-wise change in coefficients
                    }

                    std::size_t count = 0;
                    _sum = std::accumulate(series.cbegin(), series.cend(), _sum, [&](const TT& sum, const TT& coeff) { return sum + (count++ % 2 ? -coeff : coeff); });
                };

                // STL interface

                size_type   size()                  const { return _v.size(); }
                value_type  operator[](size_type i) const { return this->at(index_type::convert(i)); }
                value_type  at(size_type i)         const { return this->at(index_type::convert(i)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(index_type i)        const { return Contraction<E1, E2, G3>(_u, _sum, _g).at(i); }

            private:

                const E1& _u;
                const E2& _v;
                const G3& _g;
                TT _sum;
            };
            */
        } // namespace SpinWeightedGaunt

    } // namespace stl

} // namespace Multipole

namespace Multipole
{
    namespace stl
    {
        namespace SpinWeightedGaunt
        {
            template <typename E1,
                typename E2,
                typename G3,
                std::size_t L = E1::l_max,
                std::size_t S = E1::s_max,
                Parity P = E1::parity,
                typename I = typename E1::index_type,
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
                auto contract(const G3& gaunt,
                              const SpinWeightedSpherical::Expression<E1, L, S, P, I, typename E1::value_type>& lhs,
                              const SpinWeightedSpherical::Expression<E2, L, S, P, I, typename E2::value_type>& rhs)
            {
                return Contraction<E1, E2, G3>(lhs, rhs, gaunt);
            };
            
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
            
            template <typename E1,
                typename E2,
                typename G3,
                typename PT,
                std::size_t L = E1::l_max,
                std::size_t S = E1::s_max,
                Parity P = E1::parity,
                typename I = typename E1::index_type,
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
                auto neumann(const G3& gaunt,
                             const SpinWeightedSpherical::Expression<E1, L, S, P, I, typename E1::value_type>& lhs,
                             const SpinWeightedSpherical::Expression<E2, L, S, P, I, typename E2::value_type>& rhs,
                             const PT& percentile)
            {
                SpinWeightedSpherical::Vector<L, S, E2::parity, typename E2::index_type::value_type, typename E2::value_type> running_x = rhs, temp;
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

                //std::cout << "Step 1 has deviation = " << deviation << std::endl;

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

                    //std::cout << "Step " << I << " has deviation = " << deviation << std::endl;
                }

                //std::cout << "sum / f_00_Y_00 = factor : " << sum << " / " << rhs.at(typename E2::index_type{ 0, 0, 0 }) << " = " << sum * f_00_Y_00_inv << std::endl;

                return lhs * (sum * f_00_Y_00_inv);
            };
            
        } // namespace SpinWeightedGaunt

    } // namespace stl

} // namespace Multipole