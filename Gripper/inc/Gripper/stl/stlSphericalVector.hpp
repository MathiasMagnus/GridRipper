#ifndef STLSPHERICALVECTOR_HPP
#define STLSPHERICALVECTOR_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlSphericalExtent.hpp>
#include <Gripper/stl/stlRuntime.hpp>           // For ability to read global runtime configuration

// GSL includes
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf.h>

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <future>                   // Needed for guided parallel processing
#include <utility>                  // Needed for std::pair
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

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

                typedef VT value_type;

                // STL typedefs

                typedef std::vector<value_type>                         container_type;
                typedef typename container_type::size_type              size_type;
                typedef typename container_type::iterator               iterator_type;
                typedef typename container_type::const_iterator         const_iterator_type;
                typedef typename container_type::reverse_iterator       reverse_iterator_type;
                typedef typename container_type::const_reverse_iterator const_reverse_iterator_type;

                // Lattice typedefs

                typedef Extent                                  extent_type;
                typedef typename extent_type::index_type        index_type;
            };


            template <typename ET, typename VT>
            class Expression : public Traits<VT>
            {
            public:

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


            template <typename VT>
            class Vector : public Expression<Vector<VT>, VT>
            {
            public:

                // Common interface

                Vector() { CALLING }
                Vector(const Vector& in) : m_data(in.m_data) { CALLING }
                Vector(Vector&& in) : m_data(std::move(in.m_data)), m_extent(std::move(in.m_extent)) { CALLING }
                ~Vector() { CALLING }

                Vector& operator=(Vector&) { m_data = in.m_data; return *this; }
                Vector& operator=(Vector&&) { m_data = std::move(in.m_data); return *this; }

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

                Vector(const Extent& ext) : m_data(static_cast<size_type>(ext.final() - ext.initial())), m_extent(ext) { CALLING }

                value_type& at(const index_type& pos) { return m_data.at(pos.i); }
                const value_type& at(const index_type& pos) const { return m_data.at(pos.i); }

                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, typename ValueType>
                Vector(Expression<VecExpr, ValueType> const& expr)
                {
                    ENTERING

                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    Gripper::stl::clog << Gripper::DEBUG << "m_extent = {" << m_extent.final().i - m_extent.initial().i << "}\texpr.extent() = {" << v.extent().final().i - v.extent().initial().i << "}";
                    Gripper::stl::clog << Gripper::DEBUG << "m_data.size() = " << m_data.size() << "\texpr.size() = " << expr.size();

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Spherical::Vector::Vector{resize()}";
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    if (Gripper::stl::Private::g_runtime.getCompilerParallelization()) for (index_type i = m_extent.initial(); i != m_extent.final(); ++i) this->at(i) = v.at(i);
                    else
                    {
                        std::vector<std::pair<index_type, index_type>> ranges(std::thread::hardware_concurrency());
                        std::vector<std::future<void>> jobs(std::thread::hardware_concurrency());
                        for (std::size_t i = 0; i < ranges.size(); ++i)
                        {
                            ranges.at(i).first = m_extent.initial() + (m_extent.final() - m_extent.initial()) * (i) / ranges.size();
                            ranges.at(i).second = m_extent.initial() + (m_extent.final() - m_extent.initial()) * (i + 1) / ranges.size();
                        }

                        for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                        {
                            for (index_type i = ranges.at(j).first; i < ranges.at(j).second; ++i) this->at(i) = v.at(i);
                        });
                        for (auto& job : jobs) job.wait();

                        LEAVING
                    }
                }

            private:

                container_type m_data;
                extent_type m_extent;
            };


            template <typename E1, typename E2>
            class Sum : public Expression<Sum<E1, E2>, decltype(std::declval<typename E1::value_type>() + std::declval<typename E2::value_type>())>
            {
            private:

                typedef decltype(std::declval<typename E1::value_type>() + std::declval<typename E2::value_type>()) return_type;

            public:

                Sum(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] + _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) + _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) + _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E1, typename E2>
            class Dif : public Expression<Dif<E1, E2>, decltype(std::declval<typename E1::value_type>() - std::declval<typename E2::value_type>())>
            {
            private:

                typedef decltype(std::declval<typename E1::value_type>() - std::declval<typename E2::value_type>()) return_type;

            public:

                Dif(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] - _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) - _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) - _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E1, typename E2>
            class Mul : public Expression<Mul<E1, E2>, decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>())>
            {
            private:

                typedef decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>()) return_type;

            public:

                Mul(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] * _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) * _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) * _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E1, typename E2>
            class Div : public Expression<Div<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
            {
            private:

                typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

            public:

                Div(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] / _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) / _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) / _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E>
            class Scale : public Expression<Scale<E>, decltype(std::declval<double>() * std::declval<typename E::value_type>())>
            {
            private:

                typedef decltype(std::declval<double>() * std::declval<typename E::value_type>()) return_type;

            public:

                Scale(double alpha, Expression<E, typename E::value_type> const& v) : _alpha(alpha), _v(v) {}

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _alpha * _v[i]; }
                return_type at(size_type i)         const { return _alpha * _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _alpha * _v.at(i); }

            private:

                double _alpha;
                E const& _v;
            };

        } // namespace Spherical

        namespace SpinWeightedSpherical
        {
            template <typename IT, typename VT>
            class Traits
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

                typedef Extent<IT>                                      extent_type;
                typedef typename extent_type::index_type                index_type;
            };


            template <typename ET, typename IT, typename VT>
            class Expression : public Traits<IT, VT>
            {
            public:

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


            template <typename IT, typename VT>
            class Vector : public Expression<Vector<IT, VT>, IT, VT>
            {
            public:

                // Common interface

                Vector() = default;
                Vector(const Vector& in) = default;
                Vector(Vector&& in) = default;
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

                Vector(const typename index_type::value_type& l_max,
                       const typename index_type::value_type& s_max) :
                    m_extent(extent_type({ 0, 0, -s_max }, { l_max, l_max, s_max })),
                    m_data(convert(m_extent.final) + 1, static_cast<value_type>(0)) { CALLING }

                value_type& at(const index_type& pos) { return m_data.at(convert(pos)); }
                const value_type& at(const index_type& pos) const { return m_data.at(convert(pos)); }

                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, typename IndexType, typename ValueType>
                Vector(Expression<VecExpr, IndexType, ValueType> const& expr) : Vector()
                {
                    ENTERING

                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    Gripper::stl::clog << Gripper::DEBUG << "m_extent = { " << m_extent.initial.l << "," << m_extent.initial.m << "," << m_extent.initial.s << ";" << m_extent.final.l << "," << m_extent.final.m << "," << m_extent.final.s << " }\texpr.extent() = { " << v.extent().initial.l << "," << v.extent().initial.m << "," << v.extent().initial.s << ";" << v.extent().final.l << "," << v.extent().final.m << "," << v.extent().final.s << " }";
                    Gripper::stl::clog << Gripper::DEBUG << "m_data.size() = " << m_data.size() << "\texpr.size() = " << expr.size();

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Spherical::Vector::Vector{resize()}";
                        m_extent = v.extent();
                        m_data.resize(v.size());
                    }

                    /*if (Gripper::stl::Private::g_runtime.getCompilerParallelization()) for (index_type i = m_extent.initial; i != m_extent.final; ++i) this->at(i) = v.at(i);
                    else
                    {
                        std::vector<std::pair<index_type, index_type>> ranges(std::thread::hardware_concurrency());
                        std::vector<std::future<void>> jobs(std::thread::hardware_concurrency());
                        for (std::size_t i = 0; i < ranges.size(); ++i)
                        {
                            ranges.at(i).first = m_extent.initial + (m_extent.final - m_extent.initial) * (i) / ranges.size();
                            ranges.at(i).second = m_extent.initial + (m_extent.final - m_extent.initial) * (i + 1) / ranges.size();
                        }

                        for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                        {
                            for (index_type i = ranges.at(j).first; i < ranges.at(j).second; ++i) this->at(i) = v.at(i);
                        });
                        for (auto& job : jobs) job.wait();
                    }*/

                    // WARNING!!!! Single-device evaluator. Will fail with multi-device.
                    if (Gripper::stl::Private::g_runtime.getCompilerParallelization()) for (size_type i = 0; i < m_data.size(); ++i) this->at(i) = v.at(convert(i));
                    else
                    {
                        std::vector<std::pair<size_type, size_type>> ranges(std::thread::hardware_concurrency());
                        std::vector<std::future<void>> jobs(std::thread::hardware_concurrency());
                        for (std::size_t i = 0; i < ranges.size(); ++i)
                        {
                            ranges.at(i).first = m_data.size() * (i) / ranges.size();
                            ranges.at(i).second = m_data.size() * (i + 1) / ranges.size();
                        }

                        for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                        {
                            for (size_type i = ranges.at(j).first; i < ranges.at(j).second; ++i) this->at(convert(i)) = v.at(convert(i));
                        });
                        for (auto& job : jobs) job.wait();
                    }

                    LEAVING
                }

                // Utility function (should not exist!)

                index_type convert(const size_type& i) const 
                {
                    IT i_corr = static_cast<IT>(i / (2 * m_extent.final.s + 1)); // truncation on division is intensional (omit std::floor and double up-cast)
                    IT l = static_cast<IT>(std::round((-1. + std::sqrt(1. + 4 * i_corr)) / 2.));
                    IT m = i_corr - l*(l + 1);
                    IT s = (i % (2 * m_extent.final.s + 1)) - m_extent.final.s;

                    return index_type{ l, m, s };
                }

                size_type convert(const index_type& i) const
                {
                    return (i.l * (i.l + 1) + i.m) * (2 * m_extent.final.s + 1) + i.s + m_extent.final.s;
                }

                static index_type convert(typename const IT& max_s, const size_type& i)
                {
                    IT i_corr = static_cast<IT>(i / (2 * max_s + 1)); // truncation on division is intensional (omit std::floor and double up-cast)
                    IT l = static_cast<IT>(std::round((-1. + std::sqrt(1. + 4 * i_corr)) / 2.));
                    IT m = i_corr - l*(l + 1);
                    IT s = (i % (2 * max_s + 1)) - max_s;

                    return index_type{ l, m, s };
                }

            private:

                extent_type m_extent;
                container_type m_data;
            };

            template <typename E1,
                typename E2,
                typename I1,
                typename I2,
                typename IT = decltype(std::declval<typename I1::value_type>() + std::declval<typename I2::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() + std::declval<typename E2::value_type>())
            >
            class Sum : public Expression<Sum<E1, E2, I1, I2>, IT, VT>
            {
            private:

                typedef Extent<IT> extent_type;
                typedef VT         return_type;

            public:

                Sum(Expression<E1, I1, typename E1::value_type> const& u, Expression<E2, I2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] + _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) + _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) + _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E1,
                typename E2,
                typename I1,
                typename I2,
                typename IT = decltype(std::declval<typename I1::value_type>() - std::declval<typename I2::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() - std::declval<typename E2::value_type>())
            >
            class Dif : public Expression<Dif<E1, E2, I1, I2>, IT, VT>
            {
            private:

                typedef Extent<IT> extent_type;
                typedef VT         return_type;

            public:

                Dif(Expression<E1, I1, typename E1::value_type> const& u, Expression<E2, I2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] - _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) - _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) - _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E1,
                typename E2,
                typename I1,
                typename I2,
                typename IT = decltype(std::declval<typename I1::value_type>() * std::declval<typename I2::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>())
            >
            class Mul : public Expression<Mul<E1, E2, I1, I2>, IT, VT>
            {
            private:

                typedef Extent<IT> extent_type;
                typedef VT         return_type;

            public:

                Mul(Expression<E1, I1, typename E1::value_type> const& u, Expression<E2, I2, typename E2::value_type> const& v) : _u(u), _v(v) { assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] * _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) * _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _u.at(i) * _v.at(i); }

            private:

                E1 const& _u;
                E2 const& _v;
            };


            template <typename E,
                typename I,
                typename V,
                typename VT = decltype(std::declval<typename E::value_type>() + std::declval<V>())
            >
            class Scale : public Expression<Scale<E, I, VT>, I, VT>
            {
            private:

                typedef Extent<I> extent_type;
                typedef VT        return_type;

            public:

                Scale(const V& alpha, Expression<E, I, typename E::value_type> const& v) : _alpha(alpha), _v(v) {}

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _alpha * _v[i]; }
                return_type at(size_type i)         const { return _alpha * _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const { return _alpha * _v.at(i); }

            private:

                V _alpha;
                E const& _v;
            };

        } // namespace SpinWeightedSpherical

    } // namespace stl

} // namespace Multipole


////////////////////////////////////////////
// Spherical::Vector non-member operators //
////////////////////////////////////////////

// Unary 
template <typename E, typename T> Multipole::stl::Spherical::Scale<E> const operator+(Multipole::stl::Spherical::Expression<E, T> const& v) { return Multipole::stl::Spherical::Scale<E>(1.0, v); }
template <typename E, typename T> Multipole::stl::Spherical::Scale<E> const operator-(Multipole::stl::Spherical::Expression<E, T> const& v) { return Multipole::stl::Spherical::Scale<E>(-1.0, v); }

// Binary
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Spherical::Sum<E1, E2> const operator+(Multipole::stl::Spherical::Expression<E1, T1> const& u, Multipole::stl::Spherical::Expression<E2, T2> const& v) { return Multipole::stl::Spherical::Sum<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Spherical::Dif<E1, E2> const operator-(Multipole::stl::Spherical::Expression<E1, T1> const& u, Multipole::stl::Spherical::Expression<E2, T2> const& v) { return Multipole::stl::Spherical::Dif<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Spherical::Mul<E1, E2> const operator*(Multipole::stl::Spherical::Expression<E1, T1> const& u, Multipole::stl::Spherical::Expression<E2, T2> const& v) { return Multipole::stl::Spherical::Mul<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Spherical::Div<E1, E2> const operator/(Multipole::stl::Spherical::Expression<E1, T1> const& u, Multipole::stl::Spherical::Expression<E2, T2> const& v) { return Multipole::stl::Spherical::Div<E1, E2>(u, v); }

template <typename E, typename T> Multipole::stl::Spherical::Scale<E> const operator*(const double alpha, Multipole::stl::Spherical::Expression<E, T> const& v) { return Multipole::stl::Spherical::Scale<E>(alpha, v); }
template <typename E, typename T> Multipole::stl::Spherical::Scale<E> const operator*(Multipole::stl::Spherical::Expression<E, T> const& v, const double alpha) { return Multipole::stl::Spherical::Scale<E>(alpha, v); }
/*
namespace Multipole
{
    namespace stl
    {
        namespace Spherical
        {
            double assocLegendrePoly(const IndexPair lm, const double x) { return gsl_sf_legendre_sphPlm(lm.l.i, lm.m.i, x); }

            std::complex<double> harmonic(const IndexPair lm, double theta, double phi) { return std::pow(-1, lm.m.i) * std::sqrt((2 * lm.l.i) / (4 * M_PI) * gsl_sf_fact(lm.l.i - lm.m.i) / gsl_sf_fact(lm.l.i + lm.m.i)) * assocLegendrePoly(lm, std::cos(theta)) * std::exp(std::complex<double>(0, 1) * static_cast<double>(lm.m.i) * phi); }
        }
    }
}
*/
////////////////////////////////////////////////////////
// SpinWeightedSpherical::Vector non-member operators //
////////////////////////////////////////////////////////

// Unary 
template <typename E, typename I, typename T> Multipole::stl::SpinWeightedSpherical::Scale<E, I, T> const operator+(Multipole::stl::SpinWeightedSpherical::Expression<E, I, T> const& v) { return Multipole::stl::SpinWeightedSpherical::Scale<E, I, T>(1.0, v); }
template <typename E, typename I, typename T> Multipole::stl::SpinWeightedSpherical::Scale<E, I, T> const operator-(Multipole::stl::SpinWeightedSpherical::Expression<E, I, T> const& v) { return Multipole::stl::SpinWeightedSpherical::Scale<E, I, T>(-1.0, v); }

// Binary
template <typename E1, typename E2, typename I1, typename I2, typename T1, typename T2> Multipole::stl::SpinWeightedSpherical::Sum<E1, E2, I1, I2, T1, T2> const operator+(Multipole::stl::SpinWeightedSpherical::Expression<E1, I1, T1> const& u, Multipole::stl::SpinWeightedSpherical::Expression<E2, I2, T2> const& v) { return Multipole::stl::Spherical::Sum<E1, E2, I1, I2>(u, v); }
template <typename E1, typename E2, typename I1, typename I2, typename T1, typename T2> Multipole::stl::SpinWeightedSpherical::Dif<E1, E2, I1, I2, T1, T2> const operator-(Multipole::stl::SpinWeightedSpherical::Expression<E1, I1, T1> const& u, Multipole::stl::SpinWeightedSpherical::Expression<E2, I2, T2> const& v) { return Multipole::stl::Spherical::Dif<E1, E2, I1, I2>(u, v); }
template <typename E1, typename E2, typename I1, typename I2, typename T1, typename T2> Multipole::stl::SpinWeightedSpherical::Mul<E1, E2, I1, I2, T1, T2> const operator*(Multipole::stl::SpinWeightedSpherical::Expression<E1, I1, T1> const& u, Multipole::stl::SpinWeightedSpherical::Expression<E2, I2, T2> const& v) { return Multipole::stl::Spherical::Mul<E1, E2, I1, I2>(u, v); }

template <typename E, typename I, typename T> Multipole::stl::SpinWeightedSpherical::Scale<E, I, T> const operator*(const T& alpha, Multipole::stl::SpinWeightedSpherical::Expression<E, I, T> const& v) { return Multipole::stl::Spherical::Scale<E, I, T>(alpha, v); }
template <typename E, typename I, typename T> Multipole::stl::SpinWeightedSpherical::Scale<E, I, T> const operator*(Multipole::stl::SpinWeightedSpherical::Expression<E, I, T> const& v, const T& alpha) { return Multipole::stl::Spherical::Scale<E, I, T>(alpha, v); }

#endif // STLSPHERICALVECTOR_HPP