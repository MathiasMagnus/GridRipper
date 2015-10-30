#ifndef STLRADIALVECTOR_HPP
#define STLRADIALVECTOR_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlRadialExtent.hpp>
#include <Gripper/stl/stlRuntime.hpp>           // For ability to read global runtime configuration

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <future>                   // Needed for guided parallel processing
#include <utility>                  // Needed for std::pair, std::declval
//#include <type_traits>              // Needed for Expression typedefs


namespace Multipole
{
    namespace stl
    {
        namespace Radial
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

                //Vector& operator=(Vector&) { m_data = in.m_data; return *this; }
                //Vector& operator=(Vector&&) { m_data = std::move(in.m_data); return *this; }

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

                value_type& at(const index_type& pos) { return m_data.at(pos.r); }
                const value_type& at(const index_type& pos) const { return m_data.at(pos.r); }

                extent_type extent() const { return m_extent; }

                // Exrpession interface

                template <typename VecExpr, typename ValueType>
                Vector(Expression<VecExpr, ValueType> const& expr)
                {
                    ENTERING

                    // Extract type from encapsulating expression
                    VecExpr const& v = expr;

                    Gripper::stl::clog << Gripper::DEBUG << "m_extent = {" << m_extent.final().r - m_extent.initial().r << "}\texpr.extent() = {" << v.extent().final().r - v.extent().initial().r << "}";
                    Gripper::stl::clog << Gripper::DEBUG << "m_data.size() = " << m_data.size() << "\texpr.size() = " << expr.size();

                    // Resize if needed
                    if (m_extent != v.extent())
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Radial::Vector::Vector{resize()}";
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
                            ranges.at(i).first  = m_extent.initial() + (m_extent.final() - m_extent.initial()) * static_cast<index_type::value_type>((i    ) / ranges.size());
                            ranges.at(i).second = m_extent.initial() + (m_extent.final() - m_extent.initial()) * static_cast<index_type::value_type>((i + 1) / ranges.size());
                        }

                        for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                        {
                            for (index_type i = ranges.at(j).first; i < ranges.at(j).second; ++i) this->at(i) = v.at(i);
                        });
                        for (auto& job : jobs) job.wait();
                    }

                    LEAVING
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

        } // namespace Radial

    } // namespace stl

} // namespace Multipole


/////////////////////////////////////////
// Radial::Vector non-member operators //
/////////////////////////////////////////

// Unary 
template <typename E, typename T> Multipole::stl::Radial::Scale<E> const operator+(Multipole::stl::Radial::Expression<E, T> const& v) { return Multipole::stl::Radial::Scale<E>(1.0, v); }
template <typename E, typename T> Multipole::stl::Radial::Scale<E> const operator-(Multipole::stl::Radial::Expression<E, T> const& v) { return Multipole::stl::Radial::Scale<E>(-1.0, v); }

// Binary
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Radial::Sum<E1, E2> const operator+(Multipole::stl::Radial::Expression<E1, T1> const& u, Multipole::stl::Radial::Expression<E2, T2> const& v) { return Multipole::stl::Radial::Sum<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Radial::Dif<E1, E2> const operator-(Multipole::stl::Radial::Expression<E1, T1> const& u, Multipole::stl::Radial::Expression<E2, T2> const& v) { return Multipole::stl::Radial::Dif<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Radial::Mul<E1, E2> const operator*(Multipole::stl::Radial::Expression<E1, T1> const& u, Multipole::stl::Radial::Expression<E2, T2> const& v) { return Multipole::stl::Radial::Mul<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Radial::Div<E1, E2> const operator/(Multipole::stl::Radial::Expression<E1, T1> const& u, Multipole::stl::Radial::Expression<E2, T2> const& v) { return Multipole::stl::Radial::Div<E1, E2>(u, v); }

template <typename E, typename T> Multipole::stl::Radial::Scale<E> const operator*(const double alpha, Multipole::stl::Radial::Expression<E, T> const& v) { return Multipole::stl::Radial::Scale<E>(alpha, v); }
template <typename E, typename T> Multipole::stl::Radial::Scale<E> const operator*(Multipole::stl::Radial::Expression<E, T> const& v, const double alpha) { return Multipole::stl::Radial::Scale<E>(alpha, v); }

#endif // STLRADIALVECTOR_HPP