#ifndef STLFIELD_HPP
#define STLFIELD_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlExpansionIndex.hpp>
#include <Gripper/stl/stlExpansionExtent.hpp>
#include <Gripper/stl/stlExpansionConstantFactor.hpp>
#include <Gripper/stl/stlRuntime.hpp>           // For ability to read global runtime configuration

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <future>                   // Needed for guided parallel processing
#include <utility>                  // Needed for std::pair
#include <complex>                  // Needed for Derive operator


namespace Multipole
{
    namespace stl
    {
        namespace Expansion
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

                typedef Extent                        extent_type;
                typedef Index                         index_type;
                typedef Extent::radial_index_type     radial_index_type;
                typedef Extent::spherical_index_type  spherical_index_type;
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
                value_type  at(const index_type& i) const { return static_cast<ET const&>(*this).at(r, i); }

                // Expression interface

                operator ET&()             { return static_cast<ET&>(*this); }
                operator ET const&() const { return static_cast<const ET&>(*this); }

            protected:

                index_type sizeToIndex(const size_type& pos) const
                {
                    return index_type(static_cast<radial_index_type::value_type>(pos) / (extent().spherical().final() - extent().spherical().initial()).i,
                        static_cast<spherical_index_type::value_type>(pos) % (extent().spherical().final() - extent().spherical().initial()).i);
                }

                size_type indexToSize(const index_type& index) const { return index.radial().r * (extent().spherical().final() - extent().spherical().initial()).i + index.spherical().i; }
            };


            template <typename VT>
            class Field : public Expression<Field<VT>, VT>
            {
            public:

                // Common interface

                Field() { CALLING }
                Field(const Field& in) : m_data(in.m_data), m_extent(in.m_extent) { CALLING }
                Field(Field&& in) : m_data(std::move(in.m_data)), m_extent(std::move(in.m_extent)) { CALLING }
                ~Field() { CALLING }

                //Field& operator=(Field& in) { Gripper::stl::clog << Gripper::DEBUG << __FUNCTION__ << "(Field&) called"; m_extent = in.m_extent; m_data = in.m_data; return *this; }
                //Field& operator=(Field&& src) { Gripper::stl::clog << Gripper::DEBUG << __FUNCTION__ << "(Field&&) called"; m_extent = std::move(src.m_extent); m_data = std::move(src.m_data); return *this; }

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

                value_type&       at(size_type pos)       { return m_data.at(pos); }
                const value_type& at(size_type pos) const { return m_data.at(pos); }

                value_type&       operator[](size_type pos)       { return m_data[pos]; }
                const value_type& operator[](size_type pos) const { return m_data[pos]; }

                value_type*       data()       { return m_data.data(); }    // Returns pointer to the underlying memory
                const value_type* data() const { return m_data.data(); }    // Returns pointer to the underlying memory

                size_type size() const { return m_data.size(); }            // Returns the number of elements inside the Matrix
                void clear() { m_data.clear(); }                            // Clear the contents of the Matrix

                // Lattice interface

                Field(const extent_type& ext) : m_extent(ext), m_data((ext.radial().final() - ext.radial().initial()).r * (ext.spherical().final() - ext.spherical().initial()).i) { CALLING }

                value_type&       at(const index_type& pos)       { return m_data.at(pos.radial().r * (m_extent.spherical().final() - m_extent.spherical().initial()).i + pos.spherical().i); }
                const value_type& at(const index_type& pos) const { return m_data.at(pos.radial().r * (m_extent.spherical().final() - m_extent.spherical().initial()).i + pos.spherical().i); }

                extent_type extent() const { return m_extent; }

                // Expression interface

                template <typename FieldExpr, typename ValueType>
                Field(Expression<FieldExpr, ValueType> const& expr)
                {
                    ENTERING

                    // Extract type from encapsulating expression
                    FieldExpr const& f = expr;

                    Gripper::stl::clog << Gripper::DEBUG << "m_extent = {" << m_extent.radial().final().r - m_extent.radial().initial().r << ", " << m_extent.spherical().final().i - m_extent.spherical().initial().i << "}\texpr.extent() = {" << f.extent().radial().final().r - f.extent().radial().initial().r << ", " << f.extent().spherical().final().i - f.extent().spherical().initial().i << "}";
                    Gripper::stl::clog << Gripper::DEBUG << "m_data.size() = " << m_data.size() << "\texpr.size() = " << expr.size();

                    // Resize if needed
                    if (m_extent != f.extent())
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Field::Field{resize()}";
                        m_extent = f.extent();
                        m_data.resize(f.size());
                    }

                    if (Gripper::stl::Private::g_runtime.getCompilerParallelization())
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Field::Field{evaluating serial}";

                        for (radial_index_type r = m_extent.radial().initial(); r < m_extent.radial().final(); ++r)
                            for (spherical_index_type i = m_extent.spherical().initial(); i < m_extent.spherical().final(); ++i)
                                this->at(index_type(r, i)) = f.at(index_type(r, i));
                    }
                    else
                    {
                        Gripper::stl::clog << Gripper::DEBUG << "Field::Field{evaluating parallel}";

                        std::vector<std::pair<radial_index_type, radial_index_type>> ranges(std::thread::hardware_concurrency());
                        std::vector<std::future<void>> jobs(std::thread::hardware_concurrency());
                        for (std::size_t i = 0; i < ranges.size(); ++i)
                        {
                            ranges.at(i).first  = m_extent.radial().initial() + static_cast<radial_index_type>((i)     * static_cast<size_t>(m_extent.radial().final().r - m_extent.radial().initial().r) / ranges.size());
                            ranges.at(i).second = m_extent.radial().initial() + static_cast<radial_index_type>((i + 1) * static_cast<size_t>(m_extent.radial().final().r - m_extent.radial().initial().r) / ranges.size());
                        }

                        for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                        {
                            for (radial_index_type r = ranges.at(j).first; r < ranges.at(j).second; ++r)
                                for (spherical_index_type i = m_extent.spherical().initial(); i < m_extent.spherical().final(); ++i)
                                    this->at(index_type(r, i)) = f.at(index_type(r, i));
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

                Sum(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { CALLING assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] + _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) + _v.at(i); }

                // Lattice interface

                extent_type extent()                                        const { return _v.extent(); }
                return_type at(radial_index_type r, spherical_index_type i) const { return _u.at(r, i) + _v.at(r, i); }

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

                Dif(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { CALLING assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] - _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) - _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(const index_type& i) const { return _u.at(i) - _v.at(i); }

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

                Mul(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { CALLING assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type pos) const
                {
                    value_type result = 0;
                    index_type index = sizeToIndex(pos);

                    for (Multipole::stl::Gaunt::Matrix::iterator_type it = Gripper::stl::gaunt.marker(i.spherical()).first; it < Gripper::stl::gaunt.marker(i.spherical()).second; ++it)
                    {
                        result += _u[indexToSize(index_type(index.radial(), it->i2))] * _v[indexToSize(index_type(index.radial(), it->i3))] * it->value;
                    }

                    return result;
                }
                return_type at(size_type pos)         const
                {
                    value_type result = 0;
                    index_type index = sizeToIndex(pos);

                    for (Multipole::stl::Gaunt::Matrix::iterator_type it = Gripper::stl::gaunt.marker(i.spherical()).first; it < Gripper::stl::gaunt.marker(i.spherical()).second; ++it)
                    {
                        result += _u.at(indexToSize(index_type(index.radial(), it->i2))) * _v.at(indexToSize(index_type(index.radial(), it->i3))) * it->value;
                    }

                    return result;
                }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(const index_type& i) const
                {
                    value_type result = 0;

                    for (Multipole::stl::Gaunt::Matrix::iterator_type it = Gripper::stl::gaunt.marker(i.spherical()).first; it < Gripper::stl::gaunt.marker(i.spherical()).second; ++it)
                    {
                        result += _u.at(index_type(i.radial(), it->i2)) * _v.at(index_type(i.radial(), it->i3)) * it->value;
                    }

                    return result;
                }

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

                Div(Expression<E1, typename E1::value_type> const& u, Expression<E2, typename E2::value_type> const& v) : _u(u), _v(v) { CALLING assert(u.size() == v.size()); }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _u[i] / _v[i]; }
                return_type at(size_type i)         const { return _u.at(i) / _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(const index_type& i) const { return _u.at(i) / _v.at(i); }

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

                Scale(double alpha, Expression<E, typename E::value_type> const& v) : _alpha(alpha), _v(v) { CALLING }

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return _alpha * _v[i]; }
                return_type at(size_type i)         const { return _alpha * _v.at(i); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(const index_type& i) const { return _alpha + _v.at(i); }

            private:

                double _alpha;
                E const& _v;
            };


            template <typename E>
            class Power : public Expression<Power<E>, typename E::value_type>
            {
            public:

                Power(int power, Expression<E, typename E::value_type> const& v) : _power(power), _v(v) { CALLING }

                // STL interface

                size_type   size()                    const { return _v.size(); }
                value_type  operator[](size_type pos) const { return at(sizeToIndex(pos)); }
                value_type  at(size_type pos)         const { return at(sizeToIndex(pos)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(const index_type& i) const
                {
                    value_type result = 0;

                    for (int p = 1; p < _power; ++p)
                        for (Multipole::stl::Gaunt::Matrix::iterator_type it = Gripper::stl::gaunt.marker(i.spherical()).first; it < Gripper::stl::gaunt.marker(i.spherical()).second; ++it)
                            result += _v.at(index_type(i.radial(), it->i2)) * _v.at(index_type(i.radial(), it->i3)) * it->value;

                    return result;
                }

            private:

                int _power;
                E const& _v;
            };
            
            
            template <typename E>
            class RadialDerivate : public Expression<RadialDerivate<E>, typename E::value_type>
            {
            public:

                RadialDerivate(double dr, Expression<E, typename E::value_type> const& v) : _dr(dr), _v(v) { CALLING }

                // STL interface

                size_type   size()                    const { return _v.size(); }
                value_type  operator[](size_type pos) const { return at(sizeToIndex(pos)); }
                value_type  at(size_type pos)         const { return at(sizeToIndex(pos)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(const index_type& i) const
                {
                    if (extent_type::radial_extent_type(extent().radial().initial(), extent().radial().final() - 2).contains(i.radial())) // regular case
                    {
                        return (8 * (_v.at(i - extent_type(i.radial() + 1, i.spherical())) - _v.at(i - extent_type(i.radial() - 1, i.spherical()))) + _v.at(i - extent_type(i.radial() - 2, i.spherical())) - _v.at(i - extent_type(i.radial() + 2, i.spherical()))) / (12 * _dr);
                    }
                    else
                    {
                        if (i == extent().radial().final() - 2) // O2 case
                        {
                            return (_v.at(i - extent_type(i.radial() + 1, i.spherical())) - _v.at(i - extent_type(i.radial() - 1, i.spherical()))) / (2 * _dr);
                        }
                        
                        if (i == extent().radial().final() - 1) // O1 case
                        {
                            return (_v.at(i) - _v.at(i - extent_type(i.radial() - 1, i.spherical()))) / _dr;
                        }
                    }
                }

            private:

                double _dr;
                E const& _v;
            };
            
            
            template <typename E, int power>
            class PhiDerivate : public Expression<PhiDerivate<E, power>, std::complex<typename E::value_type>>
            {
            public:

                PhiDerivate(Expression<E, typename E::value_type> const& v) : _v(v) { CALLING }

                // STL interface

                size_type   size()                    const { return _v.size(); }
                value_type  operator[](size_type pos) const { return at(sizeToIndex(pos)); }
                value_type  at(size_type pos)         const { return at(sizeToIndex(pos)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(const index_type& i) const { return (power % 2 ? std::complex<value_type>(1, 0) : std::complex<value_type>(0, 1)) * std::pow(Spherical::IndexPair(i.spherical()).m.i, power) * _v.at(i); } // i * m * value

            private:

                E const& _v;
            };

            template <typename E, int power>
            class CPhiDerivate : public Expression<CPhiDerivate<E, power>, typename E::value_type>
            {
            public:

                CPhiDerivate(Expression<E, typename E::value_type> const& v) : _v(v) { CALLING }

                // STL interface

                size_type   size()                    const { return _v.size(); }
                value_type  operator[](size_type pos) const { return at(sizeToIndex(pos)); }
                value_type  at(size_type pos)         const { return at(sizeToIndex(pos)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                value_type  at(const index_type& i) const { return (power % 2 ? value_type(1, 0) : value_type(0, 1)) * std::pow(Spherical::IndexPair(i.spherical()).m.i, power) * _v.at(i); } // i * m * value

            private:

                E const& _v;
            };
            
            namespace Mixed
            {
                template <typename E1, typename E2>
                class MulFR : public Expression<MulFR<E1, E2>, decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>()) return_type;

                public:

                    MulFR(Expression<E1, typename E1::value_type> const& f, Radial::Expression<E2, typename E2::value_type> const& r) : _f(f), _r(r) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _f[i] * _r[i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i]; }
                    return_type at(size_type i)         const { return _f.at(i) * _r.at(i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _f.at(i) * _r.at(i.radial()); }

                private:

                    E1 const& _f;
                    E2 const& _r;
                };


                template <typename E1, typename E2>
                class DivFR : public Expression<DivFR<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivFR(Expression<E1, typename E1::value_type> const& f, Radial::Expression<E2, typename E2::value_type> const& r) : _f(f), _r(r) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _f[i] / _r[i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i]; }
                    return_type at(size_type i)         const { return _f.at(i) / _r.at(i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _f.at(i) / _r.at(i.radial()); }

                private:

                    E1 const& _f;
                    E2 const& _r;
                };


                template <typename E1, typename E2>
                class DivRF : public Expression<DivRF<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivRF(Radial::Expression<E1, typename E1::value_type> const& r, Expression<E2, typename E2::value_type> const& f) : _f(f), _r(r) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _r[i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i] / _f[i]; }
                    return_type at(size_type i)         const { return _r.at(i / (_f.extent().spherical().final() - _f.extent().spherical().initial()).i) / _f.at(i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _r.at(i.radial()) / _f.at(i); }

                private:

                    E1 const& _r;
                    E2 const& _f;
                };


                template <typename E1, typename E2>
                class MulFS : public Expression<MulFS<E1, E2>, decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>()) return_type;

                public:

                    MulFS(Expression<E1, typename E1::value_type> const& f, Spherical::Expression<E2, typename E2::value_type> const& s) : _f(f), _s(s) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _f[i] * _s[i % (_f.extent().spherical().final() - _f.extent().spherical().initial().i)]; }
                    return_type at(size_type i)         const { return _f.at(i) * _s.at(i % (_f.extent().spherical().final() - _f.extent().spherical().initial()).i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _f.at(i) * _s.at(i.spherical()); }

                private:

                    E1 const& _f;
                    E2 const& _s;
                };


                template <typename E1, typename E2>
                class DivFS : public Expression<DivFS<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivFS(Expression<E1, typename E1::value_type> const& f, Spherical::Expression<E2, typename E2::value_type> const& s) : _f(f), _s(s) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _f[i] / _s[i % (_f.extent().spherical().final() - _f.extent().spherical().initial()).i]; }
                    return_type at(size_type i)         const { return _f.at(i) / _s.at(i % (_f.extent().spherical().final() - _f.extent().spherical().initial()).i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _f.at(i) / _s.at(i.spherical()); }

                private:

                    E1 const& _f;
                    E2 const& _s;
                };


                template <typename E1, typename E2>
                class DivSF : public Expression<DivSF<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivSF(Spherical::Expression<E1, typename E1::value_type> const& s, Expression<E2, typename E2::value_type> const& f) : _f(f), _s(s) { CALLING }

                    // STL interface

                    size_type   size()                  const { return _f.size(); }
                    return_type operator[](size_type i) const { return _s[i % (_f.extent().spherical().final() - _f.extent().spherical().initial()).i] / _f[i]; }
                    return_type at(size_type i)         const { return _s.at(i % (_f.extent().spherical().final() - _f.extent().spherical().initial()).i) / _f.at(i); }

                    // Lattice interface

                    extent_type extent()                const { return _f.extent(); }
                    return_type at(const index_type& i) const { return _s.at(i.spherical()) / _f.at(i); }

                private:

                    E1 const& _s;
                    E2 const& _f;
                };


                template <typename E, typename F>
                class AngularIntegral : public Radial::Expression < AngularIntegral<E, F>, typename E::value_type >
                {
                private:

                    typedef decltype(std::declval<Multipole::stl::Radial::Vector<typename E::value_type>>()) return_type;

                public:

                    AngularIntegral(Expression<E, typename E::value_type> const& f, F func) : _f(f), _func(func) {}

                    // STL interface

                    size_type   size()                  const { return (_f.extent().radial().final() - _f.extent().radial().initial()).r; }
                    return_type operator[](size_type i) const { return 0; }
                    return_type at(size_type i)         const { return 0; }

                    // Lattice interface

                    extent_type extent()                 const { return _f.extent().radial(); }
                    return_type at(const index_type& in) const
                    {
                        return_type result = 0;

                        for (spherical_index_type i = _f.extent().spherical().initial(); _f.extent().contains(i); ++i)
                        {
                            for (double theta = -M_PI / 2; theta < M_PI / 2; theta += 0.1)
                            {
                                for (double phi = -M_PI; phi < M_PI; phi += 0.1)
                                {
                                }
                            }
                        }
                    }

                private:

                    E const& _f;
                    F const& _func;
                };

            } // namesapce Mixed

        } // namespace Expansion

    } // namespace stl

} // namespace Multipole

////////////////////////////////
// Field non-member operators //
////////////////////////////////

// Unary
template <typename E, typename T> Multipole::stl::Expansion::Scale<E> const operator+(Multipole::stl::Expansion::Expression<E, T> const& v) { return Multipole::stl::Expansion::Scale<E>(1.0, v); }
template <typename E, typename T> Multipole::stl::Expansion::Scale<E> const operator-(Multipole::stl::Expansion::Expression<E, T> const& v) { return Multipole::stl::Expansion::Scale<E>(-1.0, v); }

// Binary
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Sum<E1, E2> const operator+(Multipole::stl::Expansion::Expression<E1, T1> const& u, Multipole::stl::Expansion::Expression<E2, T2> const& v) { return Multipole::stl::Expansion::Sum<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Dif<E1, E2> const operator-(Multipole::stl::Expansion::Expression<E1, T1> const& u, Multipole::stl::Expansion::Expression<E2, T2> const& v) { return Multipole::stl::Expansion::Dif<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mul<E1, E2> const operator*(Multipole::stl::Expansion::Expression<E1, T1> const& u, Multipole::stl::Expansion::Expression<E2, T2> const& v) { return Multipole::stl::Expansion::Mul<E1, E2>(u, v); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Div<E1, E2> const operator/(Multipole::stl::Expansion::Expression<E1, T1> const& u, Multipole::stl::Expansion::Expression<E2, T2> const& v) { return Multipole::stl::Expansion::Div<E1, E2>(u, v); }

template <typename E, typename T> Multipole::stl::Expansion::Scale<E> const operator*(const double alpha, Multipole::stl::Expansion::Expression<E, T> const& v) { return Multipole::stl::Expansion::Scale<E>(alpha, v); }
template <typename E, typename T> Multipole::stl::Expansion::Scale<E> const operator*(Multipole::stl::Expansion::Expression<E, T> const& v, const double alpha) { return Multipole::stl::Expansion::Scale<E>(alpha, v); }

////////////////////////////////
// Field non-member functions //
////////////////////////////////

namespace Multipole
{
    namespace stl
    {
        namespace Expansion
        {
            template <typename T> Power<Field<T>> pown(int power, const Field<T>& base) { return Power<Field<T>>(power, base); }

            template <typename T> RadialDerivate<Field<T>> d_dr(const Field<T>& in, double dr) { return RadialDerivate<Field<T>>(dr, in); }

            template <typename T, int power> PhiDerivate<Field<T>, power> d_dphi(const Field<T>& in) { return PhiDerivate<Field<T>, power>(in); }

            template <typename T, int power> CPhiDerivate<Field<T>, power> d_dphi_complex(const Field<T>& in) { return CPhiDerivate<Field<T>, power>(in); }
        }
    }
}

/////////////////////
// Mixed operators //
/////////////////////

template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::MulFR<E1, E2> const operator*(Multipole::stl::Expansion::Expression<E1, T1> const& f, Multipole::stl::Radial::Expression<E2, T2> const& r)    { return Multipole::stl::Expansion::Mixed::MulFR<E1, E2>(f, r); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::MulFR<E1, E2> const operator*(Multipole::stl::Radial::Expression<E2, T2> const& r, Multipole::stl::Expansion::Expression<E1, T1> const& f)    { return Multipole::stl::Expansion::Mixed::MulFR<E1, E2>(f, r); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::DivFR<E1, E2> const operator/(Multipole::stl::Expansion::Expression<E1, T1> const& f, Multipole::stl::Radial::Expression<E2, T2> const& r)    { return Multipole::stl::Expansion::Mixed::DivFR<E1, E2>(f, r); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::DivRF<E1, E2> const operator/(Multipole::stl::Radial::Expression<E1, T1> const& r, Multipole::stl::Expansion::Expression<E2, T2> const& f)    { return Multipole::stl::Expansion::Mixed::DivRF<E1, E2>(r, f); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::MulFS<E1, E2> const operator*(Multipole::stl::Expansion::Expression<E1, T1> const& f, Multipole::stl::Spherical::Expression<E2, T2> const& s) { return Multipole::stl::Expansion::Mixed::MulFS<E1, E2>(f, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::MulFS<E1, E2> const operator*(Multipole::stl::Spherical::Expression<E2, T2> const& s, Multipole::stl::Expansion::Expression<E1, T1> const& f) { return Multipole::stl::Expansion::Mixed::MulFS<E1, E2>(f, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::DivFS<E1, E2> const operator/(Multipole::stl::Expansion::Expression<E1, T1> const& f, Multipole::stl::Spherical::Expression<E2, T2> const& s) { return Multipole::stl::Expansion::Mixed::DivFS<E1, E2>(f, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Mixed::DivSF<E1, E2> const operator/(Multipole::stl::Spherical::Expression<E1, T1> const& s, Multipole::stl::Expansion::Expression<E2, T2> const& f) { return Multipole::stl::Expansion::Mixed::DivSF<E1, E2>(s, f); }

#endif // STLFIELD_HPP