#ifndef STLEXPANSIONCONSTANTFACTOR_HPP
#define STLEXPANSIONCONSTANTFACTOR_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlExpansionExtent.hpp>
#include <Gripper/stl/stlRuntime.hpp>           // For ability to read global runtime configuration

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <future>                   // Needed for guided parallel processing
#include <utility>                  // Needed for std::pair


namespace Multipole
{
    namespace stl
    {
        namespace Expansion
        {
            namespace Constant
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

                    typedef Extent                                      extent_type;
                    typedef Extent::index_type                          index_type;
                    typedef typename extent_type::radial_index_type     radial_index_type;
                    typedef typename extent_type::spherical_index_type  spherical_index_type;
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
                };


                template <typename VT>
                class Field : public Expression<Field<VT>, VT>
                {
                public:

                    // Common interface

                    Field() {}
                    Field(const Field& in) : m_data(in.m_data), m_extent(in.m_extent) {}
                    Field(Field&& in) : m_data(std::move(in.m_data)), m_extent(std::move(in.m_extent)) {}
                    ~Field() {}

                    Field& operator=(Field&) { m_radial_size = in.m_radial_size; m_spherical_size = in.m_spherical_size; m_data = in.m_data; return *this; }
                    Field& operator=(Field&&) { m_radial_size = in.m_radial_size; m_spherical_size = in.m_spherical_size; m_data = std::move(in.m_data); return *this; }

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

                    value_type& at(size_type pos) { return m_data.at(pos); }
                    const value_type& at(size_type pos) const { return m_data.at(pos); }
                    value_type& operator[](size_type pos) { return m_data[pos]; }
                    const value_type& operator[](size_type pos) const { return m_data[pos]; }

                    value_type* data() { return m_data.data(); }                  // Returns pointer to the underlying memory
                    const value_type* data() const { return m_data.data(); }      // Returns pointer to the underlying memory

                    size_type size() const { return m_data.size(); }             // Returns the number of elements inside the Matrix
                    void clear() { m_data.clear(); }                             // Clear the contents of the Matrix

                    // Lattice interface

                    Field(const extent_type& ext) : m_extent(ext), m_data((ext.radial().final() - ext.radial().initial()).r * (ext.spherical().final() - ext.spherical().initial()).i) {}

                    value_type&       at(const index_type& pos)       { return m_data.at(pos.radial().r * (m_extent.spherical().final() - m_extent.spherical().initial()).i + pos.spherical().i); }
                    const value_type& at(const index_type& pos) const { return m_data.at(pos.radial().r * (m_extent.spherical().final() - m_extent.spherical().initial()).i + pos.spherical().i); }

                    extent_type extent() const { return m_extent; }

                    // Expression interface

                    template <typename FieldExpr, typename ValueType>
                    Field(Expression<FieldExpr, ValueType> const& expr)
                    {
                        // Extract type from encapsulating expression
                        FieldExpr const& f = expr;

                        // Resize if needed
                        if (m_extent != f.extent())
                        {
                            m_extent = f.extent();
                            m_data.resize(f.size());
                        }

                        if (Gripper::stl::Private::g_runtime.getCompilerParallelization())
                        {
                            for (radial_index_type r = m_extent.radial().initial(); r < m_extent.radial().final(); ++r)
                                for (spherical_index_type i = m_extent.spherical().initial(); i < m_extent.spherical().final(); ++i)
                                    this->at(index_type(r, i)) = f.at(index_type(r, i));
                        }
                        else
                        {
                            std::vector<std::pair<radial_index_type, radial_index_type>> ranges(std::thread::hardware_concurrency());
                            std::vector<std::future<void>> jobs(std::thread::hardware_concurrency());
                            for (std::size_t i = 0; i < ranges.size(); ++i)
                            {
                                ranges.at(i).first = m_extent.radial().initial() + (m_extent.radial().final() - m_extent.radial().initial()) * (i) / ranges.size();
                                ranges.at(i).second = m_extent.radial().initial() + (m_extent.radial().final() - m_extent.radial().initial()) * (i + 1) / ranges.size();
                            }

                            for (std::size_t j = 0; j < jobs.size(); ++j) jobs.at(j) = std::async(std::launch::async, [=]()
                            {
                                for (radial_index_type r = ranges.at(j).first; r < ranges.at(j).second; ++r)
                                    for (spherical_index_type i = m_extent.spherical().initial(); i < m_extent.spherical().final(); ++i)
                                        this->at(index_type(r, i)) = f.at(index_type(r, i));
                            });
                            for (auto& job : jobs) job.wait();
                        }
                    }

                private:

                    container_type m_data;
                    extent_type m_extent;
                };


                template <typename E1, typename E2>
                class MulRS : public Expression<MulRS<E1, E2>, decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>()) return_type;

                public:

                    MulRS(Radial::Expression<E1, typename E1::value_type> const& r, Spherical::Expression<E2, typename E2::value_type> const& s) : _r(r), _s(s) {}

                    // STL interface

                    size_type   size()                  const { return _r.size() * _s.size(); }
                    return_type operator[](size_type i) const { return _r[i / _s.size()] * _s[i % _s.size()]; }
                    return_type at(size_type i)         const
                    {
                        //Gripper::stl::clog << Gripper::DEBUG << "_r.at(" << i << " / " << _s.size() << ") * _s.at(" << i << " % " << _s.size() << ") = " << _r.at(i / _s.size()) * _s.at(i % _s.size());
                        return _r.at(i / _s.size()) * _s.at(i % _s.size());
                    }

                    // Lattice interface

                    extent_type extent()                const { return extent_type(_r.extent(), _s.extent()); }
                    return_type at(const index_type& i) const { return _r.at(i.radial()) * _s.at(i.spherical()); }

                private:

                    E1 const& _r;
                    E2 const& _s;
                };


                template <typename E1, typename E2>
                class DivRS : public Expression<DivRS<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivRS(Radial::Expression<E1, typename E1::value_type> const& r, Spherical::Expression<E2, typename E2::value_type> const& s) : _r(r), _s(s) {}

                    // STL interface

                    size_type   size()                  const { return _r.size() * _s.size(); }
                    return_type operator[](size_type i) const { return _r[i / _s.size()] / _s[i % _s.size()]; }
                    return_type at(size_type i)         const { return _r.at(i / _s.size()) / _s.at(i % _s.size()); }

                    // Lattice interface

                    extent_type extent()                const { return extent_type(_r.extent(), _s.extent()); }
                    return_type at(const index_type& i) const { return _r.at(i.radial()) / _s.at(i.spherical()); }

                private:

                    E1 const& _r;
                    E2 const& _s;
                };


                template <typename E1, typename E2>
                class DivSR : public Expression<DivSR<E1, E2>, decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>())>
                {
                private:

                    typedef decltype(std::declval<typename E1::value_type>() / std::declval<typename E2::value_type>()) return_type;

                public:

                    DivSR(Spherical::Expression<E1, typename E1::value_type> const& s, Radial::Expression<E2, typename E2::value_type> const& r) : _r(r), _s(s) {}

                    // STL interface

                    size_type   size()                  const { return _r.size() * _s.size(); }
                    return_type operator[](size_type i) const { return _s[i % _s.size()] / _r[i / _s.size()]; }
                    return_type at(size_type i)         const { return _s.at(i % _s.size()) / _r.at(i / _s.size()); }

                    // Lattice interface

                    extent_type extent()                const { return extent_type(_r.extent(), _s.extent()); }
                    return_type at(const index_type& i) const { return _s.at(i.spherical()) / _r.at(i.radial()); }

                private:

                    E1 const& _s;
                    E2 const& _r;
                };

            } // namespace Constant

        } // namespace Expansion

    } // namespace stl

} // namespace Multipole

/////////////////////
// Mixed operators //
/////////////////////

template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Constant::MulRS<E1, E2> const operator*(Multipole::stl::Radial::Expression<E1, T1> const& r, Multipole::stl::Spherical::Expression<E2, T2> const& s) { return Multipole::stl::Expansion::Constant::MulRS<E1, E2>(r, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Constant::MulRS<E1, E2> const operator*(Multipole::stl::Spherical::Expression<E2, T2> const& s, Multipole::stl::Radial::Expression<E1, T1> const& r) { return Multipole::stl::Expansion::Constant::MulRS<E1, E2>(r, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Constant::DivRS<E1, E2> const operator/(Multipole::stl::Radial::Expression<E1, T1> const& r, Multipole::stl::Spherical::Expression<E2, T2> const& s) { return Multipole::stl::Expansion::Constant::DivRS<E1, E2>(r, s); }
template <typename E1, typename E2, typename T1, typename T2> Multipole::stl::Expansion::Constant::DivSR<E1, E2> const operator/(Multipole::stl::Spherical::Expression<E1, T1> const& s, Multipole::stl::Radial::Expression<E2, T2> const& r) { return Multipole::stl::Expansion::Constant::DivSR<E1, E2>(s, r); }

#endif // STLEXPANSIONCONSTANTFACTOR_HPP