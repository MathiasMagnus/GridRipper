#ifndef STLGAUNTMATRIX_HPP
#define STLGAUNTMATRIX_HPP

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
// BONUS SPHERICAL INCLUDE AT THE END OF FILE
#include <Gripper/stl/stlMultipoleDefs.hpp>     // Forward declaration of Multipole types
#include <Gripper/stl/stlSphericalIndex.hpp>
#include <Gripper/stl/stlGauntIndex.hpp>
#include <Gripper/stl/stlGauntCoefficient.hpp>

#include <Gripper/stl/stlMultipoleTypes.hpp>
#include <Gripper/stl/stlSphericalExtent.hpp>
#include <Gripper/stl/stlSphericalVector.hpp>

// Standard C++ includes
#include <vector>                   // Needed for internal storage
#include <ostream>                  // Needed for output
#include <sstream>                  // Needed for output
#include <numeric>                  // std::accumulate
#include <algorithm>                // std::find_if
#include <cassert>                  // assert


namespace Multipole
{
    namespace stl
    {
        namespace Gaunt
        {
            EXPORT std::ostream& operator<<(std::ostream& os, const Matrix& mat);


            class EXPORT Matrix
            {
                friend EXPORT std::ostream& operator<<(std::ostream& os, const Matrix& mat);

            public:

                // Common typedefs

                typedef Coefficient                            value_type;

                // STL typedefs

                typedef std::vector<value_type>                container_type;
                typedef container_type::size_type              size_type;
                typedef container_type::iterator               iterator_type;
                typedef container_type::const_iterator         const_iterator_type;
                typedef container_type::reverse_iterator       reverse_iterator_type;
                typedef container_type::const_reverse_iterator const_reverse_iterator_type;

                // Private typedefs

                typedef std::pair<value_type, value_type>      value_with_error;

                // Common interface

                Matrix();
                Matrix(const Matrix& in);
                Matrix(Matrix&& in);
                ~Matrix();

                Matrix& operator=(Matrix&);
                Matrix& operator=(Matrix&&);

                // STL interface

                iterator_type begin();
                const_iterator_type begin() const;
                const_iterator_type cbegin() const;

                iterator_type end();
                const_iterator_type end() const;
                const_iterator_type cend() const;

                reverse_iterator_type rbegin();
                const_reverse_iterator_type rbegin() const;
                const_reverse_iterator_type crbegin() const;

                reverse_iterator_type rend();
                const_reverse_iterator_type rend() const;
                const_reverse_iterator_type crend() const;

                value_type at(size_type pos);                           // Returns specified element of the Matrix with bounds checking
                const value_type& at(size_type pos) const;              // Returns specified const element of the Matrix with bounds checking
                value_type operator[](size_type pos);                   // Returns specified element of the Matrix
                const value_type& operator[](size_type pos) const;      // Returns specified const element of the Matrix

                value_type* data();                                     // Returns pointer to the underlying memory
                const value_type* data() const;                         // Returns pointer to the underlying memory

                size_type size() const;                                 // Returns the number of elements inside the Matrix
                void clear();                                           // Clear the contents of the Matrix

                // Lattice interface

                Matrix(const Multipole::Gaunt::Matrix& in);

                std::pair<iterator_type, iterator_type>& marker(const Index& index);

            private:

                Spherical::Index m_L_max;                           // L_max the matrix has been created for
                container_type m_data;                              // Container to store values

                std::vector<std::pair<iterator_type, iterator_type>> m_markers;
            };

        } // namespace Gaunt

        namespace SpinWeightedGaunt
        {
            template <typename IndexInternalType = SpinWeightedSpherical::Index<>::value_type, typename ValueType = double>
            using Matrix = Multipole::SpinWeightedGaunt::Matrix<IndexInternalType, ValueType>;
            
            template <typename E1,
                typename E2,
                typename G3,
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
            class Contraction : public SpinWeightedSpherical::Expression<Contraction<E1, E2, G3>, IT, VT>
            {
            public:

                typedef SpinWeightedSpherical::Extent<IT> extent_type;
                typedef VT                                return_type;

                Contraction(const SpinWeightedSpherical::Expression<E1, typename E1::index_type::value_type, typename E1::value_type>& u,
                            const SpinWeightedSpherical::Expression<E2, typename E1::index_type::value_type, typename E2::value_type>& v,
                            const G3& g) :
                    _u(u), _v(v), _g(g) {} // { (assert(u.extent() == v.extent()) && (v.extent() == extent_type({ 0, 0, -g.getS() }, { g.getL(), g.getL(), g.getS() }))); };

                // STL interface

                size_type   size()                  const { return _v.size(); }
                return_type operator[](size_type i) const { return this->at(SpinWeightedSpherical::Vector<IT, VT>::convert(_g.getS(), i)); }
                return_type at(size_type i)         const { return this->at(SpinWeightedSpherical::Vector<IT, VT>::convert(_g.getS(), i)); }

                // Lattice interface

                extent_type extent()                const { return _v.extent(); }
                return_type at(index_type i)        const
                {
                    // Find inside the sorted sparse matrix the first element that has the same 1st index as the result index
                    auto first = std::find_if(_g.cbegin(), _g.cend(), [&](const typename G3::value_type& elem) { return elem.first.at(0) == i; });
                    // Find the first element from here that does not
                    auto last = std::find_if_not(first, _g.cend(), [&](const typename G3::value_type& elem) { return elem.first.at(0) == i; });
                    // Sum up the result
                    return std::accumulate(first, last, static_cast<return_type>(0), [&](const return_type& accum, const typename G3::value_type& elem) { return accum + _u.at(elem.first.at(1)) * _v.at(elem.first.at(2)) * elem.second; });
                    //{
                    //    if ((elem.first.at(1) == index_type{ 3, -3, 1 }) && (elem.first.at(2) == index_type{ 3, 3, -1 }))
                    //        std::cout << "Match found with coeff = " << elem.second << std::endl;
                    //
                    //    return accum + _u.at(elem.first.at(1)) * _v.at(elem.first.at(2)) * elem.second;
                    //});
                }

            private:

                const E1& _u;
                const E2& _v;
                const G3& _g;
            };
            
        } // namespace SpinWeightedGaunt

    } // namespace stl

} // namespace Multipole


namespace Gripper
{
    namespace stl
    {
        EXPORT extern Multipole::stl::Gaunt::Matrix gaunt;
    }
} // namespace Gripper

namespace Multipole
{
    namespace stl
    {
        namespace SpinWeightedGaunt
        {
            template <typename E1,
                typename E2,
                typename G3,
                typename IT = decltype(std::declval<typename E1::index_type::value_type>() + std::declval<typename E2::index_type::value_type>() + std::declval<typename G3::index_type::value_type>()),
                typename VT = decltype(std::declval<typename E1::value_type>() * std::declval<typename E2::value_type>() * std::declval<typename G3::mapped_type>())>
                auto contract(const G3& gaunt, const E1& lhs, const E2& rhs) { return Contraction<E1, E2, G3>(lhs, rhs, gaunt); };
        } // namespace SpinWeightedGaunt
    }
}

#endif // STLGAUNT_HPP