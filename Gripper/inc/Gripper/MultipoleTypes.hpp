#pragma once

// Reading settings from CMake build configuration
#include <Gripper/Gripper_Config.hpp>
#include <Gripper/Gripper_Export.hpp>

// Gripper includes
#include <Gripper/GauntMath.hpp>     // cmath with _USE_MATH_DEFINES + utlity
#include <Gripper/MultipoleDefs.hpp> // Forward declaration of Multipole types
#include <Gripper/Logger.hpp>        // Needed for logging

// GSL includes
#include <gsl/gsl_sf_coupling.h>

// Standard C++ includes
#include <vector>                   // Needed for internal storage and providing results
#include <map>                      // Needed for internal storage and providing results
#include <tuple>                    // Needed for internal storage and providing results
#include <future>                   // Needed for asynchronous operation
#include <mutex>                    // Needed for synchronization
#include <queue>                    // Needed for task scheduling
#include <fstream>                  // Needed for file I/O
#include <sstream>                  // Needed for log output
#include <ostream>                  // Needed for log STL interface
#include <chrono>                   // Needed for profiling
#include <memory>                   // Needed for unique_ptr
#include <limits>

/*
namespace Multipole
{
    namespace Spherical
    {
        template <typename ArithmeticType = std::int8_t>
        class Index
        {
        public:

            typedef ArithmeticType value_type;

            Index() = default;
            Index(const Index&) = default;
            Index(Index&&) = default;
            ~Index() = default;

            Index& operator=(const Index& in) = default;
            Index& operator=(Index&& in) = default;

            template <typename TT>
            Index(const TT l_in, const TT m_in) : l(l_in), m(m_in) {}
            Index(std::initializer_list<value_type> init)
            {
                //static_assert(init.size() == 2, "Initializer-list of Index must have 2 elements: {l, m}");
                l = *init.begin();
                m = *(init.begin() + 1);
            }

            bool operator<(const Index& rhs) const
            {
                if (l < rhs.l) return true;
                if (rhs.l < l) return false;

                if (m < rhs.m) return true;
                return false;
            }

            bool operator==(const Index& rhs) const
            {
                return (l == rhs.l) && (m == rhs.m);
            }

            bool operator<=(const Index& rhs) const
            {
                return (*this < rhs) || (*this == rhs);
            }

            value_type l;
            value_type m;
        };

        template <typename AT>
        Index<AT> next(const Index<AT>& index)
        {
            bool rewind_m_and_step_l = index.m == index.l;

            return Index<AT>
            {
                rewind_m_and_step_l ? index.l + 1 : index.l,
                rewind_m_and_step_l ? -(index.l + 1) : index.m + 1
            };
        }

    } // namespace Spherical

    namespace Gaunt
    {
        template <typename IndexInternalType, typename ValueType>
        class Matrix;

        template <typename IndexInternalType, typename ValueType>
        std::ostream& operator<<(std::ostream& os, const Matrix<IndexInternalType, ValueType>& mat)
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

        template <typename IndexInternalType = Spherical::Index<>::value_type, typename ValueType = double>
        class Matrix
        {
        public:

            typedef SpinWeightedSpherical::Index<IndexInternalType>     index_type;
            typedef ValueType                                           mapped_type;
            typedef std::array<index_type, 3>                           key_type;
            typedef std::pair<const key_type, mapped_type>              value_type;
            typedef std::vector<value_type>                             container_type;
            typedef typename container_type::size_type                  size_type;
            typedef typename container_type::iterator                   iterator_type;
            typedef typename container_type::const_iterator             const_iterator_type;
            typedef typename container_type::reverse_iterator           reverse_iterator_type;
            typedef typename container_type::const_reverse_iterator     const_reverse_iterator_type;
            typedef std::pair<const_iterator_type, const_iterator_type> marker_type;
            typedef std::vector<marker_type>                            accelerator_structure_type;

            Matrix() = default;
            Matrix(const Matrix&) = default;
            Matrix(Matrix&&) = default;
            ~Matrix() = default;

            Matrix(const IndexInternalType l_max) : m_data(), m_l_max(l_max) { compute(); }

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
            IndexInternalType getL() const { return m_l_max; }  // Returns L_max the matrix has been created for
            const marker_type& get_marker(const index_type& index) const { return m_accel.at(Spherical::Index<IndexInternalType>::convert(index)); }

        private:

            container_type m_data;                              // Container to store values
            accelerator_structure_type m_accel;                 // Container to store accelerator markers
            IndexInternalType m_l_max;                          // l_max the matrix has been created for

            void compute();                                     // Compute entire matrix

            container_type computeL(Spherical::Index,
                std::vector<container_type>*);                  // Compute coefficients with given L1 value
        };

    } // namespace Gaunt

} // namespace Multipole

///////////////////////////////////////////////////////
// SpinWeightedSpherical::Index non-member operators //
///////////////////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::SpinWeightedSpherical::Index<AT>& lhs, const Multipole::SpinWeightedSpherical::Index<AT>& rhs) { return (lhs.l == rhs.l) && (lhs.m == rhs.m) && (lhs.s == rhs.s); }
template <typename AT> bool operator!=(const Multipole::SpinWeightedSpherical::Index<AT>& lhs, const Multipole::SpinWeightedSpherical::Index<AT>& rhs) { return (lhs.l != rhs.l) || (lhs.m != rhs.m) || (lhs.s != rhs.s); }
*/