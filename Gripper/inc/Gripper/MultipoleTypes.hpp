#ifndef MULTIPOLETYPES_HPP
#define MULTIPOLETYPES_HPP

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


namespace Multipole
{
    namespace Radial
    {
        class EXPORT Index
        {
        public:
            Index();
            Index(const Index& in);
            Index(const IndexType& in);
            ~Index();

            explicit operator const IndexType&();

            Index& operator+=(const Index& rhs);
            Index& operator-=(const Index& rhs);
            Index& operator*=(const Index& rhs);
            Index& operator/=(const Index& rhs);

            Index operator++(int rhs);
            Index& operator++();
            Index operator--(int rhs);
            Index& operator--();

            IndexType r;
        };

    } // namespace Radial

    namespace Spherical
    {
        class EXPORT Index
        {
        public:
            Index();
            Index(const Index& in);
            Index(const IndexType& in);
            ~Index();

            explicit operator const IndexType&();
            
            Index& operator+=(const Index& rhs);
            Index& operator-=(const Index& rhs);
            Index& operator*=(const Index& rhs);
            Index& operator/=(const Index& rhs);

            Index operator++(int rhs);
            Index& operator++();
            Index operator--(int rhs);
            Index& operator--();

            IndexType i;
        };

        class EXPORT IndexPair
        {
        public:
            IndexPair();
            IndexPair(const IndexPair& in);
            IndexPair(const Index& l_in, const Index& m_in);
            IndexPair(const Gaunt::Index& i_in);
            ~IndexPair();

            Index l;
            Index m;
        };

    } // namespace Spherical

	namespace Spin
	{
		class EXPORT Index
		{
		public:
			Index();
			Index(const Index& in);
			Index(const IndexType& in);
			~Index();

			explicit operator const IndexType&();

			Index& operator+=(const Index& rhs);
			Index& operator-=(const Index& rhs);
			Index& operator*=(const Index& rhs);
			Index& operator/=(const Index& rhs);

			Index operator++(int rhs);
			Index& operator++();
			Index operator--(int rhs);
			Index& operator--();

			IndexType i;
		};

	} // Spin

    namespace Gaunt
    {
        class EXPORT Index
        {
        public:
            Index();
            Index(const Index& in);
            Index(const IndexType& i_in);
            Index(const Spherical::IndexPair& lm_in);
            ~Index();

            explicit operator const IndexType();

            Index& operator+=(const Index& rhs);
            Index& operator-=(const Index& rhs);
            Index& operator*=(const Index& rhs);
            Index& operator/=(const Index& rhs);

            Index operator++(int rhs);
            Index& operator++();
            Index operator--(int rhs);
            Index& operator--();

            IndexType i;
        };

        class EXPORT Coefficient
        {
        public:
            Coefficient();
            Coefficient(const Coefficient& in);
            Coefficient(const Index& i1_in, const Index& i2_in, const Index& i3_in, const ValueType& value_in);
            Coefficient(const Spherical::IndexPair& l1m1_in, const Spherical::IndexPair& l2m2_in, const Spherical::IndexPair& l3m3_in, const ValueType& value_in);
            ~Coefficient();

            Spherical::IndexPair l1m1;
            Spherical::IndexPair l2m2;
            Spherical::IndexPair l3m3;
            ValueType value;
        };

        EXPORT std::ostream& operator<<(std::ostream& os, const Matrix& mat);

        class EXPORT Matrix
        {
            friend EXPORT std::ostream& operator<<(std::ostream& os, const Matrix& mat);

        public:

            Matrix();
            Matrix(const Matrix& in);
            Matrix(Matrix&& in);
            Matrix(Spherical::Index& L_max_in);
            ~Matrix();

            Matrix& operator=(Matrix&);
            Matrix& operator=(Matrix&&);

            typedef std::vector<Coefficient>                container_type;
            typedef std::pair<ValueType, ValueType>         value_with_error;
            typedef container_type::size_type               size_type;
            typedef container_type::iterator                iterator_type;
            typedef container_type::const_iterator          const_iterator_type;
            typedef container_type::reverse_iterator        reverse_iterator_type;
            typedef container_type::const_reverse_iterator  const_reverse_iterator_type;

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

            Coefficient at(size_type pos);                      // Returns specified element of the Matrix with bounds checking
            const Coefficient& at(size_type pos) const;         // Returns specified const element of the Matrix with bounds checking
            Coefficient operator[](size_type pos);              // Returns specified element of the Matrix
            const Coefficient& operator[](size_type pos) const; // Returns specified const element of the Matrix

            Coefficient* data();                                // Returns pointer to the underlying memory
            const Coefficient* data() const;                    // Returns pointer to the underlying memory

            size_type size() const;                             // Returns the number of elements inside the Matrix
            void clear();                                       // Clear the contents of the Matrix

            Spherical::Index getL() const;                      // Returns L_max the matrix has been created for

        private:

            Spherical::Index m_L_max;                           // L_max the matrix has been created for
            container_type m_data;                              // Container to store values

            void compute();                                     // Compute entire matrix

            void computeL(Spherical::Index,
                std::vector<container_type>*);                  // Compute coefficients with given L1 value
            /*
            ValueType wigner3jSymbol(
                Spherical::IndexPair& LM1,
                Spherical::IndexPair& LM2,
                Spherical::IndexPair& LM3);                     // Compute Wigner3J-symbol of given spherical indice

            bool triangleSelectionFail(
                Spherical::IndexPair&,
                Spherical::IndexPair&,
                Spherical::IndexPair&);                         // Helper function

            bool mSelectionFail(
                Spherical::IndexPair&,
                Spherical::IndexPair&,
                Spherical::IndexPair&);                         // Helper function

            bool choose(
                Spherical::Index&,
                Spherical::Index&,
                value_with_error&);                             // Helper function

            ValueType factorial(int i);                         // Helper function
            */
        };

    } // namespace Gaunt

	namespace SpinWeightedSpherical
	{
		template <typename ArithmeticType = std::int8_t>
		struct Index
		{
			typedef ArithmeticType value_type;

			Index() = default;
			Index(const Index&) = default;
			Index(Index&&) = default;
			~Index() = default;

            Index& operator=(const Index& in) = default;

			template <typename TT>
			Index(const TT l_in, const TT m_in, const TT s_in) : l(l_in), m(m_in), s(s_in) {}
			Index(std::initializer_list<value_type> init)
			{
				//static_assert(init.size() == 3, "Initializer-list of Index must have 3 elements: {l, m, s}");
				l = *init.begin();
				m = *(init.begin() + 1);
				s = *(init.begin() + 2);
			}

			bool operator<(const Index& rhs) const
			{
				if (l < rhs.l) return true;
				if (rhs.l < l) return false;

				if (m < rhs.m) return true;
				if (rhs.m < m) return false;

				if (s < rhs.s) return true;
				return false;
			}

            bool operator==(const Index& rhs) const
            {
                return (l == rhs.l) && (m == rhs.m) && (s == rhs.s);
            }

            bool operator<=(const Index& rhs) const
            {
                return (*this < rhs) || (*this == rhs);
            }

            static Index convert(const std::size_t& i, const value_type& max_s)
            {
                value_type i_corr = static_cast<value_type>(i / (2 * max_s + 1)); // truncation on division is intensional (omit std::floor and double up-cast)
                value_type l = static_cast<value_type>(std::round((-1. + std::sqrt(1. + 4 * i_corr)) / 2.));
                value_type m = i_corr - l*(l + 1);
                value_type s = (i % (2 * max_s + 1)) - max_s;

                return Index{ l, m, s };
            }

            static std::size_t convert(const Index& i, const value_type& max_s)
            {
                return (i.l * (i.l + 1) + i.m) * (2 * max_s + 1) + i.s + max_s;
            }

			value_type l;
			value_type m;
			value_type s;
		};

        template <typename AT>
        Index<AT> next(const Index<AT>& index, const AT& max_l, const AT& max_s)
        {
            bool rewind_s_and_step_m = index.s == max_s;
            bool rewind_m_and_step_l = rewind_s_and_step_m && (index.m == index.l);

            return Index<AT>
            {
                rewind_m_and_step_l ? index.l + 1 : index.l,
                rewind_s_and_step_m ? (rewind_m_and_step_l ? -(index.l + 1) : index.m + 1) : index.m,
                rewind_s_and_step_m ? -max_s : index.s + 1
            };
        }

	} // namespace SpinWeightedSpherical

	namespace SpinWeightedGaunt
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
                    os << "{ " << index.l << " }{ " << index.m << " }{ " << index.s << " }\t";
				os << "= " << elem.second << "\n";
			}

			os.flush();

			return os;
		}

		template <typename IndexInternalType = SpinWeightedSpherical::Index<>::value_type, typename ValueType = double>
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

			Matrix(const IndexInternalType l_max, const IndexInternalType s_max) : m_data(), m_l_max(l_max), m_s_max(s_max) { compute(); }

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
			IndexInternalType getS() const { return m_s_max; }  // Returns S_max the matrix has been created for
            const marker_type& get_marker(const index_type& index) const { return m_accel.at(SpinWeightedSpherical::Index<IndexInternalType>::convert(index, m_s_max)); }

		private:

			container_type m_data;                              // Container to store values
            accelerator_structure_type m_accel;                 // Container to store accelerator markers
			IndexInternalType m_l_max;                          // l_max the matrix has been created for
			IndexInternalType m_s_max;                          // s_max the matrix has been created for			

			void compute()                                      // Compute entire matrix
			{
				ENTERING

				Gripper::clog << Gripper::LogLevel::INFO << "Computing Gaunt coefficients up to l_max,s_max = " << static_cast<int>(m_l_max) << "," << static_cast<int>(m_s_max);

				auto start = std::chrono::high_resolution_clock::now();

				std::vector<std::future<container_type>> jobs;

				for (IndexInternalType l1 = 0; l1 <= m_l_max; ++l1)
					for (IndexInternalType m1 = -m_l_max; m1 <= m_l_max; ++m1)
						jobs.push_back(std::async(std::launch::async, &Multipole::SpinWeightedGaunt::Matrix<IndexInternalType, ValueType>::computeLM, this, l1, m1));
				for (auto& job : jobs)
					job.wait();

				Gripper::clog << Gripper::LogLevel::DEBUG << "Merging results";

				for (auto& job : jobs)
				{
					auto temp = job.get();
                    std::copy(temp.cbegin(), temp.cend(), std::back_inserter(m_data));
				}
                
                for (auto i = index_type{ 0, 0, -m_s_max }; i <= index_type{ m_l_max, m_l_max, m_s_max }; i = SpinWeightedSpherical::next(i, m_l_max, m_s_max))
                {
                    auto first = std::find_if(m_data.cbegin(), m_data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                    auto last = std::find_if_not(first, m_data.cend(), [&](const value_type& elem) { return elem.first.at(0) == i; });
                    m_accel.push_back(marker_type(first, last));
                }
                
				auto end = std::chrono::high_resolution_clock::now();

				Gripper::clog << Gripper::LogLevel::TIMING << "Computation took " << std::chrono::duration_cast<std::chrono::seconds>(end.time_since_epoch() - start.time_since_epoch()).count() << " seconds.";

				LEAVING
			}

			container_type computeLM(const IndexInternalType L1, const IndexInternalType M1) // Compute coefficients with given L1 value
			{
				ENTERING

				Gripper::clog << Gripper::LogLevel::ALL << "Compute thread L1,M1 = " << L1 << "," << M1 << " (id = " << std::this_thread::get_id() << ") started.";

				container_type result;

				// Compute all coefficients with designated L1,M1 value

				for (IndexInternalType S1 = -m_s_max; S1 <= m_s_max; ++S1)
					for (IndexInternalType L2 = 0; L2 <= m_l_max; ++L2)
						for (IndexInternalType M2 = -L2; M2 <= L2; ++M2)
							for (IndexInternalType S2 = -m_s_max; S2 <= m_s_max; ++S2)
								for (IndexInternalType L3 = 0; L3 <= m_l_max; ++L3)
									for (IndexInternalType M3 = -L3; M3 <= L3; ++M3)
										for (IndexInternalType S3 = -m_s_max; S3 <= m_s_max; ++S3)
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
                                                result.push_back(value_type(key_type{ index_type{ L1, M1, S1 },
                                                                                      index_type{ L2, M2, S2 },
                                                                                      index_type{ L3, M3, S3 } },
                                                                            res));
										}					

				Gripper::clog << Gripper::LogLevel::ALL << "Compute thread L = " << L1 << "," << M1 << " finished.";

				LEAVING
				
				return result;
			}
		};

	} // namespace SpinWeightedGaunt

} // namespace Multipole

namespace Gripper
{
    EXPORT extern Multipole::Gaunt::Matrix gaunt;
} // namespace Gripper


////////////////////////////////////////
// Radial::Index non-member operators //
////////////////////////////////////////

// Unary 
Multipole::Radial::Index operator+(const Multipole::Radial::Index& rhs);
Multipole::Radial::Index operator-(const Multipole::Radial::Index& rhs);

// Binary
Multipole::Radial::Index operator+(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
Multipole::Radial::Index operator-(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
Multipole::Radial::Index operator*(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
Multipole::Radial::Index operator/(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
Multipole::Radial::Index operator%(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);

bool operator< (const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
bool operator> (const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
bool operator<=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
bool operator>=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
bool operator==(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);
bool operator!=(const Multipole::Radial::Index& lhs, const Multipole::Radial::Index& rhs);

///////////////////////////////////////////
// Spherical::Index non-member operators //
///////////////////////////////////////////

// Unary 
Multipole::Spherical::Index operator+(const Multipole::Spherical::Index& rhs);
Multipole::Spherical::Index operator-(const Multipole::Spherical::Index& rhs);

// Binary
Multipole::Spherical::Index operator+(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
Multipole::Spherical::Index operator-(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
Multipole::Spherical::Index operator*(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
Multipole::Spherical::Index operator/(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
Multipole::Spherical::Index operator%(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);

bool operator< (const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
bool operator> (const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
bool operator<=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
bool operator>=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
bool operator==(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);
bool operator!=(const Multipole::Spherical::Index& lhs, const Multipole::Spherical::Index& rhs);

///////////////////////////////////////////////
// Spherical::IndexPair non-member operators //
///////////////////////////////////////////////

// Binary
bool operator==(const Multipole::Spherical::IndexPair& lhs, const Multipole::Spherical::IndexPair& rhs);
bool operator!=(const Multipole::Spherical::IndexPair& lhs, const Multipole::Spherical::IndexPair& rhs);

///////////////////////////////////////
// Gaunt::Index non-member operators //
///////////////////////////////////////

// Binary
Multipole::Gaunt::Index operator+(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
Multipole::Gaunt::Index operator-(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
Multipole::Gaunt::Index operator*(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
Multipole::Gaunt::Index operator/(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
Multipole::Gaunt::Index operator%(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);

bool operator< (const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
bool operator> (const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
bool operator<=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
bool operator>=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
bool operator==(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);
bool operator!=(const Multipole::Gaunt::Index& lhs, const Multipole::Gaunt::Index& rhs);

/////////////////////////////////////////////
// Gaunt::Coefficient non-member operators //
/////////////////////////////////////////////

bool operator< (const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);
bool operator> (const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);
bool operator<=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);
bool operator>=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);
bool operator==(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);
bool operator!=(const Multipole::Gaunt::Coefficient& lhs, const Multipole::Gaunt::Coefficient& rhs);

///////////////////////////////////////////////////////
// SpinWeightedSpherical::Index non-member operators //
///////////////////////////////////////////////////////

// Binary
template <typename AT> bool operator==(const Multipole::SpinWeightedSpherical::Index<AT>& lhs, const Multipole::SpinWeightedSpherical::Index<AT>& rhs) { return (lhs.l == rhs.l) && (lhs.m == rhs.m) && (lhs.s == rhs.s); }
template <typename AT> bool operator!=(const Multipole::SpinWeightedSpherical::Index<AT>& lhs, const Multipole::SpinWeightedSpherical::Index<AT>& rhs) { return (lhs.l != rhs.l) || (lhs.m != rhs.m) || (lhs.s != rhs.s); }

#endif // MULTIPOLETYPES_HPP