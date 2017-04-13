#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>
#include <Gripper/stl/stlSpinWeightedSphericalIndex.hpp>    // math::sws::index

// Standard C++ includes
#include <cassert>              // assert
#include <cstddef>              // std::size_t
#include <initializer_list>     // std::initializer_list
#include <ostream>              // std::ostream
#include <iterator>             // std::iterator, std::random_access_iterator_tag


namespace math
{
    namespace sws
    {
        /// <summary>Class template representing the extent of spherical expansions.</summary>
        ///
        template <typename ArithemticType = typename index<>::value_type>
        class extent
        {
        public:

            // ConstIterator forward-declarations

            class const_bidir_iterator;
            class reverse_bidir_const_iterator;

            // Lattice type aliases

            typedef index<ArithemticType> index_type;
            typedef typename index_type::value_type index_internal_type;

            // ConstIterator aliases

            using const_iterator_type = const_bidir_iterator;
            using const_reverse_iterator_type = reverse_bidir_const_iterator;

            // ConstSTL member class

        protected:

            using iter_base = std::iterator<std::bidirectional_iterator_tag, index_type, std::ptrdiff_t, const index_type*, const index_type&>;

        public:

            /// <summary>Utility class to iterate over the index values within an extent.</summary>
            ///
            class const_bidir_iterator : public iter_base
            {
            public:

                // Iterator aliases

                using typename iter_base::iterator_category;
                using typename iter_base::value_type;
                using typename iter_base::difference_type;
                using typename iter_base::pointer;
                using typename iter_base::reference;

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                const_bidir_iterator() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                const_bidir_iterator(const const_bidir_iterator&) = default;

                /// <summary>Default move constructor.</summary>
                ///
                const_bidir_iterator(const_bidir_iterator&&) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~const_bidir_iterator() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                const_bidir_iterator& operator=(const const_bidir_iterator&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                const_bidir_iterator& operator=(const_bidir_iterator&&) = default;

                /// <summary>Index constructor.</summary>
                ///
                const_bidir_iterator(const value_type& value) : _value(value) {}

                // Iterator concept

                /// <summary>Dereference operator.</summary>
                ///
                const reference operator*() const { return _value; }

                /// <summary>Prefix increment operator.</summary>
                ///
                const_bidir_iterator& operator++()
                {
                    bool rewind_m_and_step_l = _value.m == _value.l;

                    _value.l = rewind_m_and_step_l ? _value.l + 1 : _value.l;
                    _value.m = rewind_m_and_step_l ? -_value.l : _value.m + 1;

                    return *this;
                }

                // InputIterator concept

                /// <summary>Equality operator.</summary>
                /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                ///
                inline bool operator==(const const_bidir_iterator& rhs)
                {
                    return _value == rhs._value;
                }

                /// <summary>Unequality operator.</summary>
                /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                ///
                inline bool operator!=(const const_bidir_iterator& rhs)
                {
                    return _value != rhs._value;
                }

                /// <summary>Arrow operator.</summary>
                ///
                const pointer operator->() const { return &_value; }

                /// <summary>Postfix increment operator.</summary>
                ///
                const_bidir_iterator operator++(int)
                {
                    const_bidir_iterator tmp = *this;

                    ++*this;

                    return tmp;
                }

                // BidirectionalIterator concept

                /// <summary>Prefix decrement operator.</summary>
                ///
                const_bidir_iterator& operator--()
                {
                    bool rewind_m_and_step_l = _value.m == -_value.l;

                    _value.l = rewind_m_and_step_l ? _value.l - 1 : _value.l;
                    _value.m = rewind_m_and_step_l ? _value.l : _value.m - 1;

                    return *this;
                }

                /// <summary>Postfix decrement operator.</summary>
                ///
                const_bidir_iterator operator--(int)
                {
                    const_bidir_iterator tmp = *this;

                    --*this;

                    return tmp;
                }

            private:

                value_type _value;
            };

            /// <summary>Utility class to iterate over the index values within an extent in reverse order.</summary>
            ///
            class reverse_bidir_const_iterator : public std::iterator<std::bidirectional_iterator_tag,
                                                          index_type,
                                                          std::ptrdiff_t,
                                                          const index_type*,
                                                          const index_type&>
            {
            public:

                // Iterator aliases

                using typename iter_base::iterator_category;
                using typename iter_base::value_type;
                using typename iter_base::difference_type;
                using typename iter_base::pointer;
                using typename iter_base::reference;

                // Constructors / Destructors / Assignment operators

                /// <summary>Default constructor.</summary>
                /// <remarks>Default constructed objects are in an invalid state.</remarks>
                ///
                reverse_bidir_const_iterator() = default;

                /// <summary>Default copy constructor.</summary>
                ///
                reverse_bidir_const_iterator(const reverse_bidir_const_iterator&) = default;

                /// <summary>Default move constructor.</summary>
                ///
                reverse_bidir_const_iterator(reverse_bidir_const_iterator&&) = default;

                /// <summary>Default destructor.</summary>
                ///
                ~reverse_bidir_const_iterator() = default;

                /// <summary>Default copy assignment operator.</summary>
                ///
                reverse_bidir_const_iterator& operator=(const reverse_bidir_const_iterator&) = default;

                /// <summary>Default move assignment operator.</summary>
                ///
                reverse_bidir_const_iterator& operator=(reverse_bidir_const_iterator&&) = default;

                /// <summary>Index constructor.</summary>
                ///
                reverse_bidir_const_iterator(const value_type& value) : _value(value) {}

                // Iterator concept

                /// <summary>Dereference operator.</summary>
                ///
                const reference operator*() const { return _value; }

                /// <summary>Prefix increment operator.</summary>
                ///
                reverse_bidir_const_iterator& operator++()
                {
                    bool rewind_m_and_step_l = _value.m == -_value.l;

                    _value.l = rewind_m_and_step_l ? _value.l - 1 : _value.l;
                    _value.m = rewind_m_and_step_l ? _value.l : _value.m - 1;

                    return *this;
                }

                // InputIterator concept

                /// <summary>Equality operator.</summary>
                /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                ///
                inline bool operator==(const reverse_bidir_const_iterator& rhs)
                {
                    return _value == rhs._value;
                }

                /// <summary>Unequality operator.</summary>
                /// <remarks>This operator should be non-member, but partially specializing it is no good.</remarks>
                ///
                inline bool operator!=(const reverse_bidir_const_iterator& rhs)
                {
                    return _value != rhs._value;
                }

                /// <summary>Arrow operator.</summary>
                ///
                const pointer operator->() const { return &_value; }

                /// <summary>Postfix increment operator.</summary>
                ///
                reverse_bidir_const_iterator operator++(int)
                {
                    reverse_bidir_const_iterator tmp = *this;

                    ++*this;

                    return tmp;
                }

                // BidirectionalIterator concept

                /// <summary>Prefix decrement operator.</summary>
                ///
                reverse_bidir_const_iterator& operator--()
                {
                    bool rewind_m_and_step_l = _value.m == _value.l;

                    _value.l = rewind_m_and_step_l ? _value.l + 1 : _value.l;
                    _value.m = rewind_m_and_step_l ? -_value.l : _value.m + 1;

                    return *this;
                }

                /// <summary>Postfix decrement operator.</summary>
                ///
                reverse_bidir_const_iterator operator--(int)
                {
                    reverse_bidir_const_iterator tmp = *this;

                    --*this;

                    return tmp;
                }

            private:

                value_type _value;
            };

            // Constructors / Destructors / Assignment operators

            /// <summary>Default constructor.</summary>
            /// <remarks>Default constructed objects are in an invalid state.</remarks>
            ///
            extent() = default;

            /// <summary>Default copy constructor.</summary>
            ///
            extent(const extent&) = default;

            /// <summary>Default move constructor.</summary>
            ///
            extent(extent&&) = default;

            /// <summary>Default destructor.</summary>
            ///
            ~extent() = default;

            /// <summary>Default copy assignment operator.</summary>
            ///
            extent& operator=(const extent&) = default;

            /// <summary>Default move assignment operator.</summary>
            ///
            extent& operator=(extent&&) = default;

            // ConstSTL interface

            /// <summary>Iterator to the first index inside the extent.</summary>
            ///
            const_iterator_type cbegin() { return const_iterator_type(_initial); }

            /// <summary>Iterator to one past the last index inside the extent.</summary>
            ///
            const_iterator_type cend() { return ++const_iterator_type(_final); }

            /// <summary>Iterator to the first index inside the extent.</summary>
            ///
            const_reverse_iterator_type crbegin() { return const_reverse_iterator_type(_final); }

            /// <summary>Iterator to one past the last index inside the extent.</summary>
            ///
            const_reverse_iterator_type crend() { return ++const_reverse_iterator_type(_initial); }

            // Lattice interface

            /// <summary>Constructor that validates contents.</summary>
            ///
            extent(const index_type& initial, const index_type& final) : _initial(initial), _final(final) { /*assert(_initial <= _final);*/ }

            /// <summary>Initializer list constructor that validates contents.</summary>
            ///
            extent(std::initializer_list<index_type> init) : _initial(*(init.begin())), _final(*(init.begin() + 1))
            {
                //static_assert(init.size() == 2, "Size of std::initializer_list<Index> to Multipole::stl::SWS::Extent must be 2");
            
                //assert(_initial <= _final);
            }

            /// <summary>Returns the origin of the extent.</summary>
            ///
            const index_type& initial() const { return _initial; }

            /// <summary>Returns the end of the extent.</summary>
            ///
            const index_type& final() const { return _final; }
            /*
            /// <summary>Tests whether an index is inside the extent.</summary>
            ///
            bool contains(const index_type& index) const { return (index >= _initial) && (index <= _final) ? true : false; }
            */
        private:

            index_type _initial;
            index_type _final;
        };

        /// <summary>Equality operator.</summary>
        ///
        template <typename AT> 
        bool operator==(const extent<AT>& lhs,
                        const extent<AT>& rhs)
        {
            return (lhs.final() == rhs.final()) && (lhs.initial() == rhs.initial());
        }

        /// <summary>Unequality operator.</summary>
        ///
        template <typename AT>
        bool operator!=(const extent<AT>& lhs,
                        const extent<AT>& rhs)
        {
            return (lhs.final() != rhs.final()) || (lhs.initial() != rhs.initial());
        }

    } // namespace sws

} // namespace math

namespace std
{
    /// <summary>STL stream operator overload for formatted console output.</summary>
    ///
    template <typename AT>
    ostream& operator<<(ostream& os, const math::sws::extent<AT>& ext)
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

        os << "[ " << ext.initial() << ", " << ext.final() << " ]";

        return os;
    }

} // namespace std
