#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>

// Standard C++ includes
#include <cstddef>          // std::ptrdiff_t
#include <iterator>         // std::iterator, std::random_access_iterator_tag

/*
/// <summary>Iterator for indexable-only objects.</summary>
///
template <typename T, typename IT, typename VT>
class const_index_iterator :
    public std::iterator<std::random_access_iterator_tag, VT, std::ptrdiff_t, const VT*, const VT&>
{
public:

    using indexable_type = T;
    using index_type = IT;
    using value_type = VT;

    // Constructors / Destructors / Assignment operators

    /// <summary>Default constructor.</summary>
    /// <remarks>Default constructed objects are in an invalid state.</remarks>
    ///
    const_index_iterator() = default;

    /// <summary>Default copy constructor.</summary>
    ///
    const_index_iterator(const const_index_iterator& in) = default;

    /// <summary>Default move constructor.</summary>
    ///
    const_index_iterator(const_index_iterator&& in) = default;

    /// <summary>Default destructor.</summary>
    ///
    ~const_index_iterator() = default;

    /// <summary>Default copy assignment operator.</summary>
    ///
    const_index_iterator& operator=(const const_index_iterator&) = default;

    /// <summary>Default move assignment operator.</summary>
    ///
    const_index_iterator& operator=(const_index_iterator&&) = default;

    /// <summary>Construct from indexable object and initial index.</summary>
    ///
    const_index_iterator(const indexable_type& indexable, index_type init) : _cont(&indexable), _i(init) {}

    //
    // Iterator concept
    //
    const reference operator*() const { return _value; }

    const_index_iterator& operator++()
    {
        _value += _d;

        return *this;
    }

    //
    // InputIterator concept
    //

    // These operators should be non-member, but partially specializing them is no good.
    inline bool operator==(const const_index_iterator& rhs)
    {
        return _value == rhs._value &&
            _d == rhs._d;
    }

    inline bool operator!=(const const_index_iterator& rhs)
    {
        return _value != rhs._value ||
            _d != rhs._d;
    }

    const pointer operator->() const { return &_value; }

    const_index_iterator operator++(int)
    {
        arithmetic_progression_iterator tmp = *this;

        ++*this;

        return tmp;
    }

    //
    // BidirectionalIterator concept
    //
    const_index_iterator& operator--()
    {
        _value -= _d;

        return *this;
    }

    const_index_iterator operator--(int)
    {
        arithmetic_progression_iterator tmp = *this;

        --*this;

        return tmp;
    }

    //
    // RandomAccessIterator concept
    //
    const_index_iterator& operator+=(int i)
    {
        _value += i * _d;

        return *this;
    }

    const_index_iterator& operator-=(int i)
    {
        _value -= i * _d;

        return *this;
    }

    const reference operator[](Integral i)
    {
        return _value + i * _d;
    }

    inline bool operator<(const const_index_iterator& rhs)
    {
        return _value < rhs._value;
    }

    inline bool operator>(const const_index_iterator& rhs)
    {
        return _value > rhs._value;
    }

    inline bool operator<=(const const_index_iterator& rhs)
    {
        return !(_value > rhs._value);
    }

    inline bool operator>=(const const_index_iterator& rhs)
    {
        return !(_value < rhs._value);
    }


private:

    indexable_type const* _cont;
    index_type _i;
};*/