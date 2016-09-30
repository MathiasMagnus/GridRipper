#pragma once

// Standard C++ includes
#include <cstddef>          // std::ptrdiff_t
#include <iterator>         // std::iterator, std::random_access_iterator_tag

template <typename Integral>
class arithmetic_progression_iterator :
    public std::iterator<std::random_access_iterator_tag, Integral, std::ptrdiff_t, const Integral*, const Integral&>
{
public:

    arithmetic_progression_iterator() = default;
    arithmetic_progression_iterator(const arithmetic_progression_iterator&) = default;
    arithmetic_progression_iterator(arithmetic_progression_iterator&& x) = default;
    ~arithmetic_progression_iterator() = default;

    arithmetic_progression_iterator(value_type init, value_type d) : _value(init), _d(d) {}

    value_type diff() { return _d; }

    //
    // Iterator concept
    //
    const reference operator*() const { return _value; }

    arithmetic_progression_iterator& operator++()
    {
        _value += _d;

        return *this;
    }

    //
    // InputIterator concept
    //

    // These operators should be non-member, but partially specializing them is no good.
    inline bool operator==(const arithmetic_progression_iterator& rhs)
    {
        return _value == rhs._value &&
            _d == rhs._d;
    }

    inline bool operator!=(const arithmetic_progression_iterator& rhs)
    {
        return _value != rhs._value ||
            _d != rhs._d;
    }

    const pointer operator->() const { return &_value; }

    arithmetic_progression_iterator operator++(int)
    {
        arithmetic_progression_iterator tmp = *this;
        
        ++*this;

        return tmp;
    }

    //
    // BidirectionalIterator concept
    //
    arithmetic_progression_iterator& operator--()
    {
        _value -= _d;

        return *this;
    }

    arithmetic_progression_iterator operator--(int)
    {
        arithmetic_progression_iterator tmp = *this;

        --*this;

        return tmp;
    }

    //
    // RandomAccessIterator concept
    //
    arithmetic_progression_iterator& operator+=(int i)
    {
        _value += i * _d;

        return *this;
    }

    arithmetic_progression_iterator& operator-=(int i)
    {
        _value -= i * _d;

        return *this;
    }

    const reference operator[](Integral i)
    {
        return _value + i * _d;
    }

    inline bool operator<(const arithmetic_progression_iterator& rhs)
    {
        return _value < rhs._value;
    }

    inline bool operator>(const arithmetic_progression_iterator& rhs)
    {
        return _value > rhs._value;
    }

    inline bool operator<=(const arithmetic_progression_iterator& rhs)
    {
        return !(_value > rhs._value);
    }

    inline bool operator>=(const arithmetic_progression_iterator& rhs)
    {
        return !(_value < rhs._value);
    }
    
private:

    value_type _value, _d;
};

//
//     arithmetic_progression_iterator non-member operators
//
// RandomAccessIterator concept
//

template <typename Integral>
arithmetic_progression_iterator<Integral> operator+(const arithmetic_progression_iterator<Integral>& lhs,
                                                    const typename arithmetic_progression_iterator<Integral>::difference_type count)
{
    return arithmetic_progression_iterator<Integral>(*lhs + count * lhs.diff(), lhs.diff());
}

template <typename Integral>
arithmetic_progression_iterator<Integral> operator+(const typename arithmetic_progression_iterator<Integral>::difference_type count,
                                                    const arithmetic_progression_iterator<Integral>& rhs)
{
    return arithmetic_progression_iterator<Integral>(*rhs + count * rhs.diff(), rhs.diff());
}

template <typename Integral>
arithmetic_progression_iterator<Integral> operator-(const arithmetic_progression_iterator<Integral>& lhs,
    const typename arithmetic_progression_iterator<Integral>::difference_type count)
{
    return arithmetic_progression_iterator<Integral>(*lhs - count * lhs.diff(), lhs.diff());
}

template <typename Integral>
arithmetic_progression_iterator<Integral> operator-(const typename arithmetic_progression_iterator<Integral>::difference_type count,
    const arithmetic_progression_iterator<Integral>& rhs)
{
    return arithmetic_progression_iterator<Integral>(*rhs - count * rhs.diff(), rhs.diff());
}
