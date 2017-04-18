#pragma once

// Gripper includes
#include <Gripper/stl/stlConfig.hpp>

// Standard C++ includes
#include <cstddef>          // std::ptrdiff_t
#include <iterator>         // std::iterator, std::random_access_iterator_tag


namespace stl
{
    template <typename Integral>
    class arithmetic_progression_iterator : public std::iterator<std::random_access_iterator_tag, Integral, std::ptrdiff_t, const Integral*, const Integral&>
    {
    private:

        using iter_base = std::iterator<std::random_access_iterator_tag, Integral, std::ptrdiff_t, const Integral*, const Integral&>;

    public:

        // Iterator aliases

        using typename iter_base::iterator_category;
        using typename iter_base::value_type;
        using typename iter_base::difference_type;
        using typename iter_base::pointer;
        using typename iter_base::reference;

        arithmetic_progression_iterator() = default;
        arithmetic_progression_iterator(const arithmetic_progression_iterator&) = default;
        arithmetic_progression_iterator(arithmetic_progression_iterator&&) = default;
        ~arithmetic_progression_iterator() = default;

        arithmetic_progression_iterator& operator=(const arithmetic_progression_iterator&) = default;
        arithmetic_progression_iterator& operator=(arithmetic_progression_iterator&&) = default;

        arithmetic_progression_iterator(value_type init, value_type d = 1) : _value(init), _d(d) {}

        value_type diff() const { return _d; }

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
    // arithmetic_progression_iterator non-member operators
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

    template <typename Integral>
    typename arithmetic_progression_iterator<Integral>::difference_type operator-(const arithmetic_progression_iterator<Integral>& lhs,
        const arithmetic_progression_iterator<Integral>& rhs)
    {
        if (lhs.diff() != rhs.diff())
            throw std::runtime_error{ "stl::arithmetic_progression_iterator misuse: strides do not match, distance does not exist." };

        if (*rhs - *lhs % lhs.diff() != 0)
            throw std::runtime_error{ "stl::arithmetic_progression_iterator misuse: begin, end values with designated stride does not form a range; distance does not exist." };

        return (*rhs - *lhs) / lhs.diff();
    }

} // namespace stl
