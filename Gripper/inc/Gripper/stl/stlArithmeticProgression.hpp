#pragma once

// Standard C++ includes
#include <cstddef>          // std::ptrdiff_t
#include <iterator>         // std::iterator, std::random_access_iterator_tag


namespace stl
{
    //
    // Forward declarations
    //

    template <typename Integral>
    class arithmetic_progression_iterator;

    template <typename Integral>
    arithmetic_progression_iterator<Integral> operator+(const arithmetic_progression_iterator<Integral>& lhs,
                                                        const typename arithmetic_progression_iterator<Integral>::difference_type count);
    template <typename Integral>
    arithmetic_progression_iterator<Integral> operator+(const typename arithmetic_progression_iterator<Integral>::difference_type count,
                                                        const arithmetic_progression_iterator<Integral>& rhs);
    template <typename Integral>
    arithmetic_progression_iterator<Integral> operator-(const arithmetic_progression_iterator<Integral>& lhs,
                                                        const typename arithmetic_progression_iterator<Integral>::difference_type count);
    template <typename Integral>
    arithmetic_progression_iterator<Integral> operator-(const typename arithmetic_progression_iterator<Integral>::difference_type count,
                                                        const arithmetic_progression_iterator<Integral>& rhs);
    template <typename Integral>
    typename arithmetic_progression_iterator<Integral>::difference_type operator-(const arithmetic_progression_iterator<Integral>& lhs,
                                                        const arithmetic_progression_iterator<Integral>& rhs);

    template <typename Integral>
    class arithmetic_progression_iterator : public std::iterator<std::random_access_iterator_tag, Integral, std::ptrdiff_t, const Integral*, const Integral&>
    {
    public:

        // arithmetic_progression_iterator aliases

        using iter_base = std::iterator<std::random_access_iterator_tag, Integral, std::ptrdiff_t, const Integral*, const Integral&>;

        // Iterator aliases

        using typename iter_base::iterator_category;
        using typename iter_base::value_type;
        using typename iter_base::difference_type;
        using typename iter_base::pointer;
        using typename iter_base::reference;

        // friend declarations

        friend difference_type operator-<>(const arithmetic_progression_iterator&, const arithmetic_progression_iterator&);
        friend arithmetic_progression_iterator operator+<>(const arithmetic_progression_iterator&, const difference_type);
        friend arithmetic_progression_iterator operator+<>(const difference_type, const arithmetic_progression_iterator&);
        friend arithmetic_progression_iterator operator-<>(const arithmetic_progression_iterator&, const difference_type);
        friend arithmetic_progression_iterator operator-<>(const difference_type, const arithmetic_progression_iterator&);

        // Constructors, destructors, assignment operators

        arithmetic_progression_iterator() : m_value(), m_final(), m_diff(0) {}
        arithmetic_progression_iterator(const arithmetic_progression_iterator&) = default;
        arithmetic_progression_iterator(arithmetic_progression_iterator&&) = default;
        ~arithmetic_progression_iterator() = default;

        arithmetic_progression_iterator& operator=(const arithmetic_progression_iterator&) = default;
        arithmetic_progression_iterator& operator=(arithmetic_progression_iterator&&) = default;

        arithmetic_progression_iterator(value_type init, value_type fin, value_type diff = 1) : m_value(init), m_final(fin), m_diff(diff)
        {
            if ((fin - init) % diff != 0)
                throw std::domain_error{ "stl::arithmetic_progression_iterator misuse: initial, final and stride values cannot form a range." };
        }

        // arithmetic_progression_iterator member functions

        value_type diff() const { return m_diff; }

        //
        // Iterator concept
        //
        const reference operator*() const { return m_value; }

        arithmetic_progression_iterator& operator++()
        {
            m_value += m_diff;

            return *this;
        }

        //
        // InputIterator concept
        //

        // These operators should be non-member, but partially specializing them is no good.
        inline bool operator==(const arithmetic_progression_iterator& rhs)
        {
            if (rhs.is_end())
                return m_value == m_final;
			if (this->is_end())
				return rhs.m_value == rhs.m_final;
            else
                return m_value == rhs.m_value &&
                m_final == rhs.m_final &&
                m_diff == rhs.m_diff;
        }

        inline bool operator!=(const arithmetic_progression_iterator& rhs)
        {
			return !(*this == rhs);
        }

        pointer operator->() const { return &m_value; }

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
            m_value -= m_diff;

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
            m_value += i * m_diff;

            return *this;
        }

        arithmetic_progression_iterator& operator-=(int i)
        {
            m_value -= i * m_diff;

            return *this;
        }

        const reference operator[](Integral i)
        {
            return m_value + i * m_diff;
        }

        inline bool operator<(const arithmetic_progression_iterator& rhs)
        {
            if (rhs.is_end()) return true;
			if (this->is_end()) return false;

            return m_value < rhs.m_value;
        }

        inline bool operator>(const arithmetic_progression_iterator& rhs)
        {
			if (rhs.is_end()) return false;
			if (this->is_end()) return true;

            return m_value > rhs.m_value;
        }

        inline bool operator<=(const arithmetic_progression_iterator& rhs)
        {
            return !(m_value > rhs.m_value);
        }

        inline bool operator>=(const arithmetic_progression_iterator& rhs)
        {
            return !(m_value < rhs.m_value);
        }

    private:

        bool is_end() const { return !static_cast<bool>(m_diff); }

        value_type m_value, m_final, m_diff;
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
        if (rhs.is_end())
        {
            return (*lhs - lhs.m_final) / lhs.diff();
        }
        else
        {
            if (lhs.diff() != rhs.diff())
                throw std::runtime_error{ "stl::arithmetic_progression_iterator misuse: strides do not match, distance does not exist." };

            if (lhs.m_final != rhs.m_final)
                throw std::runtime_error{ "stl::arithmetic_progression_iterator misuse: end values do not match, distance does not exist." };

            if ((*rhs - *lhs) % lhs.diff() != 0)
                throw std::runtime_error{ "stl::arithmetic_progression_iterator misuse: begin, end values with designated stride does not form a range; distance does not exist." };

            return (*rhs - *lhs) / lhs.diff();
        }
    }

} // namespace stl
