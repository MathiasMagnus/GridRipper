//////////////////////////////////////////////////////////////////////
//                                                                  //
// Unexported templated PDE solver to be used by client codebase    //
//                                                                  //
// Author: Nagy-Egri Máté Ferenc                                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#pragma once


// Standard C++ includes
#include <cassert>
#include <tuple>
#include <array>
#include <functional>

namespace PDE
{
    template <typename E>
    class Expression
    {
    public:

        // StateVector interface

        template <int N> auto get() const { return static_cast<const E&>(*this).template get<N>(); }

        // Expression interface

        operator E&()             { return static_cast<E&>(*this); }
        operator E const&() const { return static_cast<const E&>(*this); }
    };


    template <typename E1, typename E2>
    class Sum : public Expression<Sum<E1, E2>>
    {
    public:

        Sum(const Expression<E1>& u, const Expression<E2>& v) : _u(u), _v(v) { /*ASSERT MISSING*/ }

        // StateVector interface

        template <int N> auto get() const { return _u.template get<N>() + _v.template get<N>(); }

    private:

        const E1& _u;
        const E2& _v;
    };


    template <typename T, typename E>
    class Scale : public Expression<Scale<T, E>>
    {
    public:

        Scale(const T& alpha, const Expression<E>& v) : _alpha(alpha), _v(v) {}

        // StateVector interface

        template <int N> auto get() const { return _alpha * _v.template get<N>(); }

    private:

        T _alpha;
        const E& _v;
    };

    namespace impl
    {
        template <int N>
        struct Assign
        {
            template <typename ET, typename... T>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                std::get<N>(dst) = src.template get<N>();
                Assign<N - 1>::assign(dst, src);
            }
        };


        template <>
        struct Assign<0>
        {
            template <typename ET, typename... T>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                std::get<0>(dst) = src.template get<0>();
            }
        };
    }


    template <typename... T>
    class StateVector : public Expression<StateVector<T...>>
    {
    private:

        std::tuple<T...> m_values;

    public:

        // Common interface

        StateVector() = default;
        StateVector(const StateVector&) = default;
        StateVector(StateVector&&) = default;
        ~StateVector() = default;

        StateVector& operator=(const StateVector&) = default;
        StateVector& operator=(StateVector&&) = default;

        // StateVector interface

        StateVector(const T&... values) : m_values(std::make_tuple(values...)) {}
        StateVector(T&&... values) : m_values(std::make_tuple(values...)) {}

        template <int N> const auto& get() const { return std::get<N>(m_values); }
        template <int N> auto& get() { return std::get<N>(m_values); }

        // Expression interface

        template <typename StateVecExpr>
        StateVector(const Expression<StateVecExpr>& expr)
        {
            // Extract type from encapsulating expression
            const StateVecExpr& v = expr;

            impl::Assign<(sizeof...(T)) - 1>::assign(m_values, v);
        }
    };


    namespace RK4
    {
        template <typename T, typename States>
        class Solver;

        template <typename T, typename... States>
        class Solver<T, StateVector<States...>>
        {
        public:

            Solver() = default;
            Solver(const Solver&) = default;
            Solver(Solver&& src) = default;
            ~Solver() = default;

            void setLHS(StateVector<States...>&& lhs_in) { m_lhs = lhs_in; }
            const StateVector<States...>& getLHS() { return m_lhs; }
            
            void setEquationSystem(std::function<StateVector<States...>(const StateVector<States...>&)> evaluate) { m_func = evaluate; }
            std::function<StateVector<States...>(const StateVector<States...>&)> getEquationSystem() { return m_func; }

            T iterate(const T& dt)
            {
                m_k.at(k1) = m_func(m_lhs);

                tmp1 = static_cast<T>(0.5) * dt * m_k.at(k1);
                tmp2 = m_lhs + tmp1;
                m_k.at(k2) = m_func(tmp2);
                
                tmp1 = static_cast<T>(0.5) * dt * m_k.at(k2);
                tmp2 = m_lhs + tmp1;
                m_k.at(k3) = m_func(tmp2);

                tmp1 = static_cast<T>(1.0) * dt * m_k.at(k3);
                tmp2 = m_lhs + tmp1;
                m_k.at(k4) = m_func(tmp2);

                tmp1 = static_cast<T>(2.0) * m_k.at(k2);
                tmp2 = m_k.at(k1) + tmp1;
                tmp1 = static_cast<T>(2.0) * m_k.at(k3);
                tmp3 = tmp2 + tmp1;
                tmp1 = tmp3 + m_k.at(k4);
                tmp2 = (dt / static_cast<T>(6.0)) * tmp1;
                tmp3 = m_lhs + tmp2;

                m_lhs = tmp3;

                //m_lhs = m_lhs + (dt / static_cast<T>(6.0))*(m_k.at(k1) + static_cast<T>(2.0) * m_k.at(k2) + static_cast<T>(2.0) * m_k.at(k3) + m_k.at(k4));

                return dt;
            }

        private:

            enum Step
            {
                k1 = 0,
                k2 = 1,
                k3 = 2,
                k4 = 3
            };

            StateVector<States...> m_lhs, tmp1, tmp2, tmp3;
            std::array<StateVector<States...>, 4> m_k;
            std::function<StateVector<States...>(const StateVector<States...>&)> m_func;
        };

    } // namespace RK4
    
} // namespace PDE


///////////////////////////////////////////
// PDE::StateVector non-member operators //
///////////////////////////////////////////

template <typename T, typename E> auto operator*(const T& alpha, const PDE::Expression<E>& v) { return PDE::Scale<T, E>(alpha, v); }
template <typename E1, typename E2> auto operator+(const PDE::Expression<E1>& u, const PDE::Expression<E2>& v) { return PDE::Sum<E1, E2>(u, v); }