//////////////////////////////////////////////////////////////////////
//                                                                  //
// Unexported templated PDE solver to be used by client codebase    //
//                                                                  //
// Author: Nagy-Egri M�t� Ferenc                                    //
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
    template <typename T>
    class Expression
    {
    public:

        // StateVector interface

        template <int N> auto get() const { return static_cast<const T&>(*this).template get<N>(); }

        // Expression interface

        operator T&()             { return static_cast<T&>(*this); }
        operator T const&() const { return static_cast<const T&>(*this); }
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


    template <typename E>
    class Scale : public Expression<Scale<E>>
    {
    public:

        Scale(const float alpha, const Expression<E>& v) : _alpha(alpha), _v(v) {}

        // StateVector interface

        template <int N> auto get() const { return _alpha * _v.template get<N>(); }

    private:

        float _alpha;
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

        StateVector() {}
        StateVector(const StateVector& in) : m_values(in.m_values) {}
        StateVector(StateVector&& src) : m_values(std::move(src.m_values)) {}
        ~StateVector() {}

        StateVector& operator=(const StateVector& rhs) { m_values = rhs.m_values; return *this; }

        // StateVector interface

        StateVector(const T&... values) : m_values(std::make_tuple(values...)) {}
        StateVector(T&&... values) : m_values(std::make_tuple(values...)) {}

        template <int N> auto& get() const { return std::get<N>(m_values); }
        template <int N> auto& get() { return std::get<N>(m_values); }

        // Expression interface

        template <typename StateVecExpr>
        StateVector(Expression<StateVecExpr> const& expr)
        {
            // Extract type from encapsulating expression
            StateVecExpr const& v = expr;

            impl::Assign<(sizeof...(T)) - 1>::assign(m_values, v);
        }
    };


    namespace RK4
    {
        template <typename... States>
        class Solver
        {
        public:

            Solver() = default;
            Solver(const Solver& in) = delete;
            Solver(Solver&& src) = default;
            ~Solver() = default;

            void setLHS(StateVector<States...>&& lhs_in) { m_lhs = lhs_in; }
            const StateVector<States...>& getLHS() { return m_lhs; }
            
            void setEquationSystem(std::function<StateVector<States...>(const StateVector<States...>&)> evaluate) { m_func = evaluate; }
            std::function<StateVector<States...>(const StateVector<States...>&)> getEquationSystem() { return m_func; }

            float iterate(float dt)
            {
                m_k.at(k1) = m_func(m_lhs);

                m_k.at(k2) = m_func(m_lhs + 0.5f * dt * m_k.at(k1));
                
                m_k.at(k3) = m_func(m_lhs + 0.5f * dt * m_k.at(k2));

                m_k.at(k4) = m_func(m_lhs + 1.0f * dt * m_k.at(k3));

                m_lhs = m_lhs + (dt / 6.0f)*(m_k.at(k1) + 2.0f * m_k.at(k2) + 2.0f * m_k.at(k3) + m_k.at(k4));

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

            StateVector<States...> m_lhs;
            std::array<StateVector<States...>, 4> m_k;
            std::function<StateVector<States...>(const StateVector<States...>&)> m_func;
        };

    } // namespace RK4
    
} // namespace PDE


///////////////////////////////////////////
// PDE::StateVector non-member operators //
///////////////////////////////////////////

template <typename E> auto operator*(const float alpha, const PDE::Expression<E>& v) { return PDE::Scale<E>(alpha, v); }
template <typename E1, typename E2> auto operator+(const PDE::Expression<E1>& u, const PDE::Expression<E2>& v) { return PDE::Sum<E1, E2>(u, v); }