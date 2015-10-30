//////////////////////////////////////////////////////////////////////
//                                                                  //
// Unexported templated PDE solver to be used by client codebase    //
//                                                                  //
// Author: Nagy-Egri Máté Ferenc                                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#ifndef PDE_HPP
#define PDE_HPP


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

        template <int N> auto get() const { return static_cast<const T&>(*this).get<N>(); }

        // Expression interface

        operator T&()             { return static_cast<T&>(*this); }
        operator T const&() const { return static_cast<const T&>(*this); }
    };


    template <typename E1, typename E2>
    class Sum : public Expression<Sum<E1, E2>>
    {
    public:

        Sum(Expression<E1> const& u, Expression<E2> const& v) : _u(u), _v(v) { /*ASSERT MISSING*/ }

        // StateVector interface

        //template <int N, typename R = std::result_of<operator+(E1::get<N>(), E2::get<N>())>::type> const R& get() { return _u.get<N>() + _v.get<N>(); }
        template <int N> auto get() const { return _u.get<N>() + _v.get<N>(); }

    private:

        E1 const& _u;
        E2 const& _v;
    };


    template <typename E>
    class Scale : public Expression<Scale<E>>
    {
    public:

        Scale(float alpha, Expression<E> const& v) : _alpha(alpha), _v(v) { /*ASSERT MISSING*/ }

        // StateVector interface

        //template <int N, typename R = std::result_of<operator*(float, E::get<N>())>::type> const R& get() { return _alpha * _v.get<N>(); }
        template <int N> auto get() const { return _alpha * _v.get<N>(); }

    private:

        float _alpha;
        E const& _v;
    };


    template <typename E> Scale<E> const operator*(const float alpha, Expression<E> const& v) { return Scale<E>(alpha, v); }
    template <typename E1, typename E2> Sum<E1, E2> const operator+(Expression<E1> const& u, Expression<E2> const& v) { return Sum<E1, E2>(u, v); }


    template <typename... T>
    class StateVector : public Expression<StateVector<T...>>
    {
    private:

        std::tuple<T...> m_values;

        template <int N>
        struct Assign
        {
            template <typename ET>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                std::get<N>(dst) = src.get<N>();
                Assign<N - 1>::assign(dst, src);
            }
        };


        template <>
        struct Assign<0>
        {
            template <typename ET>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                std::get<0>(dst) = src.get<0>();
            }
        };

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

        //template <int N, typename R = std::tuple_element<N, std::tuple<T...>>::type> const R& get() { return std::get<N>(m_values); }
        //template <int N, typename R = std::tuple_element<N, std::tuple<T...>>::type> R& get() { return std::get<N>(m_values); }

        template <int N> auto get() const { return std::get<N>(m_values); }
        template <int N> auto get() { return std::get<N>(m_values); }

        // Expression interface

        template <typename StateVecExpr>
        StateVector(Expression<StateVecExpr> const& expr)
        {
            // Extract type from encapsulating expression
            StateVecExpr const& v = expr;

            Assign<(sizeof...(T)) - 1>::assign(m_values, v);
        }
    };


    namespace RK4
    {
        template <typename... States>
        class Solver
        {
        public:

            Solver() {}
            Solver(const Solver& in) = delete;
            Solver(Solver&& src) = delete;
            ~Solver() {}

            void setLHS(StateVector<States...>&& lhs_in) { m_lhs = lhs_in; }
            StateVector<States...>& getLHS() { return m_lhs; }
            
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

#endif // PDE_HPP