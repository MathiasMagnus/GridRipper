//
// NOTE: This example code shows integration of the initial data
//       of a spinning black hole, without the use of elliptic
//       equations. Mathematical formulae used throughout the code
//       originate from the following articles, manusciprts:
//
// BIBLIOGRAPHY: [1] Black hole initial data without elliptic equations
//                   I. R�cz and J. Winicour
//
//               [2] Black hole initial data by integrating a strongly hyperbolic system
//                   M. Babiuc-Hamilton, A. Balogh, A. L�szl�, I. R�cz and J. Winicour
//
//               [3] Syllabus
//                   I. R�cz
//
// AUTHOR: Nagy-Egri M�t� Ferenc
//

#include <STL_Example5_Schwarzschild.hpp>

// Simulation aliases
using integral = std::int16_t;
using real = double;
using complex = std::complex<real>;
using spherical_index_internal_type = std::int16_t;
using spherical_index = math::sws::index<spherical_index_internal_type>;
using spherical_extent = math::sws::extent<spherical_index_internal_type>;
template <std::size_t S, math::sws::parity P> using real_coeff_vector = math::sws::vector<S, P, spherical_index_internal_type, real>;
template <std::size_t S, math::sws::parity P> using complex_coeff_vector = math::sws::vector<S, P, spherical_index_internal_type, complex>;

// Global constants
constexpr real pi = math::constants::pi<real>;

int main()
{
    // Simulation params
    //using integral = std::int32_t;      // Type used to represent integral values
    //using real = double;                // Type used to represent real values
    //using complex = std::complex<real>; // Type used to represent complex values
    constexpr integral L_max = 5;       // Maximum multipole values for series expansion
    //constexpr integral S_max = 3;       // Maximum spin values for series expansion
    //constexpr real drho = 0.025f;       // Separation of rho coordinates
    //constexpr real rho_min = 3;         // Innermost radial lattice coordinate
    //const integral rho_max = 200;        // Coordinate # of outermost lattice site
    //const integral rho_n = 120;          // Coordinate # to start integration from
    constexpr real M = 100.f;               // Black hole mass param
    //const real neumann_length = 0.001f; // Neumann-series expansion precision

    constexpr integral max_theta_evals = 50; // Evaluation point count for chebysev integral
    constexpr integral max_phi_evals = 50;   // Evaluation point count for periodic integral
    constexpr real convergence = 1e-6;       // Threshold for adaptive convergence test

    // Constant definitions
    constexpr real one = 1, two = 2, three = 3, four = 4, eight = 8, one_per_two = one / two, three_per_two = three / two;

    // Helper definitions
    auto H = [=](const real r) { return M / r; };

    // Analytic solution (theta, phi independant, totally spherically symmetric) // based on (3.36)
    auto analytic_solution = [=](const real rho,
                                 const real, // theta
                                 const real) // phi
    {
        return -four * M / (rho * rho * std::sqrt(one + two * H(rho)));
    };

    // Binds the rho parameter of a spherical function or series expansion to a given value
    auto spherical_surface = [](const auto f, const real rho)
    {
        return [=](const real theta, const real phi) { return f(rho, theta, phi); };
    };
    auto expansion_surface = [](const auto f, const real rho)
    {
        return [=](const spherical_index idx) { return f(rho, idx); };
    };

    // Expand an analytic spherical function 'f' defined in spherical-coordinates over the spin-weighted spherical harmonics
    // Effectively changes the signature from {rho,theta,phi} to {rho,{l,m,s}}
    // Uses fixed count function evaluations (assumes worst case scenario)
    auto series_expand = [=,
                          m = max_theta_evals,
                          n = max_phi_evals](const auto f)
    {
        return [=](const real rho, const spherical_index idx)
        {
            return math::chebysev<real>(0, pi, m, [=](const real& theta)
            {
                return math::periodic<real>(0, 2 * pi, n, [=](const real& phi)
                {
#ifdef _DEBUG
                    floating::scoped_exception_enabler fpe;
#endif
                    //return std::conj(math::Y_lms(theta, phi, idx.l, idx.m, idx.s)) * f(rho, theta, phi) * std::sin(theta);
                    return std::conj(math::Y_lms(theta, phi, idx.l, idx.m, idx.s)) * f(rho, theta, phi) * std::pow(rho, 2) * std::sin(theta);
                });
            });
        };
    };

    auto adaptive_series_expand = [=](const auto f)
    {
        using counter = stl::arithmetic_progression_iterator<integral>; // Iterator used as series of numbers

        auto converge_2d = [=](auto g, counter from1, counter to1, counter from2, counter to2)
        {
            std::result_of_t<decltype(g)(const integral&, const integral&)> res = 0;
            counter it2 = to2;
            counter it1 = std::adjacent_find(from1, to1, [=, &it2, &res](const integral& n1, const integral& m1) mutable
            {
                std::cout << "converge_2d it1 at [n1,m1]=[" << n1 << "," << m1 << "] with res = " << res << std::endl;
                it2 = std::find_if(from2, to2, [=, &res](const integral& i2) mutable
                {
                    //return std::abs(g(n1, i2) - (res = g(m1, i2))) < convergence;

                    //return std::abs(g(n1, i2) - g(m1, i2)) < convergence;

                    std::cout << "converge_2d it2 at " << i2 << std::endl;
                    auto gn1i2 = g(n1, i2);
                    auto gm1i2 = g(m1, i2);
                    res = gm1i2;
                    auto diff = std::abs(gn1i2 - gm1i2);

                    if (diff > 1e15)
                        std::clog << "Something's fishy!" << std::endl;
                    
                    return diff < convergence;
                });

                if (it2 != to2)
                    return true;
                else
                    return false;
            });

            if (it1 == to1)
                //throw std::runtime_error{ "Function did not converge on the specified index ranges." };
                std::clog << "Quadrature did not converge on target with precision " << convergence << std::endl;

            if (it2 == to2)
                throw std::logic_error{ "it2 is not expected to be end iterator if it1 was not." };

            //return g(*it1, *it2);
            return res;
        };

        return [=](const real rho, const spherical_index idx)
        {
            auto expand_f = [=](integral m, integral n)
            {
                return math::chebysev<real>(0, pi, m, [=](const real& theta)
                {
                    return math::periodic<real>(0, 2 * pi, n, [=](const real& phi)
                    {
                        floating::scoped_exception_enabler fpe;

                        auto base = math::Y_lms(theta, phi, idx.l, idx.m, idx.s);
                        auto s = std::sin(theta);

                        if (std::abs(base * s) > 1e10)
                            std::clog << "Root of fishy!" << std::endl;

                        return std::conj(math::Y_lms(theta, phi, idx.l, idx.m, idx.s)) * f(rho, theta, phi) * std::pow(rho, 2) * std::sin(theta);
                    });
                });
            };

            return converge_2d(expand_f, counter{ 2, max_theta_evals, 2 }, counter{}, counter{ 2, max_phi_evals }, counter{});
        };
    };
    
    // Reduce a spin-weighted spherical expansion to the original analytic spherical function
    // Effectively changes the signature from {rho,{l,m,s}} to {rho,theta,phi}
    auto reduce_series = [=](const auto f)
    {
        return [=,
                ext = spherical_extent({ 0, 0, 0 },
                                       { L_max, L_max, 0 })](const real rho,
                                                             const real theta,
                                                             const real phi)
        {
            using F = decltype(f);
            using T = std::result_of_t<F(const real, const spherical_index)>;

            return std::real(std::accumulate(ext.cbegin(),
                                             ext.cend(),
                                             static_cast<T>(0),
                                             [=](const T& cum_sum,
                                                 const spherical_index& idx)
            {
                //std::cout << "reduce_series at " << idx << std::endl;

                return cum_sum + f(rho, idx);
            }));
        };
    };
    
    //using namespace math::sws;
    //
    //real_coeff_vector<0, parity::even> a{ L_max };
    //real_coeff_vector<0, parity::even> b{ L_max };
    //
    //real_coeff_vector<0, parity::even> c;
    //
    //c = -a;
    //c = a + b;
    //
    //complex_coeff_vector<0, parity::even> d{ L_max };

    // Transform f(rho,theta,phi) initial condition to a coefficient vector of a_{l,m,s}
    //std::transform(d.extent().cbegin(),
    //               d.extent().cend(),
    //               d.begin(),
    //               expansion_surface(series_expand(analytic_solution),
    //                                 rho_min));

    auto Y_000 = [](real, // rho
                    real theta,
                    real phi)
    {
        return math::Y_lms(theta, phi, 0, 0, 0);
    };

    auto const_one = [](real, // rho
                        real, // theta
                        real) // phi
    {
        return 1.;
    };

    auto ellipsoid = [](real x,
                        real y,
                        real z)
    {
        constexpr real a = (real)1.1;
        constexpr real b = (real)1.1;
        constexpr real c = (real)1.1;

        return std::sqrt(std::pow(a * x, 2) + std::pow(b * y, 2) + std::pow(c * z, 2));
    };

    auto descartes_to_spherical = [](const auto f)
    {
        return [=](const real rho,
                   const real theta,
                   const real phi)
        {
            return f(rho*std::sin(theta)*std::cos(phi),
                     rho*std::sin(theta)*std::sin(phi),
                     rho*std::cos(theta));
        };
    };

    auto spherical_to_descartes = [](const auto f)
    {
        return [=](const real x, const real y, const real z)
        {
            /// <cite>http://keisan.casio.com/exec/system/1359533867</cite>
            return f(std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)),
                     (real)1 / std::tan(y / x),
                     (real)1 / std::tan(std::sqrt(std::pow(x, 2) + std::pow(y, 2)) / z));
        };
    };

    //auto reconstructed_ellipsoid = spherical_to_descartes(reduce_series(adaptive_series_expand(descartes_to_spherical(ellipsoid))));

    std::ofstream out("out.dat", std::ios::ate | std::ios::app);

    for (real theta = 0; theta <= math::constants::pi<real>; theta += math::constants::pi<real> / 100)
        for (real phi = 0; phi <= 2 * math::constants::pi<real>; phi += 2 * math::constants::pi<real> / 100)
        {
            std::cout << theta << '\t' << phi << std::endl;
            constexpr real rho = 1.0;

            auto origi = descartes_to_spherical(ellipsoid);
            auto trans = reduce_series(adaptive_series_expand(origi));

            auto origi_val = origi(rho, theta, phi);
            auto trans_val = trans(rho, theta, phi);

            std::cout <<
                theta << '\t' <<
                phi << '\t' <<
                origi_val << '\t' <<
                trans_val <<
                std::endl;
        }

    out.close();

    //std::cout << Y_000(0, 0, 0) << "\t" << one_per_two * std::sqrt(one / math::pi<real>) << std::endl;
    //std::cout << adaptive_series_expand(const_one)(1, { 0, 0, 0 }) << "\t" << one_per_two * std::sqrt(one / math::pi<real>) * 4 * math::pi<real> << std::endl;
    
    //std::ofstream out( "out.dat", std::ios::ate | std::ios::app );      
    //
    //for (real theta = 0; theta <= math::pi<real>; theta += math::pi<real> / 100)
    //    for (real phi = 0; phi <= 2 * math::pi<real>; phi += 2 * math::pi<real> / 100)
    //        out <<
    //        theta << "\t"<<
    //        phi << "\t" <<
    //        Y_000(0, theta, phi).real() << "\t" <<
    //        Y_000(0, theta, phi).imag() << std::endl;
    //
    //out.close();
    
    //std::cout << Y_000(0, M_PI, M_PI) << "\t" <<
    //    0.5 * std::sqrt(1. / M_PI) << "\t" <<
    //    0.5 * std::sqrt(1. / M_PI) * 4 * M_PI << std::endl;

    //std::transform(d.extent().cbegin(),
    //               d.extent().cend(),
    //               d.begin(),
    //               expansion_surface(series_expand(const_one, 12, 30), rho_min));
    //
    //std::cout << "expansion_surface(series_expand(const_one, 12, 30),\n";
    //std::cout << "                  " << rho_min << ") = " << std::endl;
    //std::cout << "{l,m,s}\t\t(Re,Im)" << std::endl;
    //std::for_each(d.extent().cbegin(), d.extent().cend(), [it = d.cbegin()](const auto& idx) mutable
    //{
    //    std::cout << idx << "\t" << *(it++) << std::endl;
    //});
    //
    //std::cout << std::endl << std::endl;
    //
    //std::transform(d.extent().cbegin(),
    //    d.extent().cend(),
    //    d.begin(),
    //    expansion_surface(series_expand(const_one, 24, 60), rho_min));
    //
    //std::cout << "expansion_surface(series_expand(const_one, 24, 60),\n";
    //std::cout << "                  " << rho_min << ") = " << std::endl;
    //std::cout << "{l,m,s}\t\t(Re,Im)" << std::endl;
    //std::for_each(d.extent().cbegin(), d.extent().cend(), [it = d.cbegin()](const auto& idx) mutable
    //{
    //    std::cout << idx << "\t" << *(it++) << std::endl;
    //});
    //
    //std::cout << std::endl << std::endl;
    //
    //std::transform(d.extent().cbegin(),
    //    d.extent().cend(),
    //    d.begin(),
    //    expansion_surface(series_expand(const_one, 36, 120), rho_min));
    //
    //std::cout << "expansion_surface(series_expand(const_one, 12, 30),\n";
    //std::cout << "                  " << rho_min << ") = " << std::endl;
    //std::cout << "{l,m,s}\t\t(Re,Im)" << std::endl;
    //std::for_each(d.extent().cbegin(), d.extent().cend(), [it = d.cbegin()](const auto& idx) mutable
    //{
    //    std::cout << idx << "\t" << *(it++) << std::endl;
    //});

    // Simulation type aliases
    /*
    using namespace Multipole::stl;
    using radial_index = integral;
    using radial_extent = Radial::Extent<radial_index>;
    using spherical_extent = SWS::Extent<L_max, S_max, integral>;
    using spherical_index = SWS::Index<L_max, S_max, integral>;
    using gaunt_matrix = SpinWeightedGaunt::Matrix<L_max, S_max, integral, real>;
    using real_coeff_vector = SWS::Vector<L_max, S_max, Parity::Even, integral, real>;
    using complex_coeff_vector = SWS::Vector<L_max, S_max, Parity::Even, integral, complex>;
    using real_radial_coeff_vector = Radial::Vector<real, integral, real_coeff_vector>;
    using complex_radial_coeff_vector = Radial::Vector<real, integral, complex_coeff_vector>;
    using state_vector = PDE::StateVector<real_coeff_vector, complex_coeff_vector>;
    using solver = PDE::RK4::Solver<real, state_vector>;

    // Application type aliases
    using timer = std::chrono::high_resolution_clock;
    using time_point = std::chrono::high_resolution_clock::time_point;

    // Convenience enum for equation numbering
    enum Eq
    {
        K = 0,
        k = 1
    };

    // Simulation variable declarations
    constexpr real one = 1, two = 2, three = 3, four = 4, eight = 8, one_per_two = one / two, three_per_two = three / two;
    radial_extent rho_ext;
    spherical_extent lms_ext;
    real_coeff_vector N_kalap, K_kalap, kappa_null, kappa, K_initial, a, d;
    complex_coeff_vector k_initial, b;
    gaunt_matrix gaunt;*/
    /*
    solver rk4;
    real_radial_coeff_vector K_result;
    complex_radial_coeff_vector k_result;

    // Auxiliary variables
    real_coeff_vector K_sq, k_k_bar, round_braces, square_braces, K_curly_braces, K_curly_mul_N, F_K_second_term, F_K_third_term, neumann_temp1, neumann_temp2, kappa_wtf, K_curly_braces_wtf_lhs, K_curly_braces_wtf_rhs, K_final;
    complex_coeff_vector k_sq, k_bar_sq, F_K_second_term_round, F_K_second_term_round_cc, k_round_braces, k_round_braces_cc, k_square_braces, k_curly_braces, k_curly_mul_N, f_k, f_k_second_term, k_final;

    // Application variable declarations
    time_point start, end;

    // Initialize variables
    std::cout << "Initializing... "; std::cout.flush(); start = timer::now();

    rho_ext = radial_extent(rho_min, rho_max);
    lms_ext = spherical_extent({ 0, 0, 0 }, { L_max, L_max, S_max });

    N_kalap    = real_coeff_vector(lms_ext);
    K_kalap    = real_coeff_vector(lms_ext);
    kappa_null = real_coeff_vector(lms_ext);
    kappa      = real_coeff_vector(lms_ext);
    a          = real_coeff_vector(lms_ext);
    b          = complex_coeff_vector(lms_ext);
    d          = real_coeff_vector(lms_ext);
    K_initial  = real_coeff_vector(lms_ext);
    k_initial  = complex_coeff_vector(lms_ext);

    K_result = real_radial_coeff_vector(drho, rho_ext);
    k_result = complex_radial_coeff_vector(drho, rho_ext);

    gaunt = gaunt_matrix(lms_ext);
    auto contract = [&](const auto& lhs, const auto& rhs) { return SpinWeightedGaunt::contract(gaunt, lhs, rhs); };
    auto divide = [&](const auto& lhs, const auto& rhs) { return SpinWeightedGaunt::neumann(gaunt, lhs, rhs, neumann_temp1, neumann_temp2, neumann_length); };
    auto real_cast = [&](const auto& val) { return SWS::real_cast(val); };
    auto H = [&](const real r) { return M / r; };
    auto analytic = [&](const real rho) { return -four * M / (rho * rho * std::sqrt(one + two * H(rho))); }; // based on (3.36)
    auto update_constants = [contract, lms_ext, H, M, real_cast, one, two, one_per_two, eight, &N_kalap, &K_kalap, &kappa_null, &a, &b, &d](const real rho)
    {
        for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            if (i == spherical_index{ 0, 0, 0 })
            {
                N_kalap.at(i) = std::sqrt(one + two * H(rho)); // [1] eq. (3.30)
                K_kalap.at(i) = two / (rho * std::sqrt(one + two * H(rho))); // [1] eq. (3.41) and (3.25)
                kappa_null.at(i) = (eight * M * M) / (rho * rho * std::pow(two * M + rho, 2)); // [3] eq. (1.6) and discussion
                a.at(i) = rho * rho;
                b.at(i) = 0;
            }
            else
            {
                N_kalap.at(i) = 0;
                K_kalap.at(i) = 0;
                kappa_null.at(i) = 0;
                a.at(i) = 0;
                b.at(i) = 0;
            }
        }

        N_kalap = N_kalap / (one_per_two * std::sqrt(one / pi<real>));
        K_kalap = K_kalap / (one_per_two * std::sqrt(one / pi<real>));
        kappa_null = kappa_null /  (one_per_two * std::sqrt(one / pi<real>));
        a = a / (one_per_two * std::sqrt(one / pi<real>));
        b = b / (one_per_two * std::sqrt(one / pi<real>));
        
        d = contract(a, a) - real_cast(contract(b, conjugate(b)));
    };

    // Initialize states
    for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
    {
        if (i == spherical_index{ 0, 0, 0 })
        {
            K_initial.at(i) = analytic(rho_n * drho); // based on (3.36)
            k_initial.at(i) = 0; // based on (3.34)
        }
        else
        {
            K_initial.at(i) = 0;
            k_initial.at(i) = 0;
        }
    }

    // Allocate results
    for (radial_index r = rho_ext.initial(); rho_ext.contains(r); ++r)
    {
        K_result.at(r) = real_coeff_vector(lms_ext);
        k_result.at(r) = complex_coeff_vector(lms_ext);
    }
    K_result.at(rho_n) = K_initial;
    k_result.at(rho_n) = k_initial;

    // Initialize solver
    rk4.lhs() = state_vector(K_initial, k_initial);
    rk4.equation() = [&](state_vector& result, const state_vector& state)
    {
        //
        // NOTE: temporary values used to speed-up calculations. Expression Templating benefits from writing the results
        //       of costly operations if their elements are accessed multiple times. Contractions are both costly and
        //       one element will be accessed multiple times by subsequent contractions.
        //
        using namespace SWS; // For edth/edth_bar operators

        // Calculate kappa
        k_sq = contract(state.get<k>(), state.get<k>());
        K_sq = contract(state.get<K>(), state.get<K>());
        k_k_bar = real_cast(contract(state.get<k>(), conjugate(state.get<k>())));
        k_bar_sq = contract(conjugate(state.get<k>()), conjugate(state.get<k>()));
        round_braces = two * contract(a, k_k_bar) - real_cast(contract(b, k_bar_sq) + contract(conjugate(b), k_sq));
        std::cout << "d{0, 0, 0} = " << d.at({ 0, 0, 0 }) << "\t" << "divide(round_braces, d) = ";
        square_braces = divide(round_braces, d) - one_per_two * K_sq - kappa_null;

        //
        // NOTE: workaround of MSVC internal compiler error. "Simplify code around line..."
        //
        // ORIGINAL: kappa = divide(square_braces, two * state.get<K>());
        //
        kappa_wtf = two * state.get<K>();
        kappa = divide(square_braces, kappa_wtf);

        // Calculate K
        K_curly_braces_wtf_lhs = contract(a, real_cast(edth(conjugate(state.get<k>())) + edth_bar(state.get<k>())));
        K_curly_braces_wtf_rhs = real_cast(contract(b, edth_bar(conjugate(state.get<k>()))) + contract(conjugate(b), edth(state.get<k>())));
        K_curly_braces = K_curly_braces_wtf_lhs - K_curly_braces_wtf_rhs;
        K_curly_mul_N = contract(one_per_two * N_kalap, K_curly_braces);

        F_K_second_term_round = contract(a, conjugate(state.get<k>())) - contract(conjugate(b), state.get<k>());
        F_K_second_term_round_cc = contract(a, state.get<k>()) - contract(b, conjugate(state.get<k>()));
        F_K_second_term = real_cast(contract(F_K_second_term_round, edth(N_kalap)) + contract(F_K_second_term_round_cc, edth_bar(conjugate(N_kalap))));
        F_K_third_term = contract(kappa - one_per_two * state.get<K>(), K_kalap);

        K_final = divide(K_curly_mul_N, d) + divide(F_K_second_term, d) + contract(N_kalap, F_K_third_term);

        // Calculate k
        k_round_braces = contract(a, state.get<k>()) - contract(b, conjugate(state.get<k>()));
        k_round_braces_cc = contract(a, conjugate(state.get<k>())) - contract(conjugate(b), state.get<k>());
        k_square_braces = contract(k_round_braces, edth(state.get<k>())) + contract(k_round_braces_cc, edth(state.get<k>()));
        k_curly_braces = contract(kappa, edth(state.get<K>())) - divide(k_square_braces, d);
        k_curly_mul_N = contract(N_kalap, k_curly_braces);

        f_k_second_term = divide(one_per_two * edth(kappa_null), state.get<K>()) + contract(K_kalap, state.get<k>());
        f_k = -contract(kappa - one_per_two * state.get<K>(), edth(N_kalap)) + contract(N_kalap, f_k_second_term);

        k_final = -divide(k_curly_mul_N, state.get<K>()) - f_k;

        result = PDE::make_equation(
            four * four * K_final,
            k_final
            );
    };

    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // Integrate
    std::cout << "Integrating... "; std::cout.flush(); start = timer::now();
    
    // Outbound integration
    for (radial_index r = rho_n + 1; rho_ext.contains(r); ++r)
    {
        std::cout << "rho = " << r * drho << std::endl;
        update_constants(r * drho);
        rk4.iterate(+drho);
        //{
        //    std::stringstream file_name;
        //    
        //    file_name << "rk4_K_" << r << ".dat";
        //    std::ofstream data_file(file_name.str().c_str(), std::ios::ate);
        //
        //    std::size_t lms_counter = 0;
        //    for (auto lms = lms_ext.initial(); lms_ext.contains(lms); ++lms, ++lms_counter)
        //    {
        //        data_file <<
        //            lms_counter << "\t" <<
        //            rk4.m_k.at(0).get<K>().at(lms) << "\t" <<
        //            rk4.m_k.at(1).get<K>().at(lms) << "\t" <<
        //            rk4.m_k.at(2).get<K>().at(lms) << "\t" <<
        //            rk4.m_k.at(3).get<K>().at(lms) << std::endl;
        //    }
        //}

        K_result.at(r) = rk4.lhs().get<K>();
        k_result.at(r) = rk4.lhs().get<k>();
    }
    
    // Inbound integration
    rk4.lhs() = state_vector(K_initial, k_initial);
    for (radial_index r = rho_n - 1; rho_ext.contains(r); --r)
    {
        std::cout << "rho = " << r * drho << std::endl;
        update_constants(r * drho);
        rk4.iterate(-drho);

        K_result.at(r) = rk4.lhs().get<K>();
        k_result.at(r) = rk4.lhs().get<k>();
    }
    
    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    
    // Export data for plot
    std::cout << "Exporting... "; std::cout.flush(); start = timer::now();*/
    /*
    {
        std::stringstream file_name;

        file_name << "K_" << rho_n << ".dat";
        std::ofstream data_file(file_name.str().c_str(), std::ios::ate);

        std::size_t lms_counter = 0;
        for (auto lms = lms_ext.initial(); lms_ext.contains(lms); ++lms, ++lms_counter)
        {
            data_file << lms_counter << "\t" << K_initial.at(lms) << "\t" << k_initial.at(lms) << std::endl;
        }
    }

    for (radial_index r = rho_ext.initial(); rho_ext.contains(r); ++r)
    {
        std::stringstream file_name;

        file_name << "K_" << r << ".dat";
        std::ofstream data_file(file_name.str().c_str(), std::ios::ate);

        std::size_t lms_counter = 0;
        for (auto lms = lms_ext.initial(); lms_ext.contains(lms); ++lms, ++lms_counter)
        {
            data_file << lms_counter << "\t" << K_result.at(r).at(lms) << "\t" << k_result.at(r).at(lms) << std::endl;
        }
    }
    */
    /*
    std::ofstream data_file("output.dat", std::ios::ate);
    data_file.imbue(std::locale(""));

    for (radial_index r = rho_ext.initial(); rho_ext.contains(r); ++r)
    {    
        data_file <<
            r * drho << "\t" <<
            K_result.at(r).at(spherical_index{ 0, 0, 0 }) << "\t" <<
            //k_result.at(r).at(spherical_index{ 0, 0, 0 }) << " \t" <<
            analytic(r * drho) << std::endl;
    }
    
    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    */
    return EXIT_SUCCESS;
}
