#pragma warning(disable : 4503)

#include <STL_Example5_Schwarzschild.hpp>

int main()
{
    // Simulation params
    using integral = std::int32_t;      // Type used to represent integral values
    using real = double;                // Type used to represent real values
    using complex = std::complex<real>; // Type used to represent complex values
    constexpr integral L_max = 3;       // Maximum multipole values for series expansion
    constexpr integral S_max = 3;       // Maximum spin values for series expansion
    const real drho = 0.25f;            // Separation of rho coordinates
    const integral rho_min = 32;        // Coordinate # of innermost lattice site
    const integral rho_max = 40;        // Coordinate # of outermost lattice site
    const integral rho_n = 32;          // Coordinate # to start integration from
    const real M = 100.f;               // Black hole mass param
    const real neumann_length = 0.001f; // Neumann-series expansion precision

    // Simulation type aliases
    using namespace Multipole::stl;
    using radial_index = integral;
    using radial_extent = Radial::Extent<radial_index>;
    using spherical_extent = SpinWeightedSpherical::Extent<L_max, S_max, integral>;
    using spherical_index = SpinWeightedSpherical::Index<L_max, S_max, integral>;
    using gaunt_matrix = SpinWeightedGaunt::Matrix<L_max, S_max, integral, real>;
    using real_coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, real>;
    using complex_coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, complex>;
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
    gaunt_matrix gaunt;

    real_coeff_vector K_sq, k_k_bar, round_braces, square_braces, K_curly_braces, K_curly_mul_N, F_K_second_term, F_K_third_term;
    complex_coeff_vector k_sq, k_bar_sq, F_K_second_term_round, F_K_second_term_round_cc, k_round_braces, k_round_braces_cc, k_square_braces, k_curly_braces, k_curly_mul_N, f_k, f_k_second_term;

    solver rk4;

    real_radial_coeff_vector K_result;
    complex_radial_coeff_vector k_result;

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
    auto contract = [&](const auto& lhs,const auto& rhs) { return SpinWeightedGaunt::contract(gaunt, lhs, rhs); };
    auto divide = [&](const auto& lhs, const auto& rhs) { return SpinWeightedGaunt::neumann(gaunt, lhs, rhs, neumann_length); };
    auto real_cast = [&](const auto& val) { return SpinWeightedSpherical::real_cast(val); };
    auto H = [&](const real& r) { return M / r; };
    auto update_constants = [&](const real rho)
    {
        for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            if (i == spherical_index{ 0, 0, 0 })
            {
                N_kalap.at(i) = std::sqrt(one + two * H(rho)); // based on (3.30)
                K_kalap.at(i) = two / (rho * std::sqrt(one + two * H(rho))); // based on (3.41) and (3.25)
                kappa_null.at(i) = (eight * M * M) / (rho * rho * std::pow(two * M + rho, 2));
                a.at(i) = rho * rho;
                b.at(i) = 0;
                d = contract(a, a) - real_cast(contract(b, conjugate(b)));
            }
            else
            {
                N_kalap.at(i) = 0;
                K_kalap.at(i) = 0;
                kappa_null.at(i) = 0;
                a.at(i) = 0;
                b.at(i) = 0;
                d = contract(a, a) - real_cast(contract(b, conjugate(b)));
            }
        }    
    };
    auto rho = rho_n * drho;

    // Initialize states
    for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
    {
        if (i == spherical_index{ 0, 0, 0 })
        {
            K_initial.at(i) = -four * M / (rho * rho * std::pow(one + two * H(rho), three_per_two)); // based on (3.36)
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

    // Initialize solver
    rk4.lhs() = state_vector(K_initial, k_initial);
    rk4.equation() = [&](state_vector& result, const state_vector& state)
    {
        std::cout << "update_constants" << std::endl;
        update_constants(rho);

        //
        // NOTE: temporary values used to speed-up calculations. Expression Templating benefits from writing the results
        //       of costly operations if their elements are accessed multiple times. Contractions are both costly and
        //       one element will be accessed multiple times by subsequent contractions.
        //
        std::cout << "Calculate kappa" << std::endl;
        // Calculate kappa
        k_sq = contract(state.get<k>(), state.get<k>());
        K_sq = contract(state.get<K>(), state.get<K>());
        k_k_bar = real_cast(contract(state.get<k>(), conjugate(state.get<k>())));
        k_bar_sq = contract(conjugate(state.get<k>()), conjugate(state.get<k>()));
        round_braces = two * contract(a, k_k_bar) - real_cast(contract(b, k_bar_sq) + contract(conjugate(b), k_sq));
        square_braces = divide(round_braces, d) - one_per_two * K_sq - kappa_null;

        kappa = divide(square_braces, two * state.get<K>());
        std::cout << "Calculate K" << std::endl;
        // Calculate K
        K_curly_braces = contract(a, real_cast(SpinWeightedSpherical::edth(conjugate(state.get<k>())) + SpinWeightedSpherical::edth_bar(state.get<k>()))) - real_cast(contract(b, SpinWeightedSpherical::edth_bar(conjugate(state.get<k>()))) + contract(conjugate(b), SpinWeightedSpherical::edth(state.get<k>())));
        K_curly_mul_N = contract(one_per_two * N_kalap, K_curly_braces);

        F_K_second_term_round = contract(a, conjugate(state.get<k>())) - contract(conjugate(b), state.get<k>());
        F_K_second_term_round_cc = contract(a, state.get<k>()) - contract(b, conjugate(state.get<k>()));
        F_K_second_term = real_cast(contract(F_K_second_term_round, SpinWeightedSpherical::edth(N_kalap)) + contract(F_K_second_term_round_cc, SpinWeightedSpherical::edth_bar(conjugate(N_kalap))));
        F_K_third_term = contract(kappa - one_per_two * state.get<K>(), K_kalap);
        std::cout << "Calculate k" << std::endl;
        // Calculate k
        k_round_braces = contract(a, state.get<k>()) - contract(b, conjugate(state.get<k>()));
        k_round_braces_cc = contract(a, conjugate(state.get<k>())) - contract(conjugate(b), state.get<k>());
        k_square_braces = contract(k_round_braces, SpinWeightedSpherical::edth(state.get<k>())) + contract(k_round_braces_cc, SpinWeightedSpherical::edth(state.get<k>()));
        k_curly_braces = contract(kappa, SpinWeightedSpherical::edth(state.get<K>())) - divide(k_square_braces, d);
        k_curly_mul_N = contract(N_kalap, k_curly_braces);

        f_k_second_term = divide(one_per_two * SpinWeightedSpherical::edth(kappa_null), state.get<K>()) + contract(K_kalap, state.get<k>());
        f_k = -contract(kappa - one_per_two * state.get<K>(), SpinWeightedSpherical::edth(N_kalap)) + contract(N_kalap, f_k_second_term);
        std::cout << "Calculate result" << std::endl;

        result = PDE::make_equation(
            divide(K_curly_mul_N, d) + divide(F_K_second_term, d) + contract(N_kalap, F_K_third_term),
            -divide(k_curly_mul_N, state.get<K>()) - f_k
            );
    };

    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // Integrate
    std::cout << "Integrating... "; std::cout.flush(); start = timer::now();

    // Outbound integration
    for (radial_index r = rho_n; rho_ext.contains(r); ++r)
    {
        std::cout << "rho = " << rho << std::endl;
        rk4.iterate(+drho);
        rho += drho;

        K_result.at(r) = rk4.lhs().get<K>();
        k_result.at(r) = rk4.lhs().get<k>();
    }
    /*
    // Inbound integration
    rk4.setLHS(state_vector(K_initial, k_initial));
    for (radial_index r = rho_n; rho_ext.contains(r); --r)
    {
        rk4.iterate(-drho);
        rho -= drho;

        K_result.at(r) = rk4.getLHS().get<K>();
        k_result.at(r) = rk4.getLHS().get<k>();
    }
    */
    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // Export data for plot
    std::cout << "Exporting... "; std::cout.flush(); start = timer::now();

    std::ofstream data_file("output.dat", std::ios::ate);

    for (radial_index r = rho_ext.initial(); rho_ext.contains(r); ++r)
    {
        data_file <<
            r * drho << "\t" <<
            K_result.at(r).at(spherical_index{ 0, 0, 0 }) << "\t" <<
            k_result.at(r).at(spherical_index{ 0, 0, 0 }) << " \t" <<
            -four * M / (std::pow(r * drho, 2) * std::sqrt(one + 2 * H(r * drho))) << std::endl;
    }

    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    return EXIT_SUCCESS;
}
