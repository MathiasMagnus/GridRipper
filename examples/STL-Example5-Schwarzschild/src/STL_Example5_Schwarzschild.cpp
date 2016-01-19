#include <STL_Example5_Schwarzschild.hpp>

int main(int argc, char* argv[])
{
    // Simulation params
    using integral = std::int32_t;      // Type used to represent integral values
    using real = float;                 // Type used to represent real values
    using complex = std::complex<real>; // Type used to represent complex values
    constexpr integral L_max = 7;       // Maximum multipole values for series expansion
    constexpr integral S_max = 3;       // Maximum spin values for series expansion
    const real drho = 0.01f;            // Separation of rho coordinates
    const integral rho_min = 1;         // Coordinate # of innermost lattice site
    const integral rho_max = 3;         // Coordinate # of outermost lattice site
    const integral rho_n = 2;           // Coordinate # to start integration from
    const real M = 100;                 // Black hole mass param

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
    real one_per_two, three_per_two;
    radial_extent rho_ext;
    spherical_extent lms_ext;
    real_coeff_vector N_kalap, K_kalap, kappa_null, kappa_var, K_var, temp1, temp2;
    complex_coeff_vector k_var, a, b, d;
    gaunt_matrix gaunt;

    solver rk4;

    real_radial_coeff_vector K_result;
    complex_radial_coeff_vector k_result;

    // Application variable declarations
    time_point start, end;

    // Initialize variables
    std::cout << "Initializing... "; std::cout.flush(); start = timer::now();

    one_per_two = 0.5f;
    three_per_two = 0.5f;
    rho_ext = radial_extent(rho_min, rho_max);
    lms_ext = spherical_extent({ 0, 0, 0 }, { L_max, L_max, S_max });

    N_kalap    = real_coeff_vector(lms_ext);
    K_kalap    = real_coeff_vector(lms_ext);
    kappa_null = real_coeff_vector(lms_ext);
    kappa_var  = real_coeff_vector(lms_ext);
    a          = complex_coeff_vector(lms_ext);
    b          = complex_coeff_vector(lms_ext); // NOTE: In reality, this is the only complex: a,d are real
    d          = complex_coeff_vector(lms_ext);
    K_var      = real_coeff_vector(lms_ext);
    k_var      = complex_coeff_vector(lms_ext);

    K_result = real_radial_coeff_vector(drho, rho_ext);
    k_result = complex_radial_coeff_vector(drho, rho_ext);

    gaunt.calculate();
    auto contract = [&gaunt](const auto& lhs,const auto& rhs) { return SpinWeightedGaunt::contract(gaunt, lhs, rhs); };
    auto H = [&](const real& r) { return M / r; };
    auto rho = rho_n * drho;

    for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
    {
        if (i == spherical_index{ 0, 0, 0 })
        {
            N_kalap.at(i) = std::sqrt(1 + 2 * H(rho)); // based on (3.30)
            K_kalap.at(i) = 2.f / (rho * std::sqrt(1 + 2 * H(rho))); // based on (3.41) and (3.25)
            K_var.at(i) = -4 * M / (std::pow((rho), 2) * std::pow(1 + 2 * H(rho), three_per_two)); // based on (3.36)
            k_var.at(i) = 0; // based on (3.34)
        }
        else
        {
            N_kalap.at(i) = 0;
            K_kalap.at(i) = 0;
            K_var.at(i) = 0;
            k_var.at(i) = 0;
        }
    }

    // Initialize solver
    rk4.setLHS(state_vector(K_var, k_var));
    rk4.setEquationSystem([&](const state_vector& state)
    {
        for (auto i = lms_ext.initial(); lms_ext.contains(i); ++i)
        {
            N_kalap.at(i) = std::sqrt(1 + 2 * H(rho)); // based on (3.30)
            K_kalap.at(i) = 2.f / (rho * std::sqrt(1 + 2 * H(rho))); // based on (3.41) and (3.25)
            kappa_null.at(i) = (8 * M * M) / (rho * rho * std::pow(2 * M + rho, 2));
            a.at(i) = rho * rho;
            b.at(i) = 0;
        }
        d = contract(a, a);// - contract(b, complex(1, -1) * b);
        kappa_var = ((2.f * contract(a , contract(state.get<k>(), complex(1, -1) * state.get<k>())) - contract(b , contract(complex(1, -1) * state.get<k>(), complex(1, -1) * state.get<k>())) - contract(complex(1, -1) * b, contract(state.get<k>(), state.get<k>()))) / d - one_per_two * contract(state.get<K>(), state.get<K>() - kappa_null)) / ( 2.f * state.get<K>() );
        
        //
        // NOTE: temporary values used to speed-up calculations. Expression Templating benefits from writing the results
        //       of costly operations if their elements are accessed multiple times. Contractions are both costly and
        //       one element will be accessed multiple times by subsequent contractions.
        //
        temp1 = one_per_two * N_kalap * (a * (SpinWeightedSpherical::eth(complex(1, -1) * state.get<k>()) + SpinWeightedSpherical::eth_bar(state.get<k>())) - contract(b, SpinWeightedSpherical::eth_bar(complex(1, -1) * state.get<k>())) - contract(complex(1, -1) * b, SpinWeightedSpherical::eth(state.get<k>()))) / d + contract(N_kalap, (kappa_var - one_per_two * state.get<K>()));
        temp2 = -N_kalap * (contract(kappa_var, SpinWeightedSpherical::eth(state.get<K>())) - (contract(contract(a, state.get<k>()) - contract(b, state.get<k>()), SpinWeightedSpherical::eth(complex(1, -1) * state.get<k>())) + contract(contract(a, complex(1, -1) *  state.get<k>()) - contract(complex(1, -1) * b, state.get<k>()), SpinWeightedSpherical::eth(state.get<k>()))) / d) / state.get<K>() -
            contract(N_kalap, K_kalap);

        rho += drho;

        return state_vector(
            contract(temp1, K_kalap),
            contract(temp2, state.get<k>()));
    });

    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    // Integrate
    std::cout << "Integrating... "; std::cout.flush(); start = timer::now();

    // Outbound integration
    for (radial_index r = rho_n; rho_ext.contains(r); ++r)
    {
        std::cout << "r = " << r << std::endl;

        rk4.iterate(drho);

        K_result.at(r) = rk4.getLHS().get<K>();
        k_result.at(r) = rk4.getLHS().get<k>();
    }
    /*
    // Inbound integration
    for (radial_index r = rho; rho_ext.contains(r); --r)
    {
        rk4.iterate(drho);

        K_result.at(rho) = rk4.getLHS().get<K>();
        k_result.at(rho) = rk4.getLHS().get<k>();
    }
    */
    end = timer::now(); std::cout << "done. " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    return EXIT_SUCCESS;
}