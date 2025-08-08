#include <random>
#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>

enum class OptionType { CALL, PUT };

struct Option {
    double S;       // Spot price
    double K;       // Strike price
    double T;       // Time to maturity (in years)
    double r;       // Risk-free rate
    double sigma;   // Volatility
    OptionType type;
};

double monte_carlo_antithetic(const Option& opt, int n_sims, std::mt19937& gen) {
    std::normal_distribution<> dist(0.0, 1.0);
    double payoff_sum = 0.0;

    for (int i = 0; i < n_sims / 2; ++i) {
        double z = dist(gen);
        for (double zi : {z, -z}) {
            double ST = opt.S * std::exp((opt.r - 0.5 * opt.sigma * opt.sigma) * opt.T +
                                         opt.sigma * std::sqrt(opt.T) * zi);
            double payoff = (opt.type == OptionType::CALL)
                            ? std::max(ST - opt.K, 0.0)
                            : std::max(opt.K - ST, 0.0);
            payoff_sum += payoff;
        }
    }

    return std::exp(-opt.r * opt.T) * (payoff_sum / n_sims);
}

double monte_carlo_control_variate(const Option& opt, int n_sims, std::mt19937& gen) {
    std::normal_distribution<> dist(0.0, 1.0);
    std::vector<double> payoffs(n_sims), controls(n_sims);

    for (int i = 0; i < n_sims; ++i) {
        double z = dist(gen);
        double ST = opt.S * std::exp((opt.r - 0.5 * opt.sigma * opt.sigma) * opt.T +
                                     opt.sigma * std::sqrt(opt.T) * z);
        payoffs[i] = (opt.type == OptionType::CALL)
                     ? std::max(ST - opt.K, 0.0)
                     : std::max(opt.K - ST, 0.0);
        controls[i] = ST;
    }

    double mean_payoff = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / n_sims;
    double mean_control = std::accumulate(controls.begin(), controls.end(), 0.0) / n_sims;

    double cov = 0.0, var_control = 0.0;
    for (int i = 0; i < n_sims; ++i) {
        cov += (payoffs[i] - mean_payoff) * (controls[i] - mean_control);
        var_control += std::pow(controls[i] - mean_control, 2);
    }
    cov /= n_sims;
    var_control /= n_sims;

    double b = cov / var_control;
    double expected_control = opt.S * std::exp(opt.r * opt.T);
    double adjusted_payoff = mean_payoff - b * (mean_control - expected_control);

    return std::exp(-opt.r * opt.T) * adjusted_payoff;
}

double finite_diff_delta(const Option& opt, double bump, int n_sims, std::mt19937& gen) {
    Option up = opt; up.S += bump;
    Option down = opt; down.S -= bump;
    return (monte_carlo_antithetic(up, n_sims, gen) - monte_carlo_antithetic(down, n_sims, gen)) / (2.0 * bump);
}

double finite_diff_vega(const Option& opt, double bump, int n_sims, std::mt19937& gen) {
    Option up = opt; up.sigma += bump;
    Option down = opt; down.sigma -= bump;
    return (monte_carlo_antithetic(up, n_sims, gen) - monte_carlo_antithetic(down, n_sims, gen)) / (2.0 * bump);
}

int main() {
    Option opt;
    int n_sims;
    double bump;
    int opt_type_input;

    // Prompt user for input
    std::cout << "Enter Spot Price (S): ";
    std::cin >> opt.S;

    std::cout << "Enter Strike Price (K): ";
    std::cin >> opt.K;

    std::cout << "Enter Time to Maturity (T in years): ";
    std::cin >> opt.T;

    std::cout << "Enter Risk-Free Rate (r): ";
    std::cin >> opt.r;

    std::cout << "Enter Volatility (sigma): ";
    std::cin >> opt.sigma;

    std::cout << "Option Type (0 for CALL, 1 for PUT): ";
    std::cin >> opt_type_input;
    opt.type = (opt_type_input == 0) ? OptionType::CALL : OptionType::PUT;

    std::cout << "Enter Number of Simulations: ";
    std::cin >> n_sims;

    std::cout << "Enter Bump Size for Greeks (e.g., 0.01): ";
    std::cin >> bump;

    std::mt19937 gen(42);  // Seed once and reuse

    double price_antithetic = monte_carlo_antithetic(opt, n_sims, gen);

    gen.seed(42);
    double price_control = monte_carlo_control_variate(opt, n_sims, gen);

    gen.seed(42);
    double delta = finite_diff_delta(opt, bump, n_sims, gen);

    gen.seed(42);
    double vega = finite_diff_vega(opt, bump, n_sims, gen);

    std::cout << "\n--- Results ---\n";
    std::cout << "Monte Carlo (Antithetic):      " << price_antithetic << "\n";
    std::cout << "Monte Carlo (Control Variate): " << price_control << "\n";
    std::cout << "Delta (Finite Difference):     " << delta << "\n";
    std::cout << "Vega (Finite Difference):      " << vega << "\n";

    return 0;
}

