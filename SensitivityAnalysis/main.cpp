#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <string>
using namespace std;

// Share option structure (should be in common header)
struct Option {
    double S;      //Spot price
    double K;      //Strike price
    double T;      //Time to maturity(years)
    double r;      //Risk-free rate
    double sigma;  //Volatility
    std::string type;//"call" or "put"

    Option(double S = 100, double K = 105, double T = 1.0,
           double r = 0.05, double sigma = 0.2, std::string type = "call")
        : S(S), K(K), T(T), r(r), sigma(sigma), type(type) {}
};

// Sensitivity analysis class
class SensitivityAnalyzer {
public:
    using PricingFunction = std::function<double(const Option&)>;

    struct ParameterRange {
        std::string name;
        std::vector<double> values;
    };

    struct AnalysisResult {
        std::map<std::string, std::vector<double>> parameter_values;
        std::map<std::string, std::vector<double>> prices;
        std::map<std::string, std::vector<double>> deltas;
    };

    SensitivityAnalyzer(PricingFunction pricing_func)
        : pricing_func_(pricing_func) {}

    AnalysisResult analyze(const Option& base_option,
                   const std::vector<ParameterRange>& ranges) {
        AnalysisResult result;
        const double base_price = pricing_func_(base_option);

        for (const auto& range : ranges) {
            std::vector<double> prices;
            std::vector<double> deltas;

            for (double value : range.values) {
                Option modified = base_option;
                if (range.name == "sigma") modified.sigma = value;
                else if (range.name == "r") modified.r = value;
                else if (range.name == "T") modified.T = value;
                else if (range.name == "S") modified.S = value;
                else if (range.name == "K") modified.K = value;

                double price = pricing_func_(modified);
                prices.push_back(price);
                deltas.push_back(price - base_price);
            }

            result.parameter_values[range.name] = range.values;
            result.prices[range.name] = prices;
            result.deltas[range.name] = deltas;
        }
        return result;
    }

    void export_csv(const AnalysisResult& result, const std::string& filename) {
        std::ofstream file(filename);
        file << "Parameter,Value,Price,Delta\n";

        for (const auto& [param, values] : result.parameter_values) {
            for (size_t i = 0; i < values.size(); ++i) {
                file << param << ","
                     << values[i] << ","
                     << result.prices.at(param)[i] << ","
                     << result.deltas.at(param)[i] << "\n";
            }
        }
    }

    void print_dominant_factor(const AnalysisResult& result) {
        std::string dominant_param;
        double max_range = 0.0;

        for (const auto& [param, deltas] : result.deltas) {
            auto [min_it, max_it] = std::minmax_element(deltas.begin(), deltas.end());
            double range = *max_it - *min_it;

            if (range > max_range) {
                max_range = range;
                dominant_param = param;
            }
        }

        std::cout << "\nDominant Factor Analysis:\n";
        std::cout << "Most sensitive parameter: " << dominant_param << "\n";
        std::cout << "Price variation range: " << max_range << "\n";
    }

private:
    PricingFunction pricing_func_;
};

//Example pricing function (Black-Scholes)
double black_scholes(const Option& opt) {
    double d1 = (log(opt.S / opt.K) + (opt.r + 0.5 * pow(opt.sigma, 2)) * opt.T)
                / (opt.sigma * sqrt(opt.T));
    double d2 = d1 - opt.sigma * sqrt(opt.T);

    if (opt.type == "call") {
        return opt.S * 0.5 * (1 + erf(d1 / sqrt(2)))
               - opt.K * exp(-opt.r * opt.T) * 0.5 * (1 + erf(d2 / sqrt(2)));
    } else {
        return opt.K * exp(-opt.r * opt.T) * 0.5 * (1 + erf(-d2 / sqrt(2)))
               - opt.S * 0.5 * (1 + erf(-d1 / sqrt(2)));
    }
}

//Helper function to generate value ranges
std::vector<double> linspace(double start, double end, size_t points) {
    std::vector<double> result(points);
    double step = (end - start) / (points - 1);

    for (size_t i = 0; i < points; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

int main() {
    //Create base option
    Option base_option(100, 105, 1.0, 0.05, 0.2, "call");

    //Set up sensitivity analyzer with pricing function
    SensitivityAnalyzer analyzer(black_scholes);

    //Define parameter ranges to test
    std::vector<SensitivityAnalyzer::ParameterRange> ranges = {
        {"sigma", linspace(0.1, 0.5, 20)},
        {"r", linspace(0.01, 0.1, 10)},
        {"T", linspace(0.1, 2.0, 15)}
    };

    //Run analysis
    auto results = analyzer.analyze(base_option, ranges);

    //Export results
    analyzer.export_csv(results, "sensitivity_analysis.csv");
    analyzer.print_dominant_factor(results);

    std::cout << "Analysis complete. Results saved to sensitivity_analysis.csv\n";

    return 0;
}