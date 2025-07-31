#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

class BlackScholesModel {
private:
    double S;  // Current stock price
    double K;  // Strike price
    double T;  // Time to expiration (in years)
    double r;  // Risk-free interest rate
    double sigma; // Volatility

    // Standard normal cumulative distribution function
    double normalCDF(double x) const {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    }

    // Standard normal probability density function
    double normalPDF(double x) const {
        return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
    }

    // Calculate d1 parameter
    double calculateD1() const {
        return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    }

    // Calculate d2 parameter
    double calculateD2() const {
        return calculateD1() - sigma * sqrt(T);
    }

public:
    // Constructor
    BlackScholesModel(double stockPrice, double strikePrice, double timeToExpiry,
                     double riskFreeRate, double volatility)
        : S(stockPrice), K(strikePrice), T(timeToExpiry), r(riskFreeRate), sigma(volatility) {

        if (S <= 0 || K <= 0 || T <= 0 || sigma <= 0) {
            throw std::invalid_argument("All parameters must be positive");
        }
    }

    // European Call Option Price
    double callPrice() const {
        if (T == 0) {
            return std::max(S - K, 0.0);
        }

        double d1 = calculateD1();
        double d2 = calculateD2();

        return S * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
    }

    // European Put Option Price
    double putPrice() const {
        if (T == 0) {
            return std::max(K - S, 0.0);
        }

        double d1 = calculateD1();
        double d2 = calculateD2();

        return K * exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
    }

    // Call Option Delta (sensitivity to stock price)
    double callDelta() const {
        if (T == 0) {
            return (S > K) ? 1.0 : 0.0;
        }
        return normalCDF(calculateD1());
    }

    // Put Option Delta
    double putDelta() const {
        if (T == 0) {
            return (S < K) ? -1.0 : 0.0;
        }
        return normalCDF(calculateD1()) - 1.0;
    }

    // Gamma (sensitivity of delta to stock price) - same for call and put
    double gamma() const {
        if (T == 0) return 0.0;
        double d1 = calculateD1();
        return normalPDF(d1) / (S * sigma * sqrt(T));
    }

    // Theta (time decay) for call option
    double callTheta() const {
        if (T == 0) return 0.0;
        double d1 = calculateD1();
        double d2 = calculateD2();

        return (-S * normalPDF(d1) * sigma / (2.0 * sqrt(T))
                - r * K * exp(-r * T) * normalCDF(d2)) / 365.0; // Per day
    }

    // Theta (time decay) for put option
    double putTheta() const {
        if (T == 0) return 0.0;
        double d1 = calculateD1();
        double d2 = calculateD2();

        return (-S * normalPDF(d1) * sigma / (2.0 * sqrt(T))
                + r * K * exp(-r * T) * normalCDF(-d2)) / 365.0; // Per day
    }

    // Vega (sensitivity to volatility) - same for call and put
    double vega() const {
        if (T == 0) return 0.0;
        double d1 = calculateD1();
        return S * normalPDF(d1) * sqrt(T) / 100.0; // Per 1% change in volatility
    }

    // Rho for call option (sensitivity to interest rate)
    double callRho() const {
        if (T == 0) return 0.0;
        double d2 = calculateD2();
        return K * T * exp(-r * T) * normalCDF(d2) / 100.0; // Per 1% change in rate
    }

    // Rho for put option
    double putRho() const {
        if (T == 0) return 0.0;
        double d2 = calculateD2();
        return -K * T * exp(-r * T) * normalCDF(-d2) / 100.0; // Per 1% change in rate
    }

    // Implied volatility calculation using Newton-Raphson method
    double impliedVolatility(double marketPrice, bool isCall, double tolerance = 1e-6, int maxIterations = 100) const {
        double vol = 0.3; // Initial guess

        for (int i = 0; i < maxIterations; ++i) {
            BlackScholesModel tempModel(S, K, T, r, vol);
            double price = isCall ? tempModel.callPrice() : tempModel.putPrice();
            double vega = tempModel.vega() * 100.0; // Convert back from percentage

            if (std::abs(vega) < tolerance) break;

            double diff = price - marketPrice;
            if (std::abs(diff) < tolerance) return vol;

            vol = vol - diff / vega;

            if (vol <= 0) vol = 0.01; // Ensure positive volatility
        }

        return vol;
    }

    // Display all option information
    void displayOptionInfo() const {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n=== BLACK-SCHOLES OPTION PRICING MODEL ===" << std::endl;
        std::cout << "Stock Price (S): $" << S << std::endl;
        std::cout << "Strike Price (K): $" << K << std::endl;
        std::cout << "Time to Expiry (T): " << T << " years" << std::endl;
        std::cout << "Risk-free Rate (r): " << r * 100 << "%" << std::endl;
        std::cout << "Volatility (σ): " << sigma * 100 << "%" << std::endl;

        std::cout << "\n--- OPTION PRICES ---" << std::endl;
        std::cout << "Call Option Price: $" << callPrice() << std::endl;
        std::cout << "Put Option Price: $" << putPrice() << std::endl;

        std::cout << "\n--- THE GREEKS ---" << std::endl;
        std::cout << "Call Delta: " << callDelta() << std::endl;
        std::cout << "Put Delta: " << putDelta() << std::endl;
        std::cout << "Gamma: " << gamma() << std::endl;
        std::cout << "Call Theta: $" << callTheta() << " per day" << std::endl;
        std::cout << "Put Theta: $" << putTheta() << " per day" << std::endl;
        std::cout << "Vega: $" << vega() << " per 1% vol change" << std::endl;
        std::cout << "Call Rho: $" << callRho() << " per 1% rate change" << std::endl;
        std::cout << "Put Rho: $" << putRho() << " per 1% rate change" << std::endl;
    }

    // Verify put-call parity: C - P = S - K*e^(-rT)
    bool verifyPutCallParity(double tolerance = 1e-6) const {
        double leftSide = callPrice() - putPrice();
        double rightSide = S - K * exp(-r * T);
        return std::abs(leftSide - rightSide) < tolerance;
    }
};

// Example usage and testing
int main() {
    try {
        // Example: Tesla stock option (hypothetical parameters)
        double stockPrice = 250.0;      // Current Tesla stock price
        double strikePrice = 260.0;     // Strike price
        double timeToExpiry = 0.25;     // 3 months to expiration
        double riskFreeRate = 0.05;     // 5% risk-free rate
        double volatility = 0.35;       // 35% volatility

        BlackScholesModel tesla(stockPrice, strikePrice, timeToExpiry, riskFreeRate, volatility);

        // Display comprehensive option information
        tesla.displayOptionInfo();

        // Verify put-call parity
        std::cout << "\n--- VALIDATION ---" << std::endl;
        std::cout << "Put-Call Parity Check: " << (tesla.verifyPutCallParity() ? "PASSED" : "FAILED") << std::endl;

        // Example of implied volatility calculation
        double marketCallPrice = 15.0; // Hypothetical market price
        double impliedVol = tesla.impliedVolatility(marketCallPrice, true);
        std::cout << "Implied Volatility from market price $" << marketCallPrice << ": "
                  << impliedVol * 100 << "%" << std::endl;

        // Scenario analysis
        std::cout << "\n--- SCENARIO ANALYSIS ---" << std::endl;
        std::vector<double> spotPrices = {220, 240, 260, 280, 300};

        std::cout << "Spot\tCall Price\tPut Price\tCall Delta\tPut Delta" << std::endl;
        std::cout << "----\t---------\t--------\t----------\t---------" << std::endl;

        for (double spot : spotPrices) {
            BlackScholesModel scenario(spot, strikePrice, timeToExpiry, riskFreeRate, volatility);
            std::cout << spot << "\t" << std::setprecision(2) << scenario.callPrice() << "\t\t"
                      << scenario.putPrice() << "\t\t" << std::setprecision(3)
                      << scenario.callDelta() << "\t\t" << scenario.putDelta() << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

