#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

class BlackScholesCalculator {
private:
    // Normal cumulative distribution function
    double normalCDF(double x) {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    }

    // Normal probability density function
    double normalPDF(double x) {
        return exp(-0.5 * x * x) / sqrt(2.0 * M_PI);
    }

public:
    struct OptionResults {
        double callPrice;
        double putPrice;
        double callDelta;
        double putDelta;
        double gamma;
        double callTheta;
        double putTheta;
        double vega;
        double callRho;
        double putRho;
    };

    OptionResults calculate(double S, double K, double T, double r, double sigma, double q = 0.0) {
        OptionResults results;

        // Check for invalid inputs
        if (T <= 0 || sigma <= 0) {
            results = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            return results;
        }

        // Calculate d1 and d2
        double d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);

        // Calculate cumulative normal distributions
        double Nd1 = normalCDF(d1);
        double Nd2 = normalCDF(d2);
        double NminusD1 = normalCDF(-d1);
        double NminusD2 = normalCDF(-d2);

        // Calculate option prices
        results.callPrice = S * exp(-q * T) * Nd1 - K * exp(-r * T) * Nd2;
        results.putPrice = K * exp(-r * T) * NminusD2 - S * exp(-q * T) * NminusD1;

        // Calculate Greeks
        double phi_d1 = normalPDF(d1);

        // Delta
        results.callDelta = exp(-q * T) * Nd1;
        results.putDelta = exp(-q * T) * (Nd1 - 1.0);

        // Gamma (same for calls and puts)
        results.gamma = (exp(-q * T) * phi_d1) / (S * sigma * sqrt(T));

        // Theta (daily)
        results.callTheta = ((-S * phi_d1 * sigma * exp(-q * T)) / (2.0 * sqrt(T))
                           - r * K * exp(-r * T) * Nd2
                           + q * S * exp(-q * T) * Nd1) / 365.0;

        results.putTheta = ((-S * phi_d1 * sigma * exp(-q * T)) / (2.0 * sqrt(T))
                          + r * K * exp(-r * T) * NminusD2
                          - q * S * exp(-q * T) * NminusD1) / 365.0;

        // Vega (same for calls and puts, per 1% volatility change)
        results.vega = (S * exp(-q * T) * phi_d1 * sqrt(T)) / 100.0;

        // Rho (per 1% interest rate change)
        results.callRho = (K * T * exp(-r * T) * Nd2) / 100.0;
        results.putRho = (-K * T * exp(-r * T) * NminusD2) / 100.0;

        return results;
    }
};

void printHeader() {
    cout << "\n========================================\n";
    cout << "  Black-Scholes Options Calculator\n";
    cout << "         with Greeks\n";
    cout << "========================================\n\n";
}

void printResults(const BlackScholesCalculator::OptionResults& results) {
    cout << "\n" << string(50, '=') << "\n";
    cout << "                 RESULTS\n";
    cout << string(50, '=') << "\n\n";

    // Option Prices
    cout << "OPTION PRICES:\n";
    cout << string(20, '-') << "\n";
    cout << fixed << setprecision(4);
    cout << "Call Price:  $" << setw(10) << results.callPrice << "\n";
    cout << "Put Price:   $" << setw(10) << results.putPrice << "\n\n";

    // Greeks
    cout << "GREEKS:\n";
    cout << string(20, '-') << "\n";

    cout << "Delta (price sensitivity):\n";
    cout << "  Call Delta: " << setw(10) << results.callDelta << "\n";
    cout << "  Put Delta:  " << setw(10) << results.putDelta << "\n\n";

    cout << "Gamma (delta sensitivity): " << setw(8) << results.gamma << "\n\n";

    cout << "Theta (daily time decay):\n";
    cout << "  Call Theta: " << setw(10) << results.callTheta << "\n";
    cout << "  Put Theta:  " << setw(10) << results.putTheta << "\n\n";

    cout << "Vega (volatility sensitivity): " << setw(6) << results.vega << "\n\n";

    cout << "Rho (interest rate sensitivity):\n";
    cout << "  Call Rho:   " << setw(10) << results.callRho << "\n";
    cout << "  Put Rho:    " << setw(10) << results.putRho << "\n";

    cout << "\n" << string(50, '=') << "\n";
}

void printExplanation() {
    cout << "\nGREEKS EXPLANATION:\n";
    cout << string(30, '-') << "\n";
    cout << "Delta:  How much option price changes per $1 change in stock price\n";
    cout << "Gamma:  How much Delta changes per $1 change in stock price\n";
    cout << "Theta:  How much option price decreases per day (time decay)\n";
    cout << "Vega:   How much option price changes per 1% change in volatility\n";
    cout << "Rho:    How much option price changes per 1% change in interest rate\n\n";
}

int main() {
    BlackScholesCalculator calculator;
    double S, K, T, r, sigma, q;
    char choice;

    printHeader();

    do {
        cout << "Enter the following parameters:\n\n";

        cout << "Spot Price (current stock price): $";
        cin >> S;

        cout << "Strike Price: $";
        cin >> K;

        cout << "Time to Expiry (in years, e.g., 0.25 for 3 months): ";
        cin >> T;

        cout << "Risk-Free Rate (as decimal, e.g., 0.05 for 5%): ";
        cin >> r;

        cout << "Volatility (as decimal, e.g., 0.20 for 20%): ";
        cin >> sigma;

        cout << "Dividend Yield (as decimal, 0 if none): ";
        cin >> q;

        // Input validation
        if (T <= 0) {
            cout << "\nError: Time to expiry must be positive!\n";
            continue;
        }
        if (sigma <= 0) {
            cout << "\nError: Volatility must be positive!\n";
            continue;
        }
        if (S <= 0 || K <= 0) {
            cout << "\nError: Prices must be positive!\n";
            continue;
        }

        // Calculate and display results
        BlackScholesCalculator::OptionResults results = calculator.calculate(S, K, T, r, sigma, q);
        printResults(results);
        printExplanation();

        cout << "Do you want to calculate another option? (y/n): ";
        cin >> choice;
        cout << "\n";

    } while (choice == 'y' || choice == 'Y');

    cout << "Thank you for using the Black-Scholes Calculator!\n";
    return 0;
}

// Compilation instructions:
// g++ -o blackscholes blackscholes.cpp -std=c++11