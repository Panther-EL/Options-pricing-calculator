#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
using namespace std;

class BinomialTreeCalculator {
private:
    int steps;
    double dt;
    double u, d, p;  // up factor, down factor, probability
    vector<vector<double>> stockPrices;
    vector<vector<double>> optionValues;

    void buildStockTree(double S, double r, double sigma, double T, double q) {
        dt = T / steps;
        u = exp(sigma * sqrt(dt));
        d = 1.0 / u;
        p = (exp((r - q) * dt) - d) / (u - d);

        // Initialize stock price tree
        stockPrices.assign(steps + 1, vector<double>(steps + 1, 0.0));

        // Fill the stock price tree
        for (int i = 0; i <= steps; i++) {
            for (int j = 0; j <= i; j++) {
                stockPrices[i][j] = S * pow(u, j) * pow(d, i - j);
            }
        }
    }

    double calculatePayoff(double S, double K, bool isCall) {
        if (isCall) {
            return max(S - K, 0.0);
        } else {
            return max(K - S, 0.0);
        }
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

    BinomialTreeCalculator(int n = 50) : steps(n) {}

    void setSteps(int n) {
        steps = n;
    }

    double calculateOption(double S, double K, double T, double r, double sigma,
                          double q, bool isCall, bool isAmerican = false) {
        buildStockTree(S, r, sigma, T, q);

        // Initialize option value tree
        optionValues.assign(steps + 1, vector<double>(steps + 1, 0.0));

        // Calculate option values at expiration
        for (int j = 0; j <= steps; j++) {
            optionValues[steps][j] = calculatePayoff(stockPrices[steps][j], K, isCall);
        }

        // Work backwards through the tree
        double discountFactor = exp(-r * dt);

        for (int i = steps - 1; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                // European option value
                double europeanValue = discountFactor * (p * optionValues[i + 1][j + 1] +
                                                       (1 - p) * optionValues[i + 1][j]);

                if (isAmerican) {
                    // American option - compare with early exercise value
                    double exerciseValue = calculatePayoff(stockPrices[i][j], K, isCall);
                    optionValues[i][j] = max(europeanValue, exerciseValue);
                } else {
                    optionValues[i][j] = europeanValue;
                }
            }
        }

        return optionValues[0][0];
    }

    OptionResults calculateGreeks(double S, double K, double T, double r, double sigma,
                                 double q, bool isAmerican = false) {
        OptionResults results;

        // Base case calculations
        results.callPrice = calculateOption(S, K, T, r, sigma, q, true, isAmerican);
        results.putPrice = calculateOption(S, K, T, r, sigma, q, false, isAmerican);

        // Delta calculation (finite difference)
        double deltaShift = S * 0.01; // 1% shift
        double callPriceUp = calculateOption(S + deltaShift, K, T, r, sigma, q, true, isAmerican);
        double callPriceDown = calculateOption(S - deltaShift, K, T, r, sigma, q, true, isAmerican);
        double putPriceUp = calculateOption(S + deltaShift, K, T, r, sigma, q, false, isAmerican);
        double putPriceDown = calculateOption(S - deltaShift, K, T, r, sigma, q, false, isAmerican);

        results.callDelta = (callPriceUp - callPriceDown) / (2 * deltaShift);
        results.putDelta = (putPriceUp - putPriceDown) / (2 * deltaShift);

        // Gamma calculation
        results.gamma = (callPriceUp - 2 * results.callPrice + callPriceDown) / (deltaShift * deltaShift);

        // Theta calculation (finite difference with time)
        if (T > 1.0/365.0) { // At least 1 day remaining
            double timeShift = 1.0/365.0; // 1 day
            double callPriceTheta = calculateOption(S, K, T - timeShift, r, sigma, q, true, isAmerican);
            double putPriceTheta = calculateOption(S, K, T - timeShift, r, sigma, q, false, isAmerican);

            results.callTheta = callPriceTheta - results.callPrice;
            results.putTheta = putPriceTheta - results.putPrice;
        } else {
            results.callTheta = 0;
            results.putTheta = 0;
        }

        // Vega calculation
        double vegaShift = 0.01; // 1% volatility shift
        double callPriceVegaUp = calculateOption(S, K, T, r, sigma + vegaShift, q, true, isAmerican);
        double putPriceVegaUp = calculateOption(S, K, T, r, sigma + vegaShift, q, false, isAmerican);

        results.vega = callPriceVegaUp - results.callPrice;

        // Rho calculation
        double rhoShift = 0.01; // 1% rate shift
        double callPriceRhoUp = calculateOption(S, K, T, r + rhoShift, sigma, q, true, isAmerican);
        double putPriceRhoUp = calculateOption(S, K, T, r + rhoShift, sigma, q, false, isAmerican);

        results.callRho = callPriceRhoUp - results.callPrice;
        results.putRho = putPriceRhoUp - results.putPrice;

        return results;
    }

    void printTreeInfo(double S, double r, double sigma, double T, double q) {
        buildStockTree(S, r, sigma, T, q);
        cout << "\nTREE PARAMETERS:\n";
        cout << string(20, '-') << "\n";
        cout << fixed << setprecision(4);
        cout << "Number of steps:     " << steps << "\n";
        cout << "Time step (dt):      " << dt << "\n";
        cout << "Up factor (u):       " << u << "\n";
        cout << "Down factor (d):     " << d << "\n";
        cout << "Risk-neutral prob:   " << p << "\n\n";
    }
};

void printHeader() {
    cout << "\n========================================\n";
    cout << "    Binomial Tree Options Calculator\n";
    cout << "         with Greeks\n";
    cout << "========================================\n\n";
}

void printResults(const BinomialTreeCalculator::OptionResults& results, bool isAmerican = false) {
    cout << "\n" << string(50, '=') << "\n";
    cout << "                 RESULTS";
    if (isAmerican) cout << " (AMERICAN)";
    else cout << " (EUROPEAN)";
    cout << "\n" << string(50, '=') << "\n\n";

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
    cout << "\nBINOMIAL TREE MODEL EXPLANATION:\n";
    cout << string(40, '-') << "\n";
    cout << "The binomial tree model discretizes stock price movements into\n";
    cout << "up and down moves over discrete time periods.\n\n";

    cout << "GREEKS EXPLANATION:\n";
    cout << string(30, '-') << "\n";
    cout << "Delta:  How much option price changes per $1 change in stock price\n";
    cout << "Gamma:  How much Delta changes per $1 change in stock price\n";
    cout << "Theta:  How much option price decreases per day (time decay)\n";
    cout << "Vega:   How much option price changes per 1% change in volatility\n";
    cout << "Rho:    How much option price changes per 1% change in interest rate\n\n";

    cout << "ADVANTAGES OF BINOMIAL MODEL:\n";
    cout << string(30, '-') << "\n";
    cout << "• Can handle American options (early exercise)\n";
    cout << "• Can handle dividend payments\n";
    cout << "• More intuitive than Black-Scholes\n";
    cout << "• Converges to Black-Scholes as steps increase\n\n";
}

int main() {
    BinomialTreeCalculator calculator;
    double S, K, T, r, sigma, q;
    int steps;
    char choice, optionType, exerciseType;

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

        cout << "Number of time steps (recommended: 50-200): ";
        cin >> steps;

        cout << "Exercise style - (E)uropean or (A)merican: ";
        cin >> exerciseType;

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
        if (steps <= 0 || steps > 1000) {
            cout << "\nError: Number of steps must be between 1 and 1000!\n";
            continue;
        }

        // Set parameters
        calculator.setSteps(steps);
        bool isAmerican = (exerciseType == 'A' || exerciseType == 'a');

        // Display tree parameters
        calculator.printTreeInfo(S, r, sigma, T, q);

        // Calculate and display results
        BinomialTreeCalculator::OptionResults results = calculator.calculateGreeks(S, K, T, r, sigma, q, isAmerican);
        printResults(results, isAmerican);
        printExplanation();

        cout << "Do you want to calculate another option? (y/n): ";
        cin >> choice;
        cout << "\n";

    } while (choice == 'y' || choice == 'Y');

    cout << "Thank you for using the Binomial Tree Calculator!\n";
    return 0;
}

// Compilation instructions:
// g++ -o binomialtree binomialtree.cpp -std=c++11