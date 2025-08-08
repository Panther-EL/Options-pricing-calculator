#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <map>

class SensitivityAnalyzer {
private:
    struct Variable {
        std::string name;
        double baseValue;
        double minRange;
        double maxRange;
        std::vector<double> testValues;
        std::vector<double> outputValues;
        double sensitivity;

        Variable(const std::string& n, double base, double min, double max)
            : name(n), baseValue(base), minRange(min), maxRange(max), sensitivity(0.0) {}
    };

    std::vector<Variable> variables;
    std::function<double(const std::vector<double>&)> targetFunction;
    double baseOutput;
    std::string functionDescription;

public:
    void setFunction(std::function<double(const std::vector<double>&)> func, const std::string& desc) {
        targetFunction = func;
        functionDescription = desc;
    }

    void addVariable(const std::string& name, double baseValue, double minRange, double maxRange) {
        variables.emplace_back(name, baseValue, minRange, maxRange);
        std::cout << "Added variable: " << name << " (base: " << baseValue
                  << ", range: [" << minRange << ", " << maxRange << "])\n";
    }

    void clearVariables() {
        variables.clear();
        std::cout << "All variables cleared.\n";
    }

    void displayVariables() const {
        if (variables.empty()) {
            std::cout << "No variables defined.\n";
            return;
        }

        std::cout << "\n=== DEFINED VARIABLES ===\n";
        std::cout << std::setw(15) << "Variable" << std::setw(12) << "Base Value"
                  << std::setw(15) << "Min Range" << std::setw(15) << "Max Range\n";
        std::cout << std::string(60, '-') << "\n";

        for (const auto& var : variables) {
            std::cout << std::setw(15) << var.name
                      << std::setw(12) << std::fixed << std::setprecision(4) << var.baseValue
                      << std::setw(15) << var.minRange
                      << std::setw(15) << var.maxRange << "\n";
        }
        std::cout << "\n";
    }

    void calculateBaseOutput() {
        if (variables.empty() || !targetFunction) {
            std::cout << "Error: No variables or function defined.\n";
            return;
        }

        std::vector<double> baseValues;
        for (const auto& var : variables) {
            baseValues.push_back(var.baseValue);
        }

        baseOutput = targetFunction(baseValues);
        std::cout << "Base output calculated: " << std::fixed << std::setprecision(6)
                  << baseOutput << "\n";
    }

    void performOneAtATimeAnalysis(int numPoints = 10) {
        if (variables.empty() || !targetFunction) {
            std::cout << "Error: No variables or function defined.\n";
            return;
        }

        calculateBaseOutput();
        std::cout << "\nPerforming One-at-a-Time Sensitivity Analysis...\n";
        std::cout << "Testing " << numPoints << " points per variable\n\n";

        for (size_t i = 0; i < variables.size(); i++) {
            std::cout << "Analyzing variable: " << variables[i].name << "... ";

            variables[i].testValues.clear();
            variables[i].outputValues.clear();

            // Generate test values across the range
            double step = (variables[i].maxRange - variables[i].minRange) / (numPoints - 1);

            for (int j = 0; j < numPoints; j++) {
                double testValue = variables[i].minRange + j * step;
                variables[i].testValues.push_back(testValue);

                // Create input vector with all base values except the current variable
                std::vector<double> inputs;
                for (size_t k = 0; k < variables.size(); k++) {
                    if (k == i) {
                        inputs.push_back(testValue);
                    } else {
                        inputs.push_back(variables[k].baseValue);
                    }
                }

                double output = targetFunction(inputs);
                variables[i].outputValues.push_back(output);
            }

            // Calculate sensitivity (max output change / max input change)
            auto minMaxOutput = std::minmax_element(variables[i].outputValues.begin(),
                                                   variables[i].outputValues.end());
            double outputRange = *minMaxOutput.second - *minMaxOutput.first;
            double inputRange = variables[i].maxRange - variables[i].minRange;

            variables[i].sensitivity = (inputRange != 0) ? outputRange / inputRange : 0.0;

            std::cout << "Done!\n";
        }
    }

    void displaySensitivityResults() const {
        if (variables.empty()) {
            std::cout << "No analysis results available.\n";
            return;
        }

        std::cout << "\n=== SENSITIVITY ANALYSIS RESULTS ===\n";
        std::cout << "Function: " << functionDescription << "\n";
        std::cout << "Base output: " << std::fixed << std::setprecision(6) << baseOutput << "\n\n";

        std::cout << std::setw(15) << "Variable" << std::setw(15) << "Sensitivity"
                  << std::setw(15) << "Min Output" << std::setw(15) << "Max Output"
                  << std::setw(15) << "Output Range\n";
        std::cout << std::string(75, '-') << "\n";

        // Sort variables by sensitivity for ranking
        std::vector<size_t> indices(variables.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
                  [this](size_t i, size_t j) {
                      return variables[i].sensitivity > variables[j].sensitivity;
                  });

        for (size_t idx : indices) {
            const auto& var = variables[idx];
            if (!var.outputValues.empty()) {
                auto minMax = std::minmax_element(var.outputValues.begin(), var.outputValues.end());
                double minOutput = *minMax.first;
                double maxOutput = *minMax.second;
                double range = maxOutput - minOutput;

                std::cout << std::setw(15) << var.name
                          << std::setw(15) << std::setprecision(6) << var.sensitivity
                          << std::setw(15) << minOutput
                          << std::setw(15) << maxOutput
                          << std::setw(15) << range << "\n";
            }
        }
        std::cout << "\n";
    }

    void displayDetailedResults(const std::string& variableName) const {
        auto it = std::find_if(variables.begin(), variables.end(),
                               [&variableName](const Variable& v) { return v.name == variableName; });

        if (it == variables.end()) {
            std::cout << "Variable '" << variableName << "' not found.\n";
            return;
        }

        const Variable& var = *it;
        if (var.testValues.empty()) {
            std::cout << "No analysis data available for variable '" << variableName << "'.\n";
            return;
        }

        std::cout << "\n=== DETAILED RESULTS FOR: " << variableName << " ===\n";
        std::cout << std::setw(15) << "Input Value" << std::setw(15) << "Output Value"
                  << std::setw(20) << "% Change from Base\n";
        std::cout << std::string(50, '-') << "\n";

        for (size_t i = 0; i < var.testValues.size(); i++) {
            double percentChange = ((var.outputValues[i] - baseOutput) / baseOutput) * 100.0;
            std::cout << std::setw(15) << std::fixed << std::setprecision(4) << var.testValues[i]
                      << std::setw(15) << std::setprecision(6) << var.outputValues[i]
                      << std::setw(19) << std::setprecision(2) << percentChange << "%\n";
        }
        std::cout << "\n";
    }

    void performTornadoAnalysis() const {
        if (variables.empty()) {
            std::cout << "No analysis results available.\n";
            return;
        }

        std::cout << "\n=== TORNADO DIAGRAM (Variable Importance Ranking) ===\n";

        // Create sorted list by sensitivity
        std::vector<std::pair<double, std::string>> sensitivities;
        for (const auto& var : variables) {
            sensitivities.emplace_back(var.sensitivity, var.name);
        }

        std::sort(sensitivities.rbegin(), sensitivities.rend());

        // Find max sensitivity for scaling
        double maxSens = sensitivities.empty() ? 1.0 : sensitivities[0].first;
        if (maxSens == 0) maxSens = 1.0;

        for (const auto& pair : sensitivities) {
            double normalized = (pair.first / maxSens) * 50; // Scale to 50 chars max
            std::cout << std::setw(15) << pair.second << " |";
            for (int i = 0; i < static_cast<int>(normalized); i++) {
                std::cout << "█";
            }
            std::cout << " (" << std::fixed << std::setprecision(4) << pair.first << ")\n";
        }
        std::cout << "\n";
    }
};

// Predefined mathematical functions for testing
namespace TestFunctions {
    // Linear function: a*x1 + b*x2 + c*x3
    double linear(const std::vector<double>& vars) {
        if (vars.size() < 3) return 0;
        return 2.5 * vars[0] + 1.8 * vars[1] + 0.9 * vars[2];
    }

    // Quadratic function: a*x1^2 + b*x2^2 + c*x1*x2
    double quadratic(const std::vector<double>& vars) {
        if (vars.size() < 2) return 0;
        return 0.5 * vars[0] * vars[0] + 0.3 * vars[1] * vars[1] + 0.7 * vars[0] * vars[1];
    }

    // Financial NPV function: sum of discounted cash flows
    double npv(const std::vector<double>& vars) {
        if (vars.size() < 4) return 0;
        double initialInvestment = vars[0];
        double cashFlow = vars[1];
        double discountRate = vars[2];
        double periods = vars[3];

        double npvValue = -initialInvestment;
        for (int i = 1; i <= static_cast<int>(periods); i++) {
            npvValue += cashFlow / std::pow(1 + discountRate, i);
        }
        return npvValue;
    }

    // Physics projectile motion: range calculation
    double projectileRange(const std::vector<double>& vars) {
        if (vars.size() < 3) return 0;
        double velocity = vars[0];      // m/s
        double angle = vars[1];         // degrees
        double gravity = vars[2];       // m/s^2

        double angleRad = angle * M_PI / 180.0;
        return (velocity * velocity * std::sin(2 * angleRad)) / gravity;
    }
}

void displayMainMenu() {
    std::cout << "\n=== SENSITIVITY ANALYSIS TOOL ===\n";
    std::cout << "1. Define custom variables\n";
    std::cout << "2. Use predefined test functions\n";
    std::cout << "3. Run sensitivity analysis\n";
    std::cout << "4. Display results summary\n";
    std::cout << "5. View detailed variable results\n";
    std::cout << "6. Show tornado diagram\n";
    std::cout << "7. Display current variables\n";
    std::cout << "8. Clear all variables\n";
    std::cout << "9. About sensitivity analysis\n";
    std::cout << "10. Exit\n";
    std::cout << "Choose option (1-10): ";
}

void displayTestFunctions() {
    std::cout << "\n=== PREDEFINED TEST FUNCTIONS ===\n";
    std::cout << "1. Linear Function: 2.5*x1 + 1.8*x2 + 0.9*x3\n";
    std::cout << "2. Quadratic Function: 0.5*x1² + 0.3*x2² + 0.7*x1*x2\n";
    std::cout << "3. Financial NPV: Net Present Value calculation\n";
    std::cout << "4. Physics: Projectile motion range\n";
    std::cout << "Choose function (1-4): ";
}

int getValidatedInt(const std::string& prompt, int min, int max) {
    int value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value && value >= min && value <= max) {
            return value;
        }
        std::cout << "Invalid input. Please enter a number between " << min << " and " << max << ".\n";
        std::cin.clear();
        std::cin.ignore(10000, '\n');
    }
}

double getValidatedDouble(const std::string& prompt) {
    double value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            return value;
        }
        std::cout << "Invalid input. Please enter a valid number.\n";
        std::cin.clear();
        std::cin.ignore(10000, '\n');
    }
}

std::string getValidatedString(const std::string& prompt) {
    std::string value;
    std::cout << prompt;
    std::cin.ignore();
    std::getline(std::cin, value);
    return value;
}

void setupPredefinedFunction(SensitivityAnalyzer& analyzer, int choice) {
    analyzer.clearVariables();

    switch (choice) {
        case 1: // Linear function
            analyzer.setFunction(TestFunctions::linear, "Linear: 2.5*x1 + 1.8*x2 + 0.9*x3");
            analyzer.addVariable("x1", 10.0, 5.0, 15.0);
            analyzer.addVariable("x2", 8.0, 4.0, 12.0);
            analyzer.addVariable("x3", 6.0, 3.0, 9.0);
            break;

        case 2: // Quadratic function
            analyzer.setFunction(TestFunctions::quadratic, "Quadratic: 0.5*x1² + 0.3*x2² + 0.7*x1*x2");
            analyzer.addVariable("x1", 5.0, 2.0, 8.0);
            analyzer.addVariable("x2", 4.0, 1.0, 7.0);
            break;

        case 3: // Financial NPV
            analyzer.setFunction(TestFunctions::npv, "NPV: Net Present Value");
            analyzer.addVariable("Initial_Investment", 100000, 80000, 120000);
            analyzer.addVariable("Annual_Cash_Flow", 25000, 20000, 30000);
            analyzer.addVariable("Discount_Rate", 0.1, 0.05, 0.15);
            analyzer.addVariable("Periods", 10, 5, 15);
            break;

        case 4: // Projectile motion
            analyzer.setFunction(TestFunctions::projectileRange, "Physics: Projectile Range");
            analyzer.addVariable("Velocity_m/s", 50.0, 30.0, 70.0);
            analyzer.addVariable("Angle_degrees", 45.0, 15.0, 75.0);
            analyzer.addVariable("Gravity_m/s²", 9.81, 9.0, 10.5);
            break;
    }
}

void displayAbout() {
    std::cout << "\n=== ABOUT SENSITIVITY ANALYSIS ===\n";
    std::cout << "Sensitivity analysis studies how uncertainty in output can be attributed\n";
    std::cout << "to different sources of uncertainty in input variables.\n\n";

    std::cout << "This tool performs:\n";
    std::cout << "• One-at-a-time analysis: Varies one parameter while keeping others constant\n";
    std::cout << "• Sensitivity ranking: Orders variables by their impact on output\n";
    std::cout << "• Tornado diagrams: Visual representation of variable importance\n\n";

    std::cout << "Applications:\n";
    std::cout << "• Financial modeling (NPV, risk assessment)\n";
    std::cout << "• Engineering design optimization\n";
    std::cout << "• Scientific model validation\n";
    std::cout << "• Decision making under uncertainty\n\n";
}

int main() {
    SensitivityAnalyzer analyzer;

    std::cout << "Welcome to the Sensitivity Analysis Tool!\n";

    while (true) {
        displayMainMenu();
        int choice = getValidatedInt("", 1, 10);

        switch (choice) {
            case 1: { // Define custom variables
                int numVars = getValidatedInt("Enter number of variables to define: ", 1, 10);

                std::cout << "\nNote: You'll need to implement your own function.\n";
                std::cout << "For now, using a simple sum function for demonstration.\n\n";

                analyzer.setFunction([](const std::vector<double>& vars) {
                    double sum = 0;
                    for (double v : vars) sum += v;
                    return sum;
                }, "Custom: Sum of all variables");

                for (int i = 0; i < numVars; i++) {
                    std::string name = getValidatedString("Enter variable " + std::to_string(i+1) + " name: ");
                    double base = getValidatedDouble("Enter base value: ");
                    double minVal = getValidatedDouble("Enter minimum range value: ");
                    double maxVal = getValidatedDouble("Enter maximum range value: ");

                    if (minVal >= maxVal) {
                        std::cout << "Warning: Minimum should be less than maximum. Swapping values.\n";
                        std::swap(minVal, maxVal);
                    }

                    analyzer.addVariable(name, base, minVal, maxVal);
                }
                break;
            }

            case 2: { // Use predefined functions
                displayTestFunctions();
                int funcChoice = getValidatedInt("", 1, 4);
                setupPredefinedFunction(analyzer, funcChoice);
                std::cout << "Predefined function setup complete!\n";
                break;
            }

            case 3: { // Run analysis
                int numPoints = getValidatedInt("Enter number of test points per variable (5-50): ", 5, 50);
                analyzer.performOneAtATimeAnalysis(numPoints);
                std::cout << "Sensitivity analysis completed!\n";
                break;
            }

            case 4: // Display results summary
                analyzer.displaySensitivityResults();
                break;

            case 5: { // View detailed results
                std::string varName = getValidatedString("Enter variable name for detailed results: ");
                analyzer.displayDetailedResults(varName);
                break;
            }

            case 6: // Tornado diagram
                analyzer.performTornadoAnalysis();
                break;

            case 7: // Display variables
                analyzer.displayVariables();
                break;

            case 8: // Clear variables
                analyzer.clearVariables();
                break;

            case 9: // About
                displayAbout();
                break;

            case 10: // Exit
                std::cout << "Thank you for using the Sensitivity Analysis Tool!\n";
                return 0;
        }

        std::cout << "\nPress Enter to continue...";
        std::cin.ignore();
        std::cin.get();
    }

    return 0;
}