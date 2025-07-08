// C++ adapation of identifying catenary linear and directed-cycle compartmental models computationally
// Author: Brett Hajdaj
// Date: 2025 06 20

#include <iostream>
#include <vector>
#include <cln/cln.h>
#include <ginac/ginac.h>
#include <unordered_set>
#include <chrono>
#include <sstream>
#include <gmpxx.h>

using namespace std;
using namespace GiNaC;

/*
 * Struct to hold user input for the model
 */
struct UserInput {
    int n;
    int in;
    int p;
    int d;
    bool is_directed = false;
    vector<ex> parameters;
    vector<int> all_compartments;
    vector<int> P;
    unordered_set<int> leak_compartments;
    map<string, ex> symbol_map;

    // Helper function to get or create symbols
    ex get_symbol(const string& name) {
        auto it = symbol_map.find(name);
        if (it != symbol_map.end()) {
            return it->second;
        }

        symbol s(name);
        symbol_map[name] = s;
        return s;
    }
};

/**
 * Creates symbols for the edge parameters and leak parameters based on user input
 * @param input the user input containing the number of compartments and leak locations
 */
void create_symbols(UserInput& input) {
    // Create edge parameters
    for (int i = 0; i < input.n - 1; i++) {
        string forward = "k" + to_string(i + 2) + to_string(i + 1);
        string backward = "k" + to_string(i + 1) + to_string(i + 2);
        
        input.parameters.push_back(input.get_symbol(forward));
        input.parameters.push_back(input.get_symbol(backward));
    }

    // Create leak parameters
    for (int comp : input.leak_compartments) {
        string leak_name = "k0" + to_string(comp);
        input.parameters.push_back(input.get_symbol(leak_name));
    }
}

/**
 * Computes the outgoing sums for each compartment
 * @param input the user input containing the number of compartments and leak locations
 * @param compartment_set the set of compartments to compute the outgoing sums for
 * @return a vector of expressions representing the outgoing sums for each compartment
 */
vector<ex> out_Ml(UserInput& input, const vector<int>& compartment_set) {
    vector<ex> outgoing_sums;
    outgoing_sums.reserve(compartment_set.size());

    for (int elem : compartment_set) {
        bool is_leak = input.leak_compartments.count(elem);

        if (elem == 1 && input.n == 1) {
            // Special case for single compartment
            ex sum = 0;
            if (is_leak)
                sum += input.get_symbol("k01");
            outgoing_sums.push_back(sum);
            continue;
        }
    
        // If n > 1,
        if (elem == 1) {
            ex sum = input.get_symbol("k21");
            if (is_leak)
                sum += input.get_symbol("k01");
            outgoing_sums.push_back(sum);
        }
        else if (elem == input.n) {
            string kn = "k" + to_string(input.n - 1) + to_string(input.n);
            ex sum = input.get_symbol(kn);
            if (is_leak)
                sum += input.get_symbol("k0" + to_string(input.n));
            outgoing_sums.push_back(sum);
        }
        else {
            string k_prev = "k" + to_string(elem - 1) + to_string(elem);
            string k_next = "k" + to_string(elem + 1) + to_string(elem);
            ex sum = input.get_symbol(k_prev) + input.get_symbol(k_next);
            if (is_leak)
                sum += input.get_symbol("k0" + to_string(elem));
            outgoing_sums.push_back(sum);
        }
    }
    return outgoing_sums;
}

/**
 * Computes the product of edge parameters between two compartments
 * @param i the starting compartment
 * @param j the ending compartment
 * @param input the user input containing the parameters
 * @return the product of all edge parameters between compartments i and j
 */
ex kappa_ij(UserInput& input, int i, int j) {
    ex product = 1;
    for (int k = i; k < j; k++) {
        string param = "k" + to_string(k + 1) + to_string(k);
        product *= input.get_symbol(param);
    }
    return product;
}

/**
 * Computes the product of edge parameters along a subset
 * @param subset the vector of integers representing the subset of compartments
 * @param input the user input containing the parameters
 * @return the product of all edge parameters over the subset
 */
ex kappa(UserInput& input, const vector<int>& subset) {
    ex product = 1;
    for (int curr : subset) {
        string forward = "k" + to_string(curr) + to_string(curr + 1);
        string backward = "k" + to_string(curr + 1) + to_string(curr);
        product *= input.get_symbol(forward) * input.get_symbol(backward);
    }
    return product;
}

/**
 * Helper function to generate subsets of [n-1] without consecutive numbers
 * @param n the number of compartments
 * @param index the current index in the recursion
 * @param current_subset the current subset being built
 * @param gamma_set the set of all valid subsets found so far
 */
void gamma_helper(int n, int index, vector<int>& current_subset, vector<vector<int>>& gamma_set) {
    // Base case index > n
    if (index > n) {
        if (!current_subset.empty())
            gamma_set.push_back(current_subset);
        return;
    }

    // Path 1: Exclude the current index
    gamma_helper(n, index + 1, current_subset, gamma_set);

    // Path 2: Include the current index if it does not create consecutive numbers
    if (current_subset.empty() || current_subset.back() + 1 != index) {
        current_subset.push_back(index);
        gamma_helper(n, index + 2, current_subset, gamma_set);
        current_subset.pop_back();
    }
}

/**
 * Generates the set of nonempty subsets of [n-1] without consecutive numbers
 * @param input the user input containing the number of compartments
 * @return a vector of vectors representing the subsets
 */
vector<vector<int>> gamma(UserInput& input) {
    vector<vector<int>> gamma_set;
    vector<int> current_subset;
    gamma_helper(input.n - 1, 1, current_subset, gamma_set);
    return gamma_set;
}

/**
 * Takes the set of nonempty subsets of [n-1] without consecutive numbers,
 * and returns a new set where each element of the subset is incremented by 1.
 * @param gamma_set the set of subsets to be processed
 * @return a new set of subsets with each element incremented by 1
 */
vector<vector<int>> gamma_plus(vector<vector<int>>& gamma_set) {
    vector<vector<int>> gamma_plus_set;
    for (const auto& subset : gamma_set) {
        vector<int> temp_subset;
        for (const auto& elem : subset) {
            // Increment each element of the subset by 1
            temp_subset.push_back(elem + 1);
        }
        // Combine and sort subset with the incremented subset
        temp_subset.insert(temp_subset.end(), subset.begin(), subset.end());
        sort(temp_subset.begin(), temp_subset.end());
        
        // Add the combined subset to the gamma_plus_set
        gamma_plus_set.push_back(temp_subset);
    }
    return gamma_plus_set;
}


/**
 * Computes the intersection of two sets
 * @param set1 the first set
 * @param set2 the second set
 * @return a vector containing the elements that are in both sets
 */
vector<int> intersection_of_sets(const vector<int>& set1, const vector<int>& set2) {
    vector<int> intersection = {};
    for (const auto& elem : set1) {
        if (find(set2.begin(), set2.end(), elem) != set2.end()) {
            intersection.push_back(elem);
        }
    }
    return intersection;
}

/**
 * Computes the elementary symmetric polynomials for a given set of variables
 * @param variables the vector of variables for which to compute the polynomials
 * @param k the degree of the polynomial
 * @return a vector containing the elementary symmetric polynomials up to degree k
 */
vector<ex> elementary_symmetric_polynomials(const vector<ex>& variables, int k) {
    int n = variables.size();
    // Invalid case
    if (k == 0) 
        return {ex(1)};

    vector<ex> e(k + 1, 0);
    e[0] = 1;

    // Compute elementary symmetric polynomials
    for (int i = 0; i < n; ++i) {
        int limit = min(i + 1, k);
        for (int j = limit; j >= 1; --j) {
            e[j] = e[j] + e[j - 1] * variables[i];
        }
    }
    return e;
}

/**
 * Computes a new set of compartments excluding specified elements
 * @param existing_set the original set of compartments
 * @param excluded_elems the elements to be excluded from the new set
 * @return a new vector containing the elements of existing_set excluding excluded_elems
 */
vector<int> compartment_set(vector<int> existing_set, vector<int> excluded_elems = {}) {
    vector<int> new_set;

    // If no exclusions, return existing set
    if (excluded_elems.empty())
        return existing_set;

    for (const auto& elem : existing_set) {
        if (find(excluded_elems.begin(), excluded_elems.end(), elem) == excluded_elems.end()) {
            new_set.push_back(elem);
        }
    }
    return new_set;
}

/**
 * Computes the union of two sets
 * @param set1 the first set
 * @param set2 the second set
 * @return a vector containing the union of the two sets
 */
vector<int> union_sets(const vector<int>& set1, const vector<int>& set2) {
    vector<int> result = set1;
    for (const auto& elem : set2) {
        if (find(result.begin(), result.end(), elem) == result.end()) {
            result.push_back(elem);
        }
    }
    return result;
}

/**
 * Computes the Jacobian matrix of coefficients with respect to parameters
 * @param coefficients The vector of coefficient expressions
 * @param parameters The vector of parameters to differentiate against
 * @return A matrix (vector of vectors) of partial derivatives
 */
matrix compute_jacobian(const vector<ex>& coefficients, const vector<ex>& parameters) {
    size_t n = coefficients.size();
    size_t m = parameters.size();
    matrix J(n, m);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            J(i, j) = expand(coefficients[i].diff(ex_to<symbol>(parameters[j])));
        }
    }
    return J;
}

/**
 * Helper function to compute the coefficient map for the catenary model
 * @param input the user input containing the number of compartments and leak locations
 * @param subsets the set of nonempty subsets of [n-1] without consecutive numbers
 * @param subsets_plus the set of nonempty subsets of [n-1] without consecutive numbers + 1
 * @param full_kappa the full kappa expression for the model
 * @param i the index for which we compute the sum
 * @param is_tilde flag to indicate if we are computing for tilde coefficients
 * @return the computed sum as an expression
 */
ex helper_sum(UserInput& input, vector<vector<int>>& subsets, vector<vector<int>>& subsets_plus, ex full_kappa, int i, bool is_tilde = false) {
    ex sum = 0;
    int iterator = 0;
    if (!is_tilde) {
        for (const auto& subset : subsets) {
            int length = subset.size();
            int e_index = i - (2 * length);
            if (length * 2 <= i) {
                vector<ex> set_for_esp = out_Ml(input, compartment_set(input.all_compartments, subsets_plus[iterator]));
                sum += kappa(input, subset) * elementary_symmetric_polynomials(set_for_esp, e_index)[e_index];
            }
            iterator++;
        }
    }
    else {
        for (const auto& subset : subsets) {
            int length = subset.size();
            int e_index = i - input.d - (2 * length);
            vector<int> intersection = intersection_of_sets(subset, input.P);
            if (length * 2 <= i && intersection.empty()) {
                vector<ex> set_for_esp = out_Ml(input, compartment_set(input.all_compartments, union_sets(subsets_plus[iterator], input.P)));
                if (e_index < 0 || e_index > set_for_esp.size())
                    continue;
                sum += kappa(input, subset) * elementary_symmetric_polynomials(set_for_esp, e_index)[e_index];
            }
            iterator++;
        }
    }
    return sum;
}

/**
 * Computes the coefficient map for the catenary model
 * @param input the user input containing the number of compartments and leak locations
 * @param subsets the set of nonempty subsets of [n-1] without consecutive numbers
 * @param subsets_plus the set of nonempty subsets of [n-1] without consecutive numbers + 1
 * @return a vector of expressions representing the coefficients
 */
vector<ex> compute_coefficient_map(UserInput& input, vector<vector<int>>& subsets, vector<vector<int>>& subsets_plus) {
    vector<ex> coefficients;
    vector<ex> full_esp_list = elementary_symmetric_polynomials(out_Ml(input, input.all_compartments), input.n);

    vector<int> tilde_set = compartment_set(input.all_compartments, input.P);
    vector<ex> tilde_esp_list = elementary_symmetric_polynomials(out_Ml(input, tilde_set), tilde_set.size());
    
    ex full_kappa = kappa_ij(input, input.in, input.p);
    ex a_i;
    ex a_i_tilde;
    ex curr_tilde_esp;

    // First loop computes the regular coefficients a_i
    for (int i = 1; i <= input.n; i++) {
        a_i = full_esp_list[i] - helper_sum(input, subsets, subsets_plus, full_kappa, i);
        coefficients.push_back(a_i);
    }
    
    // Second loop computes the tilda coefficients a_i tilde
    for (int i = input.d; i <= input.n; i++) {
        if (i - input.d < 0 || i - input.d >= tilde_esp_list.size())
            curr_tilde_esp = 0;
        else
            curr_tilde_esp = tilde_esp_list[i - input.d];

        a_i_tilde = full_kappa * curr_tilde_esp - helper_sum(input, subsets, subsets_plus, full_kappa, i, true);
        coefficients.push_back(a_i_tilde);
    }
        
    return coefficients;
}

/**
 * Function to handle user input for the model
 * @return a UserInput struct containing the user-defined parameters
 */
UserInput user_input() {
    UserInput input;

    // User input for number of compartments
    cout << "Enter the number of compartments: ";
    cin >> input.n;

    // User input for input node
    cout << "Enter the input: ";
    cin >> input.in;

    // User input for output node
    cout << "Enter the output: ";
    cin >> input.p;

    // User input for number of leaks
    int num_leaks;
    cout << "Enter the number of leaks: ";
    cin >> num_leaks;

    input.d = input.p - input.in;

    // Store leak compartments
    for (int i = 0; i < num_leaks; ++i) {
        int leak_compartment;
        cout << "Enter leak compartment: ";
        cin >> leak_compartment;

        // If compartment number is invalid, prompt again
        if (leak_compartment < 1 || leak_compartment > input.n) {
            cout << "Invalid compartment number. Please enter a number between 1 and " << input.n << "." << endl;
            i--;
            continue;
        }

        // If leak already exists, prompt again
        if (find(input.leak_compartments.begin(), input.leak_compartments.end(), leak_compartment) != input.leak_compartments.end()) {
            cout << "Leak already exists. Please enter a different leak compartment." << endl;
            i--;
            continue;
        }

        input.leak_compartments.insert(leak_compartment);
    }

    // Define the set of all compartments from 1 to n
    for (int i = 1; i <= input.n; i++) {
        input.all_compartments.push_back(i);
    }

    // Define the set P as all vertices from input to output inclusive
    for (int i = input.in; i <= input.p; i++)
        input.P.push_back(i);

    create_symbols(input);
    return input;
}

/**
 * Converts a matrix to a Mathematica-compatible string representation
 * @param M the matrix to convert
 * @return a string representing the matrix in Mathematica format
 */
string matrixToMathematica(const matrix& M) {
    ostringstream out;
    out << "{";
    for (int i = 0; i < M.rows(); ++i) {
        out << "{";
        for (int j = 0; j < M.cols(); ++j) {
            out << M(i, j);
            if (j < M.cols() - 1) out << ", ";
        }
        out << "}";
        if (i < M.rows() - 1) out << ",\n";
    }
    out << "}";
    return out.str();
}

int main() {
    UserInput input = user_input();
    vector<ex> all_params;
    vector<vector<int>> subsets = gamma(input);
    vector<vector<int>> subsets_plus = gamma_plus(subsets);
    vector<ex> coefficients = compute_coefficient_map(input, subsets, subsets_plus);
    auto start = chrono::high_resolution_clock::now();

    cout << "\nWelcome to the Catenary Model Identifier!" << endl;
    cout << "This program computes the Jacobian matrix of the coefficients with respect to the parameters." << endl;
    cout << "You can use Mathematica to compute the rank and determinant of the models n >= 5, but there are built in functions for n < 5 if you want to use them." << endl;
    
    matrix J = compute_jacobian(coefficients, input.parameters);
    
    cout << "\nJacobian Matrix (∂coefficients/∂parameters):" << endl;
    cout << "----------------------------------------" << endl;
    for (size_t i = 0; i < J.rows(); ++i) {
        cout << "Row " << i+1 << ": [";
        for (size_t j = 0; j < J.cols(); ++j) {
            cout << J(i,j);
            if (j < J.cols()-1) cout << ", ";
        }
        cout << "]" << endl;
    }

    // FOR TOMORROW: CODE A FUNCTION TO DO NUMERICAL RANK AND CODE A FUNCTION TO FIND SUBMATRICES IF NOT SQUARE. OTHERWISE JUST TAKE DET
    // First, check coefficients in n = 4 case, see when that failsafe in the code actually runs and if it causes issues
    // New hunch, the singular locus in the case of no leaks contains the forward parameters of all compartments and the backwards parameters of the last compartment
    cout << "-----------------------------\n";
    cout << "Jacobian Matrix in Mathematica format:\n";
    cout << matrixToMathematica(J) << endl;
    
    if (input.n < 5) {
        cout << "----------------------------------------" << endl;
        if (J.rows() == J.cols()) {
            ex determinant = J.determinant();
            determinant = collect_common_factors(determinant);
            cout << "Determinant of Jacobian Matrix: " << determinant << endl;
        }
        cout << "Rank of Jacobian Matrix: " << J.rank() << endl;
        cout << "----------------------------------------" << endl;
    }

    // End timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration in milliseconds
    chrono::duration<double, milli> duration = end - start;

    // Output time taken
    cout << "-----------------------------\n";
    cout << "Computation Time: " << duration.count() << " ms" << endl;
    cout << "-----------------------------\n";
    cout << "Thank you for using the Catenary Model Identifier!" << endl;
    return 0;
}    
