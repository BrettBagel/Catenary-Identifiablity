// C++ adapation of identifying catenary linear and directed-cycle compartmental models computationally
// Author: Brett Hajdaj
// Date: 2025 06 20

#include "/Applications/Wolfram.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions/wstp.h"
#include <ginac/ginac.h>
#include <unordered_set>
#include <stdexcept>
#include <cln/cln.h>
#include <iostream>
#include <sstream>
#include <gmpxx.h>
#include <cstring>
#include <numeric>
#include <fstream>
#include <vector>
#include <chrono>
#include <thread>
#include <string>
#include <ctime>
#include <regex>
#include <map>

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

// Improved MathematicaInterface class with better error handling and debugging

class MathematicaInterface {
private:
    WSENV env;
    WSLINK rank_link;
    WSLINK factor_link;

    void clearPendingPackets(WSLINK link) {
        if (!link) return;
        while (WSReady(link)) {
            int packet = WSNextPacket(link);
            if (packet == 0) break;
            WSNewPacket(link);
        }
    }
    
    WSLINK createKernel() {
        if (!env) {
            cerr << "Environment not initialized" << endl;
            return nullptr;
        }
        
        int error;
        vector<string> launch_commands = {
            "-linkname '/Applications/Wolfram.app/Contents/MacOS/MathKernel' -linkmode launch",
            "-linkname 'math' -linkmode launch",
            "-linkname '/Applications/Mathematica.app/Contents/MacOS/MathKernel' -linkmode launch"
        };
        
        WSLINK link = nullptr;
        for (const string& cmd : launch_commands) {
            cout << "Attempting to launch kernel with: " << cmd << endl;
            link = WSOpenString(env, cmd.c_str(), &error);
            if (link) {
                cout << "Successfully created link" << endl;
                break;
            } else {
                cout << "Failed to create link with error: " << error << endl;
            }
        }
        
        if (!link) {
            cerr << "Failed to launch Mathematica kernel with any method" << endl;
            return nullptr;
        }
        
        cout << "Attempting to connect to kernel..." << endl;
        if (!WSConnect(link)) {
            cerr << "Failed to connect to Mathematica kernel" << endl;
            cerr << "Error: " << WSError(link) << endl;
            cerr << "Error message: " << WSErrorMessage(link) << endl;
            WSClose(link);
            return nullptr;
        }
        
        cout << "Connected to kernel, waiting for ready state..." << endl;
        
        int timeout = 600; // 6 seconds total
        int checks = 0;
        while (timeout > 0 && !WSReady(link)) {
            this_thread::sleep_for(chrono::milliseconds(100));
            timeout--;
            checks++;
            if (checks % 10 == 0) {
                cout << "Still waiting for kernel... (" << (600-timeout)/10.0 << "s)" << endl;
                if (WSError(link) != WSEOK) {
                    cerr << "Connection error while waiting: " << WSError(link) << endl;
                    cerr << "Error message: " << WSErrorMessage(link) << endl;
                    WSClose(link);
                    return nullptr;
                }
            }
        }

        this_thread::sleep_for(chrono::seconds(60));
        
        if (!WSReady(link)) {
            cerr << "Mathematica kernel not ready after timeout" << endl;
            cerr << "Final error state: " << WSError(link) << endl;
            cerr << "Error message: " << WSErrorMessage(link) << endl;
            WSClose(link);
            return nullptr;
        }
        
        cout << "Kernel is ready!" << endl;
        return link;
    }
    
public:
    MathematicaInterface() : env(nullptr), rank_link(nullptr), factor_link(nullptr) {}
    
    bool initialize() {
        cout << "Initializing WSTP environment..." << endl;
        env = WSInitialize(nullptr);
        if (!env) {
            cerr << "Failed to initialize WSTP environment" << endl;
            return false;
        }
        cout << "WSTP environment initialized" << endl;
        
        cout << "Creating rank calculation kernel..." << endl;
        rank_link = createKernel();
        if (!rank_link) {
            cerr << "Failed to create rank kernel" << endl;
            return false;
        }
        
        cout << "Creating determinant calculation kernel..." << endl;
        factor_link = createKernel();
        if (!factor_link) {
            cerr << "Failed to create determinant kernel" << endl;
            return false;
        }
        
        cout << "Mathematica interface initialized successfully!" << endl;
        return true;
    }
    
    void cleanup() {
        cout << "Cleaning up Mathematica interface..." << endl;
        if (rank_link) {
            WSClose(rank_link);
            rank_link = nullptr;
        }
        if (factor_link) {
            WSClose(factor_link);
            factor_link = nullptr;
        }
        if (env) {
            WSDeinitialize(env);
            env = nullptr;
        }
        cout << "Cleanup complete" << endl;
    }
    
    ~MathematicaInterface() {
        cleanup();
    }
    
    bool isReady() const {
        return env != nullptr && rank_link != nullptr && factor_link != nullptr;
    }
    
    int calc_rank(const vector<vector<string>>& matrix) {
        if (!rank_link) {
            cerr << "Rank kernel not connected" << endl;
            return -1;
        }
        
        
        clearPendingPackets(rank_link);
        WSNewPacket(rank_link);
        
        auto start_time = chrono::steady_clock::now();
        const int timeout_seconds = 60;
        
        if (!WSPutFunction(rank_link, "MatrixRank", 1)) {
            cerr << "Failed to send MatrixRank function" << endl;
            return -1;
        }
        
        if (!WSPutFunction(rank_link, "List", matrix.size())) {
            cerr << "Failed to send outer List" << endl;
            return -1;
        }
        
        for (const auto& row : matrix) {
            if (!WSPutFunction(rank_link, "List", row.size())) {
                cerr << "Failed to send row List" << endl;
                return -1;
            }
            for (const auto& element : row) {
                if (!WSPutFunction(rank_link, "ToExpression", 1)) {
                    cerr << "Failed to send ToExpression" << endl;
                    return -1;
                }
                if (!WSPutString(rank_link, element.c_str())) {
                    cerr << "Failed to send element: " << element << endl;
                    return -1;
                }
            }
        }
        
        if (!WSEndPacket(rank_link)) {
            cerr << "Failed to end packet" << endl;
            return -1;
        }
        
        int packet_type;
        int attempts = 0;
        const int max_attempts = 100;
        
        while (attempts < max_attempts) {
            auto current_time = chrono::steady_clock::now();
            if (chrono::duration_cast<chrono::seconds>(current_time - start_time).count() > timeout_seconds) {
                cerr << "Timeout waiting for rank result" << endl;
                WSNewPacket(rank_link);
                return -1;
            }
            
            packet_type = WSNextPacket(rank_link);
            
            if (packet_type == 0) {
                cerr << "Error in communication - no packet received" << endl;
                return -1;
            }
            
            if (packet_type == RETURNPKT) {
                break;
            }
            
            switch (packet_type) {
                case TEXTPKT:
                    {
                        const char* text;
                        if (WSGetString(rank_link, &text)) {
                            cout << "Kernel text: " << text << endl;
                            WSReleaseString(rank_link, text);
                        }
                    }
                    break;
                    
                case MESSAGEPKT:
                    {
                        const char* msg;
                        if (WSGetString(rank_link, &msg)) {
                            string msg_str(msg);
                            cout << "Kernel message: " << msg << endl;
                            if (msg_str.find("error") != string::npos || msg_str.find("Error") != string::npos) {
                                cerr << "Mathematica error: " << msg << endl;
                            }
                            WSReleaseString(rank_link, msg);
                        }
                    }
                    break;
                    
                case INPUTNAMEPKT:
                    {
                        const char* name;
                        if (WSGetString(rank_link, &name)) {
                            WSReleaseString(rank_link, name);
                        }
                    }
                    break;
                    
                default:
                    cout << "Received packet type: " << packet_type << endl;
                    break;
            }
            
            WSNewPacket(rank_link);
            attempts++;
        }
        
        if (attempts >= max_attempts) {
            cerr << "Max attempts reached waiting for rank result" << endl;
            WSNewPacket(rank_link);
            return -1;
        }
        
        int rank;
        if (!WSGetInteger(rank_link, &rank)) {
            cerr << "Failed to get rank result" << endl;
            WSNewPacket(rank_link);
            return -1;
        }
        
        WSNewPacket(rank_link);
        return rank;
    }

    /**
     * Factors a GiNaC expression using Mathematica's Factor function
     * @param expr The GiNaC expression to factor
     * @return The factored expression as a string in InputForm
     * @throws runtime_error if communication with Mathematica fails
     */
    string factor_expression(const ex& expr) {
        if (!factor_link) {
            throw runtime_error("Factor kernel not connected");
        }
        
        clearPendingPackets(factor_link);
        WSNewPacket(factor_link);
        
        auto start_time = chrono::steady_clock::now();
        const int timeout_seconds = 60;
        
        ostringstream expr_stream;
        expr_stream << expr;
        string expr_str = expr_stream.str();
        
        if (!WSPutFunction(factor_link, "ToString", 2)) {
            throw runtime_error("Failed to send ToString function");
        }

        if (!WSPutFunction(factor_link, "Simplify", 1)) {
            throw runtime_error("Failed to send Simplify function");
        }
        
        if (!WSPutFunction(factor_link, "Factor", 1)) {
            throw runtime_error("Failed to send Factor function");
        }
        
        if (!WSPutFunction(factor_link, "ToExpression", 1)) {
            throw runtime_error("Failed to send ToExpression function");
        }
        
        if (!WSPutString(factor_link, expr_str.c_str())) {
            throw runtime_error("Failed to send expression string: " + expr_str);
        }
        
        if (!WSPutSymbol(factor_link, "InputForm")) {
            throw runtime_error("Failed to send InputForm symbol");
        }
        
        if (!WSEndPacket(factor_link)) {
            throw runtime_error("Failed to end packet");
        }
        
        int packet_type;
        int attempts = 0;
        const int max_attempts = 100;
        
        while (attempts < max_attempts) {
            auto current_time = chrono::steady_clock::now();
            if (chrono::duration_cast<chrono::seconds>(current_time - start_time).count() > timeout_seconds) {
                WSNewPacket(factor_link);
                throw runtime_error("Timeout waiting for factor result");
            }
            
            packet_type = WSNextPacket(factor_link);
            
            if (packet_type == 0) {
                throw runtime_error("Error in communication - no packet received during factoring");
            }
            
            if (packet_type == RETURNPKT) {
                break;
            }
            
            if (packet_type == MESSAGEPKT) {
                const char* msg;
                if (WSGetString(factor_link, &msg)) {
                    string msg_str(msg);
                    if (msg_str.find("error") != string::npos || msg_str.find("Error") != string::npos) {
                        string error_msg = "Mathematica error during factoring: " + msg_str;
                        WSReleaseString(factor_link, msg);
                        WSNewPacket(factor_link);
                        throw runtime_error(error_msg);
                    }
                    WSReleaseString(factor_link, msg);
                }
            }
            
            WSNewPacket(factor_link);
            attempts++;
        }
        
        if (attempts >= max_attempts) {
            WSNewPacket(factor_link);
            throw runtime_error("Max attempts reached waiting for factor result");
        }
        
        const char* result_cstr;
        if (!WSGetString(factor_link, &result_cstr)) {
            WSNewPacket(factor_link);
            throw runtime_error("Failed to get factored result string");
        }
        
        string result(result_cstr);
        WSReleaseString(factor_link, result_cstr);
        WSNewPacket(factor_link);
        
        if (!result.empty() && result.front() == '"' && result.back() == '"') {
            result = result.substr(1, result.length() - 2);
        }

        if (result[0] == ',') {
            result[0] = ' ';
        }
        
        return result;
    }
};

/**
 * Converts a Mathematica output string to a GiNaC expression
 * Handles basic operations and power notation from Mathematica's InputForm
 * @param mathematica_str The string from Mathematica (should be in InputForm)
 * @return ex The corresponding GiNaC expression
 * @throws runtime_error if parsing fails
 */
ex mathematica_string_to_ginac(const string& mathematica_str) {
    static map<string, symbol> symbol_cache;
    
    auto get_symbol = [&](const string& name) -> symbol {
        auto it = symbol_cache.find(name);
        if (it != symbol_cache.end()) {
            return it->second;
        }
        symbol s(name);
        symbol_cache[name] = s;
        return s;
    };
    
    try {
        string processed = mathematica_str;
        
        processed = regex_replace(processed, regex(R"(\s+)"), " ");
        processed = regex_replace(processed, regex(R"(^\s+|\s+$)"), ""); // trim
        
        processed = regex_replace(processed, regex(R"((\w+)\^(\w+))"), "pow($1,$2)");
        processed = regex_replace(processed, regex(R"((\w+)\^(\d+))"), "pow($1,$2)");
        processed = regex_replace(processed, regex(R"((\d+)\^(\w+))"), "pow($1,$2)");
        
        size_t pos = 0;
        while ((pos = processed.find('^', pos)) != string::npos) {
            int base_end = pos - 1;
            int base_start = base_end;
            
            if (processed[base_end] == ')') {
                int paren_count = 1;
                base_start--;
                while (base_start >= 0 && paren_count > 0) {
                    if (processed[base_start] == ')') paren_count++;
                    else if (processed[base_start] == '(') paren_count--;
                    base_start--;
                }
                base_start++;
            } else {
                while (base_start >= 0 && (isalnum(processed[base_start]) || processed[base_start] == '_')) {
                    base_start--;
                }
                base_start++;
            }
            
            int exp_start = pos + 1;
            int exp_end = exp_start;
            
            if (processed[exp_start] == '(') {
                int paren_count = 1;
                exp_end++;
                while (exp_end < processed.length() && paren_count > 0) {
                    if (processed[exp_end] == '(') paren_count++;
                    else if (processed[exp_end] == ')') paren_count--;
                    exp_end++;
                }
                exp_end--;
            } else {
                while (exp_end < processed.length() && (isalnum(processed[exp_end]) || processed[exp_end] == '_')) {
                    exp_end++;
                }
                exp_end--;
            }
            
            string base = processed.substr(base_start, base_end - base_start + 1);
            string exponent = processed.substr(exp_start, exp_end - exp_start + 1);
            string replacement = "pow(" + base + "," + exponent + ")";
            
            processed.replace(base_start, exp_end - base_start + 1, replacement);
            pos = base_start + replacement.length();
        }

        processed = regex_replace(processed, regex(R"((\w+|\))\s+(\w+|\())"), "$1*$2");
        processed = regex_replace(processed, regex(R"((\d)\s*([a-zA-Z]))"), "$1*$2");
        
        processed = regex_replace(processed, regex(R"((\d+)/(\d+))"), "($1)/($2)");
        
        regex var_regex(R"(\b[a-zA-Z][a-zA-Z0-9]*\b)");
        sregex_iterator vars_begin(processed.begin(), processed.end(), var_regex);
        sregex_iterator vars_end;
        
        symtab table;
        for (sregex_iterator i = vars_begin; i != vars_end; ++i) {
            string var_name = i->str();
            if (var_name != "pow" && var_name != "sin" && var_name != "cos" && 
                var_name != "exp" && var_name != "log") {
                table[var_name] = get_symbol(var_name);
            }
        }

        parser reader(table);
        ex result = reader(processed);
        
        return result;
        
    } catch (const exception& e) {
        throw runtime_error("Failed to parse Mathematica expression: " + mathematica_str + 
                          " (Error: " + string(e.what()) + ")");
    }
}

/**
 * Creates symbols for the edge parameters and leak parameters based on user input
 * @param input the user input containing the number of compartments and leak locations
 */
void create_symbols(UserInput& input) {
    for (int i = 0; i < input.n - 1; i++) {
        string forward = "k" + to_string(i + 2) + to_string(i + 1);
        string backward = "k" + to_string(i + 1) + to_string(i + 2);
        
        input.parameters.push_back(input.get_symbol(forward));
        input.parameters.push_back(input.get_symbol(backward));
    }

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
            ex sum = 0;
            if (is_leak)
                sum += input.get_symbol("k01");
            outgoing_sums.push_back(sum);
            continue;
        }
    
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
    if (index > n) {
        if (!current_subset.empty())
            gamma_set.push_back(current_subset);
        return;
    }

    gamma_helper(n, index + 1, current_subset, gamma_set);

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
            temp_subset.push_back(elem + 1);
        }
        temp_subset.insert(temp_subset.end(), subset.begin(), subset.end());
        sort(temp_subset.begin(), temp_subset.end());
        
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
    if (k == 0) 
        return {ex(1)};

    vector<ex> e(k + 1, 0);
    e[0] = 1;

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
                vector<int> exclusion_set = compartment_set(input.all_compartments, subsets_plus[iterator]);
                vector<ex> set_for_esp = out_Ml(input, exclusion_set);
                sum += kappa(input, subset) * elementary_symmetric_polynomials(set_for_esp, e_index)[e_index];
            }
            iterator++;
        }
    }
    else {
        int iterator = 0;
        for (const auto& subset : subsets) {
            int length = subset.size();
            int e_index = i - input.d - (2 * length);
            vector<int> intersection = intersection_of_sets(subset, input.P);
            if (length * 2 <= i && intersection.empty()) {
                vector<int> unionize = union_sets(subsets_plus[iterator], input.P);
                vector<int> exclusion_set = compartment_set(input.all_compartments, unionize);
                vector<ex> set_for_esp = out_Ml(input, exclusion_set);
                if (e_index < 0 || e_index > set_for_esp.size()) {
                    iterator++;
                    continue;
                }
                sum += kappa(input, subset) * elementary_symmetric_polynomials(set_for_esp, e_index)[e_index];
            }
            iterator += 1;
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

    for (int i = 1; i <= input.n; i++) {
        a_i = full_esp_list[i] - helper_sum(input, subsets, subsets_plus, full_kappa, i);
        coefficients.push_back(a_i);
    }
    
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

    cout << "Enter the number of compartments: ";
    cin >> input.n;

    cout << "Enter the input: ";
    cin >> input.in;

    cout << "Enter the output: ";
    cin >> input.p;

    int num_leaks;
    cout << "Enter the number of leaks: ";
    cin >> num_leaks;

    input.d = input.p - input.in;

    for (int i = 0; i < num_leaks; ++i) {
        int leak_compartment;
        cout << "Enter leak compartment: ";
        cin >> leak_compartment;

        if (leak_compartment < 1 || leak_compartment > input.n) {
            cout << "Invalid compartment number. Please enter a number between 1 and " << input.n << "." << endl;
            i--;
            continue;
        }

        if (find(input.leak_compartments.begin(), input.leak_compartments.end(), leak_compartment) != input.leak_compartments.end()) {
            cout << "Leak already exists. Please enter a different leak compartment." << endl;
            i--;
            continue;
        }

        input.leak_compartments.insert(leak_compartment);
    }

    for (int i = 1; i <= input.n; i++) {
        input.all_compartments.push_back(i);
    }

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

vector<vector<string>> matrix_to_str_arr(const matrix& m) {
    vector<vector<string>> result;
    result.assign(m.rows(), vector<string>(m.cols(), ""));
    for (size_t i = 0; i < m.rows(); i++) {
        for (size_t j = 0; j < m.cols(); j++) {
            ostringstream oss;
            oss << m(i, j);
            result[i][j] = oss.str();
        }
    }
    return result;
}


/**
 * Computes the GCD of all r×r submatrix determinants of a non-square Jacobian matrix
 * @param J the Jacobian matrix
 * @return the symbolic GCD of all r×r submatrix determinants
 */
ex gcd_of_maximal_minors(const matrix& J, MathematicaInterface& math) {
    int r = math.calc_rank(matrix_to_str_arr(J));
    int rows = J.rows();
    int cols = J.cols();

    if (r == 0)
        return 0;

    vector<ex> minors;

    vector<vector<int>> row_indices;
    vector<vector<int>> col_indices;

    auto generate_combinations = [](int n, int k) {
        vector<vector<int>> result;
        vector<int> comb(k);
        iota(comb.begin(), comb.end(), 0);
        while (true) {
            result.push_back(comb);
            int i = k - 1;
            while (i >= 0 && comb[i] == n - k + i) --i;
            if (i < 0) break;
            ++comb[i];
            for (int j = i + 1; j < k; ++j)
                comb[j] = comb[j - 1] + 1;
        }
        return result;
    };

    row_indices = generate_combinations(rows, r);
    col_indices = generate_combinations(cols, r);
    for (const auto& row_set : row_indices) {
        for (const auto& col_set : col_indices) {
            matrix submat(r, r);
            for (int i = 0; i < r; ++i) {
                for (int j = 0; j < r; ++j) {
                    submat(i, j) = J(row_set[i], col_set[j]);
                }
            }
            ex det = mathematica_string_to_ginac(math.factor_expression(submat.determinant()));
            if (!det.is_zero()) {
                minors.push_back(det);
            }
        }
    }

    if (minors.empty())
        return 0;

    ex result = minors[0];
    for (size_t i = 1; i < minors.size(); ++i) {
        result = gcd(result, minors[i]);
    }
    return mathematica_string_to_ginac(math.factor_expression(result));
}

/**
 * @brief Generates all valid combinations of k elements from a set of n elements
 * @param n 
 * @param k
 * @return 2d array of combinations
 */
vector<vector<int>> generate_combinations(int n, int k) {
    vector<vector<int>> result;
    vector<int> comb(k);
    iota(comb.begin(), comb.end(), 1);
    
    while (true) {
        result.push_back(comb);
        int i = k - 1;
        while (i >= 0 && comb[i] == n - k + i + 1) --i;
        if (i < 0) break;
        ++comb[i];
        for (int j = i + 1; j < k; ++j)
            comb[j] = comb[j - 1] + 1;
    }
    return result;
}

/**
 * @brief Converts a array of leak compartments to a comma-separated string
 * @param leaks array of leaks
 */
string leaks_to_string(const vector<int>& leaks) {
    if (leaks.empty()) return "";
    
    string result = to_string(leaks[0]);
    for (size_t i = 1; i < leaks.size(); ++i) {
        result += "," + to_string(leaks[i]);
    }
    return result;
}

/**
 * @brief Returns just the factors of an expression
 * 
 * @param expr 
 * @return ex 
 */
ex extract_factors(const ex& expr) {
    if (!is_a<mul>(expr)) {
        if (is_a<add>(expr)) {
            return 1;
        }
        return expr;
    }
    
    ex result = 1;
    
    for (size_t i = 0; i < expr.nops(); i++) {
        ex factor = expr.op(i);
        
        if (!is_a<add>(factor)) {
            result *= factor;
        }
    }
    
    return result;
}

/**
 * @brief Runs all combinations for a given number of compartments and outputs 
 * to CSV
 * @param n number of nodes
 * @param csv_file spreadsheet out stream
 */
void run_all_combinations(int n, ofstream& csv_file, MathematicaInterface& math) {
    int count = 0;
    ex num_cases = GiNaC::pow(2, n - 1) * n * (n + 1);
    for (int num_leaks = 0; num_leaks <= n; ++num_leaks) {
        vector<vector<int>> leak_combinations;
        
        if (num_leaks == 0) {
            leak_combinations.push_back(vector<int>{});
        } else {
            leak_combinations = generate_combinations(n, num_leaks);
        }
        
        for (const auto& leak_combo : leak_combinations) {
            for (int input = 1; input <= n; ++input) {
                for (int output = input; output <= n; ++output) {
                    UserInput input_config;
                    input_config.n = n;
                    input_config.in = input;
                    input_config.p = output;
                    input_config.d = output - input;
                    
                    for (int leak : leak_combo) {
                        input_config.leak_compartments.insert(leak);
                    }                
                    
                    for (int i = 1; i <= n; i++) {
                        input_config.all_compartments.push_back(i);
                    }
                    for (int i = input; i <= output; i++) {
                        input_config.P.push_back(i);
                    }
                    
                    create_symbols(input_config);
                                        
                    vector<vector<int>> subsets = gamma(input_config);
                    vector<vector<int>> subsets_plus = gamma_plus(subsets);
                    vector<ex> coefficients = compute_coefficient_map(input_config, subsets, subsets_plus);
                    matrix J = compute_jacobian(coefficients, input_config.parameters);
                    
                    int rank = math.calc_rank(matrix_to_str_arr(J));
                    int max_rank = J.cols();
                    bool identifiable = (rank == max_rank);
                    bool is_square = J.rows() == J.cols();
                    
                    ex locus = 0;
                    string str_locus = "";
                    bool calc_locus = (n < 6);
                    if ((J.rows() == J.cols()) && calc_locus) {
                        locus = mathematica_string_to_ginac(math.factor_expression(J.determinant()));
                    } else if (identifiable && calc_locus) {
                        locus = gcd_of_maximal_minors(J, math);
                    } else if (calc_locus) {
                        locus = 0;
                    }

                    
                    csv_file << n << "," 
                             << num_leaks << "," 
                             << "\"" << leaks_to_string(leak_combo) << "\"," 
                             << input << "," 
                             << output << "," 
                             << (identifiable ? "True" : "False") << ","
                             << (is_square ? "True" : "False") << "," 
                             << rank << "/" << max_rank << ","
                             << locus << endl; 
                    count++;
                    cout << count << "/" << num_cases << endl;
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    bool batch_mode = false;
    MathematicaInterface math;
    if (!math.initialize()) {
        cerr << "Failed to initialize Mathematica interface" << endl;
        return 1;
    }
    if (argc > 1 && strcmp(argv[1], "-s") == 0) {
        batch_mode = true;
    }
    
    if (batch_mode) {
        int max_nodes = 0;
        cout << "Max graph size: " << endl;
        cin >> max_nodes;
        cout << "Running in batch mode - generating CSV for all combinations of n=" << max_nodes << "..." << endl;
        string filename = "test.csv";
        ofstream csv_file(filename);
        csv_file << "# compartments,# leaks,Leak locations,Input,Output,Identifiable,Square,Rank,Singular Locus\n";
        for (int i = 1; i <= max_nodes; i++) {
            run_all_combinations(i , csv_file, math);
        }
        csv_file.close();
        cout << "CSV file 'singular_loci.csv' has been generated." << endl;
        return 0;
    }
    
    UserInput input = user_input();
    vector<ex> all_params;
    vector<vector<int>> subsets = gamma(input);
    vector<vector<int>> subsets_plus = gamma_plus(subsets);
    vector<ex> coefficients = compute_coefficient_map(input, subsets, subsets_plus);
    ofstream J_output_file("jacobian.txt");
    bool is_full_rank = false;

    cout << "\nWelcome to the Catenary Model Identifier!" << endl;
    cout << "This program computes the Jacobian matrix of the coefficients with respect to the parameters." << endl;
    cout << "You can use Mathematica to compute the rank and determinant of the models n >= 5, but there are built in functions for n < 5 if you want to use them." << endl;
    
    matrix J = compute_jacobian(coefficients, input.parameters);

    cout << "\nCoefficients (a_i and a_i tilde):" << endl;
    cout << "----------------------------------------" << endl;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        cout << "Coefficient a_" << (i + 1) << ": " << coefficients[i] << endl;
    }
    cout << "----------------------------------------" << endl;
    
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
    J_output_file << matrixToMathematica(J) << flush;
    J_output_file.close();
    
    if (input.n < 8) {
        int wolf_rank = math.calc_rank(matrix_to_str_arr(J));
        cout << "Wolfram calculated rank of Jacobian Matrix: " << wolf_rank << endl;

        if (wolf_rank == J.cols()) {
            cout << "The model is identifiable & full rank; rank = " << wolf_rank << "." << endl;
            is_full_rank = true;
        } 
        else
            cout << "The model is not identifiable & not full rank; rank = " << wolf_rank << "." << endl;
        cout << "----------------------------------------" << endl;
        
        ex locus = 0;
        if (J.rows() == J.cols()) {
            ex determinant = J.determinant();
            determinant = collect_common_factors(determinant);
            locus = extract_factors(determinant);
            cout << "Determinant of square Jacobian Matrix: " << determinant << endl;
            cout << "Singular locus of square Jacobian Matrix: " << locus << endl;
        } 
        else if (is_full_rank){
            cout << "Jacobian Matrix is not square but is full rank, computing GCD of maximal minors instead." << endl;
            ex gcd_result = gcd_of_maximal_minors(J, math);
            cout << "GCD of maximal minors: " << gcd_result << endl;
        }
        else
            cout << "The determinant is 0.";

        
        cout << "----------------------------------------" << endl;
    }

    cout << "The user-specified parameters are as follows:\n";
    cout << "Number of Compartments: " << input.n << endl;
    cout << "Input Node: " << input.in << endl;
    cout << "Output Node: " << input.p << endl;
    cout << "Leak Compartments: ";
    for (const auto& leak : input.leak_compartments) {
        cout << leak << " ";
    }
    cout << endl;
    cout << "-----------------------------\n";
    cout << "Thank you for using the Catenary Model Identifier!" << endl;
    return 0;
}
