#include "linear_system_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

void LinearEquation::normalize() { sort(variables.begin(), variables.end()); }

bool LinearEquation::operator<(const LinearEquation &other) const {
    if (variables != other.variables) {
        return variables < other.variables;
    }
    return target_sum < other.target_sum;
}

std::vector<std::vector<int>> create_augmented_matrix(const std::set<LinearEquation> &equations, int num_variables) {
    int num_equations = equations.size();
    std::vector<std::vector<int>> augmented_matrix(num_equations, std::vector<int>(num_variables + 1, 0));

    int row = 0;
    for (const auto &equation : equations) {
        for (int var : equation.variables) {
            augmented_matrix[row][var] = 1; // Set the coefficient for the variable
        }
        augmented_matrix[row][num_variables] = equation.target_sum; // Set the target sum in the augmented column
        ++row;
    }

    return augmented_matrix;
}

std::vector<std::vector<int>> zero_sum_expansion(const std::vector<std::vector<int>> &matrix) {
    std::vector<std::vector<int>> expanded_matrix;

    // Iterate over each row in the input matrix
    for (const auto &row : matrix) {
        int rightmost_value = row.back(); // Get the rightmost value of the row
        int count_of_ones = 0;

        // Count the number of 1's in the row (ignoring the rightmost value)
        for (size_t i = 0; i < row.size() - 1; ++i) {
            if (row[i] == 1) {
                count_of_ones++;
            }
        }

        // If the count of 1's matches the rightmost value, we split the row
        if (count_of_ones == rightmost_value && rightmost_value >= 1) {
            std::vector<size_t> one_positions;

            // Find the positions of all 1's in the row
            for (size_t i = 0; i < row.size() - 1; ++i) {
                if (row[i] == 1) {
                    one_positions.push_back(i);
                }
            }

            // Create a new row for each position of 1
            for (size_t pos : one_positions) {
                std::vector<int> new_row(row.size(), 0); // Initialize with zeros
                new_row[pos] = 1;                        // Place a 1 at the correct position
                new_row.back() = 1;                      // Set the rightmost value to 1
                expanded_matrix.push_back(new_row);
            }
        } else if (count_of_ones > 0 && rightmost_value == 0) {
            // If the rightmost value is 0, expand the row and set rightmost value to 0
            std::vector<size_t> one_positions;

            // Find the positions of all 1's in the row
            for (size_t i = 0; i < row.size() - 1; ++i) {
                if (row[i] == 1) {
                    one_positions.push_back(i);
                }
            }

            // Create a new row for each position of 1
            for (size_t pos : one_positions) {
                std::vector<int> new_row(row.size(), 0); // Initialize with zeros
                new_row[pos] = 1;                        // Place a 1 at the correct position
                new_row.back() = 0;                      // Set the rightmost value to 0
                expanded_matrix.push_back(new_row);
            }
        } else {
            // If no split needed, just add the original row
            expanded_matrix.push_back(row);
        }
    }

    return expanded_matrix;
}

std::vector<std::vector<int>> equal_sum_expansion(const std::vector<std::vector<int>> &matrix) {
    std::vector<std::vector<int>> expanded_matrix;

    // Iterate over each row in the input matrix
    for (const auto &row : matrix) {
        int rightmost_value = row.back(); // Get the rightmost value of the row
        int count_of_ones = 0;

        // Count the number of 1's in the row (ignoring the rightmost value)
        for (size_t i = 0; i < row.size() - 1; ++i) {
            if (row[i] == 1) {
                count_of_ones++;
            }
        }

        // If the count of 1's matches the rightmost value, we split the row
        if (count_of_ones == rightmost_value and rightmost_value >= 1) {
            std::vector<size_t> one_positions;

            // Find the positions of all 1's in the row
            for (size_t i = 0; i < row.size() - 1; ++i) {
                if (row[i] == 1) {
                    one_positions.push_back(i);
                }
            }

            // Create a new row for each position of 1
            for (size_t pos : one_positions) {
                std::vector<int> new_row(row.size(), 0); // Initialize with zeros
                new_row[pos] = 1;                        // Place a 1 at the correct position
                new_row.back() = 1;                      // Set the rightmost value to 1
                expanded_matrix.push_back(new_row);
            }
        } else {
            // If no split needed, just add the original row
            expanded_matrix.push_back(row);
        }
    }

    return expanded_matrix;
}

void gaussian_elimination(std::vector<std::vector<int>> &augmented_matrix, int num_variables, bool enable_logging) {
    // Convert the integer matrix to a floating-point matrix
    std::vector<std::vector<float>> A(augmented_matrix.size(), std::vector<float>(augmented_matrix[0].size()));

    for (size_t i = 0; i < augmented_matrix.size(); ++i) {
        for (size_t j = 0; j < augmented_matrix[i].size(); ++j) {
            A[i][j] = static_cast<float>(augmented_matrix[i][j]);
        }
    }

    int row_count = A.size();       // number of rows
    int column_count = A[0].size(); // number of columns
    int lead = 0;                   // This will track the leading column for the pivot element

    for (int r = 0; r < row_count; ++r) {
        // Find the pivot column
        while (lead < column_count && A[r][lead] == 0) {
            bool found_non_zero = false;
            for (int i = r; i < row_count; ++i) {
                if (A[i][lead] != 0) {
                    std::swap(A[i], A[r]);
                    found_non_zero = true;
                    break;
                }
            }
            if (!found_non_zero) {
                lead++;
            }
        }

        if (lead >= column_count)
            break;

        // Normalize the pivot row
        float lead_value = A[r][lead];
        if (lead_value != 0) {
            for (int c = 0; c < column_count; ++c) {
                A[r][c] /= lead_value;
            }
        }

        // Ensure all rows above and below this row have 0 in the current column
        for (int i = 0; i < row_count; ++i) {
            if (i != r) {
                float mult = A[i][lead];
                for (int c = 0; c < column_count; ++c) {
                    A[i][c] -= mult * A[r][c];
                }
            }
        }

        lead++;

        if (enable_logging) {
            std::cout << "After processing row " << r + 1 << ":" << std::endl;
            print_matrix(A);
        }
    }

    // Copy the results back to the original matrix, truncating the float to int
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < column_count; ++j) {
            augmented_matrix[i][j] = static_cast<int>(A[i][j]); // Truncate the floating-point value
        }
    }
}

void print_vector(const std::vector<int> &vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << " ";
        }
    }
}

void brute_force_deduction(const std::vector<int> &row, int constant,
                           std::unordered_map<int, int> &var_name_to_deduced_value, bool enable_logging) {

    if (enable_logging) {
        std::cout << " === starting brute force === " << std::endl;
    }

    std::vector<int> non_zero_indices;
    int num_variables_ignoring_last_entry = row.size() - 1;

    // modify the constant by considering already deduced variable values
    for (int j = 0; j < num_variables_ignoring_last_entry; j++) {
        if (row[j] != 0 && var_name_to_deduced_value.find(j) != var_name_to_deduced_value.end()) {
            // subtract the contribution of the deduced variable
            constant -= row[j] * var_name_to_deduced_value[j];
        }
    }
    if (enable_logging) {
        std::cout << "modified target sum: " << constant << std::endl;
    }

    // Identify variables with non-zero coefficients that haven't been deduced yet
    for (int j = 0; j < num_variables_ignoring_last_entry; j++) {
        if (row[j] != 0 && var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end()) {
            non_zero_indices.push_back(j);
        }
    }

    if (enable_logging) {
        std::cout << "varibles to be deduced with non-zero coef: ";
        print_vector(non_zero_indices);
        std::cout << std::endl;
    }

    int num_non_zero_vars = non_zero_indices.size();

    // Brute-force search over all combinations of values for non-zero variables
    int num_combinations = std::pow(2, num_non_zero_vars); // 2^num_non_zero_vars possibilities
    std::unordered_map<int, int> unique_solution;          // Store the unique solution if found
    bool solution_found = false;

    if (enable_logging) {
        std::cout << "about to start trying combinations of the values {0, 1} for the mentioned variables" << std::endl;
    }

    for (int combination = 0; combination < num_combinations; combination++) {
        std::unordered_map<int, int> variable_values; // Store current combination of variable values
        int sum = 0;

        // Assign binary values based on the current combination
        for (int i = 0; i < num_non_zero_vars; i++) {
            int var_index = non_zero_indices[i];
            int var_value = (combination >> i) & 1; // Get the i-th bit of combination (0 or 1)
            variable_values[var_index] = var_value;

            if (enable_logging) {
                std::cout << "v" << var_index << " = " << var_value << ", ";
            }

            // Compute the contribution to the sum
            sum += row[var_index] * var_value;
        }

        if (enable_logging) {
            std::cout << "sum of vars: " << sum << " target sum is " << constant << std::endl;
        }

        // Check if the current combination satisfies the modified equation
        if (sum == constant) {
            if (!solution_found) {
                // If this is the first solution, store it as the unique solution candidate
                unique_solution = variable_values;
                solution_found = true;
            } else {
                // if a second solution is found, we know the solution is not unique, so exit
                if (enable_logging) {
                    std::cout << "Multiple valid solutions found, no unique solution (brute force failed)."
                              << std::endl;
                }
                return;
            }
        }
    }

    // if only one solution is found, deduce the variables
    if (solution_found) {
        if (enable_logging) {
            std::cout << "Unique solution found: ";
        }

        // Assign the found solution to var_name_to_deduced_value
        for (const auto &pair : unique_solution) {
            bool value_not_yet_deduced = var_name_to_deduced_value.find(pair.first) == var_name_to_deduced_value.end();
            if (value_not_yet_deduced) {
                var_name_to_deduced_value[pair.first] = pair.second; // Deduce the value
                if (enable_logging) {
                    std::cout << "v" << pair.first << " = " << pair.second << " ";
                }
            }
        }
        if (enable_logging) {
            std::cout << std::endl;
        }
    } else {
        // If no solution is found
        if (enable_logging) {
            std::cout << "No solution found." << std::endl;
        }
    }
}

std::unordered_map<int, int> deduce_variables(const std::vector<std::vector<int>> &augmented_matrix, int num_variables,
                                              bool enable_logging) {

    std::unordered_map<int, int> var_name_to_deduced_value;

    int num_equations = augmented_matrix.size();

    // Logging to indicate the start of back substitution if enabled.
    if (enable_logging) {
        std::cout << "Starting back substitution:" << std::endl;
    }

    // Iterate over each equation starting from the last one (back substitution).
    for (int i = num_equations - 1; i >= 0; i--) {

        auto row = augmented_matrix[i];

        // Start with the sum on the right-hand side of the equation.
        int sum = row[num_variables];
        if (enable_logging) {
            std::cout << "working on " << i << "th row of augmented matrix" << std::endl;
            print_vector(row);
            std::cout << " sum = " << sum << std::endl;
        }

        // To track which variable might be deduced in this row.
        int variable_index = -1;

        // To count the number of non-zero coefficients in the current equation.
        int non_zero_coeff_count = 0;

        if (enable_logging) {
            std::cout << "analysing variables in row" << std::endl;
        }
        for (int j = 0; j < num_variables; j++) {
            // If the coefficient is non-zero.
            int coefficient = row[j];
            bool coefficient_is_nonzero = coefficient != 0;
            if (coefficient_is_nonzero) {
                non_zero_coeff_count++;
                if (enable_logging) {
                    std::cout << "Non-zero coefficient found at index " << j << ": " << coefficient << std::endl;
                }

                // check if the current variable has not been deduced yet.
                bool variable_not_yet_deduced = var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end();
                if (variable_not_yet_deduced) {
                    variable_index = j; // Possible variable to deduce.
                    if (enable_logging) {
                        std::cout << "Variable v" << j << " will attempt to be deduced" << std::endl;
                    }
                } else {
                    auto deduced_value = var_name_to_deduced_value[j];
                    if (deduced_value != 0) {
                        // Subtracting from left and right side of the equation
                        sum -= coefficient * deduced_value;
                        if (enable_logging) {
                            std::cout << "Variable v's value is known, subtracting " << coefficient * deduced_value
                                      << " from sum, new sum: " << sum << std::endl;
                        }
                    }
                }
            }
        }

        // removed because coeffs not all non-negative directly
        /*// Because coefficients are non-negative a zero sum means all the rest are*/
        /*// zero*/
        /*if (sum == 0 && non_zero_coeff_count >= 1) {*/
        /*    for (int j = 0; j < num_variables; j++) {*/
        /*        bool is_unassigned = var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end();*/
        /*        if (row[j] != 0 && is_unassigned) {*/
        /*            var_name_to_deduced_value[j] = 0;*/
        /*            if (enable_logging) {*/
        /*                std::cout << "Deduced variable v" << j << " = 0 by zero sum" << std::endl;*/
        /*            }*/
        /*        }*/
        /*    }*/
        /*}*/

        bool deduced_something = false;

        // Special case: if the sum equals the number of non-zero coefficients, all
        // undetermined variables must be 1.
        if (sum == non_zero_coeff_count && non_zero_coeff_count >= 1) {
            deduced_something = true;
            for (int j = 0; j < num_variables; j++) {
                bool value_not_yet_deduced = var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end();
                if (row[j] != 0 && value_not_yet_deduced) {
                    var_name_to_deduced_value[j] = 1;
                    if (enable_logging) {
                        std::cout << "Deduced variable v" << j << " = 1 by equal sum" << std::endl;
                    }
                }
            }
        }

        // If there's exactly one undetermined variable in this equation.
        if (non_zero_coeff_count == 1 && variable_index != -1) {
            deduced_something = true;
            // calculate the value of the undetermined variable by dividing the sum by
            // the coefficient.
            int deduced_value = sum / augmented_matrix[i][variable_index];
            if (enable_logging) {
                std::cout << "Deduced value for variable v" << variable_index << " = " << deduced_value
                          << " by row reduction isolation" << std::endl;
            }

            if (deduced_value == 0 || deduced_value == 1) {
                // Store the deduced value in the map.
                var_name_to_deduced_value[variable_index] = deduced_value;
                if (enable_logging) {
                    std::cout << "Deduced value for variable v" << variable_index << " = " << deduced_value
                              << " by row reduction isolation" << std::endl;
                }
            } else if (enable_logging) {
                std::cout << "Deduced value for variable v" << variable_index
                          << " is not binary, it has value: " << deduced_value << " skipping" << std::endl;
            }
        }

        brute_force_deduction(row, sum, var_name_to_deduced_value, enable_logging);
    }

    // Return the map containing the deduced variable values.
    return var_name_to_deduced_value;
}
