#include "linear_system_solver.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

using namespace std;

void LinearEquation::normalize() { sort(variables.begin(), variables.end()); }

bool LinearEquation::operator<(const LinearEquation &other) const {
    if (variables != other.variables) {
        return variables < other.variables;
    }
    return target_sum < other.target_sum;
}

vector<vector<int>> create_augmented_matrix(const set<LinearEquation> &equations, int num_variables) {
    int num_equations = equations.size();
    vector<vector<int>> augmented_matrix(num_equations, vector<int>(num_variables + 1, 0));

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

void gaussian_elimination(vector<vector<int>> &augmented_matrix, int num_variables, bool enable_logging) {
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

unordered_map<int, int> deduce_variables(const vector<vector<int>> &augmented_matrix, int num_variables,
                                         bool enable_logging) {

    unordered_map<int, int> var_name_to_deduced_value;

    int num_equations = augmented_matrix.size();

    // Logging to indicate the start of back substitution if enabled.
    if (enable_logging) {
        cout << "Starting back substitution:" << endl;
    }

    // Iterate over each equation starting from the last one (back substitution).
    for (int i = num_equations - 1; i >= 0; i--) {

        auto row = augmented_matrix[i];

        // Start with the sum on the right-hand side of the equation.
        int sum = row[num_variables];
        if (enable_logging) {
            cout << "working on " << i << "th row of augmented matrix" << endl;
            print_vector(row);
            cout << " sum = " << sum << endl;
        }

        // To track which variable might be deduced in this row.
        int variable_index = -1;

        // To count the number of non-zero coefficients in the current equation.
        int non_zero_coeff_count = 0;

        // Iterate through all the variables in the current equation.
        for (int j = 0; j < num_variables; j++) {
            // If the coefficient is non-zero.
            int coefficient = row[j];
            bool coefficient_is_nonzero = coefficient != 0;
            if (coefficient_is_nonzero) {
                non_zero_coeff_count++;
                if (enable_logging) {
                    cout << "Non-zero coefficient found at index " << j << ": " << coefficient << endl;
                }

                // Check if the current variable has not been deduced yet.
                bool not_yet_deduced = var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end();
                if (not_yet_deduced) {
                    variable_index = j; // Possible variable to deduce.
                    if (enable_logging) {
                        cout << "Variable v" << j << " might be deduced" << endl;
                    }
                } else {
                    auto deduced_value = var_name_to_deduced_value[j];
                    if (deduced_value != 0) {
                        // Subtracting from left and right side of the equation
                        sum -= coefficient * deduced_value;
                        if (enable_logging) {
                            cout << "Subtracting " << coefficient * deduced_value << " from sum, new sum: " << sum
                                 << endl;
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
        /*                cout << "Deduced variable v" << j << " = 0 by zero sum" << endl;*/
        /*            }*/
        /*        }*/
        /*    }*/
        /*}*/

        // Special case: if the sum equals the number of non-zero coefficients, all
        // undetermined variables must be 1.
        if (sum == non_zero_coeff_count && non_zero_coeff_count >= 1) {
            for (int j = 0; j < num_variables; j++) {
                if (row[j] != 0 && var_name_to_deduced_value.find(j) == var_name_to_deduced_value.end()) {
                    var_name_to_deduced_value[j] = 1;
                    if (enable_logging) {
                        cout << "Deduced variable v" << j << " = 1 by equal sum" << endl;
                    }
                }
            }
        }

        // If there's exactly one undetermined variable in this equation.
        if (non_zero_coeff_count == 1 && variable_index != -1) {
            // Calculate the value of the undetermined variable by dividing the sum by
            // the coefficient.
            int deduced_value = sum / augmented_matrix[i][variable_index];
            if (enable_logging) {
                cout << "Deduced value for variable v" << variable_index << " = " << deduced_value
                     << " by row reduction isolation" << endl;
            }

            // Ensure the deduced value is binary (0 or 1).
            if (deduced_value == 0 || deduced_value == 1) {
                // Store the deduced value in the map.
                var_name_to_deduced_value[variable_index] = deduced_value;
                if (enable_logging) {
                    cout << "Deduced value for variable v" << variable_index << " = " << deduced_value
                         << " by row reduction isolation" << endl;
                }
            } else if (enable_logging) {
                cout << "Deduced value for variable v" << variable_index
                     << " is not binary, it has value: " << deduced_value << " skipping" << endl;
            }
        }
    }

    // Return the map containing the deduced variable values.
    return var_name_to_deduced_value;
}
