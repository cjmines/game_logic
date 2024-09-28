#ifndef LINEAR_SYSTEM_SOLVER_HPP
#define LINEAR_SYSTEM_SOLVER_HPP

#include <vector>
#include <set>
#include <unordered_map>

/**
 * @brief Represents a linear equation with binary variables.
 */
struct LinearEquation {
    std::vector<int> variables; ///< Collection of binary variables in the equation.
    int target_sum; ///< The target sum of the linear equation.

    /**
     * @brief Normalizes the equation by sorting the variables.
     */
    void normalize();

    /**
     * @brief Custom comparison operator to allow storing in a set.
     * 
     * @param other The other linear equation to compare with.
     * @return True if this equation is less than the other, false otherwise.
     */
    bool operator<(const LinearEquation& other) const;
};

/**
 * @brief Creates an augmented matrix from a set of linear equations.
 * 
 * @param equations The set of linear equations.
 * @param num_variables The number of binary variables.
 * @return The augmented matrix representing the system of equations.
 */
std::vector<std::vector<int>> create_augmented_matrix(
    const std::set<LinearEquation>& equations, int num_variables);

/**
 * @brief Performs Gaussian elimination on the augmented matrix.
 * 
 * @param augmented_matrix The augmented matrix to be reduced.
 * @param num_variables The number of binary variables.
 * @param enable_logging Boolean to enable or disable logging of the steps.
 */
void gaussian_elimination(std::vector<std::vector<int>>& augmented_matrix, int num_variables, bool enable_logging = false);

/**
 * @brief Deduces the values of binary variables from the reduced augmented matrix.
 * 
 * @param augmented_matrix The reduced augmented matrix.
 * @param num_variables The number of binary variables.
 * @param enable_logging Boolean to enable or disable logging of the back substitution steps.
 * @return A map of variable indices to their deduced values.
 */
std::unordered_map<int, int> deduce_variables(const std::vector<std::vector<int>>& augmented_matrix, int num_variables, bool enable_logging = false);

/**
 * @brief Prints the matrix to the console.
 * 
 * @param matrix The matrix to be printed.
 */
template <typename T>
void print_matrix(const std::vector<std::vector<T>>& matrix);

#include "linear_system_solver.tpp"


#endif // LINEAR_SYSTEM_SOLVER_HPP
