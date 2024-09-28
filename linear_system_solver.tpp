#include "linear_system_solver.hpp"
#include <vector>
#include <iostream>
#include <iomanip> // For std::setprecision

template <typename T>
void print_matrix(const std::vector<std::vector<T>>& matrix) {
    for (const auto& row : matrix) {
        for (const T& val : row) {
            if constexpr (std::is_floating_point_v<T>) {
                std::cout << std::fixed << std::setprecision(0) << val << " ";
            } else {
                std::cout << val << " ";
            }
        }
        std::cout << std::endl;
    }
}
