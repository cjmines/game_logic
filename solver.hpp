#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "linear_system_solver.hpp"
#include "game_logic.hpp"
#include <optional>
#include <vector>
#include <set>
#include <unordered_map>
#include <random>

void print_board(std::vector<std::vector<Cell>> &board);

class Solver {
public:

    /**
     * @brief Starts the solving process.
     * 
     * @return True if the board is no-guess solvable, false otherwise.
     */
    std::optional<std::pair<int, int>> solve(std::vector<std::vector<Cell>>& board, int mine_count);


private:

    /**
     * @brief Generates linear equations for the current board state.
     * 
     * @return Set of linear equations representing the system.
     */

    std::set<LinearEquation> generate_linear_equations(std::vector<std::vector<Cell>> &board, int mine_count);

    /**
     * @brief Updates the board based on the deduced variable values.
     * 
     * @param deduced_vars Map of variable indices to their deduced values.
     */
    void update_board(std::vector<std::vector<Cell>> &board, const std::unordered_map<int, int> &deduced_vars);

    /**
     * @brief Checks if the current state of the board is solved.
     * 
     * @return True if the board is solved, false otherwise.
     */
    bool is_board_solved(std::vector<std::vector<Cell>> &board);

};

#endif // SOLVER_HPP
