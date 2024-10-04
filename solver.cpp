#include "solver.hpp"
#include "game_logic.hpp"
#include <algorithm>
#include <iostream>

std::set<LinearEquation> Solver::generate_linear_equations(std::vector<std::vector<Cell>> &board, int mine_count) {
    std::set<LinearEquation> equations;
    // TODO generate the minecount contstraint as well

    // Iterate over all cells to create linear equations
    for (int r = 0; r < board.size(); ++r) {
        for (int c = 0; c < board[r].size(); ++c) {
            auto cell = board[r][c];
            /* right now we even generate eqn's on zeros as sometimes they can help:
             * 1100
             * 211#
             * F11#
             */
            if (cell.is_revealed && !cell.is_flagged && cell.adjacent_mines >= 1) {
                LinearEquation eq;
                eq.target_sum = board[r][c].adjacent_mines;

                // Gather variables (indices of adjacent cells)
                for (int dr = -1; dr <= 1; ++dr) {
                    for (int dc = -1; dc <= 1; ++dc) {
                        int nr = r + dr;
                        int nc = c + dc;
                        if (nr >= 0 && nr < board.size() && nc >= 0 && nc < board[r].size() && (dr != 0 || dc != 0)) {
                            if (board[nr][nc].is_flagged) {
                                eq.target_sum -= 1;
                            } else if (!board[nr][nc].is_revealed) {
                                eq.variables.push_back(nr * board[0].size() + nc);
                            }
                        }
                    }
                }

                eq.normalize();
                equations.insert(eq);
            }
        }
    }

    return equations;
}

void Solver::update_board(std::vector<std::vector<Cell>> &board, const std::unordered_map<int, int> &deduced_vars) {
    for (const auto &[var, value] : deduced_vars) {

        int row = var / board[0].size();
        int col = var % board[0].size();

        if (value == 1) {
            board[row][col].is_flagged = true;
        } else if (value == 0) {
            reveal_cell(board, row, col);
            /*board[row][col].is_revealed = true;*/
        }
    }
}

bool Solver::is_board_solved(std::vector<std::vector<Cell>> &board) {
    // Check if all cells are revealed or flagged correctly
    for (const auto &row : board) {
        for (const auto &cell : row) {
            if ((cell.is_mine && !cell.is_flagged) || (!cell.is_mine && !cell.is_revealed)) {
                return false;
            }
        }
    }
    return true;
}

void island_marker(std::vector<std::vector<Cell>> &board, int i, int j) {
    // Check boundaries and whether the cell is valid (not a mine and not visited)
    if (i < 0 || i >= board.size() || j < 0 || j >= board[0].size() || board[i][j].is_mine ||
        board[i][j].adjacent_mines != 0) {
        return;
    }

    board[i][j].adjacent_mines = -1;

    island_marker(board, i + 1, j);
    island_marker(board, i - 1, j);
    island_marker(board, i, j + 1);
    island_marker(board, i, j - 1);
}

std::vector<std::pair<int, int>> find_zero_islands(std::vector<std::vector<Cell>> board_copy) {
    std::vector<std::pair<int, int>> islands;
    int rows = board_copy.size();
    if (rows == 0)
        return islands;
    int cols = board_copy[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (not board_copy[i][j].is_mine and board_copy[i][j].adjacent_mines == 0) {
                islands.emplace_back(i, j);      // Store the index
                island_marker(board_copy, i, j); // Mark the whole island
            }
        }
    }
    return islands;
}

std::optional<std::pair<int, int>> Solver::solve(std::vector<std::vector<Cell>> &board, int mine_count) {

    std::vector<std::pair<int, int>> cell_list = find_zero_islands(board);

    std::cout << "Number of islands: " << cell_list.size() << std::endl;
    std::cout << "Cells in each zero island:" << std::endl;
    for (const auto &index : cell_list) {
        std::cout << "Index inside an island: (" << index.first << ", " << index.second << ")" << std::endl;
    }

    print_board_x_ray(board);

    // Randomly shuffle the list of cells
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(cell_list.begin(), cell_list.end(), g);

    int cell_idx = 0;

    bool found_ng_solvable_position = false;
    auto original_board = board;
    auto work_board = board;

    int row, col;

    while (not found_ng_solvable_position) {

        if (cell_idx == cell_list.size()) {
            std::cout << "tried all cells, nothing was found" << std::endl;
            return std::nullopt;
        }

        // reset the board and try again
        work_board = original_board;

        bool found_safe_initial_cell = false;

        while (not found_safe_initial_cell and cell_idx < cell_list.size()) {
            // TODO about to iterate through cell list till we find one, if we can't
            // stop, also if we ge to the end of cell list we stop as well.
            row = cell_list[cell_idx].first;
            col = cell_list[cell_idx].second;

            auto cell = work_board[row][col];

            // temporarily force the starting cell to be zero
            // this has to be true to be ngs in most cases because when you open a
            // non-zero number it's usually the case that nothing else opens?
            if (not cell.is_mine && cell.adjacent_mines == 0) {
                found_safe_initial_cell = true;
            } else {
                cell_idx += 1;
            }
        }

        if (cell_idx == cell_list.size()) {
            std::cout << "tried all cells, nothing was found" << std::endl;
            return std::nullopt;
        }

        // Start with the selected cell
        reveal_cell(work_board, row, col);

        std::cout << "found a starting position = (" << row << ", " << col << "),  opening there" << std::endl;

        print_board(work_board);

        bool cant_make_progress = false;

        bool logging = false;

        while (not cant_make_progress) {
            auto equations = generate_linear_equations(work_board, mine_count);
            auto augmented_matrix = create_augmented_matrix(equations, work_board.size() * work_board[0].size());
            auto expanded_augmented_matrix = equal_sum_expansion(augmented_matrix);
            expanded_augmented_matrix = zero_sum_expansion(expanded_augmented_matrix);
            /*auto zeroed_out_matrix = zero_out_columns_based_on_rightmost(augmented_matrix);*/
            /*auto equal_sum_expanded_matrix = equal_sum_expansion(augmented_matrix);*/
            /*auto equal_sum_expanded_matrix = augmented_matrix;*/

            if (logging) {
                std::cout << "augmented matrix" << std::endl;
                print_matrix(augmented_matrix);
                std::cout << "expanded matrix" << std::endl;
                print_matrix(expanded_augmented_matrix);
            }

            gaussian_elimination(expanded_augmented_matrix, work_board.size() * work_board[0].size());

            if (logging) {
                std::cout << "after gaussian elimination" << std::endl;
                print_matrix(expanded_augmented_matrix);
            }

            auto deduced_vars =
                deduce_variables(expanded_augmented_matrix, work_board.size() * work_board[0].size(), logging);

            if (deduced_vars.empty()) {
                std::cerr << "No variables could be deduced. Trying a new cell..." << std::endl;
                cant_make_progress = true;
                cell_idx += 1;
            } else {
                update_board(work_board, deduced_vars);
                std::cout << "just updated board with information" << std::endl;
                print_board(work_board);
            }
        }
        // if we can't make any more progress then the we are either not ngsolvable
        // and we go stuck or the board is completely solved
        found_ng_solvable_position = is_board_solved(work_board);
    }

    std::cout << "ng solvable from (row, col) = (" << row << ", " << col << ")" << std::endl;

    return std::pair(row, col);
}

void print_board(std::vector<std::vector<Cell>> &board) {
    for (const auto &row : board) {
        for (const auto &cell : row) {
            if (cell.is_revealed) {
                if (cell.is_mine) {
                    std::cout << 'M';
                } else {
                    std::cout << cell.adjacent_mines;
                }
            } else if (cell.is_flagged) {
                std::cout << '@';
            } else {
                std::cout << '#';
            }
            std::cout << ' ';
        }
        std::cout << std::endl;
    }
}

void print_board_x_ray(std::vector<std::vector<Cell>> &board) {
    for (const auto &row : board) {
        for (const auto &cell : row) {
            if (cell.is_mine) {
                std::cout << 'M';
            } else {
                std::cout << cell.adjacent_mines;
            }
            std::cout << ' ';
        }
        std::cout << std::endl;
    }
}
