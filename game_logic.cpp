#include <cstdlib>
#include "game_logic.hpp"
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

// Function to read the Minesweeper board from a file
std::pair<Board, int> read_board_from_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    std::vector<std::vector<Cell>> board;
    std::string line;

    int mine_count = 0;

    while (std::getline(file, line)) {
        std::vector<Cell> row;
        std::istringstream iss(line);
        std::string value;

        while (iss >> value) {
            Cell cell;
            if (value == "M") {
                cell.is_mine = true;
                mine_count += 1;
            } else {
                cell.adjacent_mines = std::stoi(value);
            }
            row.push_back(cell);
        }
        board.push_back(row);
    }

    file.close();
    return {board, mine_count};
}

Board generate_board(int mines_count, int height, int width) {

    std::random_device rd;  // Random device to seed the engine
    std::mt19937 gen(rd()); // Mersenne Twister random number engine
    std::uniform_int_distribution<> distrib_row(0, height - 1);
    std::uniform_int_distribution<> distrib_col(0, width - 1);

    std::vector<std::vector<Cell>> board(height, std::vector<Cell>(width));
    int mines_placed = 0;

    while (mines_placed < mines_count) {
        int row = distrib_row(gen);
        int col = distrib_col(gen);
        if (!board[row][col].is_mine) {
            board[row][col].is_mine = true;
            mines_placed++;
        }
    }

    for (int row = 0; row < board.size(); row++) {
        for (int col = 0; col < board[0].size(); col++) {
            if (!board[row][col].is_mine) {
                int count = 0;
                for (int i = -1; i <= 1; i++) {
                    for (int j = -1; j <= 1; j++) {
                        int r = row + i;
                        int c = col + j;
                        if (r >= 0 && r < board.size() && c >= 0 && c < board[0].size() && board[r][c].is_mine) {
                            count++;
                        }
                    }
                }
                board[row][col].adjacent_mines = count;
            }
        }
    }
    return board;
}

bool reveal_cell(std::vector<std::vector<Cell>> &board, int row, int col) {
    if (row < 0 || row >= board.size() || col < 0 || col >= board.at(0).size() || board[row][col].is_revealed ||
        board[row][col].is_flagged) {
        return true;
    }

    board[row][col].is_revealed = true;

    if (board[row][col].is_mine) {
        return false;
    }

    int flagged_neighbors = 0;
    int adjacent_count = board[row][col].adjacent_mines;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            int r = row + i;
            int c = col + j;
            if (r >= 0 && r < board.size() && c >= 0 && c < board[0].size() && board[r][c].is_flagged) {
                flagged_neighbors++;
            }
        }
    }

    if (flagged_neighbors == adjacent_count) {
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                int r = row + i;
                int c = col + j;
                if (r >= 0 && r < board.size() && c >= 0 && c < board[0].size() && !board[r][c].is_flagged) {
                    reveal_cell(board, r, c);
                }
            }
        }
    }

    if (board[row][col].adjacent_mines == 0) {
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                reveal_cell(board, row + i, col + j);
            }
        }
    }

    return true;
}

bool reveal_adjacent_cells(std::vector<std::vector<Cell>> &board, int row, int col) {
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            int r = row + i;
            int c = col + j;
            if (!reveal_cell(board, r, c)) {
                return false;
            }
        }
    }
    return true;
}

void toggle_flag_cell(std::vector<std::vector<Cell>> &board, int row, int col) {
    if (row >= 0 && row < board.size() && col >= 0 && col < board[0].size() && !board[row][col].is_revealed) {
        board[row][col].is_flagged = !board[row][col].is_flagged;
    }
}

void flag_cell(std::vector<std::vector<Cell>> &board, int row, int col) {
    if (row >= 0 && row < board.size() && col >= 0 && col < board[0].size() && !board[row][col].is_revealed) {
        board[row][col].is_flagged = true;
    }
}

void flag_adjacent_cells(std::vector<std::vector<Cell>> &board, int row, int col) {
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            int r = row + i;
            int c = col + j;
            flag_cell(board, r, c);
        }
    }
}

bool field_clear(std::vector<std::vector<Cell>> &board) {
    int remaining_mines = 0;
    for (const auto &row : board) {
        for (const auto &cell : row) {
            bool mine_is_flagged = cell.is_mine && cell.is_flagged;
            bool non_mine_is_revealed = !cell.is_mine && cell.is_revealed;
            bool cell_correct = mine_is_flagged or non_mine_is_revealed;
            if (not cell_correct) {
                return false;
            }
        }
    }
    // all cells are correct
    return true;
}
