#include <cstdlib>
#include "game_logic.hpp"
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include "stb_image.h"

Board construct_board_from_mine_bools(const std::vector<std::vector<bool>> &mines) {
    Board board(mines.size(), std::vector<Cell>(mines[0].size()));

    // Populate the board based on mine positions
    for (size_t row = 0; row < mines.size(); ++row) {
        for (size_t col = 0; col < mines[0].size(); ++col) {
            board[row][col].is_mine = mines[row][col];
        }
    }

    // Now compute the adjacent mine counts for each cell
    int rows = board.size();
    int cols = board[0].size();
    std::vector<std::pair<int, int>> directions = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1},
                                                   {0, 1},   {1, -1}, {1, 0},  {1, 1}};

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            if (!board[row][col].is_mine) {
                int adjacent_mines = 0;
                for (const auto &[dx, dy] : directions) {
                    int new_row = row + dx;
                    int new_col = col + dy;
                    if (new_row >= 0 && new_row < rows && new_col >= 0 && new_col < cols) {
                        if (board[new_row][new_col].is_mine) {
                            adjacent_mines++;
                        }
                    }
                }
                board[row][col].adjacent_mines = adjacent_mines;
            }
        }
    }

    return board;
}

// Function to read the board from a text file
std::pair<Board, int> read_board_from_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    std::vector<std::vector<bool>> mines;
    std::string line;
    int mine_count = 0;

    // Read the board layout from the file
    while (std::getline(file, line)) {
        std::vector<bool> row;
        for (char c : line) {
            // Skip spaces
            if (c == ' ') {
                continue;
            }

            bool is_mine = (c == 'M');
            row.push_back(is_mine);
            if (is_mine) {
                mine_count += 1;
            }
        }

        // Only add non-empty rows to the mines matrix
        if (!row.empty()) {
            mines.push_back(row);
        }
    }

    file.close();
    // Construct the board from the mines matrix
    Board board = construct_board_from_mine_bools(mines);
    return {board, mine_count};
}

// Function to read the board from a PNG image file
std::pair<Board, int> read_board_from_image_file(const std::string &filename) {
    // Load the image
    int width, height, channels;
    unsigned char *image = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    if (!image) {
        throw std::runtime_error("Failed to load image");
    }

    std::vector<std::vector<bool>> mines(height, std::vector<bool>(width));
    int mine_count = 0;
    int black_threshold = 50;

    // Process the image to determine mine positions
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Assuming the image is in RGB format
            int pixel_index = (y * width + x) * channels;
            unsigned char r = image[pixel_index];
            unsigned char g = image[pixel_index + 1];
            unsigned char b = image[pixel_index + 2];

            if (r == 0) {
                int grabme = 0;
            }

            // Check if the pixel is black (mine) or white (no mine)
            if (r < black_threshold && g < black_threshold && b < black_threshold) {
                mines[y][x] = true; // It's a mine
                mine_count++;
            } else {
                mines[y][x] = false; // No mine
            }
        }
    }

    stbi_image_free(image);
    // Construct the board from the mines matrix
    Board board = construct_board_from_mine_bools(mines);
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

void set_cell_flag(std::vector<std::vector<Cell>> &board, int row, int col, bool flag_value) {
    if (row >= 0 && row < board.size() && col >= 0 && col < board[0].size() && !board[row][col].is_revealed) {
        board[row][col].is_flagged = flag_value;
    }
}

void set_adjacent_cells_flags(std::vector<std::vector<Cell>> &board, int row, int col, bool flag_value) {
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            int r = row + i;
            int c = col + j;
            set_cell_flag(board, r, c, flag_value);
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
