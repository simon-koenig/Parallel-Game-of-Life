#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <chrono>

// Function to clear the console screen
void clearScreen() {
    std::cout << "\033[2J\033[H";
}

// Function to initialize the grid with random values
void initializeGrid(std::vector<std::vector<int>>& grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            grid[i][j] = rand() % 2;
        }
    }
}

// Function to print the current state of the grid
void printGrid(const std::vector<std::vector<int>>& grid) {
    clearScreen(); // Clear the console screen
    for (const auto& row : grid) {
        for (int cell : row) {
            std::cout << (cell == 1 ? "■" : "□") << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Function to get the number of live neighbors for a given cell
int getLiveNeighborCount(const std::vector<std::vector<int>>& grid, int row, int col) {
    int liveNeighbors = 0;
    int rows = grid.size();
    int cols = grid[0].size();

    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            if (i == 0 && j == 0) {
                continue;
            }

            int neighborRow = row + i;
            int neighborCol = col + j;

            if (neighborRow >= 0 && neighborRow < rows && neighborCol >= 0 && neighborCol < cols) {
                liveNeighbors += grid[neighborRow][neighborCol];
            }
        }
    }

    return liveNeighbors;
}

// Function to update the grid based on the rules of the Game of Life
void evolveGrid(std::vector<std::vector<int>>& grid) {
    int rows = grid.size();
    int cols = grid[0].size();
    std::vector<std::vector<int>> newGrid(rows, std::vector<int>(cols, 0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int liveNeighbors = getLiveNeighborCount(grid, i, j);

            // Apply rules
            if (grid[i][j] == 1) {
                if (liveNeighbors < 2 || liveNeighbors > 3) {
                    newGrid[i][j] = 0;
                } else {
                    newGrid[i][j] = 1;
                }
            } else {
                if (liveNeighbors == 3) {
                    newGrid[i][j] = 1;
                }
            }
        }
    }

    grid = newGrid;
}

// Function to run the Game of Life for a specified number of generations
void gameOfLife(int rows, int cols, int generations) {
    std::vector<std::vector<int>> grid(rows, std::vector<int>(cols, 0));
    initializeGrid(grid, rows, cols);

    for (int generation = 0; generation < generations; ++generation) {
        printGrid(grid);
        evolveGrid(grid);
        std::this_thread::sleep_for(std::chrono::milliseconds(500));  // Add a delay to make the evolution visible
    }
}

int main() {
    std::srand(static_cast<unsigned>(std::time(nullptr))); // Seed for random number generation

    // Example: Run the Game of Life with a 10x10 grid for 20 generations
    gameOfLife(15, 15, 40);

    return 0;
}
