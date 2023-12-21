#include <iostream>
#include <vector>
#include <algorithm>

// Function to print the current state of the grid
void printGrid(const std::vector<int>& grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << (grid[i * cols + j] == 1 ? '*' : ' ') << ' ';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

// Function to determine the next state of the grid based on Game of Life rules
void updateGrid(std::vector<int>& grid, int rows, int cols) {
    std::vector<int> newGrid(grid.size(), 0);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Calculate indices of neighboring cells with periodic boundary conditions
            int iPrev = (i - 1 + rows) % rows;
            int iNext = (i + 1) % rows;
            int jPrev = (j - 1 + cols) % cols;
            int jNext = (j + 1) % cols;

            // Count live neighbors
            int liveNeighbors = grid[iPrev * cols + jPrev] + grid[iPrev * cols + j] + grid[iPrev * cols + jNext]
                              + grid[i * cols + jPrev] + grid[i * cols + jNext]
                              + grid[iNext * cols + jPrev] + grid[iNext * cols + j] + grid[iNext * cols + jNext];

            // Apply Game of Life rules
            if (grid[i * cols + j] == 1) {
                newGrid[i * cols + j] = (liveNeighbors == 2 || liveNeighbors == 3) ? 1 : 0;
            } else {
                newGrid[i * cols + j] = (liveNeighbors == 3) ? 1 : 0;
            }
        }
    }

    // Update the original grid
    std::copy(newGrid.begin(), newGrid.end(), grid.begin());
}

int main() {
    const int rows = 5;
    const int cols = 5;

    // Initial configuration (example)
    std::vector<int> grid = {
        0, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        1, 1, 1, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 
    };

    // Number of generations to simulate
    const int numGenerations = 5;

    // Simulate the Game of Life for a few generations
    for (int generation = 0; generation < numGenerations; ++generation) {
        std::cout << "Generation " << generation + 1 << ":\n";
        printGrid(grid, rows, cols);
        updateGrid(grid, rows, cols);
    }

    return 0;
}
