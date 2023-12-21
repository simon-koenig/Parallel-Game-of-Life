#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

void printGrid(const std::vector<int>& grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << (grid[i * cols + j] == 1 ? "■" : "□") << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void updateGrid(std::vector<int>& grid, int rows, int cols) {
    std::vector<int> new_grid(grid.size(), 0);

    for (int i = 0; i < rows; ++i) {

        for (int j = 0; j < cols; ++j) {
            // periodic BC
            int i_prev = (i - 1 + rows) % rows;
            int i_next = (i + 1) % rows;
            int j_prev = (j - 1 + cols) % cols;
            int j_next = (j + 1) % cols;

            int alive_neighbors = grid[i_prev * cols + j_prev]
                              + grid[i_prev * cols + j]
                              + grid[i_prev * cols + j_next]
                              + grid[i_next * cols + j_prev] 
                              + grid[i_next * cols + j] 
                              + grid[i_next * cols + j_next]
                              + grid[i * cols + j_prev] 
                              + grid[i * cols + j_next];

            if (grid[i * cols + j] == 1) {
                new_grid[i * cols + j] = (alive_neighbors == 2 || alive_neighbors == 3) ? 1 : 0;
            } else{
                new_grid[i * cols + j] = (alive_neighbors == 3) ? 1 : 0;
            }
        }
    }
    // update grid
    std::copy(new_grid.begin(), new_grid.end(), grid.begin());
}

void run_sequential(std::vector<int>& grid, int rows, int cols, int iterations){
    //printGrid(grid, rows, cols);
    for (int i = 0; i <= iterations; ++i){
        updateGrid(grid, rows, cols);
    }
    printGrid(grid, rows, cols);
}
