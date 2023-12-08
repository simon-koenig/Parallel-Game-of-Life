#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

// Function to initialize the local portion of the grid with random values (0 or 1)
void initializeLocalGrid(std::vector<std::vector<int>> &localGrid, int n_loc_r, int m_loc_c)
{
    for (int i = 0; i < n_loc_r; ++i)
    {
        for (int j = 0; j < m_loc_c; ++j)
        {
            localGrid[i][j] = rand() % 2; // Randomly initialize with 0 or 1
        }
    }
}

// Function to update the local portion of the grid based on the given condition
void updateLocalGrid(std::vector<std::vector<int>> &localGrid, int n_loc_r, int m_loc_c)
{
    std::vector<std::vector<int>> updatedLocalGrid = localGrid;

    for (int i = 0; i < n_loc_r; ++i)
    {
        for (int j = 0; j < m_loc_c; ++j)
        {
            int count = 0;

            // Check 8 neighbors
            for (int ni = i - 1; ni <= i + 1; ++ni)
            {
                for (int nj = j - 1; nj <= j + 1; ++nj)
                {
                    if (ni >= 0 && ni < n_loc_r && nj >= 0 && nj < m_loc_c && !(ni == i && nj == j))
                    {
                        count += localGrid[ni][nj];
                    }
                }
            }

            // Update the value based on the condition
            if (count >= 3)
            {
                updatedLocalGrid[i][j] = 1;
            }
            else
            {
                updatedLocalGrid[i][j] = 0;
            }
        }
    }

    // Update the original local grid with the new values
    localGrid = updatedLocalGrid;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the size of the global grid
    int n = 10; // Number of rows
    int m = 10; // Number of columns

    // Create a 2D vector to represent the global grid
    std::vector<std::vector<int>> globalGrid(n, std::vector<int>(m));

    // Create a Cartesian communicator

    MPI_Comm cartComm;
    int dims[2] = {0, 0};    // 2D grid
    int periods[2] = {0, 0}; // Enable periodic boundary conditions
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cartComm);

    // Get the coordinates and dimensions of the grid
    int coords[2];
    int dimsGrid[2];
    MPI_Cart_get(cartComm, 2, dimsGrid, periods, coords);

    // Calculate the local size of the grid for each process
    int n_loc_r = n / dimsGrid[0];
    int m_loc_c = m / dimsGrid[1];

    // Create a 2D vector to represent the local portion of the grid
    std::vector<std::vector<int>> localGrid(n_loc_r, std::vector<int>(m_loc_c));

    // Initialize the local grid with random values
    initializeLocalGrid(localGrid, n_loc_r, m_loc_c);

    // Print the initial state of the local grid for each process
    std::cout << "Process " << rank << " - Initial state:\n";
    for (int i = 0; i < n_loc_r; ++i)
    {
        for (int j = 0; j < m_loc_c; ++j)
        {
            std::cout << localGrid[i][j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';

    // Number of time steps
    int numTimeSteps = 3;

    // Perform time-stepping
    for (int step = 1; step <= numTimeSteps; ++step)
    {
        // Update the local grid based on the condition
        updateLocalGrid(localGrid, n_loc_r, m_loc_c);

        // Print the state of the local grid after each time step
        std::cout << "Process " << rank << " - State after time step " << step << ":\n";
        for (int i = 0; i < n_loc_r; ++i)
        {
            for (int j = 0; j < m_loc_c; ++j)
            {
                std::cout << localGrid[i][j] << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }

    MPI_Finalize();

    return 0;
}
