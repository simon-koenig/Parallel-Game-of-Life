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

// Function to exchange border values between neighboring processes
void exchangeBorderValues(std::vector<std::vector<int>> &localGrid, int n_loc_r, int m_loc_c, MPI_Comm cartComm)
{
    MPI_Status status;

    int west, east, north, south;

    // Get ranks of neighboring processes
    MPI_Cart_shift(cartComm, 0, 1, &north, &south);
    MPI_Cart_shift(cartComm, 1, 1, &west, &east);

    // Send and receive top border values
    MPI_Sendrecv(&localGrid[0][0], m_loc_c, MPI_INT, north, 0,
                 &localGrid[n_loc_r][0], m_loc_c, MPI_INT, south, 0, cartComm, &status);

    // Send and receive bottom border values
    MPI_Sendrecv(&localGrid[n_loc_r - 1][0], m_loc_c, MPI_INT, south, 1,
                 &localGrid[-1][0], m_loc_c, MPI_INT, north, 1, cartComm, &status);

    // Send and receive west border values
    MPI_Sendrecv(&localGrid[0][0], n_loc_r, MPI_INT, west, 2,
                 &localGrid[0][m_loc_c], n_loc_r, MPI_INT, east, 2, cartComm, &status);

    // Send and receive east border values
    MPI_Sendrecv(&localGrid[0][m_loc_c - 1], n_loc_r, MPI_INT, east, 3,
                 &localGrid[0][-1], n_loc_r, MPI_INT, west, 3, cartComm, &status);
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
    MPI_Dims_create(size, 2, dims);

    if (rank == 0)
    {
        printf("dims: [%d, %d]\n", dims[0], dims[1]);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cartComm);

    // Get the coordinates and dimensions of the grid
    int coords[2];
    // int dimsGrid[2];
    //  MPI_Cart_get(cartComm, 2, dimsGrid, periods, coords);
    MPI_Cart_coords(cartComm, rank, 2, coords);

    // Calculate the local size of the grid for each process
    int n_loc_r = n / dims[0];
    int m_loc_c = m / dims[1];

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
        // Exhcange value at borders
        exchangeBorderValues(localGrid, n_loc_r, m_loc_c, cartComm);
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
