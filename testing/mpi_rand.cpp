#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <getopt.h>
#include <assert.h>

uint8_t get_random_value(const int row_id, const int col_id, const int n,
                         const int seed)
{
  uint8_t r = 0;
  int my_seed = seed + row_id * n + col_id;
  // printf("my_seed: %d\n", my_seed);
  srand(my_seed);
  // rand();
  // printf("rand=%d\n", rand());
  rand();
  r = rand() % 2;
  return r;
}

int main(int argc, char *argv[])
{
  int rank, size; // always

  // random number generator stuff
  int seed_offset = 10;
  int n = 10; // default values
  int opt;

  int n_loc_r, n_loc_c;
  int nprows, npcols;
  int prow_idx, pcol_idx;

  // Cartesian communicator stuff
  int dims[2] = {0, 0};
  int pers[2] = {0, 0};
  int coords[2];
  MPI_Comm cartcomm;

  MPI_Init(&argc, &argv);

  static struct option long_options[] = {
      {"number", required_argument, 0, 'n'},
      {"key", required_argument, 0, 'k'},
      {0, 0, 0, 0}};

  while (1)
  {
    int option_index = 0;
    opt = getopt_long(argc, argv, "n:k:", long_options, &option_index);

    if (opt == -1)
      break;

    switch (opt)
    {
    case 'n':
      n = atoi(optarg);
      break;
    case 'k':
      seed_offset = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s -n <n> -k <k>\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  if (rank == 0)
  {
    printf("n: %d\n", n);
    printf("k: %d\n", seed_offset);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Dims_create(size, 2, dims);

  if (rank == 0)
  {
    printf("dims: [%d, %d]\n", dims[0], dims[1]);
  }
  nprows = dims[0];
  npcols = dims[1];

  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pers, reorder, &cartcomm);

  n_loc_r = n / nprows;
  n_loc_c = n / npcols;

  prow_idx = rank / npcols;
  pcol_idx = rank % npcols;

  MPI_Cart_coords(cartcomm, rank, 2, coords);
  assert(prow_idx == coords[0]);
  assert(pcol_idx == coords[1]);

  if (rank == 0)
  {
    printf("n_loc_r: %d n_loc_c: %d\n", n_loc_r, n_loc_c);
  }

  printf("%d: prow_idx: %d pcol_idx: %d\n", rank, prow_idx, pcol_idx);

  // matrix with npcols; this declaration is enough to find matrix[i][j]
  // u_int8_t(*matrix)[n_loc_c];

  // matrix = (uint8_t(*)[n_loc_c])malloc(n_loc_r * n_loc_c * sizeof(u_int8_t));

  uint8_t **matrix = new uint8_t *[n_loc_r];
  for (int i = 0; i < n_loc_r; ++i)
  {
    matrix[i] = new uint8_t[n_loc_c];
  }

  int m_offset_r = prow_idx * n_loc_r;
  int m_offset_c = pcol_idx * n_loc_c;
  printf("%d: m_offset_r: %d m_offset_c: %d\n", rank, m_offset_r, m_offset_c);

  for (int i = 0; i < n_loc_r; i++)
  {
    for (int j = 0; j < n_loc_c; j++)
    {
      matrix[i][j] =
          get_random_value(m_offset_r + i, m_offset_c + j, n, seed_offset);
    }
  }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < size; i++)
  {
    if (rank == i)
    {
      printf("%d: matrix\n", rank);
      for (int i = 0; i < n_loc_r; i++)
      {
        for (int j = 0; j < n_loc_c; j++)
        {
          printf("%d ", matrix[i][j]);
        }
        printf("\n");
      }
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Comm_free(&cartcomm);

  free(matrix);

  MPI_Finalize();

  return 0;
}
