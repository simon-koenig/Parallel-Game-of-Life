
mpi_rand: mpi_rand.c
	mpic++  -Wall  -o $@ $<

clean:
	rm -f mpi_rand	
