/*  YOUR_FIRST_NAME
 *  YOUR_LAST_NAME
 *  YOUR_UBIT_NAME
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <mpi.h>
#include <unistd.h>


int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    // ...
    int rank;
    int size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    int col = rank % q;
    int row = rank / q;
    std::vector<signed char> P; //declaring parent sub-matrix
    int b = n / q;
    // printf("size is %d \n",size);
    // printf("q is %d \n",q);
    P.resize(b * b, 0);
    sleep(rank);
    // sequential processing here for each processor
    for (int i = 0; i<b; i++)
    {
    	for (int j = 0; j<b; j++)
    	{
    		if (A[i * b + j] == 1){
    			P[i * b + j] = row*b+i;//assigning the highest row-wise vertex
    			if(P[0 + j]<=P[i * b + j]){
    				P[0 + j] = P[i * b + j];
    			}
     		}
    			
    		else
    			P[i * b + j] = 0;
    	}
    }
    for (int i = 0; i<b; i++)//replicating locally found max vertex in each process
    {
    	for (int j = 0; j<b; j++)
    	{
    		P[i * b + j] = P[0 + j];
    	}
    }

    //Printing the Parent Matrix
    std::cout << "Parent Matrix P (" << row << "," << col << ")" << std::endl;
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < b; ++j) std::cout << static_cast<int>(P[i * b + j]) << " ";
        std::cout << std::endl;
    }
    MPI_Comm col_comm;
	MPI_Comm_split(MPI_COMM_WORLD, rank%q, rank, &col_comm);
	int col_rank, col_size;
	MPI_Comm_rank(col_comm, &col_rank);
	MPI_Comm_size(col_comm, &col_size);
	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",rank, size, col_rank, col_size);

 //    if (rank == 0)
	// {	
	// 	//do allreduce along columnwise, that is col = rank%q
	// 	printf("With n = %d trapezoids, our estimate\n", n);
	// 	printf("of the integral from %f to %f = %f\n",
	// 	a, b, total);
	// }

    return -1;
} // connected_components

#endif // A1_HPP
