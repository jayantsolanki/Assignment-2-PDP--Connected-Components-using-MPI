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
    std::vector<signed char> Ptemp; //declaring temporary parent sub-matrix
    std::vector<signed char> M; //declaring temp sub-matrix, helper matrix
    std::vector<signed char> Mtemp; //declaring temp sub-matrix
    std::vector<signed char> Q; //declaring temp sub-matrix
    std::vector<signed char> P; //declaring parent sub-matrix
    std::vector<signed char> PPrime; //declaring parent sub-matrix
    int b = n / q;
    int max = 0;
    // printf("size is %d \n",size);
    // printf("q is %d \n",q);
    Ptemp.resize(b * b, 0);
    P.resize(b * b, 0);
    Q.resize(b * b, 0);
    M.resize(b * b, 0);
    Mtemp.resize(b * b, 0);
    PPrime.resize(b * b, 0);
    // sequential processing here for each processor
    for (int i = 0; i<b; i++)
    {
    	for (int j = 0; j<b; j++)
    	{
    		if (A[i * b + j] == 1){
    			Ptemp[i * b + j] = row*b+i;//assigning the highest row-wise vertex
    			if(Ptemp[0 + j]<=Ptemp[i * b + j]){//storing the local columnwise vertex at the start of the column
    				Ptemp[0 + j] = Ptemp[i * b + j];
    			}
     		}
    			
    		else
    			Ptemp[i * b + j] = 0;
    		// P[i * b + j] = 0;//piggybacking, initialising other matrices with zero values
    		// M[i * b + j] = 0;//piggybacking
    	}
    }
    // for (int i = 0; i<b; i++)//replicating locally found max vertex in each process
    // {
    // 	for (int j = 0; j<b; j++)
    // 	{
    // 		Ptemp[i * b + j] = Ptemp[0 + j];
    // 	}
    // }
    //creating communicators for the performing ALLREDUCE columnwise
    MPI_Comm col_comm;
    // int color = rank%q;
    MPI_Comm_split(comm, (rank%q), rank, &col_comm);
    int col_rank, col_size;
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);
    // printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",rank, size, col_rank, col_size);
    for (int j = 0; j<b; j++)
    	{	
    		max = Ptemp[j];
    		// printf("max %d", max);
    		MPI_Allreduce(&max, &P[j], 1, MPI_INT, MPI_MAX, col_comm);
    	}
    MPI_Comm_free(&col_comm);
    for (int i = 0; i<b; i++)//replicating globally found max vertex in each process, along column
    {
    	for (int j = 0; j<b; j++)
    	{
    		P[i * b + j] = P[0 + j];
    	}
    }

    for (int i = 0; i<b; i++)//creating the helper matrix
    {
    	for (int j = 0; j<b; j++)
    	{
    		if (A[i * b + j] == 1){
    			M[i * b + j] = P[i * b + j];//assigning the respective max vertex
     		}
    	}
    }
    //now applying global allreduce row-wise on M
    //first finding the max vertex at block level and storing into Mtemp
    Mtemp=M;
    for (int i = 0; i<b; i++)
    {
    	for (int j = 0; j<b; j++)
    	{
			if(Mtemp[i * b]<=Mtemp[i * b + j]){//storing the block level rowwise max vertex at the start of the row
				Mtemp[i * b] = Mtemp[i * b + j];
			}
    	}
    }
  
    //creating communicators for the performing ALLREDUCE row-wise
    MPI_Comm row_comm;
    // int color = rank%q;
    MPI_Comm_split(comm, (rank/q), rank, &row_comm);
    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    // printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",rank, size, row_rank, row_size);
    for (int j = 0; j<b; j++)
	{	
		max = Mtemp[j*b];
		// printf("max %d", max);
		MPI_Allreduce(&max, &Q[j*b], 1, MPI_INT, MPI_MAX, row_comm);//storing resultant value in matrix Q
	}
	for (int i = 0; i<b; i++)//replicating Q along row-wise
    {
    	for (int j = 0; j<b; j++)
    	{
    		Q[i * b + j] = Q[i*b];
    	}
    }
    MPI_Comm_free(&row_comm);
    for (int i = 0; i<b; i++)
    {
    	for (int j = 0; j<b; j++)
    	{
    		if (Q[i * b + j] == j+col*b){//adding col*b, for adjusting the actual column number
    			M[i * b + j] = P[i * b + j];//assigning the respective max vertex
     		}
    	}
    }
    //now applying global allreduce row-wise again on M
    //first finding the max vertex at block level
    Mtemp=M;
    for (int i = 0; i<b; i++)
    {
    	for (int j = 0; j<b; j++)
    	{
    		if(Mtemp[i * b]<=Mtemp[i * b + j]){//storing the block level rowwise max vertex at the start of the row
				Mtemp[i * b] = Mtemp[i * b + j];
			}
    	}
    }
     //creating communicators for the performing ALLREDUCE row-wise on helper Matrix Mtemp
    // int color = rank%q;
    MPI_Comm_split(comm, (rank/q), rank, &row_comm);
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    // printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",rank, size, row_rank, row_size);
    for (int j = 0; j<b; j++)
	{	
		max = Mtemp[j*b];
		// printf("max %d", max);
		MPI_Allreduce(&max, &PPrime[j*b], 1, MPI_INT, MPI_MAX, row_comm);//storing resultant value in matrix Q
	}
	for (int i = 0; i<b; i++)//replicating PPrime along row-wise
    {
    	for (int j = 0; j<b; j++)
    	{
    		PPrime[i * b + j] = PPrime[i*b];
    	}
    }
    MPI_Comm_free(&row_comm);

    //Printing the Parent Matrix
    sleep(rank);
    // printf("Col %d",col*b);
    std::cout << "Parent Matrix P (" << row << "," << col << ")" << std::endl;
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < b; ++j) std::cout << static_cast<int>(PPrime[i * b + j]) << " ";
        std::cout << std::endl;
    }
    return PPrime[0];
} // connected_components

#endif // A1_HPP
