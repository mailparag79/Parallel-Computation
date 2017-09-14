#include "graph.h"

void Exchange_boundary()
{
	/*Send neighbour list to each of the neighbour, using BUFFER_TAG, and send receive the same from each sender.
	  Send local data to each*/
	int i, j;
	MPI_Status stat;
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//int nnodes = proc_nodes[proc_nodes.proc_rank].num_nodes;
	for(i=0;i<proc_nodes[rank].receiver_count;++i)
	{
		if (proc_nodes[rank].recv_list[i] >= 0 && proc_nodes[rank].recv_list[i] < size)
		{
			if (EDEBUG)
			{
				fprintf(stdout, "%d receving data from %d \n", rank, proc_nodes[rank].recv_list[i]);
			}
			MPI_Irecv((void *)(&proc_nodes[rank].recvvector), proc_nodes[rank].data_size, MPI_INT, 
				proc_nodes[rank].recv_list[i], BUFFER_TAG, MPI_COMM_WORLD, rreq+i);
			if (EDEBUG)
			{
				fprintf(stdout, "%d received data from %d \n", rank, proc_nodes[rank].recv_list[i]);
			}
		}
	
	}
	for(i=0;i<proc_nodes[rank].sender_count;i++)
	{
		
		if (proc_nodes[rank].send_list[i] >= 0 && proc_nodes[rank].send_list[i] < size)
		{
			if (EDEBUG)
			{
				fprintf(stdout, "%d sending data to %d \n", rank, proc_nodes[rank].recv_list[i]);
			}
			MPI_Isend(proc_nodes[rank].data, proc_nodes[rank].data_size, MPI_INT, proc_nodes[rank].send_list[i],BUFFER_TAG, MPI_COMM_WORLD, sreq+i);
			if (EDEBUG)
			{
				fprintf(stdout, "%d sent data to %d \n", rank, proc_nodes[rank].recv_list[i]);
			}
		}
	}
	/* OK, now get into the loop waiting for the operations to finish */

	MPI_Waitall(proc_nodes[rank].receiver_count, rreq, statuses);
	MPI_Waitall(proc_nodes[rank].sender_count, sreq, statuses);
	fprintf(stdout, "%d exchange boundary completed \n", rank);
}
