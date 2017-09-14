#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <metis.h>
#include <proto.h>
#define BUFFER_TAG 100
#define EXCHANGE_TAG 101
#define DEBUG 0
#define EDEBUG 1
#define MAX_PES			8192
/*Data Structure for Distributed Graph*/
int buffer_size;
static int comm_size;
typedef enum _node_type_t
{
	INTERNAL=0,
	BOUNDARY=1,
	GHOST=2
}_node_type;

typedef struct node_t
{
	//int boundary;
	int node_id;
	int visited;
	int global_node_id;
	int local_node_id;
	int proc_id;
	_node_type node_type;
}node;

typedef struct proc_node_t
{
	void *data;
	int num_bytes;
	int data_size;
	node *nodes;
	int num_nodes;
	int *send_list;
	int *recv_list;
	int **node_matrix;
	int proc_rank;
	int *recvvector;
	int *recvptr;
	int sender_count;
	int receiver_count;
}proc_node;

//typedef int idxtype;
//typedef double realtype;

MPI_Request sreq[MAX_PES], 
              rreq[MAX_PES];    /* MPI send and receive requests */
MPI_Status statuses[MAX_PES];
MPI_Status status;

node *nodes;
proc_node *proc_nodes;

typedef void (*DataCB)(void* data, int *length, int rank, int *data_size) ;
DataCB CB_SetData;
/* This program use adjacency matrix */

/************* list of all functions **************/

void SetDataCallback(DataCB function);

void SetDataFromFile(void **dataArray, int *length, int rank, int *data_size);

/* a fatal error has occurred; print mesg and exit */
int panic(char *mesg, int n);

/*Graph Initialization*/
void Graph_init(int argc, char **argv, char *meshFile, int global_buffer_size, int *rank);

/*Read local data its a CALLBACK function*/

/*Graph Initialization*/
void InitLibrary(int nnodes, idxtype *index, idxtype *edges, int *process_id, int global_buffer_size, int rank, int size);

/* insert a new row in the middle of array ba */
int arrayshift(int ba[100][20],int i,int nodenum);

/*Loads the data from the file into individual node*/
void load_data(node*, char*);

/*Pack the strcuture of graph into an array for sending*/
void pack_array(int*);

/*Unpack the data and store it into strcuture of graph after Receiving*/
void unpack_array(int*, proc_node *, node *);

/*Find the boundaryy nodes in spanning tree and mark it as boundary*/
int find_boundary(proc_node *, node *, int);

/*Exchanges the data between the boundary nodes of each processor*/
void Exchange_boundary();

/* add an edge, m=start vertice, n=end vertice */
void graphadd(int **A,int width, int m,int n, int w);

/* allocate memory for the matrix */
int **graphalloc(int nodenum);

/* free memory */
void graphfree(int **A,int nodenum);

/*Check if Process already added*/

int processAdded(int *list, int proc_id, int n);

void printConnections(int **A, int numnodes);

void InitArrayList(int **, int );
int ReadGraphPartitionFile(char *meshfile, int numberOfParts, int **processor_id, idxtype **xadj, idxtype **adjncy);

int GetNumberOfNodes(char *meshfile, idxtype **elmnts, idxtype *ne, int *etype);
