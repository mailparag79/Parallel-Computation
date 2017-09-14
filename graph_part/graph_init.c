#include "graph.h"

void InitLibrary(int nnodes, idxtype *index, idxtype *edges, int *process_id, int global_buffer_size, int rank, int size)
{

  int nodenum,edgenum,i,j,jj, k, **A,mybegin,myend, ii,gId,
  p_src_id, p_dst_id, boundary_src, boundary_dst, local_numnodes, weight, no_processor;
  int c, c1;
	  
  int **send_buffer, **recv_buffer, *data;

  if ((nodes = (node*)malloc(nnodes*sizeof(node))) == NULL)
  {
	  fprintf(stderr, "%d Unable to allocate memory \n", rank);
	  exit(0);
  }
  if ((proc_nodes = (proc_node*)malloc(size*sizeof(proc_node))) == NULL)
  {
	  fprintf(stderr, "%d Unable to allocate memory \n", rank);
	  exit(0);
  }
  
  buffer_size = global_buffer_size;
  jj = 0;
  A = (int**)malloc(nnodes*sizeof(int*));
  for(i=0;i<nnodes;i++)
	  A[i] = (int*)calloc(nnodes, sizeof(int));

  j=0;
  int local_node_count = 0;
  if (DEBUG)
  {
	  fprintf(stdout, "%d nnodes is %d \n", rank, nnodes);
	  for (j=0;j<nnodes; j++)
	  {
		  fprintf(stderr, "nodes %d is in %d processor \n", j, process_id[j]);
	  }
	  fprintf(stdout, "%d Building Node structure \n", rank);
  }
  for (i=0;i<nnodes;i++)
  {
	 /*for each entry in the index array,
		scan the edges and add the edge for each node*/
	j=0;
	c=0;
	nodes[i].proc_id = process_id[i];
	nodes[i].visited = 0;
	nodes[i].global_node_id = i;
	if(nodes[i].proc_id == rank)
	{
		proc_nodes[nodes[i].proc_id].num_nodes++;
		nodes[i].local_node_id = local_node_count++;
	}
	int start = index[i];
	int end = index[i+1];
	if (DEBUG)
		fprintf(stderr, "%d: node=%d start=%d, end=%d\n", rank, i, start, end);
	for (j=start;j<end;j++)
	{
		if (DEBUG)
			fprintf(stderr, "%d: adding onnection from %d %d \n", rank, i, edges[j]);
		graphadd(A, nnodes, i, edges[j], 1);
	}
	fprintf(stderr, "%d: done %d node processing\n", rank, i);
  	MPI_Barrier(MPI_COMM_WORLD);
  }
  fprintf(stderr, "%d: done initializing library \n", rank);
  if (DEBUG)
  {
	  for(i=0;i<nnodes;i++)
	  {
		if (nodes[i].proc_id == rank)
		{
			fprintf(stderr, "INITIAL_NODE_LIST: %d, %d %d \n", rank, nodes[i].global_node_id, nodes[i].local_node_id);
		}
	  }
  }
  //printConnections(A, nnodes);

 
  //make local tree based on global_buffer_size
  int inode_type = 0;
  fprintf(stdout, "%d initializing node type \n", rank);
  proc_nodes[rank].nodes = (node*)malloc(nnodes*sizeof(node));
  for(i=0;i<nnodes;i++){

	  proc_nodes[rank].nodes[i].global_node_id = nodes[i].global_node_id;
	  if (nodes[i].proc_id == rank){
		  proc_nodes[nodes[i].proc_id].nodes[i].local_node_id = nodes[i].local_node_id;
	  }
	  proc_nodes[rank].nodes[i].proc_id = nodes[i].proc_id;
	  
	  for(j=0;j<nnodes;j++){
		  if(A[i][j]){
			  if(nodes[j].proc_id != nodes[i].proc_id && (buffer_size !=0)){
				nodes[i].node_type = inode_type | BOUNDARY;
				buffer_size = buffer_size - 1;
			  }
			  else if(nodes[j].proc_id != nodes[i].proc_id && (buffer_size ==0)){
				  nodes[i].node_type = inode_type | GHOST;
			  }
			  else{
				nodes[i].node_type = inode_type | INTERNAL;
			  }
		  }
		  proc_nodes[rank].nodes[i].node_type = nodes[i].node_type;
	  }
	  buffer_size = global_buffer_size;
  }//end of for
  buffer_size = global_buffer_size;
  fprintf(stdout, "%d initializing send recv list \n", rank);
  //create send/receive list
  if(DEBUG)
  {
	  for(i=0;i<nnodes;i++){
		  if(proc_nodes[rank].nodes[i].proc_id == rank)
			fprintf(stderr, "NODE_ID: %d, %d \n", rank, proc_nodes[rank].nodes[i].global_node_id);
	  }
  }
  proc_nodes[rank].send_list = (int*)malloc(size*sizeof(int));
  proc_nodes[rank].recv_list = (int*)malloc(size*sizeof(int));
  c = 0; c1 = 0;
  //traverse local tree based on process_id
  //initialize send and receive list
  fprintf(stdout, "%d Initializing send array of size %d \n", rank, size);
  InitArrayList(&proc_nodes[rank].send_list, size);
  fprintf(stdout, "%d Initializing recv array of size %d \n", rank, size);
  InitArrayList(&proc_nodes[rank].recv_list, size);
  for(i=0;i<nnodes;i++){
	  if (proc_nodes[rank].nodes[i].proc_id != rank)
		  continue;
	  for(j=0;j<nnodes;j++){

		  //These nodes are directly connected, so it exist in both send and receive list
		  if((A[i][j] >= 1) && (proc_nodes[rank].nodes[i].proc_id != proc_nodes[rank].nodes[j].proc_id))
		  {

			  //TODO: Make more space efficiant using Key value pair instead of pigeon hole.
			  if(!processAdded(proc_nodes[rank].send_list, proc_nodes[rank].nodes[j].proc_id, c))
			  {
				proc_nodes[rank].send_list[c] = proc_nodes[rank].nodes[j].proc_id;
				++c;
				++proc_nodes[rank].sender_count;
			  }
			  if(!processAdded(proc_nodes[rank].recv_list, proc_nodes[rank].nodes[j].proc_id, c1))
			  {
			    proc_nodes[rank].recv_list[c1] = proc_nodes[rank].nodes[j].proc_id;
				++c1;
				++proc_nodes[rank].receiver_count;
			  }

		  }
		  //if not directly connected, and if not aready added then just add it to the recv list
		  else if((A[i][j] <= 0) && (proc_nodes[rank].nodes[i].proc_id != proc_nodes[rank].nodes[j].proc_id))
		  {
			  if(!processAdded(proc_nodes[rank].recv_list, proc_nodes[rank].nodes[j].proc_id, c1))
			  {
				  proc_nodes[rank].recv_list[c1] = proc_nodes[rank].nodes[j].proc_id;
				  ++proc_nodes[rank].receiver_count;
				  ++c1;
			  }
		  }
		  
	  }
  }


  if(EDEBUG)
  {
	  fprintf(stderr, "%d: Number of receivers: %d\n", rank, proc_nodes[rank].receiver_count);
	  for(i=0;i<size;i++){
		fprintf(stderr, "SEND: %d, %d \n", rank, proc_nodes[rank].send_list[i]);
	  }
	  for(i=0;i<size;i++){
		fprintf(stderr, "RECV: %d, %d \n", rank, proc_nodes[rank].recv_list[i]);
	  }
  }
}

void SetDataCallback(DataCB function){
	CB_SetData = function;
}

int ReadGraphPartitionFile(char *meshfile, int numberOfParts, int **processor_id, idxtype **xadj, idxtype **adjncy)
{
  int i, j, numflag=0;

  int edgecut;
  idxtype *elmnts;
  idxtype ne;
  int etype;
  int nn, k;
  char etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};
  int wgtflag = 0;
  int options[5];
  elmnts = ReadMesh(meshfile, &ne, &nn, &etype);

  *xadj = (idxtype*)malloc(nn+1*sizeof(idxtype));
  *adjncy = (idxtype*)malloc(20*nn*sizeof(idxtype));
  idxtype* wxadj = (idxtype*)calloc(nn+1, sizeof(idxtype));
  idxtype* wadjncy = (idxtype*)calloc(ne, sizeof(idxtype));

  fprintf(stdout, "Converting mesh with %d node and %d elements to nodal graph of type %d\n", nn, ne, etype);
  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, *xadj, *adjncy);

  idxtype* parts = (idxtype*)malloc(nn*sizeof(idxtype));
  *processor_id = (int*)malloc(nn*sizeof(int));
  
 /*
 * Metis partitioning Parameters
 * n: The number of vertices in the graph.
 * xadj, adjncy: The adjacency structure of the graph.
 * vwgt, adjwgt: Information about the weights of the vertices and edges
 * wgtflag: Used to indicate if the graph is weighted. wgtflag can take the following values:
	0 No weights (vwgts and adjwgt are NULL)
	1 Weights on the edges only (vwgts = NULL)
	2 Weights on the vertices only (adjwgt = NULL)
	3 Weights both on vertices and edges.
 * numflag: Used to indicate which numbering scheme is used for the adjacency structure of the graph. numflag
	0 C-style numbering is assumed that starts from 0
	1 Fortran-style numbering is assumed that starts from 1
 * nparts: The number of parts to partition the graph.
 * options: This is an array of 5 integers that is used to pass parameters for the various phases of the algorithm.
	If options[0]=0 then default values are used. If options[0]=1, then the remaining four elements of
	options[1] Determines the matching type. Possible values are:
		1 Random Matching (RM)
		2 Heavy-Edge Matching (HEM)
		3 Sorted Heavy-Edge Matching (SHEM) (Default)
	options[2] Determines the algorithm used during initial partitioning. Possible values are:
		1 Multilevel recursive bisection (Default)
	options[3] Determines the algorithm used for refinement. Possible values are:
		1 Random boundary refinement
		2 Greedy boundary refinement
		3 Random boundary refinement that also minimizes the connectivity among the subdomains (Default)
	options[4] Used for debugging purposes. Always set it to 0 (Default).
 * edgecut: Upon successful completion, this variable stores the number of edges that are cut by the partition.
 * part: This is a vector of size n that upon successful completion stores the partition vector of the graph. 
 	 The numbering of this vector starts from either 0 or 1, depending on the value of numflag.
 */ 
  options[0] = 0;
  if (numberOfParts > 8)
  {
  	fprintf(stderr, "partitioning graph of nodes (%d) and elements (%d) using metis k way method in %d parts\n", nn, ne, numberOfParts);
  	METIS_PartGraphKway(&nn, *xadj, *adjncy, NULL, NULL, &wgtflag /*0*/, &numflag /*0*/, &numberOfParts, options, &edgecut, parts);
  }
  else
  {
  	fprintf(stderr, "partitioning graph of nodes (%d) and elements (%d) using metis recursive method in %d parts\n", nn, ne, numberOfParts);
 	METIS_PartGraphRecursive(&nn, *xadj, *adjncy, wxadj, wadjncy, &wgtflag /*0*/, &numflag /*0*/, &numberOfParts, options, &edgecut, parts);
  }
  memcpy(*processor_id, parts, nn*sizeof(int));
  
  /*
  FILE *fpin;
  if ((fpin = fopen(pDistFile, "r")) == NULL) {
	printf("Failed to open file %s\n", pDistFile);
	exit(0);
  }
  i = 0;
  while (!feof(fpin)) 
  {
	fscanf(fpin, "%d\n", ((*processor_id)+i));
	i++;
  }
  fclose(fpin);*/
  free(elmnts);
  free(parts);
  return nn;
}

void Graph_init(int argc, char **argv, char *meshFile, int global_buffer_size, int *rank)
{

  int nodenum,edgenum,i,j,jj, k, **A,mybegin,myend, size, ii,
  p_src_id, p_dst_id, boundary_src, boundary_dst, local_numnodes, weight, no_processor;
  int **send_buffer, **recv_buffer, *data;
  int nnodes, *processor_id, etype;
  idxtype *index, *edges, *elmnts, ne;


  //node *nodes = (node*)malloc(nnodes*sizeof(node));
  //proc_node *proc_nodes = (proc_node*)malloc(size*sizeof(proc_node));
  MPI_Status stat;
  /*Initialize MPI*/
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  comm_size = size;
  
  if(DEBUG)
  {
  	  fprintf(stderr, "%d: MPI initailised\n", *rank);
	  fprintf(stderr, "size = %d, rank=%d \n", size, *rank);
	  fprintf(stderr, "meshfile=%s\n", meshFile);
  }
  nnodes = ReadGraphPartitionFile(meshFile, size, &processor_id, &index, &edges);
  if (DEBUG)
  {
	  fprintf(stderr, "%d nodes list (%d) \n", *rank, nnodes);
	  for(j=0;j<nnodes;j++)
	  {
		  fprintf(stderr, "%d ", index[j]);
	  }
	  fprintf (stderr, "\n %d adjan list (%d) \n", *rank, index[nnodes]);
	  for(j=0;j<index[nnodes];j++)
	  {
		  fprintf(stderr, "%d ", edges[j]);
	  }
  }
  if(DEBUG)
  {
	  fprintf (stderr, "\n%d processor list \n", *rank);
	  for (j=0;j<nnodes; j++)
	  {
		  fprintf(stderr, "%d initializing nodes %d is in %d processor \n", *rank, j, processor_id[j]);
	  }
	  fprintf(stderr, "\n");
  }
  InitLibrary(nnodes, index, edges, processor_id, global_buffer_size, *rank, size);
  //Read the local data and store it in the structure using CALLBACK function
  //
  proc_nodes[*rank].data_size = 1;
  CB_SetData(&proc_nodes[*rank].data, &proc_nodes[*rank].num_bytes, *rank, &proc_nodes[*rank].data_size);
  
  proc_nodes[*rank].recvvector = (int*)malloc(proc_nodes[*rank].data_size*proc_nodes[*rank].receiver_count*sizeof(int));

  /*Find the Minimum Spanning tree for each processor*/
  //Kruskal_MST(proc_nodes[rank]->node_matrix, proc_nodes[rank]->num_nodes, proc_nodes[rank]->num_nodes-1);

  //distribute it based on processor id

}

void InitArrayList(int **arrayList, int size)
{
	int i;
	for (i=0;i<size;i++)
	{
		*((*arrayList)+i) = -1;
		//fprintf(stdout, "Initializing array element %d to %d\n", i, *((*arrayList)+i));
	}
}

void graphadd(int **A,int width, int m,int n, int w)
{	//A[m][n]=w;
	A[m][n] = w;
}

int processAdded(int *list, int proc_id, int n)
{
	int i;
	for (i = 0;i<n;i++)
	{
		if (list[i] == proc_id)
		{
			return 1;
		}
	}
	return 0;
}

/* a fatal error has occurred; print mesg and exit */
int panic(char *mesg, int n){
  fprintf(stderr,"****** Node %d failed with error \"%s\"; exiting. ******\n",
          n, mesg);
  exit(0);
}
