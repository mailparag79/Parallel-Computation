#include <stdio.h>
#include "graph.h"

void ParseArgs(int argc, char **argv, char *meshFile, char *pDistFile);


int main(int argc,char *argv[])
{	

	FILE *fp;char newline[10],*tempch=" ";
	char filename[80];

	int rank;
	const int nnodes = 6;
	int global_buffer_size = 3;
	char meshFile[80], pDistFile[80];

	ParseArgs(argc, argv, meshFile, pDistFile);
	SetDataCallback((DataCB)&SetDataFromFile);

	/* MPI necessities */
	fprintf(stderr, "initialising graph\n");
	//Graph_init(argc, argv, nnodes, index, edges, process_id, global_buffer_size, &rank);
	Graph_init(argc, argv, meshFile, global_buffer_size, &rank);
	/*Exchange boundary*/
	Exchange_boundary();

    
	//assert(no_processor != size);

	/* freeing the matrix */
	//graphfree(A,nodenum);
	return 0;
}

void SetDataFromFile(void **dataArray, int *length, int rank, int *data_size)
{
	FILE *fp;
	char filename[80];
	//int lFileLen;
	
	sprintf(filename, "data_%d.dat", rank);
	fprintf(stderr, "%d: Copying data from %s of size %d \n", rank, filename, *data_size);
	if((fp=fopen(filename, "r"))==NULL)
		panic("SetData failed", rank);
	*data_size = ftell(fp); 
	*length = sizeof(int);
	*dataArray = (int*)calloc(*data_size, sizeof(int));
	rewind(fp);
	fread(*dataArray, *data_size, sizeof(int), fp);
	fprintf (stdout, "%d data loaded \n", rank);
	return;
}

void ParseArgs(int argc, char **argv, char *meshFile, char *pDistFile)

{
	
	if (argc < 2)
	{	
		printf("Usage: GraphDist meshFile (metis format)\n");
		exit(0);
	}

	
	strcpy(meshFile, argv[1]);
	//strcpy(pDistFile, argv[2]);
	//meshFile = argv[1];
	
	//pDistFile = argv[2];

}

