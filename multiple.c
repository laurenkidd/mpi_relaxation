/*
 * Test this example with 
 * 	mpirun -n $X ./scatter_gather
 * where for X in {2,4,8,16}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define COLS (5)
#define ROWS (8)

extern double get_change(double *x, double *y, int n);
extern void relax(double *dest,double *srce, int cols, int rows);
extern void init_grid(double **,double **,int cols, int rows);
extern void init_boundaries(double *,int,int rows);
extern void print_buffer(double *,int cols,int rows);

int main(int argc, char **argv) {

	//initialize MPI and stuff
	
   	int        rank,  size, j;
	MPI_Init(&argc, &argv);

  	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   	MPI_Comm_size( MPI_COMM_WORLD, &size );
	int subRows = ROWS/size;

    
   	MPI_Status status;
	//initialize grids
	double  *p,*p_new, *local, *local_new;
	init_grid(&local,&local_new,COLS,subRows+2);
	if (0 == rank) {
		init_grid(&p,&p_new,COLS,ROWS+2);
		init_boundaries(p,COLS,ROWS+2);
		memmove(p_new,p, COLS*(ROWS+2) * sizeof(double) );
		print_buffer(p_new,COLS,ROWS+2);
	}
	//wait for errbody
	MPI_Barrier(MPI_COMM_WORLD);

	//distribute local grid to workers
	if (rank ==0){	
		for (j=1; j<size; j++){
			MPI_Send(p+j*subRows*COLS,(subRows+2)*COLS, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Recv(local, (subRows+2)*COLS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	

	//now do the solving stuff
	double *first_buf, *last_buf;
	int count  = 1000;
	init_grid(&first_buf, &last_buf, COLS, 1);
	while(count-- > 0){
		relax(local_new, local, COLS, subRows+2);
		//make come copies
		memmove(last_buf, local_new, COLS*sizeof(double) );
		memmove(first_buf, local_new+(subRows+1)*COLS, COLS*sizeof (double) );
		int up_number = (rank == size-1)? MPI_PROC_NULL: rank+1;
		int down_number = (rank==0)? MPI_PROC_NULL : rank-1; 
		//send single row from down to up and from up to down unless you are the top or bottom row; This is sending the missing row from each
		//MPI_SendRecv(thing to send, size, type, dest, tag, thing received, size, type, source, tag, comm, status)
		MPI_Sendrecv(local_new+(1)*COLS, COLS, MPI_DOUBLE, down_number,0,  first_buf, COLS, MPI_DOUBLE, up_number, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(local_new+(subRows)*COLS, COLS, MPI_DOUBLE, up_number,0,  last_buf, COLS, MPI_DOUBLE, down_number, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		//wait for everybody to join the party
		MPI_Barrier(MPI_COMM_WORLD);

		memmove(local_new ,last_buf,COLS*sizeof(double));
		memmove(local_new + (subRows+1)*COLS,first_buf,COLS*sizeof(double));
		
		// copy updated into current grid
		memmove(local,local_new,COLS*(subRows+2)*sizeof(double));
	}


	//Send local grids back to thread 0
	MPI_Send(local_new+1*COLS,(subRows)*COLS,MPI_DOUBLE,0,(subRows+2),MPI_COMM_WORLD);
	//have thread 0 receive locals
	if ( 0 == rank ){
		for ( int i = 0; i < size; i++) {
			MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Recv(p_new+(status.MPI_SOURCE*subRows+1)*COLS,(subRows)*COLS,MPI_DOUBLE,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}
        }

	// have thread 0 print resulting matrix
	if ( 0 == rank ){
		print_buffer(p_new,COLS,ROWS+2);
		if ( 0 == rank ) {
			free(p_new);
			free(p);
		}
	}
	free(local);
	free(local_new);
	MPI_Finalize();




	/*ORIGINAL CODE
	// Do some calculations.
	int count = 0;
	while (count++ < 300) {
		relax(p_new,p,COLS,ROWS+2); //the +2 to rows: b/c c stores data row-wise & want to break up data 
                //make rows integer powers of 2, so you can make 2, or 4 or etc processes
                // +2 -> adding buffer row to top and bottom (for b.c)
		if (0.000001< get_change(p_new,p,COLS*(ROWS+2)))
			printf("Change/Iteration %f/%d\n",get_change(p_new,p,COLS*(ROWS+2)),count); //change at each iteration
		memmove(p,p_new, COLS*(ROWS+2) * sizeof(double) );
	}*/
	return 0;




}

void relax(double *new,double *old, int cols, int rows){
	for ( int j = 1 ; j < (rows)-1; j++) {
		for ( int i = 1 ; i < cols-1; i++) {
                        //assuming only charge at borders of mesh; assuming points along the border will be held const
			new[i+j*cols] = 0.25*(old[i-1+j*cols]+old[i+1+j*cols]+old[i+(j-1)*cols]+old[i+(j+1)*cols]);//relaxtion func
//			printf("  %d   \n",i+(j-1)*COLS);
//			printf("%d  %d  %d\n",i-1+j*COLS,j+i*COLS,j+1+i*COLS);
//			printf("  %d   \n",i+(j+1)*COLS);
		}
//		puts(" ");
//		puts(" ");
	}


}

void init_boundaries(double *l_p,int cols,int rows){
/* 
	for ( int i = 0 ; i < cols; i++) {
		l_p[i] = -1;
		l_p[i + cols*(rows -1)] = 1;
	}

	for ( int i = 0 ; i < rows; i++) {
		l_p[i*cols] = -1;
		l_p[(i+1)*cols -1] = 1; //column
		// l_p[i*(cols-1)+rows-1] = -1; //diagonal
	}
*/
	for (int i = 0 ; i < cols*rows; i++ )
		l_p[i] = 0;
	l_p[cols/2] = 1;
}


void init_grid(double **p, double **p_new,int cols, int rows){
	if (NULL == (*p = malloc(  cols*rows * sizeof(double) )  )  ) {
		puts("Allocation Error.");
		exit(99);
	}
	if (NULL == (*p_new = malloc(cols*rows * sizeof(double)))) {
		puts("Allocation Error.");
		exit(99);
	}
	for ( int i = 0 ; i < cols*rows; i++){

		(*p)[i] = 0. ;
		(*p_new)[i] = 0. ;
	}
}

void print_buffer(double *p, int cols, int rows){
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%f\t",p[i*cols+j]);
		}
		puts(" ");
	}


}

double get_change(double *x, double *y, int n){
	double result=0.;
	for (int i=0; i<n;i++)
		result+=x[i]*x[i]-y[i]*y[i];
	return result;
}

