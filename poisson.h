#include <stdio.h>
#include "poisson.h"
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <cmath>

#define PI 3.14

void poisson_setup (MPI_Comm &comm, int n, MPI_Comm &grid_comm, std::vector<point_t> &x)
{
	//Do MPI_Cart_create
	//MPI_Dims_create

	int taskid;
	int numtasks;
	///int dims[] = {3,3,3};
	int dims[] = {0,0,0};
	int periods[] = {0,0,0};
	

	int num_procs = std::pow(numtasks, 1/3);
    printf("NUM PROCS %d\n", num_procs);
	MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
	MPI_Dims_create (numtasks, num_procs, dims);
	MPI_Cart_create(comm, num_procs, dims, periods, 0, &grid_comm);	
		
    MPI_Comm_rank (MPI_COMM_WORLD, &taskid);
	
	int coord_proc[3];
	//int coord_proc[2];

	printf("NUM PROCS %d NUM TASKS %d\n", num_procs, numtasks);
	MPI_Cart_coords(grid_comm, taskid, num_procs, coord_proc);		

	printf("PROCESSOR COORDINATES %d %d %d\n", coord_proc[0], coord_proc[1], coord_proc[2]);
	real_t m = n/num_procs; // (number of nodes in a single dimension/number of processors in each dimension)
	
	/*
	 * Generate the points and store it in the vector
	 */
	
	int m1  = (int)m;

	int i, j, k;
	point_t t;
	
	real_t h = 1/(real_t)(n-1);
	printf("THE VALUE of H is %lf N is %d\n", h, n);
	for (i=0; i<m1; i++)
	{
		for (j=0;j<m1;j++)
		{
			for (k=0; k<m1; k++)
			{
				//t.x = (coord_proc[0] * m1 + i) * h;
				//t.y = (coord_proc[1] * m1 + j) * h;
				//t.z = (coord_proc[2] * m1 + k) * h;
				//printf("T(x) %lf T(y) %lf T(z) %lf\n", t.x, t.y, t.z);
				//x.push_back(t);
			}
		}
	}


		
	//printf("TASKID %d\n", taskid);
}


void poisson_residual(MPI_Comm &grid_comm, std::vector<real_t> &a, std::vector<real_t> &f, std::vector<real_t> &v, real_t &res)
{
	
}

void displayPoints (std::vector<point_t> x)
{
	printf("Following are the points in its domain\n");
	for(std::vector<point_t>::iterator t=x.begin(); t!=x.end(); t++)
	{
		//printf("(%lf, %lf)\n", t.coord[0], t.coord[1]);	
		//printf("(%lf, %lf)\n", (*t).x, (*t).y);	
		printf("(%lf, %lf, %lf)\n", (*t).x, (*t).y, (*t).z);	
	}
}

void populateF (std::vector<real_t> *f, std::vector<point_t> x)
{
	for(std::vector<point_t>::iterator t=x.begin(); t!=x.end(); t++)
	{
#if 0
		real_t val = sin(4 * PI * (*t).x) * sin(10 * PI * (*t).y) * sin(14 * PI * (*t).z) * (12 + (312 * PI * PI));
		f->push_back (val);
#endif
	}
}

void populateA (std::vector<real_t> *a, std::vector<point_t> x)
{
	for(std::vector<point_t>::iterator t=x.begin(); t!=x.end(); t++)
	{
		//real_t val = sin(4 * PI * (*t).x) * sin(10 * PI * (*t).y) * sin(14 * PI * (*t).z) * (12 + (312 * PI * PI));
		a->push_back (12.0);
	}
}

void populateV (std::vector<real_t> *v, std::vector<point_t> x)
{
	for(std::vector<point_t>::iterator t=x.begin(); t!=x.end(); t++)
	{
		real_t val = sin(4 * PI * (*t).x) * sin(10 * PI * (*t).y) * sin(14 * PI * (*t).z) ;
		v->push_back (val);
	}
}

int main(int argc, char *argv[])
{
	int numtasks, taskid, len;
	char hostname[MPI_MAX_PROCESSOR_NAME];


	MPI_Status st;
	MPI_Comm ncomm = NULL;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank (MPI_COMM_WORLD, &taskid);
    
    MPI_Get_processor_name (hostname, &len);
	std::vector<point_t> x;

	MPI_Comm bcomm = MPI_COMM_WORLD; 

	int n = atoi(argv[1]);	
	/*
	 * Poisson Setup
	 */
	//poisson_setup (bcomm, 3, ncomm, x);	
	//poisson_setup (bcomm, 4, ncomm, x);	
	poisson_setup (bcomm, n, ncomm, x);	

#if 0
	printf("For Processor with rank %d\n", taskid);
	displayPoints (x);
#endif

	std::vector<real_t> a;
	std::vector<real_t> f;
	std::vector<real_t> v;
	real_t res;

	/*
	 * Need to figure out what 'v' is
	 */

	populateF (&f, x);
	populateA (&f, x);
	populateV (&v, x);
	
	/*
	 * Poisson Residual
	 */

	int src, dest;

	//printf("CART SHIFT FOR THE RANK %d\n", taskid);

	
	//TOP and BOTTOM	
	MPI_Cart_shift (ncomm, 0, 1, &src, &dest);
	//LEFT AND RIGHT
	MPI_Cart_shift (ncomm, 1, 1, &src, &dest);
	//ABOVE AND BELOW
	MPI_Cart_shift (ncomm, 2, 1, &src, &dest);
	
	//printf("SRC %d DEST %d\n", src, dest);
#if 0	
	if (src == MPI_PROC_NULL)
		printf("HAIIIII\n");	
#endif
	poisson_residual (ncomm, a, f, v, res);

#if 0
	int ndims = 3;
	int dims[3] = {3,3,3};
	int periods[3] = {0,0,0};

	MPI_Comm ncomm;

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &ncomm);

	int coord[3];
	//if (taskid == 4)
	{
		MPI_Cart_coords(ncomm, taskid, 3, coord);
		printf("Coordinates of processor ranked %d is (%d, %d, %d)\n", taskid, coord[0], coord[1], coord[2]);
	}
#endif
	MPI_Finalize();
	return 0;	
}
