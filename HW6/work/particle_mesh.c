# include "mpi.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


double L = 1.0;			
int N = 32;			

double *u, *u_new;

double *gradX, *gradY;		

/* macro to convert 2D Array to a 1D array */
#define INDEX(i,j) ((N+2)*(i)+(j))

int my_rank;		

int *proc;			
int *i_min, *i_max;		
int *left_proc, *right_proc;	

int main ( int argc, char *argv[] );
void allocate_arrays ( );
void jacobi ( int num_procs, double f[] );
void make_domains ( int num_procs );
double *make_rhs ( );

/* Function Definitions */
void Form_Gradient(int num_procs, double f[]) {
  double h;
  int i;
  int j;
  MPI_Request request[4];
  int requests;
  MPI_Status status[4];
/*
  H is the lattice spacing.
*/
  h = L / ( double ) ( N + 1 );
/* 
  Update overlays using non-blocking send/receive 
*/
  requests = 0;

  if ( left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(i_min[my_rank] - 1, 1), N, MPI_DOUBLE,
      left_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(i_min[my_rank], 1), N, MPI_DOUBLE,
      left_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );
  }

  if ( right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(i_max[my_rank] + 1, 1), N, MPI_DOUBLE,
      right_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(i_max[my_rank], 1), N, MPI_DOUBLE,
      right_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );
  }
/* 
  Gradient calculations for internal vertices in my domain.
*/
  for ( i = i_min[my_rank] + 1; i <= i_max[my_rank] - 1; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      gradX[INDEX(i,j)] =
         ( (u[INDEX(i+1,j)] - u[INDEX(i-1,j)]) / h);

      gradY[INDEX(i,j)] =
         ( (u[INDEX(i,j+1)] - u[INDEX(i,j-1)]) / h);

    }
  }

  /* 
  Wait for all non-blocking communications to complete.
*/
  MPI_Waitall ( requests, request, status );
/* 
  Update gradient calculations for boundary vertices in my domain.
*/
  i = i_min[my_rank];
  for ( j = 1; j <= N; j++ )
  {
    gradX[INDEX(i,j)] =
          ( (u[INDEX(i+1,j)] - u[INDEX(i-1,j)]) / h);

    gradY[INDEX(i,j)] =
          ( (u[INDEX(i,j+1)] - u[INDEX(i,j-1)]) / h);

  }

  i = i_max[my_rank];
  if (i != i_min[my_rank])
  {
    for (j = 1; j <= N; j++)
    {
      gradX[INDEX(i,j)] =
         ( (u[INDEX(i+1,j)] - u[INDEX(i-1,j)]) / h);

      gradY[INDEX(i,j)] =
         ( (u[INDEX(i,j+1)] - u[INDEX(i,j-1)]) / h);

    }
  }

  return;

  
}

void allocate_arrays ( ) 
{
  int i;
  int ndof;

  ndof = ( N + 2 ) * ( N + 2 );

  gradX = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++)
  {
    gradX[i] = 0.0;
  }

  gradY = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++ )
  {
    gradY[i] = 0.0;
  }

  u = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++)
  {
    u[i] = 0.0;
  }

  u_new = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++ )
  {
    u_new[i] = 0.0;
  }

  return;
}

/* Set up all of the RHS border and the inner using function provided */
double *make_rhs ( ) 
{
  double *f;
  int i;
  int j;
  int k;
  double q;
  double phi0 = -10.0;
  double phi1 = 10.0;

  f = ( double * ) malloc ( ( N + 2 ) * ( N + 2 ) * sizeof ( double ) );

  for ( i = 0; i < ( N + 2 ) * ( N + 2 ); i++ )
  {
    f[i] = 0.0;
  }
 
  for ( j = 0; j < N; j++ )
  {
    
    for ( i = 0; i < N; i++ )
    {
      /* Left boundary */
      if ( i == 0 )
      {
        f[INDEX(i,j)] = phi0;
      }
      /* Right boundary */
      if ( i == N - 1 )
      {
        f[INDEX(i,j)] = phi1;
      }
      /* top and bottom */
      if ( j == 1 || j == N - 1 )
      {
        f[INDEX(i,j)] = (((1 - (double)(i)/N)*phi0) + ((((double)(i)/N)) * phi1) ) ;
      }
      else
      {
        
        f[INDEX(i,j)]=0.0;
      }
    }
  }

  return f;
}

void jacobi ( int num_procs, double f[] ) 
{
  double h;
  int i;
  int j;
  MPI_Request request[4];
  int requests;
  MPI_Status status[4];
/*
  H is the lattice spacing.
*/
  h = L / ( double ) ( N + 1 );
/* 
  Update ghost layers using non-blocking send/receive 
*/
  requests = 0;

  if ( left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(i_min[my_rank] - 1, 1), N, MPI_DOUBLE,
      left_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(i_min[my_rank], 1), N, MPI_DOUBLE,
      left_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );
  }

  if ( right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(i_max[my_rank] + 1, 1), N, MPI_DOUBLE,
      right_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(i_max[my_rank], 1), N, MPI_DOUBLE,
      right_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );
  }
/* 
  Jacobi update for internal vertices in my domain.
*/
  for ( i = i_min[my_rank] + 1; i <= i_max[my_rank] - 1; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      u_new[INDEX(i,j)] =
        0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
                 u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
                 h * h * f[INDEX(i,j)] );
    }
  }
/* 
  Wait for all non-blocking communications to complete.
*/
  MPI_Waitall ( requests, request, status );
/* 
  Jacobi update for boundary vertices in my domain.
*/
  i = i_min[my_rank];
  for ( j = 1; j <= N; j++ )
  {
    u_new[INDEX(i,j)] =
      0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
               u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
               h * h * f[INDEX(i,j)] );
  }

  i = i_max[my_rank];
  if (i != i_min[my_rank])
  {
    for (j = 1; j <= N; j++)
    {
      u_new[INDEX(i,j)] =
        0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
                 u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
                 h * h * f[INDEX(i,j)] );
    }
  }

  return;
}

void make_domains ( int num_procs ) 
{
  double d;
  double eps;
  int i;
  int p;
  double x_max;
  double x_min;
/* 
  Allocate arrays for process information.
*/
  proc = ( int * ) malloc ( ( N + 2 ) * sizeof ( int ) );
  i_min = ( int * ) malloc ( num_procs * sizeof ( int ) );
  i_max = ( int * ) malloc ( num_procs * sizeof ( int ) );
  left_proc = ( int * ) malloc ( num_procs * sizeof ( int ) );
  right_proc = ( int * ) malloc ( num_procs * sizeof ( int ) );
/* 
  Divide the range [(-eps+1)..(N+eps)] evenly among the processes.
*/
  eps = 0.0001;
  d = ( N - 1.0 + 2.0 * eps ) / ( double ) num_procs;

  for ( p = 0; p < num_procs; p++ )
  {
/* 
  Process the domain.
*/
    x_min = - eps + 1.0 + ( double ) ( p * d );
    x_max = x_min + d;
/* 
  Identify vertices belonging to this process.
*/
    for ( i = 1; i <= N; i++ )
    {
      if ( x_min <= i && i < x_max )
      {
        proc[i] = p;
      }
    }
  }

  for ( p = 0; p < num_procs; p++ )
  {
/* 
  Find the smallest vertex index in the domain.
*/
    for ( i = 1; i <= N; i++ )
    {
      if ( proc[i] == p )
      {
        break;
      }
    }
    i_min[p] = i;
/* 
  Find the largest vertex index in the domain.
*/
    for ( i = N; 1 <= i; i-- )
    {
      if ( proc[i] == p )
      {
        break;
      }
    }
    i_max[p] = i;
/* 
  Find processes to left and right. 
*/
    left_proc[p] = right_proc[p] = -1;

    if ( proc[p] != -1 ) 
    {
      if ( 1 < i_min[p] && i_min[p] <= N )
      {
        left_proc[p] = proc[i_min[p] - 1];
      }
      if ( 0 < i_max[p] && i_max[p] < N )
      {
        right_proc[p] = proc[i_max[p] + 1];
      }
    }
  }

  return;
}

void Save_Data(char* name,double grad[],int myTaskId) {

   FILE * fpGradient;
      
   int i,j;
   i=0;
   j=0;

   char buffer[16];
     
   sprintf(buffer, "%s%d.out", name,myTaskId);
   

   if((fpGradient=fopen(buffer, "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   for ( i = i_min[my_rank]; i <= i_max[my_rank]; i++ )
    {
      for ( j = 1; j <= N; j++ )
      {   
         fprintf(fpGradient,"%d, %d, %f\n",i,j,grad[INDEX(i,j)]);
      }
   }
   
   fclose(fpGradient);
}

double error_per_id ( int N,int my_rank, double a[] )

{
  int i;
  int j;
  double err;

  err = 0.0;

  int low = i_min[my_rank];
  int high = i_max[my_rank];

  printf("low is %d high is %d",low,high);

  for ( i = low; i <= high; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      err = err + a[INDEX(i,j)] * a[INDEX(i,j)];
      //printf("err %f a %f\n",err,a[INDEX(i,j)]);
    }
  }
  
  return err;
}



int main ( int argc, char *argv[] ) 

{
  double change;
  double epsilon = 1.0E-03;
  double *f;
  char file_name[100];
  int i;
  int j;
  int num_procs;
  int step;
  double *swap;
  double my_change;
  int my_n;
  int n;
  
/*
  MPI initialization.
*/
  MPI_Init ( &argc, &argv );

  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

  MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

  allocate_arrays ( );
  f = make_rhs ( );
  make_domains ( num_procs );

  step = 0;

/*
  Begin iteration.
*/
  do 
  {
    jacobi ( num_procs, f );
    ++step;
/* 
  Estimate the error 
*/
    change = 0.0;
    n = 0;

    my_change = 0.0;
    my_n = 0;

    for ( i = i_min[my_rank]; i <= i_max[my_rank]; i++ )
    {
      for ( j = 1; j <= N; j++ )
      {
        if ( u_new[INDEX(i,j)] != 0.0 ) 
        {
          my_change = my_change 
            + fabs ( 1.0 - u[INDEX(i,j)] / u_new[INDEX(i,j)] );

          my_n = my_n + 1;
        }
      }
    }
    MPI_Allreduce ( &my_change, &change, 1, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD );

    MPI_Allreduce ( &my_n, &n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    if ( n != 0 )
    {
      change = change / n;
    }
    if ( my_rank == 0 && ( step % 10 ) == 0 ) 
    {
      printf ( "  N = %d, n = %d, my_n = %d, Step %4d  Error = %g\n", N, n, my_n, step, change );
    }
/* 
  Interchange U and U_NEW.
*/
    swap = u;
    u = u_new;
    u_new = swap;
  } while ( epsilon < change );

  Form_Gradient(num_procs, f);

/* 
  Each process writes out a file the solution and gradient
*/
  Save_Data("solution",u_new,my_rank);

  Save_Data("gradientX",gradX,my_rank);

  Save_Data("gradientY",gradY,my_rank);

/*
  Terminate MPI.
*/
  MPI_Finalize ( );
/*
  Free memory.
*/
  free ( f );
 
  return 0;
}