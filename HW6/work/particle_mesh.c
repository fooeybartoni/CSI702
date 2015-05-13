# include "mpi.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#define  ARRAYSIZE   100
#define  MASTER      0
#define XYMAX 33.0
#define XYMIN 0.0
#define A 1
#define B 1
#define NUMPROC 4
#define MINCHG 0.0
#define MAXCHG 1.0

double  dataX[ARRAYSIZE];
double  dataY[ARRAYSIZE];

double gridPtChg[ARRAYSIZE];

double pListX[4][ARRAYSIZE];
double pListY[4][ARRAYSIZE];

struct partStruct {
   double x,y,fX,fY,charge;
};

typedef struct partStruct Particle;



Particle parts[NUMPROC][ARRAYSIZE];

// These may need to be not global
double max_value;
double time_step;
double delta_x,delta_y;


double L = XYMAX - XYMIN;			
int N = 32;			

double *u, *u_new;

double *eFieldX, *eFieldY;

double *rho;	

#define INDEX(i,j) ((N+2)*(i)+(j))

int my_rank;    

int *proc;      
int *i_min, *i_max;   
int *left_proc, *right_proc;  

int main ( int argc, char *argv[] );
void allocate_arrays ( );
void jacobi ( int num_procs, double rho[] );
void make_domains ( int num_procs );
double *make_rho ( );
void Save_Particle_Data(char* name,Particle p[][ARRAYSIZE],int myTaskId);
void Save_Data(char* name,double grad[],int myTaskId);
void make_particles(int num_procs,int my_rank);
void calc_grid_charges(int num_procs, int my_rank);
void calc_forces(int num_procs);
void Form_E_Field(int num_procs, double u_new[]);


/*
    Functions
*/

void make_particles(int num_procs,int my_rank) {
  int i=0;
  
/*
  H is the lattice spacing.
*/
  //h = L / ( double ) ( N + 1 ); 

printf("made it into make_particles\n");

  srand(time(NULL) * rand() * (rand()/RAND_MAX+my_rank));

  double range = (XYMAX - 1 - XYMIN); 
  double div = (RAND_MAX / range);
  double skip = range/num_procs;


  for (i=0; i<ARRAYSIZE; i++){
    
    parts[my_rank][i].x = XYMIN + (my_rank * skip) + (rand() / div )/num_procs;
    parts[my_rank][i].y = XYMIN + (rand() / div );
    parts[my_rank][i].fX = 0.1;
    parts[my_rank][i].fY = 0.1;
    parts[my_rank][i].charge = MINCHG + (rand() / (RAND_MAX + 1.0));
    //if (my_rank == 0 && i%10 == 0) {
    {  printf("%d --- rank %d %f, %f, %f, %f, %f\n",i,my_rank,
        parts[my_rank][i].x,
        parts[my_rank][i].y,
        parts[my_rank][i].fX,
        parts[my_rank][i].fY,
        parts[my_rank][i].charge);
    }

  }

  printf("just before return in make_particles\n");
  return;
}


void calc_grid_charges(int num_procs, int my_rank) {
  
  printf("in the calc grid method\n");

  double dblpx, dblpy,h=0.0;
  int px,py=0;
  int i,j,k;
  
  MPI_Request request[4];
  int requests;
  MPI_Status status[4];
/*
  H is the lattice spacing.
*/
  h = L / ( double ) ( N + 1 );

/* 
  Charge calculations for internal vertices in my domain.
*/
  int count = 0;
  //printf("before the loop in calc_grid_charges\n");
  
  // Find all of the charged particles around the point
  // iterate over the array
  //   - calculate partial charge using distance from 
  //    ratio = h - dist from grid pt / grid pt diagonal dist
  //    * charge
  // 

  for (k=0; k<ARRAYSIZE; k++) {
    dblpx = (parts[my_rank][k].x/h)+1;
    dblpy = (parts[my_rank][k].y/h)+1;
    px = (int)dblpx;
    py = (int)dblpy;
    rho[INDEX(px,py)] += parts[my_rank][k].charge/(h*h);
    //printf("h is %f, dbl of x is %f, dbl of y is %f\n",h,dblpx,dblpy);
    //printf("x is %d, y is %d, rho is %f\n",px,py,rho[INDEX(px,py)]);
   
  }
 
  return;
}

void calc_forces(int num_procs) {
  

Save_Data("eFieldX_test",eFieldX,my_rank);

  printf("in the calc forces method\n");
  /* Calculate the force given the charge and the Electric field
         F = qE
  */

//int my_rank =0;

  double dblpx, dblpy,h=0.0;
  int px,py=0;
  int i,j,k;
  double q, eX, eY = 0.0;
  
  MPI_Request request[4];
  int requests;
  MPI_Status status[4];
/*
  H is the lattice spacing.
*/
  h = L / ( double ) ( N + 1 );

/* 
  Charge calculations for internal vertices in my domain.
*/
  int count = 0;
  printf("before the loop in calc_forces\n");
  
  // Find all of the charged particles around the point
  // iterate over the array
  //   - calculate partial charge using distance from 
  //    ratio = h - dist from grid pt / grid pt diagonal dist
  //    * charge
  // 

  for (k=0; k<ARRAYSIZE; k++) {
    dblpx = ((parts[my_rank][k].x)/h)+1;
    dblpy = ((parts[my_rank][k].y)/h)+1;
    px = (int)dblpx;
    py = (int)dblpy;
    
    q = parts[my_rank][k].charge;

    eX = eFieldX[INDEX(px,py)];
    eY = eFieldY[INDEX(px,py)];
    printf("my calc_forces rank is %d",my_rank);
    printf("q is %f, eX is %f, eY is %f\n",q, eX, eY);

    // This can be expanded to interpolation of four corner grid points 
    parts[my_rank][k].fX = q * eX;
    parts[my_rank][k].fY = q * eY;

    
    printf("rank is %d, x is %d, y is %d, fX is %f, fY is %f\n",my_rank,px,py,
      parts[my_rank][k].fX,parts[my_rank][k].fY );
  }
 
  return;
}

/* Function Definitions */
void Form_E_Field(int num_procs, double u_new[]) {
  /* The Electric Field is the negative of the Potential which is the
     Gradient from Homework #5
  */

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
    MPI_Irecv ( u_new + INDEX(i_min[my_rank] - 1, 1), N, MPI_DOUBLE,
      left_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u_new + INDEX(i_min[my_rank], 1), N, MPI_DOUBLE,
      left_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );
  }

  if ( right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u_new + INDEX(i_max[my_rank] + 1, 1), N, MPI_DOUBLE,
      right_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u_new + INDEX(i_max[my_rank], 1), N, MPI_DOUBLE,
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
      eFieldX[INDEX(i,j)] = -1.0 *
         ( (u_new[INDEX(i+1,j)] - u_new[INDEX(i-1,j)]) / h);

      eFieldY[INDEX(i,j)] = -1.0 *
         ( (u_new[INDEX(i,j+1)] - u_new[INDEX(i,j-1)]) / h);

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
    eFieldX[INDEX(i,j)] = -1.0 *
          ( (u_new[INDEX(i+1,j)] - u_new[INDEX(i-1,j)]) / h);

    eFieldY[INDEX(i,j)] = -1.0 *
          ( (u_new[INDEX(i,j+1)] - u_new[INDEX(i,j-1)]) / h);

  }

  i = i_max[my_rank];
  if (i != i_min[my_rank])
  {
    for (j = 1; j <= N; j++)
    {
      eFieldX[INDEX(i,j)] = -1.0 *
         ( (u_new[INDEX(i+1,j)] - u_new[INDEX(i-1,j)]) / h);

      eFieldY[INDEX(i,j)] = -1.0 *
         ( (u_new[INDEX(i,j+1)] - u_new[INDEX(i,j-1)]) / h);

    }
  }

  return;

  
}

void allocate_arrays ( ) 
{
  int i;
  int ndof;

  ndof = ( N + 2 ) * ( N + 2 );

  //rho = ( double * ) malloc ( ndof * sizeof ( double ) );
  //for ( i = 0; i < ndof; i++)
  //{
  //  rho[i] = 0.0;
  //}

  eFieldX = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++)
  {
    eFieldX[i] = 0.0;
  }

  eFieldY = ( double * ) malloc ( ndof * sizeof ( double ) );
  for ( i = 0; i < ndof; i++ )
  {
    eFieldY[i] = 0.0;
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
double *make_rho ( ) 
{

  // MADE ALL ZEROS to start off with
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
        //f[INDEX(i,j)] = 0.0;
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

void jacobi ( int num_procs, double rho[] ) 
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
                 h * h * rho[INDEX(i,j)] );
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
               h * h * rho[INDEX(i,j)] );
  }

  i = i_max[my_rank];
  if (i != i_min[my_rank])
  {
    for (j = 1; j <= N; j++)
    {
      u_new[INDEX(i,j)] =
        0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
                 u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
                 h * h * rho[INDEX(i,j)] );
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

void Save_Particle_Data(char* name,Particle p[][ARRAYSIZE],int myTaskId) {
  FILE * partOut;
      
  int i;
  i=0;
   
  char buffer[16];
     
  sprintf(buffer, "%s%d.out", name,myTaskId);
   

  if((partOut=fopen(buffer, "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
  }

  for ( i = 0; i < ARRAYSIZE; i++ )
  {
    fprintf(partOut,"%f, %f, %f, %f, %f\n",
      p[myTaskId][i].x,
      p[myTaskId][i].y,
      p[myTaskId][i].fX,
      p[myTaskId][i].fY,
      p[myTaskId][i].charge);
  }
   
  fclose(partOut);
}


MPI_Datatype create_particle_datatype()
   {
      MPI_Datatype particle_type;
      MPI_Type_contiguous (5,MPI_DOUBLE,&particle_type);
      MPI_Type_commit (&particle_type);
      return particle_type;
   }



int main ( int argc, char *argv[] ) 

{
  double change;
  double epsilon = 1.0E-03;
  //double *f;
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
  rho = make_rho ( );
  make_domains ( num_procs );


//printf("I got past prior inits\n");
  step = 0;  //what is this?

  //if (my_rank == 0) {
    printf("my rank is %d\n", my_rank);
    make_particles(num_procs,my_rank);
    printf("made it past make particles\n");
    calc_grid_charges(num_procs, my_rank);
    //exit(0);
  //}
  //make_particles(my_rank);
 // calc_grid_charges(num_procs,my_rank);

/*
  Begin iteration.
*/
  do 
  {
    jacobi ( num_procs, rho );
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

  Form_E_Field(num_procs,u_new);



  printf("just before calc_forces\n");

  calc_forces(num_procs);

/* 
  Each process writes out a file the solution and gradient
*/
  Save_Data("solution",u_new,my_rank);

  Save_Data("eFieldX",eFieldX,my_rank);

  Save_Data("eFieldY",eFieldY,my_rank);

  Save_Data("voltages",rho,my_rank);

  Save_Particle_Data("particles",parts,my_rank);
  

/*
  Terminate MPI.
*/
  MPI_Finalize ( );
/*
  Free memory.
*/
  //free ( f );
  
 
  return 0;
}