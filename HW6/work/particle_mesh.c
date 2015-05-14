# include "mpi.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#define  ARRAYSIZE   400
#define  MASTER      0
#define XYMAX 33.0
#define XYMIN 0.0
#define A 1
#define B 1
#define NUMPROC 4
#define MINCHG 0.0
#define MAXCHG 1.0
#define DUMMY 999.0

double  dataX[ARRAYSIZE];
double  dataY[ARRAYSIZE];

double gridPtChg[ARRAYSIZE];


struct partStruct {
   double x,y,fX,fY,charge,velX,velY;
};

typedef struct partStruct Particle;

MPI_Datatype particle_type;

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
void Save_Particle_Data(char* name,Particle p[][ARRAYSIZE],int myTaskId,int loop);
void Save_Data(char* name,double grad[],int myTaskId);
void make_particles(int num_procs,int my_rank);
void calc_grid_charges(int num_procs, int my_rank);
void calc_forces(int num_procs);
void Form_E_Field(int num_procs, double u_new[]);
void find_velocity(int my_rank,int proc_num,int loop);
void do_Poissons(int num_procs);

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
    
    if (i<ARRAYSIZE/4) {
      parts[my_rank][i].x = XYMIN + (my_rank * skip) + (rand() / div )/num_procs;
      parts[my_rank][i].y = XYMIN + (rand() / div );
      parts[my_rank][i].fX = 0.1;
      parts[my_rank][i].fY = 0.1;
      parts[my_rank][i].charge = MINCHG + (rand() / (RAND_MAX + 1.0));
      parts[my_rank][i].velX = 0.1;
      parts[my_rank][i].velY = 0.1;
      //if (my_rank == 0 && i%10 == 0) {
      {  printf("%d --- rank %d %f, %f, %f, %f, %f\n",i,my_rank,
          parts[my_rank][i].x,
          parts[my_rank][i].y,
          parts[my_rank][i].fX,
          parts[my_rank][i].fY,
          parts[my_rank][i].charge);
      }
    }
    else {
      parts[my_rank][i].x = DUMMY;
      parts[my_rank][i].y = DUMMY;
      parts[my_rank][i].fX = 0.1;
      parts[my_rank][i].fY = 0.1;
      parts[my_rank][i].charge = MINCHG + (rand() / (RAND_MAX + 1.0));
      parts[my_rank][i].velX = 0.1;
      parts[my_rank][i].velY = 0.1;
      //if (my_rank == 0 && i%10 == 0) {
      {  printf("%d --- rank %d %f, %f, %f, %f, %f\n",i,my_rank,
          parts[my_rank][i].x,
          parts[my_rank][i].y,
          parts[my_rank][i].fX,
          parts[my_rank][i].fY,
          parts[my_rank][i].charge);
      }



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
    if (parts[my_rank][k].x < DUMMY) {
      dblpx = (parts[my_rank][k].x/h)+1;
      dblpy = (parts[my_rank][k].y/h)+1;
      px = (int)dblpx;
      py = (int)dblpy;
      rho[INDEX(px,py)] += parts[my_rank][k].charge/(h*h);
      //printf("h is %f, dbl of x is %f, dbl of y is %f\n",h,dblpx,dblpy);
      //printf("x is %d, y is %d, rho is %f\n",px,py,rho[INDEX(px,py)]);
    }
  }
 
  return;
}

void calc_forces(int num_procs) {
  

//Save_Data("eFieldX_test",eFieldX,my_rank);

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
    if (parts[my_rank][k].x < DUMMY) {
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
  }
 
  return;
}

void find_velocity(int my_rank,int proc_num,int loop) {
/* Given that the particle array is filled with forces we can use Verlets to 
    determine the velocity.
        px(:) = px(:) + vx(:)*dt + 0.5*Fx(:)*dt*dt
        py(:) = py(:) + vy(:)*dt + 0.5*Fy(:)*dt*dt
        vx(:) = vx(:) + 0.5*Fx(:)*dt
        vy(:) = vy(:) + 0.5*Fy(:)*dt
*/

  int i=0;
  double dt = 0.1; 

  for (i=0; i<ARRAYSIZE;i++) {
    if (parts[my_rank][i].x < DUMMY) {
      parts[my_rank][i].x = parts[my_rank][i].x + parts[my_rank][i].velX * dt + 
                                  0.5 * (parts[my_rank][i].fX * dt * dt);

      parts[my_rank][i].y = parts[my_rank][i].y + parts[my_rank][i].velX * dt + 
                                  0.5 * (parts[my_rank][i].fY * dt * dt);  



      parts[my_rank][i].velX = parts[my_rank][i].velX  + 
                                  0.5 * (parts[my_rank][i].fX * dt * dt);

      parts[my_rank][i].velY = parts[my_rank][i].velX  + 
                                  0.5 * (parts[my_rank][i].fY * dt * dt); 
    } 

  }

  Save_Particle_Data("nextParts",parts,my_rank,loop);

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

void Save_Particle_Data(char* name,Particle p[][ARRAYSIZE],int myTaskId, int loop) {
  FILE * partOut;
      
  int i;
  i=0;
   
  char buffer[16];
     
  sprintf(buffer, "loop%d%s%d.out", loop,name,myTaskId);
   

  if((partOut=fopen(buffer, "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
  }

  for ( i = 0; i < ARRAYSIZE; i++ )
  {
    // add in after debugging
    if (parts[my_rank][i].x < DUMMY) {
      fprintf(partOut,"%f, %f, %f, %f, %f, %f, %f\n",
        p[myTaskId][i].x,
        p[myTaskId][i].y,
        p[myTaskId][i].fX,
        p[myTaskId][i].fY,
        p[myTaskId][i].charge,
        p[myTaskId][i].velX,
        p[myTaskId][i].velY);
    }
  }
   
  fclose(partOut);
}

void do_Poissons(int num_procs) {

  double change;
  double epsilon = 1.0E-03;
  //double *f;
  char file_name[100];
  int i;
  int j;
  
  int step;
  double *swap;
  double my_change;
  int my_n;
  int n;

  step = 0;
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


}

MPI_Datatype create_particle_datatype()
   {
      MPI_Datatype array_of_types[1];
    int array_of_blocklengths[1];
    MPI_Aint array_of_displacements[1];

    /* Create MPI Datatype for Particle struct*/
    array_of_types[0] = MPI_DOUBLE;
    array_of_blocklengths[0] = 7;
    array_of_displacements[0] = 0;
    MPI_Type_create_struct(1, array_of_blocklengths, array_of_displacements,
        array_of_types, &particle_type);
    MPI_Type_commit(&particle_type);
   }

void move_particles(int num_procs, int my_rank) {

double pListX[ARRAYSIZE];
double pListY[ARRAYSIZE];
double pListVelX[ARRAYSIZE];
double pListVelY[ARRAYSIZE];
double pListCharge[ARRAYSIZE];

double rListX[ARRAYSIZE];
double rListY[ARRAYSIZE];
double rListVelX[ARRAYSIZE];
double rListVelY[ARRAYSIZE];
double rListCharge[ARRAYSIZE];




  int i=0;
  int foo = 0;
  int bar = 0;
  double range = (XYMAX - 1 - XYMIN); 
  double skip = range/num_procs;
  double px, rpx = 0.0;

  int lFlag = 0;
  int rFlag = 0;

  int count = 1;
  int rcount = 1;
  int tag1 = 1;
  int tag2 = 123;
  int source = 0;
  int left;
  int right;
  
  MPI_Status status;
  MPI_Request request[16];
  int requests;
  //MPI_Status status[16];

  // clean and initialize transfer array
  for (i=0;i<ARRAYSIZE * 2;i++) {
    pListX[i] = 999.0;
    pListY[i] = 999.0;
  }
//printf("in part move rank %d",my_rank);
  //if (my_rank ==0) {
    for (i=0;i<ARRAYSIZE;i++) {
      px = parts[my_rank][i].x;
      //first check if it not out of bounds
      if (px >= XYMIN && px <= XYMAX) {
        printf ("particle in bounds x is %f\n",px);
        if ( my_rank < num_procs-1 &&  px > ((my_rank + 1) * skip)) {
          // move particle to the right
          //rDest = my_rank + 1;
          rFlag = 1;
          printf("Inside move particle to the RIGHT rank %d, px %f\n",my_rank,px);
          pListX[i] = parts[my_rank][i].x;
          pListY[i] = parts[my_rank][i].y;
          pListVelX[i] = parts[my_rank][i].velX;
          pListVelY[i] = parts[my_rank][i].velY;
          pListCharge[i] = parts[my_rank][i].charge; 


          //Remove from local
          parts[my_rank][i].x = DUMMY;
          parts[my_rank][i].y = DUMMY;
        }
        else if (my_rank > 0 &&  px < (XYMIN + (my_rank * skip))) {
          //move particle to the left
          printf("Inside move particle to the LEFT\n");
          //lDest = my_rank - 1;
          printf("Inside move particle to the RIGHT");
          pListX[i] = parts[my_rank][i].x;
          pListY[i] = parts[my_rank][i].y;
          pListVelX[i] = parts[my_rank][i].velX;
          pListVelY[i] = parts[my_rank][i].velY;
          pListCharge[i] = parts[my_rank][i].charge; 


          //Remove from local
          parts[my_rank][i].x = DUMMY;
          parts[my_rank][i].y = DUMMY;
        }
        else {
          // Value is good take no action
          // Nothing to do as array is initialized with dummy data          
        }
        
      }
      else {
        // Particle is out of bounds so just erase it
        parts[my_rank][i].x = DUMMY;
        parts[my_rank][i].y = DUMMY;
         
      }

    }
    
      // Send the ARRAY of particles over to other domains
      pListX[0] = 567.88;
 


      left = my_rank - 1;    /* left neighbour */
      right = my_rank + 1;   /* right neighbour */

      /*
       * here we make a decision about what to do at the ends of the chain
       * do not send let it escape
       */

      if (left < 0) left = MPI_PROC_NULL;
      if (right > num_procs - 1) right = MPI_PROC_NULL;
   
   
              //printf("receiving Shit from LEFT rank %d\n",my_rank);
        //MPI_Sendrecv(pListX, ARRAYSIZE, MPI_DOUBLE, right, tag1,
         //            rListX, ARRAYSIZE, MPI_DOUBLE, left, tag2, 
        //  MPI_COMM_WORLD, &status);

      //MPI_Sendrecv(&lFlag, 1, MPI_INT, right, tag1,
       //              &rFlag, 1, MPI_INT, left, tag2, 
       //   MPI_COMM_WORLD, &status);
    
        printf ("Process %d received %f from %d and sent %f to %d\n", 
          my_rank, lFlag, left, rFlag, right);
        


        //printf("sending Shit RIGHT rank %d\n",my_rank);
        //MPI_Isend(&pListX[0], ARRAYSIZE, MPI_DOUBLE, rDest, tag1, MPI_COMM_WORLD,request + requests++);
          

      
      // Send to the LEFT
      //if (my_rank > 0 ) {
        //printf("receiving Shit from RIGHT rank %d\n",my_rank);
        
        //MPI_Irecv(&pListX[100], ARRAYSIZE, MPI_DOUBLE, my_rank+1, tag2, 
          //MPI_COMM_WORLD, request + requests++);
        
       // printf("sending Shit LEFT rank %d\n",my_rank);
        //MPI_Isend(&pListX[0], ARRAYSIZE, MPI_DOUBLE, lDest, tag2, MPI_COMM_WORLD,request + requests++);
      //}
  //}

 
    printf("requests %d",requests);
    //MPI_Waitall ( requests, request, status );
    //MPI_Wait(&request,&status);
//MPI_Barrier(MPI_COMM_WORLD);
    //printf("Holy Shit rank %d I got a particle x is: %f\n",my_rank,pListX[100]);
  //}
}

int main ( int argc, char *argv[] ) 

{
  
  int num_procs;

  int i;
  
/*
  MPI initialization.
*/
  MPI_Init ( &argc, &argv );

  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

  MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

  

  allocate_arrays ( );
  rho = make_rho ( );
  make_domains ( num_procs );


  printf("my rank is %d\n", my_rank);
  make_particles(num_procs,my_rank);
  printf("made it past make particles\n");
  
  Save_Particle_Data("particles",parts,my_rank,0);

  for (i=1;i<12;i++) {
    calc_grid_charges(num_procs, my_rank);
    
    do_Poissons(num_procs);

    Form_E_Field(num_procs,u_new);

    calc_forces(num_procs);

    

    find_velocity(my_rank,num_procs,i);

  //move_particles(num_procs,my_rank);
  }

/* 
  Each process writes out a file the solution and gradient
*/
  Save_Data("solution",u_new,my_rank);

  Save_Data("eFieldX",eFieldX,my_rank);

  Save_Data("eFieldY",eFieldY,my_rank);

  Save_Data("voltages",rho,my_rank);

  
  

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