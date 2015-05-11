#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define  ARRAYSIZE   100
#define  MASTER      0
#define XYMAX 2.0
#define XYMIN -2.0
#define A 1
#define B 1
#define NUMGRID 4

double  dataX[ARRAYSIZE];
double  dataY[ARRAYSIZE];

double pListX[4][ARRAYSIZE];
double pListY[4][ARRAYSIZE];



struct partStruct {
   double x,y,velX,velY;
};

typedef struct partStruct Particle;

Particle answer[ARRAYSIZE];

// These may need to be not global
double max_value;
double time_step;
double delta_x,delta_y;

   //function implementation
double hwFunc(double x, double y) {
   double z;
   z = A * x  + B * (x / (pow(x,2) + pow(y,2)));
   
   return z;
}

double Partial_Derivative_X(double x, double y) {
   double H = (XYMAX - XYMIN)/ARRAYSIZE;
   double z;
   z = (hwFunc(x + H, y) - hwFunc(x-H,y)) / (2.0 * H); 
   return z;
}

double Partial_Derivative_Y(double x, double y) {
   double K = (XYMAX - XYMIN)/ARRAYSIZE;
   double z;
   z = (hwFunc(x, y + K) - hwFunc(x,y-K)) / (2.0 * K); 
   return z;
}

void Form_Gradient(Particle xycoord[], size_t size) {

   int i,j,k;
   i=0;
   j=0;
   k=0;
   double xGrad,yGrad,small_i,small_j;
   
   max_value=0.0;
   double magnitude=0.0;

   for (i=0; i<ARRAYSIZE; i++) {
      
      double px = xycoord[i].x;
      double py = xycoord[i].y;
      xGrad=Partial_Derivative_X(px,py);
      xycoord[i].velX = xGrad;
      yGrad=Partial_Derivative_Y(px,py);
      xycoord[i].velY = yGrad;
      magnitude = sqrt(pow(xGrad,2)+pow(yGrad,2));
      //printf("X: %f Y: %f GradX: %f  GradY: %f  Magnitude: %f\n",px,py,xGrad,yGrad,magnitude);
         
      if (fabs(magnitude) > max_value)
      {
         max_value = fabs(magnitude);
      }
      delta_x = delta_y = ((XYMAX - XYMIN)/NUMGRID);  // This is the 100 x 100 grid from HW #3
      int CFL =1;
      time_step = delta_x/max_value; 
   }

   printf("Maximum Vector Magnitude is: %f\n",max_value);
   printf("Delta_x is: %f, Time step is: %f\n",delta_x,time_step);

  
}

void Save_Data(Particle xycoord[],int myTaskId) {

   FILE * fpVelocity;
   FILE * fpTimeStep;
   
   int i,j;
   i=0;
   j=0;

   char buffer[16];
   
   
   sprintf(buffer, "velocity%d.out", myTaskId);
   
   

   if((fpVelocity=fopen(buffer, "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   if((fpTimeStep=fopen("time_step.out", "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   fprintf(fpTimeStep,"%f",time_step);

   for (i=0; i<ARRAYSIZE; i++) {
      if (xycoord[i].x < 90) {   
         fprintf(fpVelocity,"%f,",xycoord[i].x);
         fprintf(fpVelocity,"%f,",xycoord[i].y);
         fprintf(fpVelocity,"%f,",xycoord[i].velX);
         fprintf(fpVelocity,"%f",xycoord[i].velY);
         
         if (i != ARRAYSIZE-1) {
            fprintf(fpVelocity,"\n");         
         }
      }
   }
   
   fclose(fpVelocity);
   fclose(fpTimeStep);
}


// I would GATHER back using this struct MPI Datatype
MPI_Datatype create_particle_datatype()
   {
      MPI_Datatype particle_type;
       MPI_Datatype array_of_types[1];
       int array_of_blocklengths[1];
       MPI_Aint array_of_displacements[1],lb,extent;
       MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);

       /* Create MPI Datatype for Particle struct*/
       array_of_types[0] = MPI_DOUBLE;
       array_of_blocklengths[0] = 4;
       array_of_displacements[0] = 0;
       
       MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
           array_of_types, &particle_type);
       MPI_Type_commit(&particle_type);
       return particle_type;
   }

/* ########################################################   

      Below is the MAIN section including the MPI Portion

   ########################################################
*/

   int main(int argc, char *argv[]) {

   int   numtasks, taskid, rc, dest, offset, i, j, tag1,
         tag2, tag3, source, chunksize; 
   
   MPI_Status status;

   /***** Initializations *****/
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   if (numtasks % 4 != 0) {
      printf("Quitting. Number of MPI tasks must be divisible by 4.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(0);
   }
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   printf ("MPI task %d has started...\n", taskid);
   chunksize = (ARRAYSIZE / numtasks);
   tag3 = 1;
   tag2 = 2;
   tag1 = 3;
   /***** Master task only ******/
   if (taskid == MASTER){

      srand(time(NULL));

      double range = (XYMAX - XYMIN); 
      double div = RAND_MAX / range;
      
      for (i=0; i<ARRAYSIZE; i++){
         for (j=0; j<ARRAYSIZE; j++){
            dataX[i] = XYMIN + (rand() / div);
            dataY[i] = XYMIN + (rand() / div);
         }
      }
      
      for (i=0; i<4; i++) {
         for (j=0; j<ARRAYSIZE; j++) {
            pListX[i][j] = 99.0;
            pListY[i][j] = 99.0;
         }
      }

      int inCount1, inCount2, inCount3, inCount0 = 0;

      double x,y;
      for (i=0; i<ARRAYSIZE; i++) {
         x = dataX[i]; y=dataY[i];
         if((x>=0) &&(y>=0))
         {
            //printf("\n The co-ordinate (%f,%f) is in first Quardrant\n",x,y);
            pListX[0][i]=x;
            pListY[0][i]=y;
         }
         else if((x<0) &&(y<0))
         {
            //printf("\n The co-ordinate (%f,%f) is in Third Quardrant\n",x,y);
            pListX[2][i]=x;
            pListY[2][i]=y;
         }
         else if((x>0) &&(y<0))
         {
            //printf("\n The co-ordinate (%f,%f) is in Fourth Quardrant\n",x,y);
            pListX[3][i]=x;
            pListY[3][i]=y;
         }
         else {
            //printf("\n The co-ordinate (%f,%f) is in Second Quardrant\n",x,y);
            pListX[1][i]=x;
            pListY[1][i]=y;
         }
      }

      printf("Initialized array %f\n",pListX[0][0]);

      /* Send each task its portion of the array - master keeps 1st part */
      offset = chunksize;

      for (dest=1; dest<numtasks; dest++) {

         MPI_Send(&pListX[dest][0], ARRAYSIZE, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
         MPI_Send(&pListY[dest][0], ARRAYSIZE, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
         
         printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset);
         offset = offset + chunksize;
      }

      for (i=0; i < ARRAYSIZE; i++) {
         answer[i].x = pListX[taskid][i];
         answer[i].y = pListY[taskid][i];
         answer[i].velX = 0.0;
         answer[i].velY = 0.0;
      }

      // Now Master has to do his work
      Form_Gradient(answer,ARRAYSIZE);
      //answer = Form_Gradient(answer);

      printf("Answer x: %f, velX: %f\n", answer[0].x,answer[0].velX);

      Save_Data(answer,taskid);

      // I will save the Gather functionality for HW6 but can put it in if I 
      // have time.  I am happy with each processor calling the Gradient function
      // for each of their respective quadrants and servicing their assigned 
      // aka not '99.00' sets of coordinates.
   }

   if (taskid > MASTER) {

      int i = 0;

      /* Receive my portion of array from the master task */
      source = MASTER;
      
      MPI_Recv(&pListX[taskid][0], ARRAYSIZE, MPI_DOUBLE, source, tag2, 
      MPI_COMM_WORLD, &status);
      MPI_Recv(&pListY[taskid][0], ARRAYSIZE, MPI_DOUBLE, source, tag2, 
      MPI_COMM_WORLD, &status);
      
      //printf("hello I am in %d dataX item is %f\n",taskid,pListX[taskid][0]);
      //printf("hello I am in %d dataY item is %f\n",taskid,pListY[taskid][0]);

      // Initialize the Particle Array
      for (i=0; i < ARRAYSIZE; i++) {
         answer[i].x = pListX[taskid][i];
         answer[i].y = pListY[taskid][i];
         answer[i].velX = 0.0;
         answer[i].velY = 0.0;
      }

      Form_Gradient(answer,ARRAYSIZE);
      //answer = Form_Gradient(answer);

      printf("Answer x: %f, velX: %f\n", answer[0].x,answer[0].velX);

      Save_Data(answer,taskid);
      

      
   }

   MPI_Finalize();

}   /* end of main */