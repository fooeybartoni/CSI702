#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define  ARRAYSIZE   100
#define  MASTER      0
#define XYMAX 2.0
#define XYMIN -2.0

double  dataX[ARRAYSIZE];
double  dataY[ARRAYSIZE];

double pListX[4][ARRAYSIZE];
      double pListY[4][ARRAYSIZE];

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
/*
      for (i=0; i<ARRAYSIZE; i++) {
         if (dataX[i] < 0) { // left half of plot
            if (dataY[i] < 0) {// Quadrant 3
               // Q3
               printf("in ONE\n");
               pListX[3][inCount3] = dataX[i];
               pListY[3][inCount3] = dataY[i];
               inCount3 += 1;
            }
            else if (dataY[i] >= 0) { 
               // Q0
               printf("in TWO\n");
               pListX[0][inCount0] = dataX[i];
               pListY[0][inCount0] = dataY[i];
               inCount0 += 1;
            }
            else
               printf("Error Categorizing");
         }
         else { //right half of plot
            if (dataY[i] < 0) {// Quadrant 2
               // Q2
               printf("in Three\n");
               pListX[2][inCount2] = dataX[i];
               pListY[2][inCount2] = dataY[i];
               inCount2 += 1;
            }
            else if (dataY[i] >= 0) { //Quadrant 1
               // Q1
               printf("in FOUR\n");
               pListX[1][inCount1] = dataX[i];
               pListY[1][inCount1] = dataY[i];
               inCount1 += 1;
            }
            else
               printf("Error Categorizing");
         }
      }
*/
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
         


         
 //exit(0);
         MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
         MPI_Send(&pListX[dest][0], ARRAYSIZE, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
         //MPI_Send(&dataY[offset], chunksize, MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD);
         printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset);
         offset = offset + chunksize;

         
      }
   }

   if (taskid > MASTER) {

      /* Receive my portion of array from the master task */
      source = MASTER;
      MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      MPI_Recv(&pListX[taskid][0], ARRAYSIZE, MPI_DOUBLE, source, tag2, 
      MPI_COMM_WORLD, &status);
      //MPI_Recv(&dataY[offset], chunksize, MPI_DOUBLE, source, tag3, 
      //MPI_COMM_WORLD, &status);


      printf("hello I am in %d dataX item is %f\n",taskid,pListX[taskid][0]);
      //printf("hello I am in %d dataY item is %f\n",taskid,dataY[offset]);
   }

   MPI_Finalize();

}   /* end of main */