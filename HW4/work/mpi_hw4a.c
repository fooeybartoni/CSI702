#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 1
#define B 1
#define NUMPARTS 1000
#define XYMAX 2.0
#define XYMIN -2.0
#define NUMGRID 4
#define MASTER 0

double hwFunc(double x, double y);

double Partial_Derivative_X(double x, double y);

double Partial_Derivative_X(double x, double y);

void Form_Gradient();

struct partStruct {
   double x,y,velX,velY;
   int gridNum;
};

typedef struct partStruct Particle;

double max_value;
double time_step;
double delta_x,delta_y;

//Particle xycoord[NUMPARTS];


//-------------------------------------------
void Make_Particles(Particle *xycoord,int cnt){
   double cxy=XYMIN;
   //Particle xycoord[NUMPARTS];
   int i,j=0;
   srand(time(NULL));

   double range = (XYMAX - XYMIN); 
   double div = RAND_MAX / range;
   
   for (i=0; i<cnt; i++){
      for (j=0; j<cnt; j++){
         xycoord[i].x = XYMIN + (rand() / div);
         xycoord[i].y = XYMIN + (rand() / div);
         xycoord[i].velX = 0;
         xycoord[i].velY = 0;
         xycoord[i].gridNum = 0;
      }
   }

   //return xycoord;

}

// Sort the structs based on quadrant assignment
int SortFunc(const void* a, const void* b) {
  
   int l = ((struct partStruct *)a)->gridNum;
   int r = ((struct partStruct *)b)->gridNum;
   return (l - r);

}

//function implementation
double hwFunc(double x, double y) {
   double z;
   z = A * x  + B * (x / (pow(x,2) + pow(y,2)));
   
   return z;
}

double Partial_Derivative_X(double x, double y) {
   double H = (XYMAX - XYMIN)/NUMPARTS;
   double z;
   z = (hwFunc(x + H, y) - hwFunc(x-H,y)) / (2.0 * H); 
   return z;
}

double Partial_Derivative_Y(double x, double y) {
   double K = (XYMAX - XYMIN)/NUMPARTS;
   double z;
   z = (hwFunc(x, y + K) - hwFunc(x,y-K)) / (2.0 * K); 
   return z;
}

void Form_Gradient(Particle *xycoord,int cnt) {

   int i,j,k;
   i=0;
   j=0;
   k=0;
   double xGrad,yGrad,small_i,small_j;
   
   max_value=0.0;
   double magnitude=0.0;

   for (i=0; i<cnt; i++) {
      
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

void Save_Data(Particle *xycoord,int cnt) {

   FILE * fpVelocity;
   FILE * fpTimeStep;
   
   int i,j;
   i=0;
   j=0;

   if((fpVelocity=fopen("velocity.out", "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   if((fpTimeStep=fopen("time_step.out", "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   fprintf(fpTimeStep,"%f",time_step);

   for (i=0; i<NUMPARTS; i++) {
          
      fprintf(fpVelocity,"%f,",xycoord[i].x);
      fprintf(fpVelocity,"%f,",xycoord[i].y);
      fprintf(fpVelocity,"%f,",xycoord[i].velX);
      fprintf(fpVelocity,"%f",xycoord[i].velY);
      
      if (i != NUMPARTS-1) {
         fprintf(fpVelocity,"\n");         
      }
   }
   
   fclose(fpVelocity);
   fclose(fpTimeStep);
}

MPI_Datatype create_particle_datatype()
   {
      MPI_Datatype particle_type;
       MPI_Datatype array_of_types[2];
       int array_of_blocklengths[2];
       MPI_Aint array_of_displacements[2],lb,extent;
       MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);

       /* Create MPI Datatype for Particle struct*/
       array_of_types[0] = MPI_DOUBLE;
       array_of_blocklengths[0] = 4;
       array_of_displacements[0] = 0;
       
       array_of_types[1] = MPI_INT;
       array_of_blocklengths[1] = 1;
       array_of_displacements[1] = sizeof(Particle);

       MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
           array_of_types, &particle_type);
       MPI_Type_commit(&particle_type);
       return particle_type;
   }
double xSend[NUMPARTS], ySend[NUMPARTS];
int main(int argc, char *argv[]) {

int offset0, offset1, offset2,offset3;
   int   numtasks, taskid, rc, dest, i, j, tag1,
      tag2, tag3, source, chunksize, cnt; 
   MPI_Status status;
   
   int p0Cnt,p1Cnt,p2Cnt,p3Cnt;
   //p0Cnt=0;
   //p1Cnt=1;
   //p2Cnt=2;
   //p3Cnt=3;

   //int offset0, offset1, offset2,offset3;
   
   

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   if ( NUMPARTS % numtasks != 0) {
      printf("Quitting. Number of MPI tasks divide evenly into 1000.\n");
      MPI_Abort(MPI_COMM_WORLD, rc);
      exit(0);
   }
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   printf ("MPI task %d has started...\n", taskid);
   chunksize = (NUMPARTS / numtasks);
   tag2 = 1;
   tag1 = 2;
   tag3 = 3;

   /***** Master task only ******/
   if (taskid == MASTER){

      /* Initialize the array 
         Based on the x:-2 to 2 and y: -2 to 2 boundary
         and 4 partitions of the data n=4 processes 
         send each of the quadrants their data points to process.
         Quadrants are from upper left clockwise 0,1,2,3
      */
      MPI_Datatype particle_type = create_particle_datatype();
      //int offset0, offset1, offset2,offset3;

      Particle xycoord[NUMPARTS] = {0};
      Make_Particles(xycoord,NUMPARTS);
      i=0;
      p0Cnt=0;
      p1Cnt=0;
      p2Cnt=0;
      p3Cnt=0;
      //double xSend[NUMPARTS], ySend[NUMPARTS];
      for(i=0; i<NUMPARTS; i++) {
         xSend[i] = 0.0;
         ySend[i] = 0.0;
      }
      
      for (i=0; i<NUMPARTS; i++) {
         if (xycoord[i].x < 0) { // left half of plot
            if (xycoord[i].y < 0) {// Quadrant 3
               // Q3
               xycoord[i].gridNum = 3;
               xSend[p3Cnt] = xycoord[i].x;
               ySend[p3Cnt] = xycoord[i].y;
               p3Cnt=p3Cnt + 1;
               
            }
            else if (xycoord[i].y >= 0) { 
               // Q0
               xycoord[i].gridNum = 0;
               //xSend[0][p3Cnt] = xycoord[i].x;
               //ySend[0][p3Cnt] = xycoord[i].y;
               p0Cnt=p0Cnt + 1;
            }
            else
               printf("Error Categorizing");
         }
         else { //right half of plot
            if (xycoord[i].y < 0) {// Quadrant 2
               // Q2
               xycoord[i].gridNum = 2;
               //xSend[2][p3Cnt] = xycoord[i].x;
               //ySend[2][p3Cnt] = xycoord[i].y;
               p2Cnt=p2Cnt +1;
            }
            else if (xycoord[i].y >= 0) { //Quadrant 1
               // Q1
               xycoord[i].gridNum = 1;
               //xSend[2][p3Cnt] = xycoord[i].x;
               //ySend[2][p3Cnt] = xycoord[i].y;
               p1Cnt=p1Cnt + 1;
            }
            else
               printf("Error Categorizing");

         }
      }
      int total = p0Cnt+p1Cnt+p2Cnt+p3Cnt;
      printf("Initialized array cnts of %d,%d,%d,%d, sum: %d\n",p0Cnt,p1Cnt,p2Cnt,p3Cnt,total);
      qsort(xycoord, NUMPARTS, sizeof(xycoord[0]), SortFunc);
      
      int dest=1;
      offset0=0;
      offset1=p0Cnt;
      offset2=offset1+p1Cnt;
      offset3=offset2+p2Cnt;
      chunksize=p3Cnt;
      

      for (i=0;i<NUMPARTS;i++) {
         printf("x is %f \ty is %f\n",xSend[i],ySend[i]);



      }

printf("array size is %d\n",chunksize);


      /* Send each task its portion of the array - master keeps 1st part */
      // NEED TO REFACTOR TO CREATE IN LOOP  
      // Send P1 set
      MPI_Send(&offset1, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);

      MPI_Send(&xSend[offset1], chunksize, MPI_FLOAT, dest, tag2, MPI_COMM_WORLD);
      //MPI_Send(&ySend[offset1], chunksize, MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD);
      //MPI_Send(&xycoord[offset1], chunksize, particle_type, dest, tag2, MPI_COMM_WORLD);
      printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset1);

/*
      // Send P2 set
      dest=2;
      chunksize = p2Cnt;
      MPI_Send(&offset2, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
      //MPI_Send(&xycoord[offset2], chunksize, particle_type, dest, tag2, MPI_COMM_WORLD);
      printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset2);

      // Send P3 set
      dest=3;
      chunksize = p3Cnt;
      MPI_Send(&offset3, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
      //MPI_Send(&xycoord[offset3], chunksize, particle_type, dest, tag2, MPI_COMM_WORLD);
      printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset3);
       */         
   }


   if (taskid > 0) {

      MPI_Datatype particle_type = create_particle_datatype();
      int cntRows,chunksize;
      //double xCatch[NUMPARTS];
      //double yCatch[NUMPARTS];
      /* Receive my portion of array from the master task */
      source = MASTER;

      
      MPI_Recv(&offset1, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      //MPI_Recv(&xyCatch, 1, particle_type, source, tag2, 
      //   MPI_COMM_WORLD, &status);

      MPI_Recv(&xSend[offset1], chunksize, MPI_FLOAT, source, tag2, MPI_COMM_WORLD, &status);
      //MPI_Recv(&ySend[offset1], chunksize, MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
      
      //MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
      //printf("In task %d, cntRows is %d and have data value %f\n",taskid,chunksize,xCatch[0]);
      
   }
/*
   if (taskid == 2) {

       MPI_Datatype particle_type = create_particle_datatype();

       Particle xyCatch[NUMPARTS];
      //Receive my portion of array from the master task 
      source = MASTER;
      
      MPI_Recv(&offset2, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      //MPI_Recv(&xyCatch, 1, particle_type, source, tag2, 
      //   MPI_COMM_WORLD, &status);
      printf("offset is %d\n",offset2);
      MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
      printf("In task %d and have data value %f\n",taskid,xyCatch[0].x);
      
      
      printf("taskid is %d\n",taskid);

   }

   

   if (taskid == 3) {

       MPI_Datatype particle_type = create_particle_datatype();

      Particle xyCatch[NUMPARTS];
      // Receive my portion of array from the master task 
      source = MASTER;
      
      MPI_Recv(&offset3, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      //MPI_Recv(&xyCatch, 1, particle_type, source, tag2, 
      //   MPI_COMM_WORLD, &status);
      printf("offset is %d",offset3);
      MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
      printf("In task %d and have data value %f",taskid,xyCatch[0].x);
   } 
*/
    MPI_Finalize();      

   //Make_Particles();
   //Form_Gradient();
   //Save_Data();
}



  
