#include <stdio.h>
#include <math.h>
#include <omp.h>

#define A 1
#define B 1
#define numPts 100
#define XYMAX 2.0
#define XYMIN -2.0

double hwFunc(double x, double y);

double Partial_Derivative_X(double x, double y);

double Partial_Derivative_X(double x, double y);

void Form_Gradient();

void Form_Actual();

double gradientX[numPts][numPts];
double gradientY[numPts][numPts];
double xycoord[numPts];

double actual[numPts][numPts];

//-------------------------------------------
//function implementation
double hwFunc(double x, double y) {
   double z;
   z = A * x  + B * (x / (pow(x,2) + pow(y,2)));
   
   return z;
}

double Partial_Derivative_X(double x, double y) {
   double H = (XYMAX - XYMIN)/numPts;
   double z;
   z = (hwFunc(x + H, y) - hwFunc(x-H,y)) / (2.0 * H); 
   return z;
}

double Partial_Derivative_Y(double x, double y) {
   double K = (XYMAX - XYMIN)/numPts;
   double z;
   z = (hwFunc(x, y + K) - hwFunc(x,y-K)) / (2.0 * K); 
   return z;
}

void Form_Gradient() {

   int i,j,k;
   i=0;
   j=0;
   k=0;
   double cxy;
   double xGrad,yGrad,small_i,small_j;
   
   int target_thread_num = 4;
   omp_set_num_threads(target_thread_num);
   unsigned long times[target_thread_num];


   cxy=XYMIN;
   
   double step = (XYMAX - XYMIN)/numPts;

   for (k=0; k<numPts; k++){
      xycoord[k]=cxy;
      cxy = cxy + step;
   }

   double gX[numPts][numPts];
   double gY[numPts][numPts];

   // OMP Loop
   #pragma omp parallel for shared(gX,gY), private(i,j)
   for (i=0; i<numPts; i++) {
      for(j=0; j<numPts; j++) {

         //int thread_id = omp_get_thread_num();
         //printf("thread id: %d running",thread_id);

         xGrad=Partial_Derivative_X(xycoord[i],xycoord[j]);
         gX[i][j]=xGrad;
         yGrad=Partial_Derivative_Y(xycoord[i],xycoord[j]);
         gY[i][j]=yGrad;
      }
   }

   // Deep Copy Arrays Can be replaced by elegant pointer stuff that I am not up to right now
   for (i=0; i<numPts; i++) {
      for (j=0; j<numPts; j++) {
         gradientX[i][j] = gX[i][j];
         gradientY[i][j] = gY[i][j];
      }
   }
}

void Save_Data() {

   FILE * fpX;
   FILE * fpY;
   FILE * fpCoords;

   int i,j;
   i=0;
   j=0;

   if((fpX=fopen("ax.out", "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   if((fpY=fopen("ay.out", "w+"))==NULL) {
    printf("Cannot open file fpY.\n");
   }

   if((fpCoords=fopen("xyCoords.out", "w+"))==NULL) {
    printf("Cannot open file xyCoords.\n");
   }

   for (i=0; i<numPts; i++) {
      fprintf(fpCoords,"%f\t%f",xycoord[i],xycoord[i]);
      for (j=0; j<numPts; j++) {
         fprintf(fpX,"%f \t",gradientX[i][j]);
         fprintf(fpY,"%f \t",gradientY[i][j]);
         
      }
      if (i != numPts-1) {
         fprintf(fpX,"\n"); 
         fprintf(fpY,"\n");
         fprintf(fpCoords, "\n");
      }
   }
   fclose(fpX);
   fclose(fpY);
   fclose(fpCoords);
}

void Form_Actual() {
    int i,j;
   i=0;
   j=0;

   for (i=0; i<100; i++) {

      for (j=0; j<100; j++) {

         actual[i][j]=hwFunc(i,j);
      }
   }

}

int main() {
   Form_Gradient();
   Save_Data();
}
