#include <stdio.h>
#include <math.h>

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

typedef struct {
    double gradX;
    double gradY;
} grad;



double gradientX[numPts][numPts];
double gradientY[numPts][numPts];

grad actual[100][100];

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
   cxy=XYMIN;
   double xycoord[numPts];
   double step = (XYMAX - XYMIN)/numPts;

   for (k=0; k<numPts; k++){
      xycoord[k]=cxy;
      cxy = cxy + step;
   }

   for (i=0; i<numPts; i++) {

      for(j=0; j<numPts; j++) {
         xGrad=Partial_Derivative_X(xycoord[i],xycoord[j]);
         gradientX[i][j]=xGrad;
         yGrad=Partial_Derivative_Y(xycoord[i],xycoord[j]);
         gradientY[i][j]=yGrad;
         //printf("%f,%f,%f,%f\n",xycoord[i],xycoord[j],x,y);
      }
   }

   FILE * fpX;
   FILE * fpY;
   FILE * fpCoords;
/*
   if((fp=fopen("ax.out", "wb"))==NULL) {
    printf("Cannot open file.\n");
   }

   if(fwrite(gradientX, sizeof(double), numPts*numPts, fp) != numPts*numPts)
    printf("File write error.");
   fclose(fp);
*/

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
      fprintf(fpX,"\n"); 
      fprintf(fpY,"\n");
      fprintf(fpCoords, "\n");
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

         actual[i][j].gradX=hwFunc(i,j);
         actual[i][j].gradY=hwFunc(i,j);
      }
   }

}


int main() {
   Form_Gradient();
}
