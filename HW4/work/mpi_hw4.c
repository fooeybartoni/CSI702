#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 1
#define B 1
#define numParts 1000
#define XYMAX 2.0
#define XYMIN -2.0

double hwFunc(double x, double y);

double Partial_Derivative_X(double x, double y);

double Partial_Derivative_X(double x, double y);

void Form_Gradient();

typedef struct {
   double x,y,velX,velY;
} Particle;

//double gradientX[numParts][numParts];
//double gradientY[numParts][numParts];
Particle xycoord[numParts];

//-------------------------------------------
void Make_Particles(){
   double cxy=XYMIN;
   
   int i,j=0;
   srand(time(NULL));

   double range = (XYMAX - XYMIN); 
   double div = RAND_MAX / range;
   
   for (i=0; i<numParts; i++){
      for (j=0; j<numParts; j++){
         xycoord[k].x = XYMIN + (rand() / div);
         xycoord[k].y = XYMIN + (rand() / div);
      }
   }
}

//function implementation
double hwFunc(double x, double y) {
   double z;
   z = A * x  + B * (x / (pow(x,2) + pow(y,2)));
   
   return z;
}

double Partial_Derivative_X(double x, double y) {
   double H = (XYMAX - XYMIN)/numParts;
   double z;
   z = (hwFunc(x + H, y) - hwFunc(x-H,y)) / (2.0 * H); 
   return z;
}

double Partial_Derivative_Y(double x, double y) {
   double K = (XYMAX - XYMIN)/numParts;
   double z;
   z = (hwFunc(x, y + K) - hwFunc(x,y-K)) / (2.0 * K); 
   return z;
}

void Form_Gradient() {

   int i,j,k;
   i=0;
   j=0;
   k=0;
   double xGrad,yGrad,small_i,small_j;
   
   double max_value=0.0;
   double magnitude=0.0;

   for (i=0; i<numParts; i++) {
      
         double px = xycoord[cnt].x;
         double py = xycoord[cnt].y;
         xGrad=Partial_Derivative_X(px,py);
         xycoord[i].velX = xGrad;
         yGrad=Partial_Derivative_Y(px,py);
         xycoord[i].velY = yGrad;
         magnitude = sqrt(pow(xGrad,2)+pow(yGrad,2));
         printf("X: %f Y: %f GradX: %f  GradY: %f  Magnitude: %f\n",px,py,xGrad,yGrad,magnitude);
         
         if (magnitude > max_value)
         {
            max_value = magnitude;
         }
      }
   }

   printf("Maximum Vector Magnitude is: %f\n",max_value);

   
   
}

void Save_Data() {

   FILE * fpVelocity;
   //FILE * fpY;
   //FILE * fpCoords;

   int i,j;
   i=0;
   j=0;

   if((fpVelocity=fopen("ax.out", "w+"))==NULL) {
    printf("Cannot open file fpX.\n");
   }

   //if((fpY=fopen("ay.out", "w+"))==NULL) {
   // printf("Cannot open file fpY.\n");
   //}

   //if((fpCoords=fopen("xyCoords.out", "w+"))==NULL) {
   // printf("Cannot open file xyCoords.\n");
   //}
   
   for (i=0; i<numParts; i++) {
          
         fprintf(fpVelocity,"%f \t",coords[i]);
         fprintf(fpVelocity,"%f \t",coords[i]);
         fprintf(fpVelocity,"%f \t",coords[i]);
         fprintf(fpVelocity,"%f \t",coords[i]);
         
         
      }
      if (i != numParts-1) {
         fprintf(fpVelocity,"\n");         
      }
   }
   fclose(fpVelocity);
}

int main() {
   
   Make_Particles();
   Form_Gradient();
   Save_Data();
}
