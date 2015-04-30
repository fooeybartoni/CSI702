#include <stdio.h>
#include <math.h>

#define A 1
#define B 1
#define H 0.01
#define K 0.01

double hwFunc(double x, double y);

double Partial_Derivative_X(double x, double y);

double Partial_Derivative_X(double x, double y);

void Form_Gradient();

void Form_Actual();

typedef struct {
    double gradX;
    double gradY;
} grad;

grad gradient[100][100];

grad actual[100][100];

//-------------------------------------------
//function implementation
double hwFunc(double x, double y) {
   double z;
   z = A * x  + B * (x / (pow(x,2) + pow(y,2)));
   //printf("%f,%f,%f\n",x,y,z);
   return z;
}

double Partial_Derivative_X(double x, double y) {

   double z;
   z = (hwFunc(x + H, y) - hwFunc(x-H,y)) / (2 * H); 
   return z;
}

double Partial_Derivative_Y(double x, double y) {

   double z;
   z = (hwFunc(x, y + K) - hwFunc(x,y-K)) / (2 * K); 
   return z;
}

void Form_Gradient() {

   int i,j;
   i=0;
   j=0;
   double x,y,small_i,small_j;

   for (i=0; i<100; i++) {

      for (j=0; j<100; j++) {

         small_i = i/25.0 - 2.0;
         small_j = j/25.0 - 2.0;

         x=Partial_Derivative_X(small_i,small_j);
         gradient[i][j].gradX=x;
         y=Partial_Derivative_Y(small_i,small_j);
         gradient[i][j].gradY=y;
         printf("%f,%f,%f,%f\n",small_i,small_j,x,y);
      }
   }

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
