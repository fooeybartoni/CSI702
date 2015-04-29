#include <stdio.h>
#include <math.h>

#define A 1
#define B 1
#define H 0.1
#define K 0.1

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
   printf("x: %f y: %f z: %f\n",x,y,z);
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

   for (i=0; i<100; i++) {

      for (j=0; j<100; j++) {

         gradient[i][j].gradX=Partial_Derivative_X(i,j);
         gradient[i][j].gradY=Partial_Derivative_Y(i,j);
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
   Form_Actual();
}
