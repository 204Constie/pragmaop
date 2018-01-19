#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

using namespace std;
   const long ITER = 1000000L;
   double *t = new double [ ITER ];
   long suma = 0;
   double maxim;

void initAndCalcSuma() {
#pragma omp parallel for reduction( +:suma)
   for ( long l = 0; l < ITER; l++ ) {
     t[ l ] = 0.0;
     suma+=l;
   }
}

void initT() {
   for ( long i = 1; i < ITER/500; i++ )
#pragma omp parallel for
     for ( long j = 0; j < ITER; j++ )
        t[ j ] += 0.1+i+j;
}

void initT2() {
#pragma omp parallel for schedule( guided )
   for ( long i = 0; i < ITER/10; i++ )
     for ( long j = 0; j < i; j++ )
        t[ i ] += 0.1+i+j;
}

void findMax() {
   maxim = t[0];
   double local_max = maxim;
#pragma omp parallel firstprivate(local_max)
{
#pragma omp for
   for ( long i = 1; i < ITER; i++ )
     if ( t[i] > local_max ) local_max = t[i];

#pragma omp critical
{
     if ( local_max > maxim ) maxim = local_max;
} // critical
} // parallel
}

int main( void ) {

   initAndCalcSuma();
   initT();
   initT2();
   findMax();

   cout << "SUMA : " << suma << " max = " << maxim << endl;

} // main
