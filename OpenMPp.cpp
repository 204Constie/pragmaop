#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

using namespace std;

int main( void )
{
   int report_steps = 10000000;
   long int steps = 0;
   double angle;
   double multipl = 6.28318531 / ( RAND_MAX + 1.0 );

   double *x;
   double *y;
   int s;

#pragma omp parallel
{

#pragma omp single
 {
   s = omp_get_num_threads();
   cout << "Threads : " << s << endl;

   x = new double[ s ];
   y = new double[ s ];
   for ( int i = 0; i < s; i++ )
   {
     x[ i ] = 0.0;
     y[ i ] = 0.0;
   }
 }
} // parallel

   while (true)
   {
#pragma omp parallel for private( angle )
     for ( int th = 0; th < s; th++ )
     {
       int n = omp_get_thread_num();
       cout << "N : " << n << endl;
       for ( int i = 0; i < report_steps; i++ )
       {
         angle = multipl * random(); // liczba od 0 do 2 pi
         x[n] += sin( angle );
         y[n] += cos( angle );
       }
     }
#pragma omp single
   {
     steps += report_steps;
     double d_sum = 0.0;
     for ( int th = 0; th < s; th++ )
     {
       d_sum += sqrt( x[ th ] * x[ th ] + y[ th ] * y[ th ] );
       cout << "D_sum : " << d_sum << endl;
     }
     d_sum /= s;

     cout << "Raport po " << steps << " srednia odleglosc / sqrt( liczba krokow ) = "
          << d_sum / sqrt( steps ) << endl;
   }
  }
}
