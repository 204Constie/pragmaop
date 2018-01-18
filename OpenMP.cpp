/*

Program symuluje ruch pijanego marynarza (inaczej: bladzenie przypadkowe).
Kolejne kroki maja identyczna dlugosc, ale wykonywane sa w losowym kierunku.

Mozna pokazac, ze odleglosc marynarza od polozenia startowego (tutaj 0,0)
jest proporcjonalna do pierwiastka kwadratowego z liczby krokow.

Prosze przerobic niniejszy program tak, aby osobne procesy symulowaly
NIEZALEZNY ruch marynarzy. Jeden proces symuluje ruch jednego marynarza.

W okresowych raportach (!!!! co 10 sekund !!!! a nie co N krokow) prosze podawac srednia odleglosc
marynarzy od punktu startowego w stosunku do pierwiastka kwadratowego z liczby zasymulowanych
dla danego marynarza krokow.

*/


#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main( void )
{
   int report_steps = 10000000;
   long int steps = 0; // liczba krokow od poczatku symulacji
   double x = 0.0; // wsp. x marynarza
   double y = 0.0; // wsp. y marynarza
   double angle;   // pod tym katem marynarz dokona jednego kroku
   double multipl = 6.28318531 / ( RAND_MAX + 1.0 );

   while (true)
   {
     for ( int i = 0; i < report_steps; i++ )
     {
       angle = multipl * random(); // liczba od 0 do 2 pi
       x += sin( angle );
       y += cos( angle );
     }
     steps += report_steps;
     cout << "Raport po " << steps << " odleglosc / sqrt( liczba krokow ) = "
          << sqrt( ( x * x + y * y ) / steps ) << endl;
   }
}
