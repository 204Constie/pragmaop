// #include "Distribution.h"
//
// double Distribution::getProbability( int n, int v ) {
// 	return 1.0 / ( 5.0 + n * v  );
// }

#include "Distribution.h"

double Distribution::getProbability( int n, int ball ) {
	if ( n >= ball ) return 1.0;
	return 0;
}
