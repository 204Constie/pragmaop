#include "Distribution.h"

double Distribution::getProbability( int n, int ball ) {
	if ( n >= ball ) return 1.0;
	return 0;
}
