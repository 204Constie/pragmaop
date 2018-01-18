/*
 * Distribution.h
 *
 *  Created on: 16 maj 2016
 *      Author: oramus
 */

#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

class Distribution {
public:
	// prawdopowobienstwo wylosowania kuli o numerze v gdy jest ona losowana jako n-ta z kolei
	static double getProbability( int n, int v );
};


#endif /* DISTRIBUTION_H_ */
