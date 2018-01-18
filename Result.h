/*
 * Result.h
 *
 *  Created on: 16 maj 2016
 *      Author: oramus
 */

#ifndef RESULT_H_
#define RESULT_H_

class Result {
public:

	double maxValue;
	long maxN;
	double avgValue;
	long values;

	Result(double maxValue, long maxN,  double avgValue, long values);
	virtual ~Result();
};

#endif /* RESULT_H_ */
