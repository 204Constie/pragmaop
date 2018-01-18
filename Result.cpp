/*
 * Result.cpp
 *
 *  Created on: 16 maj 2016
 *      Author: oramus
 */

#include "Result.h"

Result::Result(double maxValue, long maxN,  double avgValue, long values) {
	this->maxValue = maxValue;
	this->maxN = maxN;
	this->avgValue = avgValue;
	this->values = values;
}

Result::~Result() {
}
