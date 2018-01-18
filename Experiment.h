/*
 * Experiment.h
 *
 *  Created on: 16 maj 2016
 *      Author: oramus
 */

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "Result.h"

class Experiment {
private:
	long *histogram; // histogram
	bool *used; // wektor pomocniczy przy losowaniu bez zwracania

	bool *usedPerThread;

	int drawsNumber; // liczba kulek
	int balls; // liczba kulek
	long hmax; // maksymalny indeks histogramu
	long hmin; // najmniejszy, rozsadny indeks histogramu

	void clearUsed(); // kasuje pamiec o wylosowanych kulach
	long singleExperimentResult();

public:
	Experiment(int balls, int drawsNumber);
	virtual ~Experiment();

	Result *calc(long experiments);

};

#endif /* EXPERIMENT_H_ */
