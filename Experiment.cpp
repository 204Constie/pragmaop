/*
 * Experiment.cpp
 *
 *  Created on: 16 maj 2016
 *      Author: oramus
 */

 // used jest kazdy wlasny
 // generator bufor jest wazny (kazdy wlasny?)
 // the code for final result is in openMP II like last one sth with ++ or stuff

#include<stdlib.h>
#include<iostream>
#include <omp.h>
#include <stdio.h>
#include <sys/time.h>

#include "Experiment.h"
#include "Distribution.h"

#define DEBUG_ON_

using namespace std;

struct drand48_data drand_Buffor;
#pragma omp threadprivate(drand_Buffor)

int flag = 0;
#pragma omp threadprivate(flag)

Experiment::Experiment(int balls, int drawsNumber) {
	this->balls = balls;
	this->drawsNumber = drawsNumber;

	hmax = 0; // wyznaczamy maksymalna sume
	hmin = 0; // i najmniejsza sume

	for (int i = 0; i < drawsNumber; i++) {
		hmax += balls - i;
		hmin += i + 1; // 1 + 2 + 3 + ... liczba losowan
	}

	cout << "Histogram min: " << hmin << " max: " << hmax << endl;

	histogram = new long[hmax + 1];

// each thread own one used array
	usedPerThread = new bool[balls];
	bool *usedPerThread = usedPerThread;

	for (long i = 0; i < hmax + 1; i++)
		histogram[i] = 0;
}

void Experiment::clearUsed() {

	for (int i = 0; i < balls; i++)
		usedPerThread[i] = false;
}

long Experiment::singleExperimentResult() {
	long sum = 0;
	int ball;
	double p;

	clearUsed();

	// struct plantSeed{
	// 	plantSeed(bool b){
	// 		int seed = (unsigned)(random() * (omp_get_thread_num()+2));
	// 		srand48_r(seed, &drand_Buffor);
	// 	}
  //
	// };
	// plantSeed ps(true);
#pragma omp parallel
{
	if(flag == 0){
		flag = 1;
		int seed = (unsigned)(random() * (omp_get_thread_num()+2));
		cout << "seed: " << seed << endl;
		srand48_r(seed, &drand_Buffor);
	}
}

	for (int i = 0; i < drawsNumber; i++) {

		double result = 0;
		drand48_r(&drand_Buffor, &result);

		ball = 1 + (int) (((double) balls * rand() ) / ( 1 + 1.0)); // rand losuje od 0 do RAND_MAX wlacznie
		cout << "ball: " << ball << endl;

		if (usedPerThread[ball - 1])
			continue;

		p = Distribution::getProbability(i + 1, ball); // pobieramy prawdopodobienstwo wylosowania tej kuli

		if ((rand() / ( RAND_MAX + 1.0)) < p) // akceptacja wyboru kuli z zadanym prawdopodobienstwem
				{
#ifdef DEBUG_ON
			cout << "Dodano kule o numerze " << ball << endl;
#endif
			usedPerThread[ball - 1] = true;
			sum += ball; // kule maja numery od 1 do balls wlacznie
			i++;
		}
	}

	// cout << "Suma = " << sum << endl;

	return sum;
}

Result * Experiment::calc(long experiments) {
	cout << "calc" << endl;

#pragma omp parallel for
	for (long l = 0; l < experiments; l++) {
		// i = singleExperimentResult() i pragma omp atomic zeby zabezpieczyc histogram
		int i = singleExperimentResult();
#pragma omp atomic
		histogram[i]++;
	}

	long maxID = 0;
	long minID = 0;
	long maxN = 0;
	long minN = experiments;
	double sum = 0.0;
	long values = 0;

#pragma omp parallel for
	for (long idx = hmin; idx <= hmax; idx++) {
		if (maxN < histogram[idx]) {
			maxN = histogram[idx];
			maxID = idx;
		}
#pragma omp atomic
		sum += idx * histogram[idx];
#pragma omp atomic
		values += histogram[idx];
	}


// indeks to wartosc, histogram -> liczba wystapien
	cout << "sum: " << sum << " values: " << values << " maxID: " << maxID << " maxN: " << maxN << endl;
	return new Result(maxID, maxN, sum / values, values);
}

Experiment::~Experiment() {
	delete[] histogram;
	// delete[] used;
	delete[] usedPerThread;
}
