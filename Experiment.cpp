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
  // cout << "1 "  << endl;
// each thread own one used array
	// used = new bool[balls];
	usedPerThread = new bool[balls];
	bool *usedPerThread = usedPerThread;
  // cout << "2 "  << endl;

	for (long i = 0; i < hmax + 1; i++)
		histogram[i] = 0;
    // cout << "3 " << endl;
}

void Experiment::clearUsed() {

	for (int i = 0; i < balls; i++)
		usedPerThread[i] = false;
		// used[i] = false;
    // cout << "4 "  << endl;
}

long Experiment::singleExperimentResult() {
	long sum = 0;
	int ball;
	double p;
  // cout << "5 "  << endl;

	clearUsed();

  // cout << "6 "  << endl;

	// int seed = (unsigned)(random() * (omp_get_thread_num()+2));

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
		// struct timeval tf;
		// gettimeofday(&tf, NULL);
		// int seed = tf.tv_sec *1000 + tf.tv_usec * 0.001 + omp_get_thread_num();
		cout << "seed: " << seed << endl;
		srand48_r(seed, &drand_Buffor);
	}
}

  // cout << "7 "  << endl;
	// srand48_r(seed, &drand_Buffor);
// #pragma omp for private(ps)
// #pragma omp for
	for (int i = 0; i < drawsNumber; i++) {
    // cout << "11 "  << endl;
		double result = 0;
		// drand48_r(drand_Buffor, &result);
		drand48_r(&drand_Buffor, &result);
    // cout << "12 "  << endl;

		// cout << "result: " << result << endl;
		double hh = ((double) balls * result) / ( RAND_MAX + 1.0);
		int jj = (int) result * 1e9;
		cout << "result: " << jj << endl;
		// cout << "balls: " << jj << endl;
		ball = 1 + (int) (((double) balls * result * 1e9) / ( RAND_MAX + 1.0)); // rand losuje od 0 do RAND_MAX wlacznie

		if (usedPerThread[ball - 1])
			continue;

    // cout << "21 "  << endl;
		p = Distribution::getProbability(i + 1, ball); // pobieramy prawdopodobienstwo wylosowania tej kuli
    // cout << "22 "  << endl;

		if ((result / ( RAND_MAX + 1.0)) < p) // akceptacja wyboru kuli z zadanym prawdopodobienstwem
				{
#ifdef DEBUG_ON
			cout << "Dodano kule o numerze " << ball << endl;
#endif
			usedPerThread[ball - 1] = true;
			sum += ball; // kule maja numery od 1 do balls wlacznie
			i++;
		}
		// cout << "6666666666666" << endl;
	}

	// cout << "Suma = " << sum << endl;

	return sum;
}

Result * Experiment::calc(long experiments) {
	cout << "calc" << endl;

// #pragma omp parallel
// {
 // for shared(experiments, histogram) private(l)\
				// reduction(+ : histogram)
#pragma omp parallel for
	for (long l = 0; l < experiments; l++) {
		// cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		// i = singleExperimentResult() i pragma omp atomic zeby zabezpieczyc histogram
		int i = singleExperimentResult();
		// cout << "i: " << i << endl;
#pragma omp atomic
		histogram[i]++;
	}
	cout << "here " << endl;
// }

	long maxID = 0;
	long minID = 0;
	long maxN = 0;
	long minN = experiments;
	double sum = 0.0;
	long values = 0;

	cout << "vars " << endl;

// for shared(maxN, maxID, hmin, hmax, histogram, sum) private(idx)\
				// reduction(+ : sum, values)
#pragma omp parallel for
	for (long idx = hmin; idx <= hmax; idx++) {
		// cout << "_-________------______------_____-----______-----_______------_________-------------" << endl;
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
	cout << "return" << endl;
	return new Result(maxID, maxN, sum / values, values);
}

Experiment::~Experiment() {
	delete[] histogram;
	// delete[] used;
	delete[] usedPerThread;
}
