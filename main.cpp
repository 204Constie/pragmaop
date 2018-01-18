#include"Result.h"
#include"Experiment.h"

#include <iostream>
#include <sys/time.h>

using namespace std;

double timeLimit;
const long CALC_TIME = 15000;
const long CHUNK_SIZE = 5000;

bool hasTimeToContinue() {
	struct timeval tf;
	gettimeofday(&tf, NULL);
	double now = tf.tv_sec * 1000 + tf.tv_usec * 0.001;

	if (now < timeLimit)
		return true;  // limit czasu nie osiagniety
	return false; // juz po czasie, pora konczyc obliczenia
}

void initializeTimeLimit(double msec) {
	struct timeval tf;
	gettimeofday(&tf, NULL);
	timeLimit = tf.tv_sec * 1000 + tf.tv_usec * 0.001 + msec; // ustawiamy czas zakonczenia obliczen
}

void show(Result *res) {
	cout << "Wynik na podstawie " << res->values << " eksperymentow " << endl
			<< " - srednia wazona to                " << res->avgValue << endl
			<< " - najczesciej pojawila sie wartosc " << res->maxValue
			<< " wystapila ona " << res->maxN << " razy" << endl;
}

int main(void) {

	Experiment *ex = new Experiment(100, 50); // 50 kul ze 100

	initializeTimeLimit(CALC_TIME);

	do {
		show(ex->calc(CHUNK_SIZE));
	} while (hasTimeToContinue());

	return 0;
}



// #include"Result.h"
// #include"Experiment.h"
//
// #include <iostream>
// #include <sys/time.h>
//
// using namespace std;
//
// double timeLimit;
// const long CALC_TIME = 15000;
// const long CHUNK_SIZE = 5000;
//
// bool hasTimeToContinue() {
// 	struct timeval tf;
// 	gettimeofday(&tf, NULL);
// 	double now = tf.tv_sec * 1000 + tf.tv_usec * 0.001;
//
// 	if (now < timeLimit)
// 		return true;  // limit czasu nie osiagniety
// 	return false; // juz po czasie, pora konczyc obliczenia
// }
//
// void initializeTimeLimit(double msec) {
// 	struct timeval tf;
// 	gettimeofday(&tf, NULL);
// 	timeLimit = tf.tv_sec * 1000 + tf.tv_usec * 0.001 + msec; // ustawiamy czas zakonczenia obliczen
// }
//
// void show(Result *res) {
// 	cout << "Wynik na podstawie " << res->values << " eksperymentow " << endl
// 			<< " - srednia wazona to                " << res->avgValue << endl
// 			<< " - najczesciej pojawila sie wartosc " << res->maxValue
// 			<< " wystapila ona " << res->maxN << " razy" << endl;
// }
//
// int main(void) {
//
// 	Experiment *ex = new Experiment(100, 25); // 25 kul ze 100
//
// 	initializeTimeLimit(CALC_TIME);
//
// 	do {
// 		show(ex->calc(CHUNK_SIZE));
// 	} while (hasTimeToContinue());
//
// 	return 0;
// }
