#ifndef SPORE_H
#define SPORE_H

#include <chrono>
#include <random>
#include "Img.h"

#define PI 3.14159265358979323846

enum Rtype{
	CAUCHY, CAUCHY_MIX  // NO means that there is no wind
};

class Sporulation{
private:
	double vonmisesvariate(double mu, double kappa);
	int **sp;
public:
	Sporulation();
	void SporeGen(Img& I, double **weather, double rate);
	void SporeSpreadDisp(Img& S_umca, Img& S_oaks, Img& I_umca, Img& I_oaks, Img& IMM, 
		Rtype rtype, double **weather, double scale1, int kappa=2, Direction wdir=NO, 
		double scale2 = 0.0,double gamma = 0.0);
	~Sporulation();
};

#endif