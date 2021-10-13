
#ifndef INTEGRATION_PROCEDURE_H
#define INTEGRATION_PROCEDURE_H
#include "Function.h"

typedef struct {
	double lowX;
	double upX;
	double lowY;
	double upY;
	double lowZ;
	double upZ;
} MyStr;

void Integration(double (*function)(double, double, double, double*), 
			int series, double *parameter, int myid, int size, int count, 
			MyStr *listINcom, double *subSumm, int *step, double *bounds);

#endif