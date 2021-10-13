#ifndef PLOTFUNCION_H 
#define PLOTFUNCION_H

// IN: 
//		- function with two parameters and an array of parameters
//		- boundaries of x-parameter
//		- step of x arising
//		- parameter y
//		- array of parameters
//		- name of the output file
// IN-OUT:
//		- value of function	
//		- string message of the result of the end of the program work
int printFunction1Dx (double *valueoffunction, double function(double, double, double*), double x1, double x2, double step, double y, double *par, char *filename, char *msg[]);

// IN: 
//		- function with two parameters and an array of parameters
//		- parameter x
//		- boundaries of y-parameter
//		- step of y arising
//		- array of parameters
//		- name of the output file
// IN-OUT:
//		- value of function	
//		- string message of the result of the end of the program work
int printFunction1Dy (double *valueoffunction, double function(double, double, double*), double x, double y1, double y2, double step,  double *par, char *filename, char *msg[]);
#endif