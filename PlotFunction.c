#include <stdio.h>

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
int printFunction1Dx (double *valueoffunction, double function(double, double, double*), double x1, double x2, double step, double y, double *par, char *filename, char *msg[]) {
	FILE *fileprint;
	
	double x;
	char *tempfilename;
	sprintf(tempfilename, "%s_x=%s..%s_y=%s", filename, x1, x2, y);

	// open file verify procedure 
	if ((&fileprint, tempfilename, "w")) {
		*msg = "Error: printFunction1D could not opened fileprint"; 
		return 1;
	} 

	// Print values of x and function
	fprintf(fileprint, "x\tfunction\n");
	for (x = x1; x <= x2; x+=step) {
		fprintf(fileprint, "%g\t%g\n", x, function(x, y, par));
	}
}

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
int printFunction1Dy (double *valueoffunction, double function(double, double, double*), double x, double y1, double y2, double step,  double *par, char *filename, char *msg[]) {
	FILE *fileprint;
	
	double y;
	char valueofY [10];

	sprintf(valueofY, "%s", x);
	// open file verify procedure 
	if ((&fileprint, filename, "w")) {
		*msg = "Error: printFunction1D could not opened fileprint"; 
		return 1;
	} 

	// Print values of x and function
	fprintf(fileprint, "x\tfunction\n");
	for (y = y1; y <= y2; y+=step) {
		fprintf(fileprint, "%g\t%g\n", x, function(x, y, par));
	}
}