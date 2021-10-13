#include <math.h>
#include "InputData.h"

double DEvariable_Z(double t_variable) {	
	return tanh((PI/2)*sinh(t_variable));
}
double DEvariable_X(double t_variable, double boundaryLow, double boundaryUp)	{	
	return DEvariable_Z(t_variable)*(boundaryUp-boundaryLow)/2 + (boundaryUp+boundaryLow)/2;
}
double DEtransformFunction3D(double (*myfunction)(double, double, double, double*),
				double Xt_variable, double Yt_variable, double Zt_variable, double *par, double *bounds) {	
	double DerivateDEX, DerivateDEY, DerivateDEZ;
	double DerivateBouX, DerivateBouY, DerivateBouZ;
	double x_variable, y_variable, z_variable;
	double LOWX, UPX;
	double LOWY, UPY;
	double LOWZ, UPZ;
	LOWX = bounds[0];
	UPX = bounds[1];
	LOWY = bounds[2];
	UPY = bounds[3];
	LOWZ = bounds[4];
	UPZ = bounds[5];
	x_variable = DEvariable_X(Xt_variable, LOWX, UPX);
	y_variable = DEvariable_X(Yt_variable, LOWY, UPY);
	z_variable = DEvariable_X(Zt_variable, LOWZ, UPZ);

	DerivateBouX = (UPX-LOWX)/2.0;
	DerivateBouY = (UPY-LOWY)/2.0;
	DerivateBouZ = (UPZ-LOWZ)/2.0;

	DerivateDEX = (PI/2)*cosh(Xt_variable)/pow(cosh((PI/2)*sinh(Xt_variable)),2);
	DerivateDEY = (PI/2)*cosh(Yt_variable)/pow(cosh((PI/2)*sinh(Yt_variable)),2);
	DerivateDEZ = (PI/2)*cosh(Zt_variable)/pow(cosh((PI/2)*sinh(Zt_variable)),2);
	return  myfunction(x_variable, y_variable, z_variable, par)*DerivateDEX*DerivateDEY*DerivateDEZ*DerivateBouX*DerivateBouY*DerivateBouZ;
}
double DEtransformFunction2D(double (*myfunction2D)(double, double, double*),
				double Xt_variable, double Yt_variable, double *par, double *bounds) {	
	double DerivateDEX, DerivateDEY;
	double DerivateBouX, DerivateBouY;
	double x_variable, y_variable;
	double LOWX, UPX;
	double LOWY, UPY;
	LOWX = bounds[0];
	UPX = bounds[1];
	LOWY = bounds[2];
	UPY = bounds[3];
	x_variable = DEvariable_X(Xt_variable, LOWX, UPX);
	y_variable = DEvariable_X(Yt_variable, LOWY, UPY);

	DerivateBouX = (UPX-LOWX)/2.0;
	DerivateBouY = (UPY-LOWY)/2.0;

	DerivateDEX = (PI/2)*cosh(Xt_variable)/pow(cosh((PI/2)*sinh(Xt_variable)),2);
	DerivateDEY = (PI/2)*cosh(Yt_variable)/pow(cosh((PI/2)*sinh(Yt_variable)),2);
	return  myfunction2D(x_variable, y_variable, par)*DerivateDEX*DerivateDEY*DerivateBouX*DerivateBouY;
}

double DEboundary() {
	double epsilon = pow(10.0,-14);
	double h = pow(2.,-10);
	double t=0;
	double x_variable;
	double k; 
	for (k = 0;k<=20*pow(2.0,10);k++) {
		t=h*k;
		x_variable = tanh((PI/2)*sinh(t));
		if (fabs(x_variable-1) < epsilon) {
			break;
		}
	}
	return t;
}
