// DoubleExponenta.h содержит определение функций, обеспечивающих 
// преобразование, типа ДВОЙНОЙ ЭКСПОНЕНТЫ переменной интегрирования 


#ifndef DOUBLEEXPONENTA_H
#define  DOUBLEEXPONENTA_H 

double DEvariable_Z(double t_variable);
double DEvariable_X(double t_variable, double boundaryLow, double boundaryUp);
double DEtransformFunction2D(double (*myfunction2D)(double, double, double*),
				double Xt_variable, double Yt_variable, double *par, double *bounds);
double DEtransformFunction3D(double (*myfunction3D)(double, double, double, double*),
				double Xt_variable, double Yt_variable, double Zt_variable, double *par, double *bounds);
double DEboundary();
#endif