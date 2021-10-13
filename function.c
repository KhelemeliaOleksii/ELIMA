#include "Function.h"
#include "InputData.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h> //malloc
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
//#include "Faddeeva.h"

// ПРИНИМАЕТ: действительный аргумент, типа double
// оценивает полученный аргумент
// после чего определяет нужное приближение
// ВОЗВРАЩАЕТ: приближенное значение комплексной функции ошибок
//				от действительного аргумента, умноженой на 
//				exp(-x^2)*erfi(x);
double DispertionFunctionApproximation (double x) {
	const double a = 1.5;
	const double b = 2.6;
	if (fabs(x) <= a) {
		const double a0 = 2./sqrt(PI);
		const double a1 = 1.;
		const double a2 = 0.33333333333333333333333333333333; //1/3
		const double a3 = 0.1; // 1/10
		const double a4 = 0.02380952380952380952380952380952; //1/42
		const double a5 = 0.00462962962962962962962962962963; //1/216
		const double a6 = 7.5757575757575757575757575757576e-4; //1/1320
		const double a7 = 1.0683760683760683760683760683761e-4; //1/9360
		const double a8 = 1.3227513227513227513227513227513e-5; //1/75600

		return exp(-x*x)*a0*x*(a1 + a2*pow(x,2) + a3*pow(x,4) + a4*pow(x,6) + 
			a5*pow(x,8) + a6*pow(x,10) + a7*pow(x,12) + a8*pow(x,14));
	} else if (fabs(x) >= b) {
		const double a0 = 2./sqrt(PI);
		const double a1 = 0.5; //1/2
		const double a2 = 0.25; //1/4
		const double a3 = 0.375; // 3/8
		const double a4 = 0.9375; //15./16
		const double a5 = 3.28125; //35*24/256 -> 840/256
		const double a6 = 14.765625; //63*120/512 -> 7560/512
		const double a7 = 81.2109375; //231*720/2048 -> 166320/2048
		const double a8 = 527.87109375; //429*5040/4096 -> 2162160/4096

		return a0*(a1 + a2/pow(x,2) + a3/pow(x,4) + a4/pow(x,6) + 
			a5/pow(x,8) + a6/pow(x,10) + a7/pow(x,12) + a8/pow(x,14))/x;
	} else if (x > 0){
		const double a0 = 2./sqrt(PI);
		const double a4 = 0.092;
		const double a1 = 3.05; //1/2
		const double a2 = 0.0121; //1/4
		const double a3 = 0.2253; // 3/8

		return a0*(a4*pow((x-a1),2) - a2*x + a3);
	} else if (x < 0) {
		const double a0 = 2./sqrt(PI);
		const double a4 = 0.092;
		const double a1 = 3.05; //1/2
		const double a2 = 0.0121; //1/4
		const double a3 = 0.2253; // 3/8

		return a0*(-a4*pow((x+a1),2) - a2*x - a3);
	} else {
		fprintf(stdout, "Error value of incoming parameter.");
		fprintf(stdout, "Prog: double DispertionFunctionApproximation (double x).");
		fprintf(stdout, "File: function.c\n");
		exit(1);
	}
}

// ПРИНИМАЕТ: 3 действительных аргумента, типа double
	// сферические координаты: phi, theta, rho
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// И.А. Ларкин. Прохождение частиц через плазму ЖЭТФ , 37(1), 264-272 (1959).
	// формула 27, обезразмеренная с востановленными единицами e=hbar=m_e=1
double Larkin27(double x, double k, double* parameter){
	double delta = parameter[2];
	double tau = parameter[3]; 
	double Wk;
	double RePhi, ImPhi;
	double UpP, LowS, PreP, PreS; 
	double Wx1, Wx2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
	double k2;
	Wk = (x) ;//+ delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
	PreP = A_Khel*sqrt(PI)/sqrt(2.*tau)/delta/k/2.;		// Palarization operator coefficient
	Wx1 = (Wk-delta*k)/sqrt(2.*tau);					// arguments of dispertion function
	Wx2 = (Wk+delta*k)/sqrt(2.*tau);
	UpP = exp(-Wx1*Wx1);
	ImPhi = PreP*UpP*(1.-exp(-2.*delta*k*Wk/tau));		// Imaginary Part of Polarization Operator
	RePhi = PreP*(DispertionFunctionApproximation(Wx2) 
		- DispertionFunctionApproximation(Wx1));		// Real Part of Polarization Operator
														
	PreS = A_Khel/sqrt(2*PI*tau)/delta;					// Energy Losses Integral Coefficient
	k2=k*k;
	LowS = 1/((k2+RePhi)*(k2+RePhi)+ImPhi*ImPhi);

	return PreS*Wk*k2*UpP*LowS;
}


// ПРИНИМАЕТ: 2 действительных аргумента, типа double
	// сферические координаты: cos(theta) = -1..1, rho = 0..infinity
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// И.А. Ларкин. Прохождение частиц через плазму ЖЭТФ , 37(1), 264-272 (1959).
	// формула 27, обезразмеренная с востановленными единицами e=hbar=m_e=1
double ELI(double x, double k, double* parameter){
	double delta = parameter[2];
	double tau = parameter[3]; 
	double Wk;
	double RePhi, ImPhi;
	double UpP, LowS, PreP, PreS; 
	double Wx1, Wx2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
	double k2;
	Wk = (x) ;//+ delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
//	Wk = (x) + delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
	PreP = A_Khel*sqrt(PI)/sqrt(2.*tau)/delta/k/2.;		// Palarization operator coefficient
	Wx1 = (Wk-delta*k)/sqrt(2.*tau);					// arguments of dispertion function
	Wx2 = (Wk+delta*k)/sqrt(2.*tau);
	UpP = exp(-Wx1*Wx1);
	ImPhi = PreP*UpP*(1.-exp(-2.*delta*k*Wk/tau));		// Imaginary Part of Polarization Operator
	//ImPhi = PreP*UpP*(1.);		// Imaginary Part of Polarization Operator
	RePhi = PreP*(DispertionFunctionApproximation(Wx2) 
		- DispertionFunctionApproximation(Wx1));		// Real Part of Polarization Operator
														
	PreS = A_Khel/sqrt(2*PI*tau)/delta;					// Energy Losses Integral Coefficient
	k2=k*k;
	LowS = 1/((k2+RePhi)*(k2+RePhi)+ImPhi*ImPhi);

	//return 2*PreS*Wk*k2*UpP*LowS;
	return PreS*Wk*k2*ImPhi*LowS/PreP; // без температурного слагаемого
}

double function2D (double x, double y, double *par) {
	//return 1.;
	return 605*y/((1+120*(1-y))*((1+120*(1-y))*(1+120*(1-y))+25*x*x*y*y));//(0..1, 0..1, 0..1) = 0.1047591113142868D  01
	//return 1/(x+y+z)/(x+y+z);//(0..1, 0..1, 0..1) = 0.8630462173553432D   00
	//return cos(x + y); // -0.4 (0..3PI)(0..3PI)
	//return exp(-(x*x+y*y+z*z))*(x+y+z+1)/((x+y+z+1)*(x+y+z+1)+1);
	//return x/pow((1-x*x), 0.5);
	//return pow(x, 0.5)*pow(z, 0.5)*pow(y, 0.5)*log(x*y*z); // -.5925925926 (0..1)
	//return pow(x, 0.5)*log(x); // -0.444444 (0..1)
	//return exp(x)*cos(x); //1.905238691 (0..PI/2)
}

double function3D (double x, double y, double z, double *par) {
	return gsl_sf_bessel_In(1, 1.);
	//return 605*y/((1+120*(1-y))*((1+120*(1-y))*(1+120*(1-y))+25*x*x*y*y));//(0..1, 0..1, 0..1) = 0.1047591113142868D  01
	//return 1/(x+y+z)/(x+y+z);//(0..1, 0..1, 0..1) = 0.8630462173553432D   00
	//return cos(x + y); // -0.4 (0..3PI)(0..3PI)
	//return exp(-(x*x+y*y+z*z))*(x+y+z+1)/((x+y+z+1)*(x+y+z+1)+1);
	//return x/pow((1-x*x), 0.5);
	//return pow(x, 0.5)*pow(z, 0.5)*pow(y, 0.5)*log(x*y*z); // -.5925925926 (0..1)
	//return pow(x, 0.5)*log(x); // -0.444444 (0..1)
	//return exp(x)*cos(x); //1.905238691 (0..PI/2)
}
// Возращает значение параметра delta_0, равного
// delta_0 = hbar*omega_P/(2*T_e,perp)*tau_perp
// pre = hbar*omega_P/2
// Plank`s const, devided on 2*Pi
	// hbar = 6.582 E-16 eV*c
// Plasma frequency - for HESR 
	// omega_P = 2.9 E8 c^(-1)	
// const value - velocity of the light
	// c = 2.99792 E10 cm/c
// const value  - energy of electron
	// Ee = m_e*c^2 =  5.11 E5 eV
//// К расчетам по Пархомчуку
//// плотность n_e = 10^7 cm^(-3)
//double delta (double Te_per) {
//	double delta_0 = 5.89089*pow(10.,-8);	// Pre подсчитан предварительно
//	return delta_0/Te_per;
//}
// К расчетам по HESR
double delta (double Te_per) {
	double delta_0 = 9.5439*pow(10.,-8);	// Pre подсчитан предварительно
	return delta_0/Te_per;
}

// Возращает значение параметра tau_i, равного
// tau_i = tau_0*(V_0/V_i)^2
// tau_0 = 1/(m_e*c^2)*c^2/V_0^2
// скорость света в вакууме
	// с = 2.99792 E10 cm/c
// энергия покоя электрона
	// Ee = m_e*c^2 =  5.11 E5 eV
// нормировочная скорость
	// V_0 = 10 E6 cm/c
double tau (double ViV0, double Te) {
	double tau_0 = 1.758811E3;	// подсчитан предварительно
	return tau_0*Te/ViV0/ViV0;
}

// Возращает значение параметра tau_per, равного
// tau_per = Te_per/me/Ve_per^2
// Поскольку предполагается, что 
// Te_per = me*Ve_per^2,
// тогда tau_perpendicular = 1 !!всегда!!
double tau_perpendicular () {
	return 1.;
}
// Возращает значение параметра tau_per, равного
// tau_par = Te_par/me/Ve_per^2
double tau_parallel (double Te_par, double Te_per) {
	return Te_par/Te_per;
}


// ПРИНИМАЕТ: 3 действительных аргумента, типа double
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// Хелемели А.В. Влияние анизотропии на энергетические потери... 
double ELIA(double gx, double gy, double gz, double* parmtr){
	// компоненты безразмерного импульса налетающей частицы
	double sx = parmtr[0]; 
	double sy = parmtr[1];
	double sz = parmtr[2];
	// абсолютное значение вектора импульса безразмерного
	double V = parmtr[3];
	// поперечная и продольная тепературы электронного газа
	double tau_per = tau_perpendicular(); 
	double tau_par = tau_parallel(parmtr[5], parmtr[4]); 
	
	double W;
	double ge, g;
	double ReKappa, ImKappa;
	double Delta = delta(parmtr[4]);
//	double UpP, LowS, PreP, PreS; 
	double Xi1, Xi2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
//	double k2;
	double PreKappa;
	double PreInt, PostInt, Up;
	g = sqrt(gx*gx + gy*gy + gz*gz);
	ge = sqrt(gx*gx + gy*gy + gz*gz*tau_par/tau_per);
	
	W = (gx*sx + gy*sy + gz*sz) - Delta*mM*g*g;
	
	if ((gx == 0) && (gy == 0) && (gz == 0)) {
		return 0;
	} else {
		// аргуметні дисперсионніх функций
		Xi1 = (W +  Delta*g*g)/sqrt(2*tau_per)/ge; // Xi+;
		Xi2 = (W -  Delta*g*g)/sqrt(2*tau_per)/ge; // Xi+;

		Up = exp(-Xi2*Xi2);
		PreKappa = A_Khel*sqrt(PI/2/tau_per)/Delta/(g*g)/ge/2;
		ImKappa = PreKappa*
			Up*(1-exp(-2*Delta*W*g*g/ge/ge/tau_per));
		ReKappa = PreKappa*
			(DispertionFunctionApproximation(Xi1) 
			- DispertionFunctionApproximation(Xi2));
		PreInt = 1/Delta/sqrt(2*PI*tau_per)/(2*PI);
		PostInt = 1/((1+ReKappa)*(1+ReKappa) + ImKappa*ImKappa);
		return PreInt*W*Up*PostInt/(g*g*g*g)/(ge); 
	}
}

// ПРИНИМАЕТ: 3 действительных аргумента, типа double
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// Хелемели А.В. Потери в замагниченном электронном газе... 

	// h=Omega_H/Omega_P is dimless parameter of magnetic field
	// beta = 2*delta_0*h/tau_perp
	
	// диэлектрическая восприимчивость представлена в виде ряда
	// по параметру a^2 = q_t^2 * Delta/h 

double ELiMA(double gx, double gy, double gz, double* parmtr){
		// it needs to obtaine this parameter for all new jobs 
	double h = 1e5;	// parameter of magnetic field
	
	// a difference between Landau levels
	double s ; 
	double N = 15;

	// components of the union velocity vector V = V_i/V_0
	double Sx = parmtr[0]; 
	double Sy = parmtr[1];
	double Sz = parmtr[2];

	// value of the union velocity vector V = V_i/V_0
	double V = parmtr[3];

	// a perpendicular and a parallel components of the temperature of the electron gas
	double tau_per = parmtr[7]; 
	double tau_par = parmtr[8]; 
	double tau;
	double W;		// dimmles frequency
	double beta;	// magnet argument of exponent
	double a2;		// parameter of Lambda-function, 
						// it is always used in form a^{2N}	
	double PreSusc;	// prefix before susseptibility
	double Delta = parmtr[6];		// parameter delta = hbar W_0 / 2m_eV_0^2 	
	double g;
	double ReKappa, ImKappa;
	double Xi1, Xi2;

//	double A_Khel = 1.;			// Polarization coefficient by Khelemelia (h=0)
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
//	double AB = 4/PI;				// Polarization coefficient by Khelemelia (h not= 0)
	double PreInt, PostInt;
	double  Up = 0, Up0;
	double Pi = 3.14159265;
	double TemperatureFactor;		// a temperature factor 1-exp(-2 w delta0/tau)
//	double theta1H, theta2H;
	double argBess, preBess;
	double ImMagnetPart = 0.;
	double ImMagnetPartTilde = 0.;
	double ReMagnetPart = 0.;
	double preScaledBess; 
	double exp_beta, exp_beta2;
	double ScaledBesselI_0, ScaledBesselI_1;
	double ScaledBesselInApprox;
	double ScaledBessel2 = 0, ScaledBessel1=0;
	////////////////////////////////////////
	// General position
	g = sqrt(gx*gx + gy*gy + gz*gz);		// a dimless wave vector	
	W = (gx*Sx + gy*Sy + gz*Sz) - Delta*mM*g*g;		// a dimless frequency, determined as hw/hw_{0} = (E_{p} - E_{p-hk})/hw_{0}
	tau = tau_par/3. + 2.0*tau_per/3.0;		// temperature

	if (h <= 0) { // procedure calculate magnetic case (h > 0)
		fprintf(stdout, "ERROR in function ELiMA. Parameter h cannot be ZERO or Less"); fflush(stdout);
		return 0;
	} 

	if ((gz == 0) || ((gx == 0) && (gy == 0) && (gz == 0))) {  // to exclude numerical divergence
		return 0;
	} else {
		PreSusc = 1/4./Delta*sqrt(Pi)/sqrt(2*tau_par);  //a prefix of susceptibility. It's no variables of integration
		PreInt = PreSusc/Pi/Pi;	

		Xi1 = (W +  Delta*gz*gz)/sqrt(2*tau_par)/gz; // variable 1 of exponenta;
		Xi2 = (W -  Delta*gz*gz)/sqrt(2*tau_par)/gz; // variable 2 of exponenta;

		a2 = (gx*gx + gy*gy)*Delta/h;	// an argument of LambdaFunction
											// the parameter is inverse to the parameter of the magnetic field h = W_H/W_0
		beta = 2.*Delta*h/tau_per;		// an argument of exp
											// the argument is proportinal to the parameter of the magnetic field h = W_H/W_0  
		exp_beta2 = exp(-beta/2.);
		exp_beta = exp(-beta);

		// for GS 
		//argBess = 2.*a2*(exp(-beta/2.)/(1.-exp(-beta)));
		argBess = 2.*a2*(exp_beta2/(1.-exp_beta));
		//preBess	= exp(-a2*(1.+exp(-beta))/(1.-exp(-beta)));	

		//preScaledBess = exp(-a2*(1.-exp(-beta/2.))/(1.+exp(-beta/2.)));	
		preScaledBess	= exp(-a2*(1.-exp_beta2)/(1.+ exp_beta2));	

		ScaledBesselI_0 = gsl_sf_bessel_I0_scaled(argBess);
		ScaledBesselI_1 = gsl_sf_bessel_I1_scaled(argBess);
		
		for (s = 0; s<=N; s++) {
			double theta1Hp, theta1Hm, theta2Hp, theta2Hm;
			double XiSH;
			double RS;
			double BesselInApprox;
			XiSH = s*h/sqrt(2.*tau_par)/gz;
			theta1Hp = Xi1 + XiSH;
			theta1Hm = Xi1 - XiSH;
			theta2Hp = Xi2 + XiSH;
			theta2Hm = Xi2 - XiSH;
			if (s == 0) {
				ScaledBesselInApprox = ScaledBesselI_0;
				ScaledBessel2 = ScaledBesselInApprox;
			} else if (s == 1){
				ScaledBesselInApprox = ScaledBesselI_1;
				ScaledBessel1 = ScaledBesselInApprox;
			} 
			else{// if (s == 2) {
				if (fabs(argBess)  > 1e-16) {
					ScaledBesselInApprox = 	ScaledBessel2 -2*(s-1)/(argBess)*ScaledBessel1;
				} else {
					ScaledBesselInApprox = 0;
				}	
				ScaledBessel2=ScaledBessel1;
				ScaledBessel1 = ScaledBesselInApprox;
			}
/*			else if (s == 2) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.31501*exp(-argBess/2.59524) + 0.66376*exp(-argBess/0.77492) + 0.01877);
				} else if (argBess > 10){
					ScaledBesselInApprox  = ScaledBesselI_0 - 
(0.00526*exp(-argBess/56.75778) + 0.07012*exp(-argBess/7.63132));
				}	
			}
			else if (s == 3) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
						(0.43971*exp(-argBess/2.87401) + 0.51964*exp(-argBess/0.64401) + 0.03649);
				} else if (argBess > 10){
					ScaledBesselInApprox = ScaledBesselI_0 - 
						(0.00965*exp(-argBess/64.4686) + 0.11764*exp(-argBess/8.68577));
				}	
			}else if (s == 4) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.39464*exp(-argBess/3.6857) + 0.55537*exp(-argBess/0.63911) + 0.04741);
				} else if (argBess > 10){
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.15563*exp(-argBess/9.78864) + 0.0142*exp(-argBess/72.4432));
				}	
			}else if (s == 5) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.59352*exp(-argBess/0.66634) + 0.34453*exp(-argBess/4.38543) + 0.05867);
				} else if (argBess > 10){
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.01778*exp(-argBess/83.15393) + 0.1749*exp(-argBess/11.27954));
				}	
			}else if (s == 6) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.61095*exp(-argBess/0.68453) + 0.31*exp(-argBess/4.55027) + 0.07473);
				} else if (argBess > 10){
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.02033*exp(-argBess/96.00131) + 0.17945*exp(-argBess/13.13505));
				}	
			}else if (s == 7) {
				if (argBess >=0 && argBess <= 10) {
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.61248*exp(-argBess/0.68843) + 0.29203*exp(-argBess/4.28334) + 0.09082);
				} else if (argBess > 10){
					ScaledBesselInApprox = ScaledBesselI_0 - 
(0.02148*exp(-argBess/112.94641) + 0.17485*exp(-argBess/15.45809));
				}	
			} else {
				fprintf(stdout, "function.c: ELiMA: BesselInApprox: A too much value of the parameter S was used");fflush(stdout);
				exit(1);
			}*/
			RS = preScaledBess * ScaledBesselInApprox;

			if ( s==0 ) {
				ImMagnetPart += - RS *
					(exp(-theta1Hm*theta1Hm)  - exp(-theta2Hp*theta2Hp));

				ReMagnetPart += RS * 
					(DispertionFunctionApproximation(theta1Hm) -
						DispertionFunctionApproximation(theta2Hp));
				if (fabs(W) <= 1e-16) {
					ImMagnetPartTilde +=  (tau/tau_par) * RS * exp(-theta2Hp*theta2Hp);
				} else {
					ImMagnetPartTilde = ImMagnetPart;
				}
			} else {
				//ImMagnetPart += - RS *
				//		(exp(-s*beta/2.)*(exp(-theta1Hm*theta1Hm)  - exp(-theta2Hp*theta2Hp)) +
				//		exp(s*beta/2.)*(exp(-theta1Hp*theta1Hp)  - exp(-theta2Hm*theta2Hm)));
				//ReMagnetPart += RS *
				//	(exp(-s*beta/2)*(DispertionFunctionApproximation(theta1Hm) -
				//				DispertionFunctionApproximation(theta2Hp)) +
				//	exp(s*beta/2)*(DispertionFunctionApproximation(theta1Hp)-
				//				DispertionFunctionApproximation(theta2Hm)));
				ImMagnetPart += - RS * exp(s*beta/2.)*
						(exp(-s*beta)*(exp(-theta1Hm*theta1Hm)  - exp(-theta2Hp*theta2Hp)) +
						(exp(-theta1Hp*theta1Hp)  - exp(-theta2Hm*theta2Hm)));
				ReMagnetPart += RS * exp(s*beta/2)*
					(exp(-s*beta)*(DispertionFunctionApproximation(theta1Hm) -
								DispertionFunctionApproximation(theta2Hp)) +
					(DispertionFunctionApproximation(theta1Hp)-
								DispertionFunctionApproximation(theta2Hm)));
				if (fabs(W) <= 1e-16) {
					ImMagnetPartTilde +=  (tau/tau_par)* RS *(
						exp(-s*beta/2.)*exp(-theta2Hm*theta2Hm)*(1+s*h/gz/gz/Delta) +
						exp(s*beta/2.)*exp(-theta2Hp*theta2Hp)*(1-s*h/gz/gz/Delta));
				} else {
					ImMagnetPartTilde = ImMagnetPart;
				}
			}

		}
		if (fabs(W) <= 1e-16) {
			TemperatureFactor = 1.;
		} else{
			TemperatureFactor = 1/(1.-exp(-2.*Delta*W/tau));
		}
		//////////////////////////////////////////
		//to get the PostInt part
		// An imaginary part of susceptibility is
		ImKappa = PreSusc*ImMagnetPart/g/g/gz;
		
		// an real part of susceptibility is
		ReKappa = PreSusc*ReMagnetPart/g/g/gz;

		PostInt = 1/((1+ReKappa)*(1+ReKappa) + ImKappa*ImKappa);
		
		//return PreInt*ImMagnetPartTilde; 
		//return PreInt*W/(g*g*g*g*gz)*TemperatureFactor*PostInt; 
		return PreInt*ImMagnetPartTilde*W/(g*g*g*g*gz)*TemperatureFactor*PostInt; 
	}
}
