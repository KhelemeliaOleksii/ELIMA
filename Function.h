// Function.h содержит определение прототипа подинтегральной функции
#ifndef FUNCTION_H
#define FUNCTION_H

// Test 2-dimentional function, to verify calculation algorithm
// x,y - variables of integration
// par[] - set of integral external parameters
double function2D (double x, double y, double *par);

// Test 3-dimentional function, to verify calculation algorithm
// x,y,z - variables of integration
// par[] - set of integral external parameters
double function3D (double x, double y, double z, double *par);

// ПРИНИМАЕТ: действительный аргумент, типа double
	// оценивает полученный аргумент
	// после чего определяет нужное приближение
// ВОЗВРАЩАЕТ: приближенное значение комплексной функции ошибок
	// от действительного аргумента, умноженой на 
	// exp(-x^2)*erfi(x);
double DispertionFunctionApproximation (double x);

// ПРИНИМАЕТ: 3 действительных аргумента, типа double
	// сферические координаты: phi, theta, rho
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// И.А. Ларкин. Прохождение частиц через плазму ЖЭТФ , 37(1), 264-272 (1959).
	// формула 27, обезразмеренная с востановленными единицами e=hbar=m_e=1
double Larkin27(double x, double k, double* parameter);

// ПРИНИМАЕТ: 2 действительных аргумента, типа double
	// сферические координаты: cos(theta) = -1..1, rho = 0..infinity
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// И.А. Ларкин. Прохождение частиц через плазму ЖЭТФ , 37(1), 264-272 (1959).
	// формула 27, обезразмеренная с востановленными единицами e=hbar=m_e=1
double ELI(double x, double k, double* parameter);

// ПРИНИМАЕТ: обратное значение нормированой скорости иона
// V_i/V_0
// Возращает значение параметра delta_0, равного
// delta_0 = hbar*omega_P/(2*m_e*c^2)*c^2/V_e^2
// pre = hbar*omega_P/(2*m_e*c^2)*c^2
// постоянная Планка, нормированная на 2*Pi
	// hbar = 6.582 E-16 eV*c
// плазменная частота, взята из эксперимента HESR 
	// omega_P = 2.9 E8 c^(-1)	
// скорость света в вакууме
	// с = 2.99792 E10 cm/c
// энергия покоя электрона
	// Ee = m_e*c^2 =  5.11 E5 eV
double delta (double Te_per);

// ПРИНИМАЕТ: значение нормированой скорости иона, V_i/V_0
	//		и значение температуры електронного газ, Te
// ВОЗРАЩАЕТ: значение параметра tau_i, равного
// tau_i = tau_0*(V_0/V_i)^2
// tau_0 = 1/(m_e*c^2)*c^2/V_0^2
// скорость света в вакууме
	// с = 2.99792 E10 cm/c
// энергия покоя электрона
	// Ee = m_e*c^2 =  5.11 E5 eV
// нормировочная скорость
	// V_0 = 10 E6 cm/c
double tau (double ViV0, double Te);

// Возращает значение параметра tau_per, равного
// tau_per = Te_per/me/Ve_per^2
double tau_perpendicular() ;

// Возращает значение параметра tau_per, равного
// tau_par = Te_par/me/Ve_per^2
double tau_parallel (double Te_par, double Te_per) ;

// ПРИНИМАЕТ: 3 действительных аргумента, типа double
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// Хелемели А.В. Влияние анизотропии на энергетические потери... 
double ELIA(double gx, double gy, double gz, double* parmtr);

// ПРИНИМАЕТ: 3 действительных аргумента, типа double
// ВОЗВРАЩАЕТ: значение функции, представленной в работе 
	// Хелемели А.В. Потери в замагниченном электронном газе... 
double ELiMA(double gx, double gy, double gz, double* parmtr);


#endif
