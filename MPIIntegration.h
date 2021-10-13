// Function.h содержит определение прототипа подинтегральной функции
#ifndef MPIINTGRATION_H
#define MPIINTGRATION_H

// ПРИНИМАЕТ:	1. IN Порядковый номер процесса
	//			2. IN Количество рабочих процессов
	//			3. IN Порядковый номер запуска подпрограммы
	//			4. IN Параметры задачи
// ВОЗВРАЩАЕТ:	5. INOUT результат интегрирования
void PreIntegration (int rank, int size, int series, double *parameters, double *result, double *worktime);


#endif
