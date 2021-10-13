// Интервалы интегрирования разбиваются на части,	//
// соответственно к количеству рабочих процессов	//
// После вызывает интегрирующий алгоритм,			//
// передавая ему как параметр масив интервалов		//
//////////////////////////////////////////////////////

/* Определяет приделы интегрирования*/
#include "Boundaries.h"  
/* Метод преобразования подинтегральной функции Double Exponenta*/
/* L.Ye NUMERICAL QUADRATURE: THEORY AND COMPUTATION*/
#include "DoubleExponenta.h" 
/*Содержит интегрируемую функцию*/
#include "Function.h"
/*Содержит некоторые необходимые параметры задачи,*/
/*значение которых опредеялет значение интеграла */
#include "InputData.h"
/*Интегрирующий алгоритм*/
#include "integration_procedure.h"
#include <stdio.h>
#include <stdlib.h> //malloc
//#include <conio.h>
#include <string.h>	//memcpy
#include <mpi.h>
#include <math.h>
#include <time.h> //

void PreIntegration (int rank, int size, int series, double *set_parmtrs, double *result, double *worktime) {
	int MPIErrorCode=1;

	double subSumm=0.0, totalSumm = 0.0;
	double intervalX, intervalY, intervalZ;
	double statisticTime = 0.0, totalstatisticTime = 0.0;
	
	//double parameter[4];	//массив, содержащий некоторые параметры интегрирования
	double bounds[6];		//массив, содержащий пределы интегирования

	MyStr *list;		//весь список интервалов интегрирования
	MyStr *sendlist;	//список, отсылаемых интервалов
	
	int *step;
	int tempstep = 0;
	int totalstep = 0;
	int i, j, k;
	double lowX, lowY, lowZ;
	int count = 0;
	double startwtime = 0.0;
	double endwtime = 0.0;

	MPI_Status status;
	MPI_Request reqSumm[1] = {MPI_REQUEST_NULL};

	////////////////////////////
	//..Создаем новый тип MPI...
	// Название типа
	MPI_Datatype strtype;
	// Данные с идентификационным номером id
	MPI_Datatype IdData;
	// содержание	
	MPI_Datatype type[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype IdType[7] = {MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

	// количество элементов
	int blocklen[6] = {1, 1, 1, 1, 1, 1};
	int IdBlocklen[7] = {1, 1, 1, 1, 1, 1, 1};
	// расстояние между элементами (в памяти)
	MPI_Aint disp[6] = {0, 0, 0, 0, 0, 0};
	MPI_Aint IdDisp[7] = {0, 0, 0, 0, 0, 0, 0};
	// конструируем тип
	MPI_Type_struct(6, blocklen, disp, type, &strtype);
	MPI_Type_struct(7, IdBlocklen, IdDisp, IdType, &IdData);
	// реестрация типа
	MPI_Type_commit(&strtype);
	MPI_Type_commit(&IdData);
	
	// засекаем время работы программы
	startwtime = MPI_Wtime();	

	list = (MyStr *)malloc(size*size*size*1000*sizeof(MyStr));//список для интервалов интегрирования
	sendlist = (MyStr *)malloc(size*size*1000*sizeof(MyStr)); // список для отсылки заданий (интервалов интегрирования)
	
	step = (int *)malloc(size*sizeof(int));			// массив счетчиков
														//, указывающих на количество итераций, 
														// проведенных каждым процессом
	for (i = 0; i <size; i++) {
		step[i] = 0;
	}
	// in ELI_La parameter[0-2] is a Direction of unit velocity vector 
	//parameter[0] =  sqrt(1.0 - wz*wz);
	//parameter[1] =  wz;

	// in ELI_La parameter[3-4] is a Temperature of electron gas, Dimentionless
	//parameter[2] =  TAU/3.;
	//parameter[3] =  TAU/3.;

	count = 0;
	if (rank == 0) {
		// преобразование координат по алгоритму двойной экспоненты
		list[count].lowX = - DEboundary();   
		list[count].upX = -list[count].lowX; 
		list[count].lowY = list[count].lowX;
		list[count].upY = -list[count].lowX;
		list[count].lowZ = list[count].lowX;
		list[count].upZ = -list[count].lowX;

		//fprintf(stdout, "%g\t%g\t%g\t%g\n", list[0].upZ, list[0].lowZ, list[0].upY, list[0].lowY);fflush(stdout);
		
		lowX = list[count].lowX ;
		lowY = list[count].lowY ;
		lowZ = list[count].lowZ ;

		// Разбиение интервала интегрирования на более мелкие части
		// шаги разбиения
		
		intervalX = fabs(list[count].upX - list[count].lowX) / size /10. ;
		intervalY = fabs(list[count].upY - list[count].lowY) / size/10. ;
		intervalZ = fabs(list[count].upZ - list[count].lowZ) / size	/10. ;
		for (j = 0; j < size*10; j++) {
			for (k = 0; k < size*10; k++) {
				for (i = 0; i < size*10; i++) {
					list[count].lowX = lowX + intervalX*k;
					list[count].upX = lowX + (intervalX)*(k+1);
					list[count].lowY = lowY + intervalY*j;
					list[count].upY = lowY + (intervalY)*(j+1);
					list[count].lowZ = lowZ + intervalZ*i;
					list[count].upZ = lowZ + (intervalZ)*(i+1);
					count++;
				}
			}
		}
		count = count/size; // количесто элементов для отправки
		// #WARNING
		// Рассылку можно заменить просто выборочным цыклом,
		// когда из цыкла каждый процесс выбирает эллемент 
		// через определенный шаг
		for (i = 1; i < size; i++) {
			// Отсылаем каждому процессу задания:
			// отсылка количества заданий
			MPI_Send(&count, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			// отсылка заданий count шт.
			memcpy(sendlist, &list[i*(count)], count*sizeof(MyStr));
			MPI_Send(sendlist, 6*count, strtype, i, 0, MPI_COMM_WORLD);
		}
	} else {
		// Получение заданий для работы:
		// получение количества присылаемых заданий
		MPI_Recv(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		// получение заданий count шт.
		//MPI_Recv(list, 4*count, strtype, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(list, 6*count, strtype, 0, 0, MPI_COMM_WORLD, &status);
		// #WARNING 
		// Получющий список может быть и меньшим, 
		// размера sendlist
	}
	bounds[0] = LOWX;
	bounds[1] = UPX;
	bounds[2] = LOWY;
	bounds[3] = UPY;
	bounds[4] = LOWZ;
	bounds[5] = UPZ;

	//if (rank == 0) {
	//	for(i = 0; i < 6; i++) {
	//	fprintf(stdout, "bound %d = %g\n", i, bounds[i]); fflush(stdout);
	//	}
	//}
	fprintf(stdout, "a"); fflush(stdout);
	/////////////////////////////////////////////////////////////
	Integration(ELiMA, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	//Integration(ELIA, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	//Integration(ELI, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	//Integration(function3D, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);

	/////////////////////////////////////////////////////////////
	//Просуммируем промежуточные значения суммы//
	// Сбор результатов от каждого процесса
	if (rank != 0) {
		MPI_Send(&tempstep, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	} else {
		totalstep += tempstep;
		for(i = 1; i < size; i++) {	
			MPI_Recv(&step[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
			totalstep +=step[i];
		}
	}
	if (rank != 0) {
		MPI_Send(&subSumm, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
	} else {
		totalSumm += subSumm;
		for(i = 1; i < size; i++) {	
			MPI_Recv(&subSumm, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
			totalSumm += subSumm;
		}
	}
	if (rank == 0) {
		step[0] = tempstep;
		for(i = 0; i < size; i++) {
			fprintf(stdout, "proc %d do %6d (%5.2f %%) from %d iterations\n", 
					i, step[i], (double)step[i]/(double)totalstep*100., totalstep);
		}
		fprintf(stdout, "summ = %e\n", totalSumm);
		fflush(stdout);
		*result = totalSumm;
		endwtime = MPI_Wtime();
		*worktime = endwtime-startwtime;
		printf("wall clock time = %f\n", *worktime);fflush(stdout);
	}


	// освобождение памяти
	free(step);
	free(list);
	free(sendlist);
	//printf("proc %d****************\n", rank);fflush(stdout);
}