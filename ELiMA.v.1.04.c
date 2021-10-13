//////////////////////////////////////////////////////////////////////
// Изменения вносимые в связи с расчетом по Пархомчуку
// 1) function.c: double delta (TePerp) поменять коєфициент для расчетов параметров HESR
///////////////////////////////////
/// Программа MPI_Larkin27 предазначена для подсчета										//
/// интегрального уравнения ф.27 из работы Ларкина											//
/// А.И. Ларкин. Прохождение частиц через плазму ЖЭТФ, 37(1), 264-272 (1959)				//
/// [A.I. Larkin, Passage of particles through plasma, Sov. Phys. JETP 10, 186-191 (1960)]	//
//////////////////////////////////////////////////////////////////////////////////////////////
/// Интегральное уравнение представляет собой двумерный интеграл
///////////////////////////////////////////////////////////////////
/// athor: O.V.Khelemelia
/// data: 21.06.2014
///////////////////////////////////////////////////////////////////
/// Для подключения "mpi.h" необходимо провести следующие манипуляции:
/// 1) project->properties->VC++ Directories: in Include Derictories add adress to mpich2/include
/// 2) project->properties->VC++ Directories: in Library Derictories add adress to mpich2/lib
/// 3) project->properties->linker->input: in Additional Dependencies add mpi.h

#include "MPIIntegration.h"
#include "Boundaries.h"
#include "PlotFunction.h"
//#include "InputData.h"
#include "Function.h"
#include "readNline.h"
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <stdlib.h> //malloc
#include <time.h>
#include <direct.h> 

// количество параметров задачи
#define COUNTPARMTRS 6	// (3-компонента скорости)
					// абсолютное значение скорости

					// ((2)поперечная и продольная температурі)
// Название программы
#define PROG_NAME "ELiMA"
// Версия программы
#define PROG_VERSION ".v.5.05"
// Формат вывода данных
#define FILE_FORMAT ".dat"
// Внешний параметр
// Поперечная температура
#define TEMPERATURAPer "_AnIso(Par)h10(5)_TePer"
// Продоьлная температура
#define TEMPERATURAPar "TePar"
// Путь к файлу
#define FILEPATH "C:\\Program Files (x86)\\MPICH2\\bin\\"
int main (int argc, char **argv) {
	// log-file
	FILE *LOGfile;

	int MPIErrorCode = 1;
	int rank, size;
	int flagMPI_Init ;
	int i;
	double worktime;
	
	MPI_Comm MPI_COMM_External;

	// Создание лог-файла запуска программы

	MPI_Init(&argc, &argv); // MPI initialization
	MPI_Initialized (&flagMPI_Init); // test MPI initialization
	
	//fprintf(stdout, "Hi\n");fflush(stdout);
	if (!flagMPI_Init) {
		// определение начала работы программы
		time_t tt;
		struct tm * ptm;
		char buf[BUFSIZ];
		if ((fopen_s(&LOGfile ,"C:\\Program Files (x86)\\MPICH2\\bin\\log.dat", "w"))) { // TRUE is equal "File was not opened"
			fprintf(stdout, "Can't open log-file\n");fflush(stdout);
			exit(1);
		}
		ptm = (struct tm*)malloc(sizeof(struct tm));
		tt = time(NULL);
		localtime_s(ptm, &tt);
		strftime(buf, BUFSIZ, "%d-%m-%Y\nSTART \n", ptm);
		fprintf(LOGfile, buf); fflush(LOGfile);

		// программа начала работу не в параллельном режиме
		// !!! завершить работу !!!
		fprintf(LOGfile, "\\------------------------------\\ \n");fflush(stdout);
		fprintf(LOGfile, "MPI was not initiated.\nProgram MPI_Larkin stoped.");fflush(stdout);
		free(ptm);
		exit(1);
	}


	//definition of number of work processors
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	//definition of inbexes of work processors
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // index number of processor

	// Создание нового MPI комуникатора
	MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_External);
	if (rank == 0) {
	/// Определение время начала работы программы
	/// дата записывается в log.dat в формате
	/// Day-Month-Year Hour(00-23):Minute(00-59):Second(00-61)
		time_t tt;
		struct tm * ptm;
		char buf[BUFSIZ];
		ptm = (struct tm*)malloc(sizeof(struct tm));
		tt = time(NULL);
		localtime_s(ptm, &tt);
		strftime(buf, BUFSIZ, "%d-%m-%Y %H:%M:%S\nSTART \n", ptm);
	
	/// Создание log-файла
		if ( (fopen_s(&LOGfile ,"C:\\Program Files (x86)\\MPICH2\\bin\\log.dat", "w")) ) { // TRUE is equal "File was not opened"
			int MPIerrorCode = 1; // File was not opened
			fprintf(stdout, "Can't open log-file\n");fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
		}
		fprintf(LOGfile, buf); fflush(LOGfile);
		free(ptm);
	}

	// Программа работает только для более чем одного процессора
	if (rank == 0) {
		if(size <= 1) {
			fprintf(stdout, "Number of processes have to be more than 1\n");fflush(stdout);
			fprintf(LOGfile, "Program was stoped. Number of processes have to be more than 1\n");fflush(LOGfile);
			MPI_Abort(MPI_COMM_WORLD, MPIErrorCode);
		}
	}

	////////////////////////////
	// Главная часть программы//
	////////////////////////////
	for (i = 2; i < 3/*выход при отсутствии параметров*/; i++) {	
		// i - порядковый номер линии с которой нужно начинать считывать даные с файла
		double *sndrecvset;
		double *set_parmtrs;
		double result = 0.;		// конечный ответ
		int series = 1;			// порядковый номер запуска программы
		double nz=0.;
		int position = 0, length; 

		// необходимые определение для созданиянового mpi-типа
		MPI_Datatype sndrcvparmtrs; // имя нового типа
		MPI_Datatype oldtypes[COUNTPARMTRS+1];	// массив под старые типы
		MPI_Aint offsets[COUNTPARMTRS+1];		// массив под растояния между элементами нового типа
		int blklens[COUNTPARMTRS+1];				// массив под количество старых типов

		// выделение памяти под массив параметров для пересілки
		length = (COUNTPARMTRS+3)*sizeof(double)+sizeof(int); //основніе параметрі(count) и обезразмеренніе (3)
		set_parmtrs = (double*)malloc(length);	
		sndrecvset = (double*)malloc(COUNTPARMTRS*sizeof(double));
		//////////////////////////////////////////////////////////////////////
		/*//	  Create a send-recv struct type like int+3*double			//
		//////////////////////////////////////////////////////////////////////
		   /*Data types*/													//
		oldtypes[0] = MPI_INT;									//
		oldtypes[1] = MPI_DOUBLE;
		oldtypes[2] = MPI_DOUBLE;
		oldtypes[3] = MPI_DOUBLE;
		oldtypes[4] = MPI_DOUBLE;
		oldtypes[5] = MPI_DOUBLE;
		oldtypes[6] = MPI_DOUBLE;
		/*Количество элементов в каждом блоке*/							//
		blklens[0] = 1;														//
		blklens[1] = 1;												//
		blklens[2] = 1;												//
		blklens[3] = 1;												//
		blklens[4] = 1;												//
		blklens[5] = 1;												//
		blklens[6] = 1;												//
			/*Смещение каждого блока, измеряемые в байтах*/					//
		offsets[0] = 0;														//
		MPI_Type_size(MPI_INT, &offsets[1]);						//
		MPI_Type_size(MPI_DOUBLE, &offsets[2]);						//
		MPI_Type_size(MPI_DOUBLE, &offsets[3]);						//
		MPI_Type_size(MPI_DOUBLE, &offsets[4]);						//
		MPI_Type_size(MPI_DOUBLE, &offsets[5]);						//
		MPI_Type_size(MPI_DOUBLE, &offsets[6]);						//
			/*Собственно создание новой структуры*/							//
		MPI_Type_struct( COUNTPARMTRS+1, blklens, offsets, oldtypes, &sndrcvparmtrs );		//
			/*Реестрация новой структуры*/									//
		MPI_Type_commit( &sndrcvparmtrs);										//

		//{
		//	int j;
		//	MPI_Type_size(sndrcvparmtrs, &j);
		//	fprintf(stdout, "%d\n", j);
		//}
		//////////////////////////////////////////////////////////////////////
		// Считывание входных параметров задачи
		if (rank == 0) {
			FILE *InputParmtrs; // файл, с которого считываются параметры задачи
			char* error_msg;
			float *par;				// массив под параметры задачи
			// Открытие файла, содержащего необходимые параметры задачи
			if (fopen_s(&InputParmtrs, "C:\\Program Files (x86)\\MPICH2\\bin\\data.dat", "r+")) { // TRUE is equal "File was not opened"
				int MPIerrorCode = 1; // File was not opened
				fprintf(LOGfile, "\\------------------------------\\ \n");fflush(LOGfile);
				fprintf(LOGfile, "Can't open input data File \n");fflush(LOGfile);
				fprintf(stdout, "Can't open input data File\n");fflush(stdout);
				free(par);
				fclose(LOGfile);
				MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
			} 
			fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
			fprintf(LOGfile, "Input data file was opened successfully\n");fflush(stdout);

			// выделение памяти под массив параметров для считівания из файла
			par = (float*)malloc(COUNTPARMTRS*sizeof(double));	
			//	считывание параметров с файла
			if(readNline(&InputParmtrs, i, par, COUNTPARMTRS, &error_msg)) { // е. TRUE, тогда подпрограмма завершилась с ошибкой
				int k, m = 1;

				// е. возможны ошибки работы во время считывания даных,
				// описание ошибки содержется в error_msg
				fprintf(LOGfile, "\\------------------------------\\ \n");fflush(LOGfile);
				fprintf(LOGfile, error_msg);fflush(LOGfile);

				// упаковка отслаемых данных: ключ и список параметров
				position = 0;	// позиция упаковки
				MPI_Pack(&m, 1, MPI_INT, set_parmtrs, length, &position, MPI_COMM_External);
				//MPI_Pack(par, count, MPI_DOUBLE, set_parmtrs, length, &position, MPI_COMM_External);
				for (k = 1; k < size; k++) {
					MPI_Send(set_parmtrs, 1, sndrcvparmtrs, k, 110, MPI_COMM_External);
				}
				break;

			}
			else { // FALSE: данные с файла были считаны удачно
				int k, m = 0;

				// упаковка отслаемых данных: ключ и список параметров
				position = 0;
				//fprintf(stdout, "proc %d m=%d\n", rank,  m);
				for (k = 0; k < COUNTPARMTRS; k++) {
					sndrecvset[k] = (double)par[k]; 
				}
				MPI_Pack(&m, 1, MPI_INT,set_parmtrs, length, &position, MPI_COMM_External);
				//MPI_Pack(sndrecvset, COUNTPARMTRS, MPI_DOUBLE, set_parmtrs, length, &position, MPI_COMM_External);
				for (k = 1; k < size; k++) {
					MPI_Send(set_parmtrs, 1, sndrcvparmtrs, k, 110, MPI_COMM_External);
				}
				for (k = 1; k < size; k++) {
						MPI_Send(sndrecvset, COUNTPARMTRS, MPI_DOUBLE, k, 111, MPI_COMM_External);
				}
				//
				//for(i = 0; i < COUNTPARMTRS; i++) {
				//	sndrecvset[i] = 0.;
				//}
				// Распаковка полученного сообщения (сигнал о завершении работы)
				//position = 0;
				//MPI_Unpack(set_parmtrs, length, &position, &m, 1, MPI_INT, MPI_COMM_External);
				//MPI_Unpack(set_parmtrs, length, &position, sndrecvset, COUNTPARMTRS, MPI_DOUBLE, MPI_COMM_External);
				//for(i = 0; i < COUNTPARMTRS; i++) {
				//	fprintf(stdout, "proc %d par[%d] = %g\n", 
				//			rank,  i, sndrecvset[i]);
				//}
			}
			free(par);
		}
		else {
			int m = 1;
			MPI_Status status;
			// получить список параметров задачи
			MPI_Recv(set_parmtrs, 1, sndrcvparmtrs, 0, 110, MPI_COMM_External, &status);
			position = 0;
			// Распаковка полученного сообщения (сигнал о завершении работы)
			MPI_Unpack(set_parmtrs, length, &position, &m, 1, MPI_INT, MPI_COMM_External);
			// е. полученный сигнал k = 1, то выход из цикла
//			fprintf(stdout, "proc %d m=%d\n", rank,  position);
			if (m) { // Некорректніе пересілаеміе данніе
				break;
			}
			else {
//				fprintf(stdout, "position is %d\n", position);fflush(stdout);
				//распаковка полученного сообщения (список параметров)
				//MPI_Unpack(set_parmtrs, length, &position, sndrecvset, 1, sndrcvparmtrs, MPI_COMM_External);
				MPI_Recv(sndrecvset, COUNTPARMTRS, MPI_DOUBLE, 0, 111, MPI_COMM_External, &status);
			}
		}
//		fprintf(stdout, "sizeofFloat = %d, sizeofInt = %d, sizeofDouble =%d\n", sizeof(float), sizeof(int), sizeof(double));
		//for(i = 0; i < COUNTPARMTRS; i++) {
		//	fprintf(stdout, "proc %d par[%d] = %g\n", 
		//			rank,  i, sndrecvset[i]);
		//}
		//fprintf(stdout, "proc %d par[0] = %g, par[1] = %g\n", rank,  par[0], par[1]);

		set_parmtrs[0] = sndrecvset[0];	// Vx
		set_parmtrs[1] = sndrecvset[1];	// Vy
		set_parmtrs[2] = sndrecvset[2];	// Vz
		set_parmtrs[3] = sndrecvset[3];	// |V|
		set_parmtrs[4] = sndrecvset[4];	// Te_per
		set_parmtrs[5] = sndrecvset[5];	// Te_par
		set_parmtrs[6] = delta(sndrecvset[4]);	// Delta
		set_parmtrs[7] = tau_perpendicular();	
		set_parmtrs[8] = tau_parallel(sndrecvset[5], sndrecvset[4]);
		
	
		//for(i = 0; i < 9; i++) {
		//	fprintf(stdout, "proc %d par[%d] = %g\n", 
		//			rank,  i, set_parmtrs[i]);
		//}
			////////////////////////////////////////////////////////////////
		//if (rank == 0) {
		//	double valueoffunction = 0;
		//	char *functionviewfile = "functionview";
		//	char *errormsg;
		//	if (printFunction1Dx(&valueoffunction, ELI, 0., 100000, 10000, 0, set_parmtrs, functionviewfile, &errormsg)) {
		//		fprintf(stdout, "ERORR");
		//	}
		//}
		//fprintf(stdout, "delta = %e\n", delta(sndrecvset[4]));
		fprintf(stdout, "Proc %d start next task\n", rank); fflush(stdout);
		//////////////////////////////////////////////////////////////////
		// Вызов функции интегрирования
		//////////////////////////////////////////////////////////////////
		PreIntegration(rank, size, series, set_parmtrs, &result, &worktime);
		//////////////////////////////////////////////////////////////////
		fprintf(stdout, "Proc %d finished the task\n", rank); fflush(stdout);

		/// Открытие файла для вывода результатов вычислений
		if (rank == 0) {
			double temperature_coeficient; 
			// создание имени файла, который будет 
			// содержать результаты вычислений
			char *file_name;			// массив под имя файла
			int length_name;			// длинна имени файла
			// Result of calculation are written to the file
			FILE *OUTResult;		// Выводится результат вычислений

			// V_e/V_o = sqrt(T_e * C^2/(m_eC^2)) = sqrt(T_e) * 41.9382450;
			// [V_e/V_o]/[T_e] = cm/eV/s
			temperature_coeficient = 41.94;
			// Длинна имени файла
			length_name = sizeof(PROG_NAME)+sizeof(PROG_VERSION)+sizeof(FILE_FORMAT)
							+sizeof(TEMPERATURAPer)+sizeof(TEMPERATURAPar)+sizeof(double)+sizeof(FILEPATH);
			// выделение памяти под массив, содержащий имя файла
			file_name = (char*)malloc(length_name);
			// образование имени файла
			sprintf_s(file_name, length_name, "%s%s%s%s%g%s%g%s", FILEPATH, PROG_NAME, PROG_VERSION, 
							 TEMPERATURAPer, set_parmtrs[4], TEMPERATURAPar, set_parmtrs[5], FILE_FORMAT);
			
			if ((fopen_s(&OUTResult , file_name, "a+"))) { // TRUE is equal "File was not opened"
															// a+ открыть с возможностью дописывать
															// начинает с конца файла
				int MPIerrorCode = 1;  // File was not opened
				fprintf(LOGfile, "\\------------------------------\\ \n");fflush(LOGfile);
				fprintf(LOGfile, "Can't open File for calculated results\n");fflush(LOGfile);
				fprintf(stdout, "Can't open File for calculated results\n");fflush(stdout);
				MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
			} 
			fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
			fprintf(LOGfile, "ResultFile was opened successfully\n");fflush(stdout);
			// Vi/V0	\t	Te \t	Delta \t	Tau \t 
			// result \t precision \t program \t 
			// bouX_ \t bouX^ \t bouK_ \t bouK^\n", 
			
			fprintf(OUTResult, "%g\t%g\t%g\t%g\t %g\t%g\t%g\t %g\t%g\t%g\t%d\t%f\t%s%s\t  %g\t%g\t%g\t%g\t%g\t%g\n", 
						set_parmtrs[0], set_parmtrs[1], set_parmtrs[2], set_parmtrs[3],
						set_parmtrs[3]*temperature_coeficient*sqrt(set_parmtrs[4]), set_parmtrs[4], set_parmtrs[5], 
						result, result/(temperature_coeficient*sqrt(set_parmtrs[4])), PRECISION, size, worktime, PROG_NAME, PROG_VERSION, 
						LOWX, UPX, LOWY, UPY, LOWZ, UPZ);fflush(OUTResult);
					
			//fprintf(OUTResult, "%g\n", result);fflush(OUTResult);
			//fprintf(OUTResult, "Vi/V0=%g\tTe=%g\tDelta=%g\tTau=%g\n", 
			//			set_parmtrs[0], set_parmtrs[1], set_parmtrs[2], set_parmtrs[3]);fflush(OUTResult);
			////fprintf(OUTResult, "Boundaries: \n phi_ = %e \t phi^ = %e\n", LOWZ, UPZ);fflush(OUTResult);
			//fprintf(OUTResult, "x_ = %e \t x^ = %e\n", LOWX, UPX);fflush(OUTResult);
			//fprintf(OUTResult, "k_ = %e \t k^ = %e\n", LOWY, UPY);fflush(OUTResult);
			//fprintf(OUTResult, "accuracy %e \n", PRECISION);fflush(OUTResult);
			if(OUTResult) {
				if(fclose(OUTResult)) {
					fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
					fprintf(LOGfile, "File was not closed\n");		
					fprintf(stdout, "File was not closed\n");		
				} else {
					fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
					fprintf(LOGfile, "ResultFile was closed successfully\n");fflush(stdout);
				}
			}
		}
		free(sndrecvset);
		free(set_parmtrs);
		MPI_Type_free(&sndrcvparmtrs);
	}
	/// Закрытие файлов для вывода (log-file, result-file)
	if(rank == 0) {
		if(LOGfile) {
			fprintf(LOGfile, "END\n");		
			fclose(LOGfile);
		}
	}

	fprintf(stdout, "proc %d finalize\n", rank);fflush(stdout);
	/*MPI-work are stoped*/
	/*noone mpi_procedure cannot be initialized after*/
	MPI_Finalize();


	return 0;
}
