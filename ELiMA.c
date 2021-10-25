//////////////////////////////////////////////////////////////////////
// A program to calculate energy losses integrall by Simpson`s method 
// with a Double Exponent Transformation of 
// the integration area algorithm
//////////////////////////////////////////////////////////////////////

// ELI - Energy Losses of Ion in uniform electron gas
// ELIA - Energy Losses of Ion in electron gas 
			// with Anisotorpic velocity distribution
// ELIMA - Energy Losses of ion in Magnetized electron gas
			// with Anisotorpic velocity distribution

// Configuration parametrs of the Microsoft Visual Studio project
// are written on ReadConfigParameters.txt file

// -- Main file: ELiMA.c
// -- Parameters (COUNTPARMTRS) to the program are imported from file:
		//"C:\\Program Files (x86)\\MPICH2\\bin\\data.dat"
// -- integration procedure start by PreIntegration() function;
// -- result is written to OUTResult_file with compiled name:
		//PROG_NAME+PROG_VERSION+
		//+MAGNETFIELDPar+MAGNETFIELDPar+
		//+TEMPERATURAPer+Value_TEMPERATURAPer+
		//+TEMPERATURAPar+Value_TEMPERATURAPar+FILE_FORMAT
// precision of calculation: InputData.h
// boundaries of calculation: Boundaries.h

// Massage Passing Interface (mpi.h) is used for parallel calculation

#include "MPIIntegration.h"
#include "Boundaries.h"
#include "PlotFunction.h"
#include "Function.h"
#include "readNline.h"
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <stdlib.h> //malloc
#include <time.h>
#include <direct.h> 

#define COUNTPARMTRS 6	// define number of import parameters
	// 3 - union velocity vector: v_x, v_y, v_z
	// absolute velocity value: V
	// 2 - electron temperature: T_perp, T_parallel
	// longitudinal magnetic field: h

#define PROG_NAME "ELiMA" //names of program:
	// ELI
	// ELIA
	// ELiMA
#define PROG_VERSION ".v.6.01"
#define FILE_FORMAT ".dat"
#define TEMPERATURAPer "TePer"
#define TEMPERATURAPar "TePar"
#define FILEPATH "C:\\Program Files (x86)\\MPICH2\\bin\\"

int main (int argc, char **argv) {
	// log-file
	FILE *LOGfile;

	int rank, size;
	int flagMPI_Init ;
	int i;
	double worktime;
	
	MPI_Comm MPI_COMM_External; // duplicate of mpi-comm-world

	MPI_Init(&argc, &argv);	// MPI initialization
	// test MPI initialization
	MPI_Initialized (&flagMPI_Init); 
	if (!flagMPI_Init) {
		//create log.dat file
		if ((fopen_s(&LOGfile ,"C:\\Program Files (x86)\\MPICH2\\bin\\log.dat", "w"))) { // TRUE is equal "File was not opened"
			fprintf(stdout, "ERROR: Can't open log-file\n");fflush(stdout);
			exit(1);
		} else {
			time_t tt;
			struct tm * ptm;
			char buf[BUFSIZ];
 			ptm = (struct tm*)malloc(sizeof(struct tm));
			tt = time(NULL);
			localtime_s(ptm, &tt);
			strftime(buf, BUFSIZ, "%d-%m-%Y\nSTART \n", ptm);
			fprintf(LOGfile, buf); fflush(LOGfile);
			free(ptm);
		}
		fprintf(LOGfile, "\\------------------------------\\ \n");fflush(stdout);
		fprintf(LOGfile, "MPI was not initiated.\nProgram MPI_Larkin stoped.");fflush(stdout);
		exit(1);
	}

	//define of number of work processors
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	//define of indexes of work processors
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

	// create a duplicate of MPI_COMM environment
	MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_External);
	
	// create(open) log.dat file

	fprintf(stdout, "a");fflush(stdout);

	if (rank == 0) {
		if ( (fopen_s(&LOGfile ,"C:\\Program Files (x86)\\MPICH2\\bin\\log.dat", "w")) ) { // TRUE is equal "File was not opened"
			int MPIerrorCode = 1; // File was not opened
			fprintf(stdout, "Can't open log-file\n");fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
		} else {
			// start time of calculation is written
			time_t tt;
			struct tm * ptm;
			char buf[BUFSIZ];
			ptm = (struct tm*)malloc(sizeof(struct tm));
			tt = time(NULL);
			localtime_s(ptm, &tt);
			strftime(buf, BUFSIZ, "%d-%m-%Y %H:%M:%S\nSTART \n", ptm);

			fprintf(LOGfile, buf); fflush(LOGfile);
			free(ptm);
		}
	}

	// Number of processes have to be more than 1
	if (rank == 0) {
		if(size <= 1) {
			int MPIErrorCode = 1;
			fprintf(stdout, "Number of processes have to be more than 1\n");fflush(stdout);
			fprintf(LOGfile, "Program was stoped. Number of processes have to be more than 1\n");fflush(LOGfile);
			MPI_Abort(MPI_COMM_WORLD, MPIErrorCode);
		}
	}
	for (i = 2; i < 3; i++) { //to create and use local variable
		double *sndrecvset;
		double *set_parmtrs;
		double result = 0.;		// result of calculation
		int series = 1;			// it`s artefact.
		double nz=0.;
		int position = 0, length; 

		// create new mpi variable to send-recieve procedure
		MPI_Datatype sndrcvparmtrs; // name of new mpi variable
		MPI_Datatype oldtypes[COUNTPARMTRS+1];	// array of old mpi variables
		MPI_Aint offsets[COUNTPARMTRS+1];		// array of distance from start for old mpi variables
		int blklens[COUNTPARMTRS+1];				// length of old mpi variables

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
				MPI_Recv(sndrecvset, COUNTPARMTRS, MPI_DOUBLE, 0, 111, MPI_COMM_External, &status);
			}
		}

		set_parmtrs[0] = sndrecvset[0];	// Vx
		set_parmtrs[1] = sndrecvset[1];	// Vy
		set_parmtrs[2] = sndrecvset[2];	// Vz
		set_parmtrs[3] = sndrecvset[3];	// |V|
		set_parmtrs[4] = sndrecvset[4];	// Te_per
		set_parmtrs[5] = sndrecvset[5];	// Te_par
		set_parmtrs[6] = delta(sndrecvset[4]);	// Delta
		set_parmtrs[7] = tau_perpendicular();	
		set_parmtrs[8] = tau_parallel(sndrecvset[5], sndrecvset[4]);

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
