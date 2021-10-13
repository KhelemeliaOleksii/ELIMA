#include "Boundaries.h"
#include "DoubleExponenta.h"
#include "Function.h"
#include "InputData.h"
//#include "utilprog.h"
#include <stdio.h>
#include <stdlib.h> //malloc
//#include <conio.h>
#include <string.h>	//memcpy
#include <mpi.h>
#include <math.h>
#include <time.h> //
#include "integration_procedure.h"
// точность вычислений


void Integration(double (*function3D)(double, double, double, double*), 
						int series, double *parameter, int myid, int size, int count, 
						MyStr *listINcom, double *subSumm, int *step, double *bounds) { 
	int KEY_exit = 1;  // ключ выхода из цикла интегрирования
	int incount= 0; //количество входящих заданий

	////////////////////////////////////////////////////////////////////////
	/// Таги ключей сообщений
	// tag 10 - количество отсылаемых(получаемых) дополнительных заданий (sendcount)
	// tag 11 - отсылаемые(получаемые) дополнительные задания
	// tag 12 - запрос на дополнительные задания
	// tag 13 - извещение об отсутсвии заданий для выполнения
	MPI_Status *status13r;
	MPI_Status status12r;
	MPI_Status status10s;
	MPI_Status status10r;
	MPI_Status status11r;
	MPI_Request req10r = MPI_REQUEST_NULL;
	MPI_Request req12s = MPI_REQUEST_NULL;
	MPI_Request req12r;
	MPI_Request *req13r;
	int flag12s = 1;
	int flag12r = 1;
	int flag10s;
	int flag10r = 0;
	int flag11r = 1;
	int *flag13r;
	int countsend;
	int k, i, j;			// iterators
	int *flagexitIn;
	int *IncreaseProcId;
	int *LeaveProcId;
	int LeaveN = 0;
	int LeaveJ = 0;
	int bye = 1;
	int InN;
	int InNTemp;
	int Inj;
	int Inprev;
	int flagExitProc = 0;
	int signexitIncom;
	int *signallIncom;
	int *signallAsk;

	MyStr *list;
	int *INcountcompare;	// количество задач в возрастающем потоке


	double subBubble;	
	double subCells;
	double stepIterationX, stepIterationY, stepIterationZ; // шаги итерации
	double func;	// промежуточное значение функции
	double lowX;	// нижняя граница по Х
	double lowY;	// нижняя граница по Y
	double lowZ;	// нижняя граница по Y
	double X, Y, Z;	// переменные интегрирования
	double t, u;	// Simpson`s coefficient for 1D
	double r, s;	// rize range of Simpson`s coefficient to 2D
	double o, p;	// rize range of Simpson`s coefficient to 3D
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MPI_Datatype sndrcvdata, sndrcvcount, sndrcvmsg;
    MPI_Datatype sndrcvcount, sndrcvmsg;
	MPI_Datatype Strtype, arhtypes, oldtypes[2];
    MPI_Aint offsets[2], arhsets;
    int blklens[2], arhblklens;

	/* sender */
    //////////////////////////////////////////////////////////////////////
	/*//       Create a strust MyStr struct type */						//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	arhtypes = MPI_DOUBLE;												//
		/*Количество элементов в каждом блоке*/							//
	arhblklens = 6;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	arhsets = 0;														//
		/*Собственно создание новой структуры*/							//
    MPI_Type_struct( 1, &arhblklens, &arhsets, &arhtypes, &Strtype);	//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &Strtype );										//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Char*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_CHAR;												//
        /*Количество элементов в каждом блоке*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*Собственно создание новой структуры*/							//
    MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvmsg );		//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &sndrcvmsg );									//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Int*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_INT;												//
        /*Количество элементов в каждом блоке*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*Собственно создание новой структуры*/							//
    MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvcount );		//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &sndrcvcount );									//
	//////////////////////////////////////////////////////////////////////
	
	///сделать ограничения на переполнение(е. 10^5 -> 10^6);
	//////////////////////////////////////
	// ..Создание динамических массивов...
	// Выделение памяти на:
	// Стек заданий
	list = (MyStr *)malloc(100000*sizeof(MyStr));
	// очереди на сообщения oб отсутствии дополнительных заданий;;
	req13r = (MPI_Request*)malloc(size*sizeof(MPI_Request));
	// статусы полученых сообщений oб отсутствии дополнительных заданий;
	status13r = (MPI_Status*)malloc(size*sizeof(MPI_Status));
	// Список флагов о прикращении работы соседних процессов
	flagexitIn = (int*)malloc(size*sizeof(int));
	// Список флагов о входящих сигналах oб отсутствии дополнительных заданий;
	flag13r = (int*)malloc(size*sizeof(int));
	// Список доступных в прошлом количеств задач у соседних процессов;
	INcountcompare = (int*)malloc(size*sizeof(int));
	// Список соседних процессов для каждого процесса;
	IncreaseProcId = (int*)malloc(size*sizeof(int));
	// Список рабочих (или завершивших работу) процессов;
	LeaveProcId = (int*)malloc(size*sizeof(int));
	// Список флагов на вхоящие задания;
	signallIncom = (int *)malloc(size*sizeof(int));
	// Список флагов запросов на дополните задания;
	signallAsk = (int *)malloc(size*sizeof(int));

	// копируем входящий список заданий
	memcpy(&list[1], listINcom, count*sizeof(MyStr));

	// начальная подсумма каждого процесса
	*subSumm = 0;

	// Изначально все соседние процессы равны
	// Предполагается, что у них нет ни одного задания
	for (k = 0; k < size; k++) {
		INcountcompare[k] = 0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	/// для каждого процесса выстраивается новое множество последовательности 
	/// соседних процессов InProcId, где данный процес является корневым(root = 0),
	/// а все последующие процессы выстаиваются в ряд по мере возрастания их порядкового номера
	/// пример: процесс 4, общее число процессов 6
	///			4(0), 5(1), 0(2), 1(3), 2(4), 3(5),
	///          где первый номер действительный, 
	///			а в скобках - номер в новом множестве
	for(i = myid; i < size; i++) {
		IncreaseProcId[i-myid] = i;
	}
	for(i = 0; i < myid; i++) {
		IncreaseProcId[size-myid+i] = i;
	}
	InN = size-1;	// количество элементов в восходящем списке
	InNTemp = InN;	// 
	Inj = 1;		// порядковый номер запрашиваемого процесса(восх),
	Inprev = 0;		// предыдущий запрашиваемый процесс
					//	 в начале им является, собственно, корневой процесс

	for (k = 1; k < size; k++) {																//
		signallIncom[k] = 1;																		//
	}

	// Обнуляем счетчики :
	// 1) сигналов о получениии сигналов о выходе
	// 2) полученных сигналов об отсутствии дополнительных заданий
	// 3) очереди на получение сигналов о выходе
	for (k = 1; k < size; k++) {																//
		flagexitIn[k] = 0;
		flag13r[k] = 0;
		req13r[k] = MPI_REQUEST_NULL;
	}
	
	//
	for (k = 1; k < size; k++) {																//
		signallAsk[k] = 1;																		//
	}																							//

	KEY_exit = 1;			// ключ выхода
//	/////////////////////////////////////////////////////////////////////////
//	//..вход в основное тело цикла,											/
//	// где проводится обработка полученных заданий count					/
//	/////////////////////////////////////////////////////////////////////////
//	fprintf(stdout, "@\n");
//	fflush(stdout);
	for(i=1; i <= InN; i++) {
		MPI_Iprobe(IncreaseProcId[i], 10, MPI_COMM_WORLD, &flag10s, &status10s);
		if(flag10s) {
			int intbuf[2];
			MPI_Recv(intbuf, 1, sndrcvcount, IncreaseProcId[i], 10, MPI_COMM_WORLD, &status10s);			
		}
		MPI_Iprobe(IncreaseProcId[i], 10, MPI_COMM_WORLD, &flag10s, &status10s);
		if(flag10s) {
			int intbuf[2];
			MPI_Recv(intbuf, 1, sndrcvcount, IncreaseProcId[i], 10, MPI_COMM_WORLD, &status10s);			
		}
	}
//	fprintf(stdout, " myid %d list[0].lowX = %f, list[0].upX = %f \n", myid, listINcom[0].lowX, listINcom[0].upX); fflush(stdout);
//	fprintf(stdout, " myid %d list[0].upX = %f, list[0].upX = %f \n", myid, UPX, UPY); fflush(stdout);
	//for (j = 0; j <=4; j++) {
	//	for (i = 0; i <=4; i++) {
	//		
	//	}
	//}


	while (KEY_exit > 0) {
		int seriestmp12;
		while(count) {
			// сделан проход по циклу
			*step+=1;
			//шаги итераций
			stepIterationX = fabs(list[count].upX - list[count].lowX) / 4.0;
			stepIterationY = fabs(list[count].upY - list[count].lowY) / 4.0;
			stepIterationZ = fabs(list[count].upZ - list[count].lowZ) / 4.0;
			//if (myid == 0) {
			//	fprintf(stdout, "z= %f, y= %f\n", list[count].lowZ, list[count].lowY); fflush(stdout);
			//}
			//промежуточные суммы
			subBubble = 0.0; subCells = 0.0;
			//double o, p;	// Simpson`s coefficient for 1D
			for (i=0; i<=4; i++) {
				Z = list[count].lowZ + i*stepIterationZ; //Z
				if ((i == 0)||(i == 4)) {
					o=1;
					p=1;
				}
				if ((i == 2)) {
					o=2;
					p=4;
				}
				if ((i == 1)||(i == 3)) {
					o=4;
					p=0;
				}
				//	double r, s;	// rize range of Simpson`s coefficient to 2D
				for (j=0; j<=4; j++) {
					Y = list[count].lowY + j*stepIterationY ; //Y
					if ((j == 0)||(j==4)) {
						r=1;
						s=1;
					}
					if ((j == 2)) {
						r=2;
						s=4;
					}
					if ((j == 1)||(j==3)) {
						r=4;
						s=0;
					}

					//	double t, u;	// rize range of Simpson`s coefficient to 3D
					for (k=0; k<=4; k++) {
						X = list[count].lowX + k*stepIterationX ;	//X
						if ((k == 0)||(k==4)) {
							t=1;
							u=1;
						}
						if ((k == 2)) {
							t=2;
							u=4;
						}
						if ((k == 1)||(k==3)) {
							t=4;
							u=0;
						}
//						func = function(X, Y, Z, parameter);						
						func = DEtransformFunction3D(function3D, X, Y, Z, parameter, bounds);						
						subBubble += func*u*s*p;
						subCells += func*t*r*o;
					}
				}
			}
			subBubble = subBubble*stepIterationX*stepIterationY*stepIterationZ*8./27.;
			subCells = subCells*stepIterationX*stepIterationY*stepIterationZ/27.;//			subCells = subCells*stepIterationX/3;
		// УСЛОВИЕ НЕОБХОДИМОЙ ТОЧНОСТИ ВЫЧИСЛЕНИЙ
			//if (fabs(subCells-subBubble) > PRECISION*fabs(subCells)) { 	
			// subCells --- это значение функции после	DE-преобразования
			// Поскольку, DE-преобразование сильно подавляет значение функции,
			// то получаем, что на некоторых участках, далеких от вершини купола
			// DE-преобразования, значение subCells -> 0
			// Как результат, получаемая точность превышает нужную
			// Вывод: условие относительной ошибки сильно замедляет работу проограммы
			if (fabs(subCells-subBubble) > PRECISION) { 			
															// не нужно ставить ">=" ,поскольку, 
															//когда subCells = 0, и subBubble = 0, 
															//то получаем замкнутый цикл
				lowX = list[count].lowX;
				lowY = list[count].lowY;
				lowZ = list[count].lowZ;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {				
						for (k = 0; k < 4; k++) {
							list[count].lowX = lowX + (stepIterationX)*(k);
							list[count].upX = lowX + (stepIterationX)*(k+1);
							list[count].lowY = lowY + stepIterationY*(j);
							list[count].upY = lowY + (stepIterationY)*(j+1);
							list[count].lowZ = lowZ + stepIterationZ*(i);
							list[count].upZ = lowZ + (stepIterationZ)*(i+1);
							count++;
						}
					}
				}
				count--;
			} else 
			{
				*subSumm += subCells;
				count -= 1;
			}
			if (count >= 5) {
				if (flag12r) {
					MPI_Irecv(&seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);
				}
				MPI_Test(&req12r, &flag12r, &status12r);
				if (flag12r) {
					//fprintf(stdout, "proc %d recv from %d ask for tasks with series%d\n", myid, status12r.MPI_SOURCE, seriestmp);
					//fflush(stdout);
					if (seriestmp12 == series) {
						MPI_Datatype sndrcvdata;
						double *bufdatasend;
						int position = 0;
						int intbuf[2];
						int mpiintsize, mpidoublesize, strtypesize;
						MyStr *sendlist;

						countsend = 0;

						if (count%2) { // непарное количество заданий
							countsend = (count-1)/2;
							count = countsend+1;
							sendlist = (MyStr*)malloc(countsend*sizeof(MyStr));
							memcpy(sendlist, &list[count+1], countsend*sizeof(MyStr));
						}
						else { // парное количество заданий
							countsend = (count)/2;																	////
							count = countsend;																		////
							sendlist = (MyStr*)malloc(countsend*sizeof(MyStr));
							memcpy(sendlist, &list[count+1], countsend*sizeof(MyStr));								////
						}
						if (countsend == 0) {
							int err = 001;
							fprintf(stdout, "ERROR_%d\n", err);fflush(stdout);
							MPI_Abort(MPI_COMM_WORLD, err);
						}

							//отсылка количества передаваемых данных с серийником
						intbuf[0] = series;			// серийник выполняемого общего задания
						intbuf[1] = countsend;		// количество частных заданий
						MPI_Send(intbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);

						//////////////////////////////////////////////////////////////////////
						/*//	  Create a send-recv struct type like MyStr*/				//
						//////////////////////////////////////////////////////////////////////
						   /*Data types*/													//
						oldtypes[0] = MPI_INT;									//
						oldtypes[1] = Strtype;												//
							/*Количество элементов в каждом блоке*/							//
						blklens[0] = 1;														//
						blklens[1] = countsend;												//
							/*Смещение каждого блока, измеряемые в байтах*/					//
						offsets[0] = 0;														//
						MPI_Type_size(MPI_INT, &offsets[1]);						//
							/*Собственно создание новой структуры*/							//
						MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
							/*Реестрация новой структуры*/									//
						MPI_Type_commit( &sndrcvdata );										//
						//////////////////////////////////////////////////////////////////////
						MPI_Type_size(MPI_INT, &mpiintsize);
						MPI_Type_size(Strtype, &strtypesize);

						bufdatasend = (double*)malloc(countsend*strtypesize+mpiintsize);
						position = 0;
						MPI_Pack(&intbuf[0], 1, MPI_INT, bufdatasend, countsend*strtypesize+mpiintsize, &position, MPI_COMM_WORLD);
						MPI_Pack(sendlist, countsend, Strtype, bufdatasend, countsend*strtypesize+mpiintsize, &position, MPI_COMM_WORLD);
						//fprintf(stdout, "q");fflush(stdout);
						MPI_Send(bufdatasend, 1, sndrcvdata, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);
						//MPI_Send(sendlist, 6*countsend, Strtype, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);
			//			fprintf(stdout, "proc %d send to proc %d %d(%d) tasks\n", myid, status12r.MPI_SOURCE, countsend, intbuf[1]);
			//			fflush(stdout);
						free(bufdatasend);
						free(sendlist);
						MPI_Type_free(&sndrcvdata);
					}
				}
			}
		}
		///////////////////////////////////////////////////////////////////
		//	if(count == 0) 	// е. у процесса закончились задания	///
		///////////////////////////////////////////////////////////////////
	//break;	
	while(!count) {
		int inbuf[2];
		int outbuf[2]; 
		int byebuf;
		int seriestmp10;
///
///получить сообщения о прекращении работы других процессов
///
		for(i = 1; i <= InNTemp; i++) {
			if(!flagexitIn[i]){ //единажды принятие выхода от других процессов
//				MPI_Irecv(&bye, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				MPI_Irecv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				flagexitIn[i] = 1;
			}
			if(!flag13r[i]) {
				MPI_Test(&req13r[i], &flag13r[i], &status13r[i]);
			}
/// е. здесь будут вноситься какие-нибудь изменения, 
/// не забыть их зделать в нижней части программы
        //...Попросить дополнительные задания
		// поставить новый флаг, связанный с flag11r
			if (flag13r[i]) {
				int seriestmp13;
				seriestmp13 = byebuf;
				if (seriestmp13 == series) {
	//				fprintf(stdout, "proc %d know proc %d has left program\n", myid, IncreaseProcId[i]);
	//				fflush(stdout);
					InNTemp --;
					// е. запрашиваемый на дополнительные задания процес вышел,
					// то запрос на дополнительные задания отменяется
	///
	///е. сообщение о прекращении работы процессом принято, тогда
	/// изымаем его из списка рабочих процессов
	///     
					if (Inj == i) {
						flag11r = 1;
						if (!flag10r) {
							MPI_Cancel(&req10r);
							flag10r = 1;
						}
					}
					if(Inj > i) {
						Inj --;
					}
					for(j = i; j <= InNTemp; j++){
						IncreaseProcId[j] = IncreaseProcId[j+1];
						signallIncom[j] = signallIncom[j+1]; 
					}
					for (j = i; j<=InNTemp; j++) {
						req13r[j] = req13r[j+1];
						flag13r[j] = flag13r[j+1] ;
					}
				} else {
					flagexitIn[i] = 0;
					flag13r[i]=0;
//					req13r[i] = 0;
				}
			}
        }
        ///
        /// Сделать вставку на удаление отосланных сообщений процессу, завершившему свою работу
        ///
///
/// е. нет работающих процессов, то завершаем работу
///
        InN = InNTemp; // новое число рабочих процессов;
		if (InN == 0) {//если нет рабочих процессов, то выход
			if (!flag12r) {
				MPI_Cancel(&req12r);
			}
			if (!flag10r) {
		        MPI_Cancel(&req10r);
			}
			//fprintf(stdout, "proc %d leaving caused by no more processes\n", myid);	fflush(stdout);
			KEY_exit  = 0;
			break;
		}
///
/// СДЕЛАТЬ и ВЫПОЛНИТЬ запрос на дополнительные задания
///
		if (flag11r) {
			int tmp = series;
			MPI_Send(&tmp, 1, MPI_INT, IncreaseProcId[Inj], 12, MPI_COMM_WORLD);
		}
		// получение запрос на дополнительные задания
		// и отсылка сообщения об отсутствии дополнительных заданий
///
/// СДЕЛАТЬ и по-возможности ПОЛУЧИТЬ запрос на дополнительные задания
///
        //..Получить запрос на дополнительные задания
        if (flag12r) {
            MPI_Irecv(&seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);										////
            flag12r = 0;																									////
        }
        MPI_Test(&req12r, &flag12r, &status12r);
 
///
/// е. запрос на дополнительные задания ПОЛУЧЕН, 
/// тогда СДЕЛАТЬ и ПОСЛАТЬ сообщение, что дополнительных заданий НЕТ
///    //..Оптравить сообщение, что дополнительных заданий нет. Обнулить счетчик данного процесса
			// е. запрашиваемый процесс вышел, то сообщение не нужно отправлять
		if (flag12r) {
			//fprintf(stdout, "proc %d has %d series and recv %d tmpser\n", myid, series, seriestmp);
			//fflush(stdout);
			if(series == seriestmp12) { // е. серийник полученного сообщения соответствует 
												// серийнику выполняемого общего задания, тогда...
				outbuf[0] = series;
				outbuf[1] = 0;
				MPI_Send(outbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
				countsend = 0;
				flag12r = 1;
				//MPI_Send(&countsend, 1, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
			} 
        }
///
/// СДЕЛАТЬ и по-возможности ПОЛУЧИТЬ сообщение 
/// о КОЛИЧЕСТВЕ запрашиваемых дополнительных заданий
///     //.. получение дополнительных заданий
        if (flag11r) {
			//MPI_Irecv(&incount, 1, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            MPI_Irecv(inbuf, 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            //MPI_Irecv(&intbuf[1], 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            flag11r = 0;
            flag10r = 0;
        }
        MPI_Test(&req10r, &flag10r, &status10r);

///
///е. ПОЛУЧЕНО сообщение с количеством дополнительных заданий, 
/// тогда получить это количество заданий
///
        if (flag10r) {
			//flag10r = 0;
			flag11r = 1;    // флаг принятия сообщения
					//			MPI_Unpack(intbuf, 2*sizeof(int), &position, &seriestmp, 1, MPI_INT, MPI_COMM_WORLD);
			seriestmp10 = inbuf[0];
			if (seriestmp10 == series) 
			{
				//MPI_Unpack(intbuf, 2*sizeof(int), &position, &incount, 1, MPI_INT, MPI_COMM_WORLD);
				incount = inbuf[1];
				//fprintf(stdout, "proc %d recv from proc %d %d tasks\n", myid, status10r.MPI_SOURCE, incount);
				//fflush(stdout);
			/// е. количество равно 0, то обнулить счетчик отправившего процесса
			   if (incount == 0) { // пришло пустое сообщение
//					flag11r = 1;    // флаг принятия сообщения
					i = 0;
					signallIncom[Inj] = 0; // обозначение отсутствия заданий процесса,
												//отправившего пустое сообщение
			/// изменить порядковый номер запрашиваемого процесса
			/// е. в предыдущем запросе было меньше заданий, пор.номер возрастает
			/// иначе пор.номер спадает
					if ((Inj <= (InN)) && (Inj > 0)) {
						if (incount >= INcountcompare[Inprev] ) {
							Inprev = Inj;
							if (Inj == (InN)) {
								Inj = 1;
							} else {
								Inj++;
							}
						} else {
							Inprev = Inj;
							if (Inj == 1) {
								Inj = (InN);
							} else {
								Inj--;
							}
						}
					}
			/// е. количество больше 0, то СДЕЛАТЬ и ПОЛУЧИТЬ сообщение с
			/// дополнительными заданиями
			/// активировать счетчик заданий процесса
				} else if (incount > 0) {
					MPI_Datatype sndrcvdata;
					int seriestmp11;
					int position = 0;
					double *bufdatarecv;
					//////////////////////////////////////////////////////////////////////
					/*//	  Create a send-recv struct type like MyStr*/				//
					//////////////////////////////////////////////////////////////////////
					   /*Data types*/													//
					oldtypes[0] = MPI_INT;									//
					oldtypes[1] = Strtype;												//
						/*Количество элементов в каждом блоке*/							//
					blklens[0] = 1;														//
					blklens[1] = incount;												//
						/*Смещение каждого блока, измеряемые в байтах*/					//
					offsets[0] = 0;														//
					MPI_Type_size(MPI_INT, &offsets[1]);						//
						/*Собственно создание новой структуры*/							//
					MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
						/*Реестрация новой структуры*/									//
					MPI_Type_commit( &sndrcvdata );										//
					//////////////////////////////////////////////////////////////////////
					bufdatarecv = (double*)malloc(incount*sizeof(MyStr)+sizeof(int));
//					fprintf(stdout, "proc %d recv from %d list of %d(%d) tasks\n", myid, status10r.MPI_SOURCE, incount, inbuf[1]);
//					fflush(stdout);

//					MPI_Recv(recvlist, 6*incount, Strtype, IncreaseProcId[Inj], 11, MPI_COMM_WORLD, &status11r);
					MPI_Recv(bufdatarecv, 1, sndrcvdata , IncreaseProcId[Inj], 11, MPI_COMM_WORLD, &status11r);

					MPI_Unpack(bufdatarecv, incount*sizeof(MyStr)+sizeof(int), &position, &seriestmp11, 1, MPI_INT, MPI_COMM_WORLD);   
					if (seriestmp11 == series) {
						MyStr *recvlist;
						recvlist =(MyStr*)malloc(incount*sizeof(MyStr));
						MPI_Unpack(bufdatarecv, incount*sizeof(MyStr)+sizeof(int), &position, recvlist, incount, Strtype, MPI_COMM_WORLD);   

//						flag11r = 1;
						memcpy(&list[1], recvlist, incount*sizeof(MyStr));										//
						signallIncom[Inj] = 1;
						count = incount;
				/// изменить порядковый номер запрашиваемого процесса
				/// е. в предыдущем запросе было меньше заданий, пор.номер возрастает
				/// иначе пор.номер спадает
						if ((Inj <= (InN)) && (Inj > 0)) {
							if (incount >= INcountcompare[Inprev] ) {
								Inprev = Inj;
								if (Inj == (InN)) {
									Inj = 1;
								} else {
									Inj++;
								}
							} else {
								Inprev = Inj;
								if (Inj == 1) {
									Inj = (InN);
								} else {
									Inj--;
								}
							}
						}
						free(recvlist);
					}
					MPI_Type_free(&sndrcvdata);
					free(bufdatarecv);
					break;
				} else if (incount < 0) {
					fprintf(stdout, "\nproc %d !!!!!SystemError!!!!\n", myid);
					fflush(stdout);
					MPI_Finalize();
				}
				INcountcompare[Inprev] = incount;	
			} 
//			else {
//				flag11r = 1;
//			}
        }
        // проверить всех ли опросил
        // е. да, то отослать сообщения о своем выходе и выйти
        // иначе продолжить опрашивать
///
/// е. все счетчики процессов обнулены, то это является достаточным условием 
/// для выхода из программы
/// СДЕЛАТЬ и ВЫПОЛНИТЬ рассылку сообщения о своем выходе
/// ОТМЕНИТЬ запросы на дополнительные задания 
/// и проверку на выход  других процессов
///
        if (!flagExitProc) {
            signexitIncom = 0;
			for (k = 1; k <= InN; k++) {		//
                signexitIncom += signallIncom[k];			//
            }
            if (!signexitIncom) {
                for (i = 1; i <= InNTemp; i++) {
					byebuf = series; 
					MPI_Send(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD);
					MPI_Cancel(&req13r[i]);
				}
				//fprintf(stdout, "proc %d leaving caused by no more tasks\n", myid);fflush(stdout);

                // сделать вставку на отмену сделланных запросов(на новые задания, на проверку действующих процессов)
				if (!flag12r) {
                    MPI_Cancel(&req12r);
				}
				if (!flag10r) {
                    MPI_Cancel(&req10r);
				}
                KEY_exit = 0; // сигнал выхода
                break;
            }
		}
	}


	}
	free(status13r);
	free(flagexitIn);
	free(req13r);
	free(flag13r);
	free(INcountcompare);
	free(IncreaseProcId);
	free(list);
	free(signallIncom);
	free(signallAsk);

	fprintf(stdout, "*\n");
	fflush(stdout);
}