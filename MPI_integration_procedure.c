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
// �������� ����������


void Integration(double (*function3D)(double, double, double, double*), 
						int series, double *parameter, int myid, int size, int count, 
						MyStr *listINcom, double *subSumm, int *step, double *bounds) { 
	int KEY_exit = 1;  // ���� ������ �� ����� ��������������
	int incount= 0; //���������� �������� �������

	////////////////////////////////////////////////////////////////////////
	/// ���� ������ ���������
	// tag 10 - ���������� ����������(����������) �������������� ������� (sendcount)
	// tag 11 - ����������(����������) �������������� �������
	// tag 12 - ������ �� �������������� �������
	// tag 13 - ��������� �� ��������� ������� ��� ����������
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
	int *INcountcompare;	// ���������� ����� � ������������ ������


	double subBubble;	
	double subCells;
	double stepIterationX, stepIterationY, stepIterationZ; // ���� ��������
	double func;	// ������������� �������� �������
	double lowX;	// ������ ������� �� �
	double lowY;	// ������ ������� �� Y
	double lowZ;	// ������ ������� �� Y
	double X, Y, Z;	// ���������� ��������������
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
		/*���������� ��������� � ������ �����*/							//
	arhblklens = 6;														//
		/*�������� ������� �����, ���������� � ������*/					//
	arhsets = 0;														//
		/*���������� �������� ����� ���������*/							//
    MPI_Type_struct( 1, &arhblklens, &arhsets, &arhtypes, &Strtype);	//
		/*���������� ����� ���������*/									//
	MPI_Type_commit( &Strtype );										//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Char*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_CHAR;												//
        /*���������� ��������� � ������ �����*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*�������� ������� �����, ���������� � ������*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*���������� �������� ����� ���������*/							//
    MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvmsg );		//
		/*���������� ����� ���������*/									//
	MPI_Type_commit( &sndrcvmsg );									//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Int*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_INT;												//
        /*���������� ��������� � ������ �����*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*�������� ������� �����, ���������� � ������*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*���������� �������� ����� ���������*/							//
    MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvcount );		//
		/*���������� ����� ���������*/									//
	MPI_Type_commit( &sndrcvcount );									//
	//////////////////////////////////////////////////////////////////////
	
	///������� ����������� �� ������������(�. 10^5 -> 10^6);
	//////////////////////////////////////
	// ..�������� ������������ ��������...
	// ��������� ������ ��:
	// ���� �������
	list = (MyStr *)malloc(100000*sizeof(MyStr));
	// ������� �� ��������� o� ���������� �������������� �������;;
	req13r = (MPI_Request*)malloc(size*sizeof(MPI_Request));
	// ������� ��������� ��������� o� ���������� �������������� �������;
	status13r = (MPI_Status*)malloc(size*sizeof(MPI_Status));
	// ������ ������ � ����������� ������ �������� ���������
	flagexitIn = (int*)malloc(size*sizeof(int));
	// ������ ������ � �������� �������� o� ���������� �������������� �������;
	flag13r = (int*)malloc(size*sizeof(int));
	// ������ ��������� � ������� ��������� ����� � �������� ���������;
	INcountcompare = (int*)malloc(size*sizeof(int));
	// ������ �������� ��������� ��� ������� ��������;
	IncreaseProcId = (int*)malloc(size*sizeof(int));
	// ������ ������� (��� ����������� ������) ���������;
	LeaveProcId = (int*)malloc(size*sizeof(int));
	// ������ ������ �� ������� �������;
	signallIncom = (int *)malloc(size*sizeof(int));
	// ������ ������ �������� �� ��������� �������;
	signallAsk = (int *)malloc(size*sizeof(int));

	// �������� �������� ������ �������
	memcpy(&list[1], listINcom, count*sizeof(MyStr));

	// ��������� �������� ������� ��������
	*subSumm = 0;

	// ���������� ��� �������� �������� �����
	// ��������������, ��� � ��� ��� �� ������ �������
	for (k = 0; k < size; k++) {
		INcountcompare[k] = 0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	/// ��� ������� �������� ������������� ����� ��������� ������������������ 
	/// �������� ��������� InProcId, ��� ������ ������ �������� ��������(root = 0),
	/// � ��� ����������� �������� ������������ � ��� �� ���� ����������� �� ����������� ������
	/// ������: ������� 4, ����� ����� ��������� 6
	///			4(0), 5(1), 0(2), 1(3), 2(4), 3(5),
	///          ��� ������ ����� ��������������, 
	///			� � ������� - ����� � ����� ���������
	for(i = myid; i < size; i++) {
		IncreaseProcId[i-myid] = i;
	}
	for(i = 0; i < myid; i++) {
		IncreaseProcId[size-myid+i] = i;
	}
	InN = size-1;	// ���������� ��������� � ���������� ������
	InNTemp = InN;	// 
	Inj = 1;		// ���������� ����� �������������� ��������(����),
	Inprev = 0;		// ���������� ������������� �������
					//	 � ������ �� ��������, ����������, �������� �������

	for (k = 1; k < size; k++) {																//
		signallIncom[k] = 1;																		//
	}

	// �������� �������� :
	// 1) �������� � ���������� �������� � ������
	// 2) ���������� �������� �� ���������� �������������� �������
	// 3) ������� �� ��������� �������� � ������
	for (k = 1; k < size; k++) {																//
		flagexitIn[k] = 0;
		flag13r[k] = 0;
		req13r[k] = MPI_REQUEST_NULL;
	}
	
	//
	for (k = 1; k < size; k++) {																//
		signallAsk[k] = 1;																		//
	}																							//

	KEY_exit = 1;			// ���� ������
//	/////////////////////////////////////////////////////////////////////////
//	//..���� � �������� ���� �����,											/
//	// ��� ���������� ��������� ���������� ������� count					/
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
			// ������ ������ �� �����
			*step+=1;
			//���� ��������
			stepIterationX = fabs(list[count].upX - list[count].lowX) / 4.0;
			stepIterationY = fabs(list[count].upY - list[count].lowY) / 4.0;
			stepIterationZ = fabs(list[count].upZ - list[count].lowZ) / 4.0;
			//if (myid == 0) {
			//	fprintf(stdout, "z= %f, y= %f\n", list[count].lowZ, list[count].lowY); fflush(stdout);
			//}
			//������������� �����
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
		// ������� ����������� �������� ����������
			//if (fabs(subCells-subBubble) > PRECISION*fabs(subCells)) { 	
			// subCells --- ��� �������� ������� �����	DE-��������������
			// ���������, DE-�������������� ������ ��������� �������� �������,
			// �� ��������, ��� �� ��������� ��������, ������� �� ������� ������
			// DE-��������������, �������� subCells -> 0
			// ��� ���������, ���������� �������� ��������� ������
			// �����: ������� ������������� ������ ������ ��������� ������ ����������
			if (fabs(subCells-subBubble) > PRECISION) { 			
															// �� ����� ������� ">=" ,���������, 
															//����� subCells = 0, � subBubble = 0, 
															//�� �������� ��������� ����
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

						if (count%2) { // �������� ���������� �������
							countsend = (count-1)/2;
							count = countsend+1;
							sendlist = (MyStr*)malloc(countsend*sizeof(MyStr));
							memcpy(sendlist, &list[count+1], countsend*sizeof(MyStr));
						}
						else { // ������ ���������� �������
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

							//������� ���������� ������������ ������ � ����������
						intbuf[0] = series;			// �������� ������������ ������ �������
						intbuf[1] = countsend;		// ���������� ������� �������
						MPI_Send(intbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);

						//////////////////////////////////////////////////////////////////////
						/*//	  Create a send-recv struct type like MyStr*/				//
						//////////////////////////////////////////////////////////////////////
						   /*Data types*/													//
						oldtypes[0] = MPI_INT;									//
						oldtypes[1] = Strtype;												//
							/*���������� ��������� � ������ �����*/							//
						blklens[0] = 1;														//
						blklens[1] = countsend;												//
							/*�������� ������� �����, ���������� � ������*/					//
						offsets[0] = 0;														//
						MPI_Type_size(MPI_INT, &offsets[1]);						//
							/*���������� �������� ����� ���������*/							//
						MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
							/*���������� ����� ���������*/									//
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
		//	if(count == 0) 	// �. � �������� ����������� �������	///
		///////////////////////////////////////////////////////////////////
	//break;	
	while(!count) {
		int inbuf[2];
		int outbuf[2]; 
		int byebuf;
		int seriestmp10;
///
///�������� ��������� � ����������� ������ ������ ���������
///
		for(i = 1; i <= InNTemp; i++) {
			if(!flagexitIn[i]){ //�������� �������� ������ �� ������ ���������
//				MPI_Irecv(&bye, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				MPI_Irecv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				flagexitIn[i] = 1;
			}
			if(!flag13r[i]) {
				MPI_Test(&req13r[i], &flag13r[i], &status13r[i]);
			}
/// �. ����� ����� ��������� �����-������ ���������, 
/// �� ������ �� ������� � ������ ����� ���������
        //...��������� �������������� �������
		// ��������� ����� ����, ��������� � flag11r
			if (flag13r[i]) {
				int seriestmp13;
				seriestmp13 = byebuf;
				if (seriestmp13 == series) {
	//				fprintf(stdout, "proc %d know proc %d has left program\n", myid, IncreaseProcId[i]);
	//				fflush(stdout);
					InNTemp --;
					// �. ������������� �� �������������� ������� ������ �����,
					// �� ������ �� �������������� ������� ����������
	///
	///�. ��������� � ����������� ������ ��������� �������, �����
	/// ������� ��� �� ������ ������� ���������
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
        /// ������� ������� �� �������� ���������� ��������� ��������, ������������ ���� ������
        ///
///
/// �. ��� ���������� ���������, �� ��������� ������
///
        InN = InNTemp; // ����� ����� ������� ���������;
		if (InN == 0) {//���� ��� ������� ���������, �� �����
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
/// ������� � ��������� ������ �� �������������� �������
///
		if (flag11r) {
			int tmp = series;
			MPI_Send(&tmp, 1, MPI_INT, IncreaseProcId[Inj], 12, MPI_COMM_WORLD);
		}
		// ��������� ������ �� �������������� �������
		// � ������� ��������� �� ���������� �������������� �������
///
/// ������� � ��-����������� �������� ������ �� �������������� �������
///
        //..�������� ������ �� �������������� �������
        if (flag12r) {
            MPI_Irecv(&seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);										////
            flag12r = 0;																									////
        }
        MPI_Test(&req12r, &flag12r, &status12r);
 
///
/// �. ������ �� �������������� ������� �������, 
/// ����� ������� � ������� ���������, ��� �������������� ������� ���
///    //..��������� ���������, ��� �������������� ������� ���. �������� ������� ������� ��������
			// �. ������������� ������� �����, �� ��������� �� ����� ����������
		if (flag12r) {
			//fprintf(stdout, "proc %d has %d series and recv %d tmpser\n", myid, series, seriestmp);
			//fflush(stdout);
			if(series == seriestmp12) { // �. �������� ����������� ��������� ������������� 
												// ��������� ������������ ������ �������, �����...
				outbuf[0] = series;
				outbuf[1] = 0;
				MPI_Send(outbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
				countsend = 0;
				flag12r = 1;
				//MPI_Send(&countsend, 1, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
			} 
        }
///
/// ������� � ��-����������� �������� ��������� 
/// � ���������� ������������� �������������� �������
///     //.. ��������� �������������� �������
        if (flag11r) {
			//MPI_Irecv(&incount, 1, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            MPI_Irecv(inbuf, 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            //MPI_Irecv(&intbuf[1], 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            flag11r = 0;
            flag10r = 0;
        }
        MPI_Test(&req10r, &flag10r, &status10r);

///
///�. �������� ��������� � ����������� �������������� �������, 
/// ����� �������� ��� ���������� �������
///
        if (flag10r) {
			//flag10r = 0;
			flag11r = 1;    // ���� �������� ���������
					//			MPI_Unpack(intbuf, 2*sizeof(int), &position, &seriestmp, 1, MPI_INT, MPI_COMM_WORLD);
			seriestmp10 = inbuf[0];
			if (seriestmp10 == series) 
			{
				//MPI_Unpack(intbuf, 2*sizeof(int), &position, &incount, 1, MPI_INT, MPI_COMM_WORLD);
				incount = inbuf[1];
				//fprintf(stdout, "proc %d recv from proc %d %d tasks\n", myid, status10r.MPI_SOURCE, incount);
				//fflush(stdout);
			/// �. ���������� ����� 0, �� �������� ������� ������������ ��������
			   if (incount == 0) { // ������ ������ ���������
//					flag11r = 1;    // ���� �������� ���������
					i = 0;
					signallIncom[Inj] = 0; // ����������� ���������� ������� ��������,
												//������������ ������ ���������
			/// �������� ���������� ����� �������������� ��������
			/// �. � ���������� ������� ���� ������ �������, ���.����� ����������
			/// ����� ���.����� �������
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
			/// �. ���������� ������ 0, �� ������� � �������� ��������� �
			/// ��������������� ���������
			/// ������������ ������� ������� ��������
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
						/*���������� ��������� � ������ �����*/							//
					blklens[0] = 1;														//
					blklens[1] = incount;												//
						/*�������� ������� �����, ���������� � ������*/					//
					offsets[0] = 0;														//
					MPI_Type_size(MPI_INT, &offsets[1]);						//
						/*���������� �������� ����� ���������*/							//
					MPI_Type_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
						/*���������� ����� ���������*/									//
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
				/// �������� ���������� ����� �������������� ��������
				/// �. � ���������� ������� ���� ������ �������, ���.����� ����������
				/// ����� ���.����� �������
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
        // ��������� ���� �� �������
        // �. ��, �� �������� ��������� � ����� ������ � �����
        // ����� ���������� ����������
///
/// �. ��� �������� ��������� ��������, �� ��� �������� ����������� �������� 
/// ��� ������ �� ���������
/// ������� � ��������� �������� ��������� � ����� ������
/// �������� ������� �� �������������� ������� 
/// � �������� �� �����  ������ ���������
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

                // ������� ������� �� ������ ���������� ��������(�� ����� �������, �� �������� ����������� ���������)
				if (!flag12r) {
                    MPI_Cancel(&req12r);
				}
				if (!flag10r) {
                    MPI_Cancel(&req10r);
				}
                KEY_exit = 0; // ������ ������
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