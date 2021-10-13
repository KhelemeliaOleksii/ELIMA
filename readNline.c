#include"readNline.h"

// функция readNline считывает 
// n-ую строку из входящего файла "in"
// возвращает массив типа double длинною count, 
// содержащийся в считаной строке
// е. в программе произошел сбой
//		возвращается значение 1
//		возращается запись ошибки msg
//	иначе:
//		возвращается значение 0
int readNline (FILE **in, unsigned int n, float *par, int count, char* msg[]) {
	FILE *tempin;
	char *filename = "tempin.dat";
	unsigned int Ntemp = 1;
	int chcount=0, chcount1 = 0, chcount2=0;
	int i;
	char d;

	// строки в файле начинают исчисляться с номера "1"
	if (n == 0) {
		*msg = "ERROR. Program \"readNline\": First line must to have number '1'(not 0)\n";
		return 1;
	}
	// ошибка создания файла
	if (fopen_s(&tempin, filename, "w+")) { 
		*msg = "ERROR. Program \"readNline\": can't create temp file\n";
		return 1;
	}
	
	// пустой файл
	// Set pointer to beginning of file:
	fseek( *in, 0L, SEEK_SET );
	if ((d = getc(*in) ) == EOF) {
		*msg = "ERROR. Program \"readNline\": incoming file is empty\n";
		// закрытие и удаление временного файла
		fclose(tempin);
		remove(filename);
		return 1;
	}

	// поиск и копирование строки под номером n
	// Set pointer to beginning of file:
	fseek( *in, 0L, SEEK_SET );
	while((d = getc(*in) ) != EOF) {
		chcount1++;
		if (Ntemp == n) {
			chcount2++;
			if (d  == '\n') {
				break;
			}
			fputc(d, tempin);
		}
		if (d  == '\n') {
			Ntemp++;
		}
	}

	// входящий файл не содержит строки
	// с указанным номером
	if (Ntemp < n) {
		*msg = "ERROR. Program \"readNline\": incoming file does not have as a line as you ask\n";
		/*fprintf(stdout, "ERROR. Incoming file does not \
						have a line with number %d:\n\
						a program\" readNline\" was tried \
						to read %d-th line in incoming file\n",
						n, n);
		*/
		// закрытие и удаление временного файла
		fclose(tempin);
		remove(filename);
		return 1;
	}

	// Read data back from file:
	// Set pointer to beginning of file:
	/* !!!! как избежать неправильного формата данных*/
	fseek( tempin, 0L, SEEK_SET );
	for (i = 0; i < count; i++) {
		if( !(fscanf_s(tempin, "%e ", &par[i])) ) {
			*msg = "ERROR. Program \"readNline\": error of reading data\n";
			// закрытие и удаление временного файла
			fclose(tempin);
			remove(filename);
			return 1;
		}
	}


	// поиск и удаление строки под номером n
	// Set pointer to beginning of file:
	fseek( *in, 0L, SEEK_SET );
	fseek( tempin, 0L, SEEK_SET );
	chcount = chcount1;
	chcount1-=chcount2;
	while(chcount1--) {   // копирование первой строки
		d = getc(*in);
		fputc(d, tempin);
	}
	fseek( *in, (chcount), SEEK_SET );
	fseek( tempin, (chcount-chcount2), SEEK_SET );
	while((d = getc(*in) ) != EOF) {   // копирование после n-ой строки
		fputc(d, tempin);
	}

	fseek( *in, 0L, SEEK_SET );
	fseek( tempin, 0L, SEEK_SET );
	while((d = getc(tempin) ) != EOF) {
		fputc(d, *in);
	}

	fclose(tempin);
	// Удаление временного файла
	if(remove(filename)) {
		*msg = "WARNING. Program \"readNline\": error of deleting the temp file\n";
		return 0;
	}

	*msg = "Program \"readNline\" finished well. Bye";
	return 0;
}