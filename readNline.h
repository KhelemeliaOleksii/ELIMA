#ifndef READNLINE_H
#define READNLINE_H

#include<conio.h>
#include<io.h>
#include <stdlib.h> // string
#include<stdio.h>

#define MSGLENGTH 500 
// функция readNline считывает 
// n-ую строку из входящего файла "in"
// возвращает массив типа double длинною count, 
// содержащийся в считаной строке
// е. в программе произошел сбой
//		возвращается значение 1
//		возращается запись ошибки msg
//	иначе:
//		возвращается значение 0
int readNline (FILE **in, unsigned int number, float *par, int count, char *msg[]);
#endif